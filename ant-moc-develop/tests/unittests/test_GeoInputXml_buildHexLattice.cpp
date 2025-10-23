/**
 * @file test_GeoInputXml_buildHexLattice.cpp
 * @brief Test Geometry reading and building process
 * @date Sep 7, 2019
 * @author An Wang, USTB (wangan@xs.ustb.edu.cn)
 */

#include "testing/geometry_testharness.h"
#include "antmoc/enum_types.h"
#include "antmoc/lattice_utils.h"

using namespace antmoc;

namespace {


// Test fixture
class test_GeoInputXml_buildHexLattice : public GeometryTestHarness {

protected:

  void SetUp() {
    // Add materials
    StringVec names = {
      "1",
    };

    for (auto &name : names)
      _mat_input->setMaterial(name, new Material(0, name.c_str()));
  }

};


/// \brief Test y-orientated HexLattice
/// \details The pitches of the lattice is (2, 2). Methods like containsPoints,
///          minSurfaceDist are tested with simple cases.
TEST_F(test_GeoInputXml_buildHexLattice, refinedYHexLatticeLayout2x2) {
  ROOT_ONLY();

  // Test input
  std::string test =
    "<lattice id='100' name='HexLattice' type='Hexagon' n_rings='2'>"
    "  <widths r='2.0' z='[1.0]*2'/>"
    "  <refines z='2'/>'"
    "  <universes>"
    "       1     "
    "     1   1   "
    "       2     "
    "     1   1   "
    "       1     "
    "       2     "
    "     2   2   "
    "       1     "
    "     2   2   "
    "       2     "
    "  </universes>"
    "  <universe id='1' name='pin-cells'>"
    "    <cell id='1' material='1' name='large' region='-1'>"
    "      <surface id='1' type='z-cylinder' coeffs='0 0 0.1'/>"
    "    </cell>"
    "  </universe>"
    "  <universe id='2' name='pin-cells'>"
    "    <cell id='2' material='1' name='small' region='+1'>"
    "      <surface id='1' type='z-cylinder' coeffs='0 0 0.2'/>"
    "    </cell>"
    "  </universe>"
    "</lattice>";

  _doc.Parse(test.c_str());

  auto lattice = buildLattice(_doc.FirstChildElement("lattice"));

  // Test each attribute of the Lattice
  IntPyVec u_ids = {1, 2};

  IntPyVec layout = { // refined
    2, 2, 2, 1, 2, 2, 2,
    2, 2, 2, 1, 2, 2, 2,
    1, 1, 1, 2, 1, 1, 1,
    1, 1, 1, 2, 1, 1, 1,
  };

  double hex_widths_r = 2.0;

  WidthVec hex_widths_z = {0.5, 0.5, 0.5, 0.5};

  int num_r = 2;

  testHexLattice(lattice, 100, "HexLattice", &layout, &hex_widths_r,
                 &hex_widths_z, &num_r, &u_ids);

  // Test if the lattice contains some points
  std::vector<Point> in_points = {
    {0, 0, 0},
    {0, 3.0, 1.0},  // top
    {0, 2.99, 1.0},  // near top
    {0, -2.99, -1.0},  // near bottom
    //{0, -3.0, -1.0},  // bottom outside the box
    {1.5*std::sqrt(3), 1.5, 0},  // upper right corner
  };

  for (auto &p : in_points)
    EXPECT_TRUE(lattice->containsPoint(&p))
      << p.toString();

  // Along z-axis
  std::vector<Point> out_points = {
    {0, 3.0, 1.01},  // out of top
    {0, -3.0, -1.01},  // out bottom
  };

  for (auto &p : out_points)
    EXPECT_FALSE(lattice->containsPoint(&p))
      << p.toString();

  // On x-y plane
  out_points = {
    {0, 3.01, 0},  // out of top
    {0, -3.01, 0},  // out bottom
    {1.5*std::sqrt(3), 1.5 + 0.01, 0},  // out of upper right corner
  };

  for (auto &p : out_points)
    EXPECT_FALSE(lattice->containsPoint(&p))
      << p.toString();

  // Test minSurfaceDist on x-y plane
  // There are 6 directions in total. The following snippet
  // calculates the distances between given points and sides
  // of the lattice cells the points reside.
  const double pi = std::acos(-1);
  DblPyVec directions = {
    pi/6, pi/2, 5*pi/6, -pi/6, -pi/2, -5*pi/6,
  };

  double cos30 = std::cos(pi/6);
  double sin30 = std::sin(pi/6);
  std::vector<Point> points = {
    // centers
    {0, 0, 0},
    {0,  2.0, 0},
    {0, -2.0, 0},
    { 2.0*cos30,  2.0*sin30, 0},
    { 2.0*cos30, -2.0*sin30, 0},
    {-2.0*cos30,  2.0*sin30, 0},
    {-2.0*cos30, -2.0*sin30, 0},
    // close to boundaries
    {0,  2.99, 0},
    {0, -2.99, 0},
    //{ 2.99*cos30,  2.99*sin30, 0},
    //{ 2.99*cos30, -2.99*sin30, 0},
    //{-2.99*cos30,  2.99*sin30, 0},
    //{-2.99*cos30, -2.99*sin30, 0},
  };


  DblPyVec2D distances = DblPyVec2D{ DblPyVec{1.0, 1.0, 1.0, 1.0, 1.0, 1.0} } * 7;
  double tan30 = std::tan(pi/6);
  double d1 = 0.01/sin30;
  double d2 = (0.01/tan30 + 1.0*tan30)*cos30*2 - d1;
  DblPyVec2D distances_b = {
    {d1, 0.01, d1, d2, 1.99, d2},
    {d2, 1.99, d2, d1, 0.01, d1},

  };
  distances += distances_b;

  for (size_t i = 0; i < points.size(); ++i)
    for (size_t a = 0; a < directions.size(); ++a) {
      double d = lattice->minSurfaceDist(&points[i], directions[a]);
      EXPECT_NEAR(distances[i][a], d, 1E-15)  // FIXME
        << points[i].toString() << '\n'
        << "  dir[" << a << "] = " << directions[a]
        << "  distance = " << d;
    }

  delete lattice;
}


/// \brief Test x-orientated HexLattice
/// \details The pitches of the lattice is (2, 2). Methods like containsPoints,
///          minSurfaceDist are tested with simple cases.
TEST_F(test_GeoInputXml_buildHexLattice, refinedXHexLatticeLayout2x2) {
  ROOT_ONLY();

  // Test input
  std::string test =
    "<lattice id='100' name='HexLattice' type='Hexagon' n_rings='2'"
    "  orientation='x'>"
    "  <widths r='2.0' z='[1.0]*2'/>"
    "  <refines z='2'/>'"
    "  <universes>"
    "     1   1   "
    "   1   2   1 "
    "     1   1   "
    "     2   2   "
    "   2   1   2 "
    "     2   2   "
    "  </universes>"
    "  <universe id='1' name='pin-cells'>"
    "    <cell id='1' material='1' name='large' region='-1'>"
    "      <surface id='1' type='z-cylinder' coeffs='0 0 0.1'/>"
    "    </cell>"
    "  </universe>"
    "  <universe id='2' name='pin-cells'>"
    "    <cell id='2' material='1' name='small' region='+1'>"
    "      <surface id='1' type='z-cylinder' coeffs='0 0 0.2'/>"
    "    </cell>"
    "  </universe>"
    "</lattice>";

  _doc.Parse(test.c_str());

  auto lattice = buildLattice(_doc.FirstChildElement("lattice"));

  // Test each attribute of the Lattice
  IntPyVec u_ids = {1, 2};

  IntPyVec layout = { // refined
    2, 2, 2, 1, 2, 2, 2,
    2, 2, 2, 1, 2, 2, 2,
    1, 1, 1, 2, 1, 1, 1,
    1, 1, 1, 2, 1, 1, 1,
  };

  double hex_widths_r = 2.0;  ///< pitch

  WidthVec hex_widths_z = {0.5, 0.5, 0.5, 0.5};   ///< nonuniform z-widths

  int num_r = 2;  ///< radial tiles

  testHexLattice(lattice, 100, "HexLattice", &layout, &hex_widths_r,
                 &hex_widths_z, &num_r, &u_ids);

  // Test if the lattice contains some points
  std::vector<Point> in_points = {
    {0, 0, 0},
    {3.0, 0, 1.0},  // rightmost
    {2.99, 0, 1.0},  // near rightmost
    {-2.99, 0, -1.0},  // near leftmost
    //{0, -3.0, -1.0},  // bottom outside the box
    {1.5, 1.5*std::sqrt(3), 0},  // upper right corner
  };

  for (auto &p : in_points)
    EXPECT_TRUE(lattice->containsPoint(&p))
      << p.toString();

  // Along z-axis
  std::vector<Point> out_points = {
    {3.0, 0, 1.01},  // out of top
    {-3.0, 0, -1.01},  // out bottom
  };

  for (auto &p : out_points)
    EXPECT_FALSE(lattice->containsPoint(&p))
      << p.toString();

  // On x-y plane
  out_points = {
    {3.01, 0, 0},  // out of top
    {-3.01, 0, 0},  // out bottom
    {1.5 + 0.01, 1.5*std::sqrt(3), 0},  // out of upper right corner
  };

  for (auto &p : out_points)
    EXPECT_FALSE(lattice->containsPoint(&p))
      << p.toString();

  // Test minSurfaceDist on x-y plane
  // There are 6 directions in total. The following snippet
  // calculates the distances between given points and sides
  // of the lattice cells the points reside.
  const double pi = std::acos(-1);
  DblPyVec directions = {
    0, pi/3, 2*pi/3, pi, -2*pi/3, -pi/3
  };

  double cos30 = std::cos(pi/6);
  double sin30 = std::sin(pi/6);

  // The first 7 points are centers of lattice cells
  std::vector<Point> points = {
    // centers
    {0, 0, 0},
    { 2.0, 0, 0},
    {-2.0, 0, 0},
    { 2.0*sin30,  2.0*cos30, 0},
    {-2.0*sin30,  2.0*cos30, 0},
    { 2.0*sin30, -2.0*cos30, 0},
    {-2.0*sin30, -2.0*cos30, 0},
    // close to boundaries
    { 2.99, 0, 0},
    {-2.99, 0, 0},
    //{ 2.99*sin30,  2.99*cos30, 0},
    //{-2.99*sin30,  2.99*cos30, 0},
    //{ 2.99*sin30, -2.99*cos30, 0},
    //{-2.99*sin30, -2.99*cos30, 0},
  };

  // Distances from centers to lattice sides
  DblPyVec2D distances = DblPyVec2D{ DblPyVec{1.0, 1.0, 1.0, 1.0, 1.0, 1.0} } * 7;
  double tan30 = std::tan(pi/6);
  double d1 = 0.01/sin30;
  double d2 = (0.01/tan30 + 1.0*tan30)*cos30*2 - d1;
  // Distances from near-boundary points to lattice sides
  DblPyVec2D distances_b = {
    {0.01, d1, d2, 1.99, d2, d1},
    {1.99, d2, d1, 0.01, d1, d2},
  };
  distances += distances_b;

  for (size_t i = 0; i < points.size(); ++i)
    for (size_t a = 0; a < directions.size(); ++a) {
      double d = lattice->minSurfaceDist(&points[i], directions[a]);
      EXPECT_NEAR(distances[i][a], d, 1E-15)  // FIXME
        << points[i].toString() << '\n'
        << "  dir[" << a << "] = " << directions[a]
        << "  distance = " << d;
    }

  delete lattice;
}
}/* namespace */
