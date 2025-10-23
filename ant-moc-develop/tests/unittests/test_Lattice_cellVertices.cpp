/// \brief Test lattice cell vertices computing
/// \date March 10, 2020
/// \author An Wang, USTB (wangan@xs.ustb.edu.cn)

#include "testing/test_utils.h"
#include "testing/geometry_testutils.h"
#include "antmoc/Lattice.h"
#include "antmoc/math_utils.h"

using namespace antmoc;

namespace {

// Test fixture
class test_Lattice_cellVertices : public testing::Test { };


/// \brief Compute lattice cell vertices for a uniform RecLattice
TEST_F(test_Lattice_cellVertices, uniformRecLattice2x2x2) {
  ROOT_ONLY();

  // Construct a 1.2*1.2*1.2 lattice
  RecLattice lattice;
  lattice.setWidth(0.6, 0.6, 0.6);
  lattice.setNumX(2);
  lattice.setNumY(2);
  lattice.setNumZ(2);
  lattice.computeSizes();

  // Lattice cells are ordered in x-y-z order
  // The first cell (0,0,0)
  const int index = 0;
  PyVector<Point> oracles = {
    {   0,    0, -0.6},
    {-0.6,    0, -0.6},
    {-0.6, -0.6, -0.6},
    {   0, -0.6, -0.6},

    {   0,    0,    0},
    {-0.6,    0,    0},
    {-0.6, -0.6,    0},
    {   0, -0.6,    0},
  };

  auto vertices = lattice.getLatticeCellVertices(index);

  for (size_t i = 0; i < oracles.size(); ++i) {
    Point &oracle = oracles[i];
    Point &vertex = vertices[i];

    EXPECT_EQ(vertex, oracle)
      << "LatticeCell: " << index << '\n'
      << "Vertex: " << vertex.toString() << '\n'
      << "Oracle: " << oracle.toString()
      << std::endl;
  }
}


/// \brief Test a y-orientated HexLattice with nr = 2 and nz = 2
TEST_F(test_Lattice_cellVertices, yOrientatedHexLattice2x2) {
  ROOT_ONLY();

  // 2x2 hexagon lattice with radial pitch 1
  HexLattice lattice;
  double width_r = 1.;
  DoubleVec widths_z = {0.3, 0.9};

  lattice.setOrientation("y");
  lattice.setWidths(width_r, widths_z);
  lattice.setNumR(2);
  lattice.setNumZ(2);
  lattice.computeSizes();

  // Lattice cells are ordered in x-y-z order
  //      cell
  // cell      cell
  //      cell
  // cell      cell
  //      cell
  double sin60 = std::sin(M_PI/3);

  // The center cell (1,1,0)
  const int index = 4;
  PyVector<Point> oracles = {
    {  1/(2*sin60),    0, -0.6 },
    {  1/(4*sin60),  0.5, -0.6 },
    { -1/(4*sin60),  0.5, -0.6 },
    { -1/(2*sin60),    0, -0.6 },
    { -1/(4*sin60), -0.5, -0.6 },
    {  1/(4*sin60), -0.5, -0.6 },

    {  1/(2*sin60),    0, -0.3 },
    {  1/(4*sin60),  0.5, -0.3 },
    { -1/(4*sin60),  0.5, -0.3 },
    { -1/(2*sin60),    0, -0.3 },
    { -1/(4*sin60), -0.5, -0.3 },
    {  1/(4*sin60), -0.5, -0.3 },
  };

  auto vertices = lattice.getLatticeCellVertices(index);

  for (size_t i = 0; i < oracles.size(); ++i) {
    Point &oracle = oracles[i];
    Point &vertex = vertices[i];

    EXPECT_EQ(vertex, oracle)
      << "LatticeCell: " << index << '\n'
      << "Vertex: " << vertex.toString() << '\n'
      << "Oracle: " << oracle.toString()
      << std::endl;
  }
}


/// \brief Test an x-orientated HexLattice with nr = 2 and nz = 2
TEST_F(test_Lattice_cellVertices, xOrientatedHexLattice2x2) {
  ROOT_ONLY();

  // 2x2 hexagon lattice with radial pitch 1
  HexLattice lattice;
  double width_r = 1.;
  DoubleVec widths_z = {0.3, 0.9};

  lattice.setOrientation("x");
  lattice.setWidths(width_r, widths_z);
  lattice.setNumR(2);
  lattice.setNumZ(2);
  lattice.computeSizes();

  // Lattice cells are ordered in x-y-z order
  //      cell      cell
  // cell      cell      cell
  //      cell      cell
  double sin60 = std::sin(M_PI/3);

  // The center cell (1,1,0)
  const int index = 4;
  PyVector<Point> oracles = {
    {  0.5,  1/(4*sin60), -0.6 },
    {    0,  1/(2*sin60), -0.6 },
    { -0.5,  1/(4*sin60), -0.6 },
    { -0.5, -1/(4*sin60), -0.6 },
    {    0, -1/(2*sin60), -0.6 },
    {  0.5, -1/(4*sin60), -0.6 },
           
    {  0.5,  1/(4*sin60), -0.3 },
    {    0,  1/(2*sin60), -0.3 },
    { -0.5,  1/(4*sin60), -0.3 },
    { -0.5, -1/(4*sin60), -0.3 },
    {    0, -1/(2*sin60), -0.3 },
    {  0.5, -1/(4*sin60), -0.3 },
  };

  auto vertices = lattice.getLatticeCellVertices(index);

  for (size_t i = 0; i < oracles.size(); ++i) {
    Point &oracle = oracles[i];
    Point &vertex = vertices[i];

    EXPECT_EQ(vertex, oracle)
      << "LatticeCell: " << index << '\n'
      << "Vertex: " << vertex.toString() << '\n'
      << "Oracle: " << oracle.toString()
      << std::endl;
  }
}

} // namespace
