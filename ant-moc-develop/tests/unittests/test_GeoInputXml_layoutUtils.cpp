/**
 * @file test_GeoInputXml_layoutUtils.cpp
 * @brief Test Geometry reading and building process
 * @date Aug 26, 2019
 * @author An Wang, USTB (wangan@xs.ustb.edu.cn)
 */

#include "testing/geometry_testharness.h"
#include "antmoc/enum_types.h"
#include "antmoc/Surface.h"
#include "antmoc/Cell.h"
#include "antmoc/Universe.h"
#include "antmoc/Lattice.h"

using namespace antmoc;

namespace {

// Test fixture
class test_GeoInputXml_layoutUtils : public GeometryTestHarness {

protected:

  void SetUp() {
    // Add materials
    StringVec names = {
      "1", "Beer", "Coffee", "Milk"
    };

    for (auto &name : names)
      _mat_input->setMaterial(name, new Material(0, name.c_str()));
  }

};


#ifndef ENABLE_MPI_

TEST_F(test_GeoInputXml_layoutUtils, illFormedLayout) {

  // Case 1: # of rows < # of lattice cells
  const std::string test1 =
    "<lattice name='RecLattice' type='Rectangle'>"
    "  <row> 1 1 </row>"
    "  <widths x='[3.50364] * 2' y='[3.50364] * 2' z='[1.0]'/>"
    "  <universe id='1'> <cell id='1' material='1'/> </universe>"
    "</lattice>";
  _doc.Parse(test1.c_str());
  EXPECT_ANY_THROW(buildLattice(_doc.FirstChildElement("lattice")));

  // Case 2: # of items in a row < # of lattice cells
  const std::string test2 =
    "<lattice name='RecLattice' type='Rectangle'>"
    "  <row> 1 1 </row>"
    "  <row> 1 </row>"
    "  <widths x='[3.50364] * 2' y='[3.50364] * 2' z='[1.0]'/>"
    "  <universe id='1'> <cell id='1' material='1'/> </universe>"
    "</lattice>";
  _doc.Parse(test2.c_str());
  EXPECT_ANY_THROW(buildLattice(_doc.FirstChildElement("lattice")));

  // Case 3: vector expressions in 'row's
  const std::string test3 =
    "<lattice name='RecLattice' type='Rectangle'>"
    "  <row> [1]*2 </row>"
    "  <row> [1]*2 </row>"
    "  <widths x='[3.50364] * 2' y='[3.50364] * 2' z='[1.0]'/>"
    "  <universe id='1'> <cell id='1' material='1'/> </universe>"
    "</lattice>";
  _doc.Parse(test3.c_str());
  EXPECT_ANY_THROW(buildLattice(_doc.FirstChildElement("lattice")));

  // Case 4: ill-formed vector expressions
  const std::string test4 =
    "<lattice name='RecLattice' type='Rectangle'>"
    "  <universes> "
    "    [1] * 2   "
    "    [1] * 2   "
    "  </universes>"
    "  <widths x='[3.50364] * 2' y='[3.50364] * 2' z='[1.0]'/>"
    "  <universe id='1'> <cell id='1' material='1'/> </universe>"
    "</lattice>";
  _doc.Parse(test4.c_str());
  EXPECT_ANY_THROW(buildLattice(_doc.FirstChildElement("lattice")));
}

#endif


TEST_F(test_GeoInputXml_layoutUtils, rowExtendedRecLattice) {
  ROOT_ONLY();

  const std::string test =
    "<global>"
    "  <surface id='1' type='z-cylinder' coeffs='0 0 0.3320'/>"
    "</global>"
    "<lattice id='101' name='RecLatticeWithRows' type='Rectangle'>"
    "  <row> 3 3 2R 2 1R </row>"
    "  <row> 2 3R 3 1 </row>"

    "  <row> 1 2R 2 2R</row>"
    "  <row> 1 1 1 2 3 3 </row>"
    "  <widths x='[3.14] * 6' y='[3.14] * 2' z='[0.15] * 2'/>"
    "  <universe id='1' name='beer'>"
    "    <cell id='1' material='Beer' name='beer_cell' region='-1'/>"
    "  </universe>"
    "  <universe id='2' name='coffee'>"
    "    <cell id='1' material='Coffee' name='coffee_cell' region='-1'/>"
    "  </universe>"
    "  <universe id='3' name='milk'>"
    "    <cell id='1' material='Milk' name='milk_cell' region='-1'/>"
    "  </universe>"
    "</lattice>";

  _doc.Parse(test.c_str());

  _geo_input->buildGlobalPrimitives(_doc.FirstChildElement("global"));
  auto lattice = buildLattice(_doc.FirstChildElement("lattice"));

  // Test each attribute of the Lattice
  IntPyVec u_ids = {1, 2, 3};

  // Test lattice layout
  IntPyVec3D layout = {
    {
     {3, 3, 3, 3, 2, 2},
     {2, 2, 2, 2, 3, 1},
    },
    {
     {1, 1, 1, 2, 2, 2},
     {1, 1, 1, 2, 3, 3},
    },
  };

  Widths3D widths_oracles = {
    WidthVec{3.14} * 6,
    WidthVec{3.14} * 2,
    WidthVec{0.15} * 2,
  };

  testRecLattice(lattice, 101, "RecLatticeWithRows", &layout, &widths_oracles, &u_ids);

  delete lattice;
}


TEST_F(test_GeoInputXml_layoutUtils, layoutWithVectorExpressions) {
  ROOT_ONLY();

  const std::string test =
    "<global>"
    "  <surface id='1' type='z-cylinder' coeffs='0 0 0.3320'/>"
    "</global>"
    "<lattice id='101' name='RecLatticeWithRows' type='Rectangle'>"
    "  <universes>"
    "   [3]*4 [2]*2"
    "   [2]*4 3 1  "

    "   [1]*3 [2]*3"
    "   [1]*3 2 [3]*2"
    "  </universes>"

    // lower priority
    "  <row> 3 </row>"
    "  <row> 2 </row>"
    "  <row> 1 </row>"
    "  <row> 1 </row>"

    "  <widths x='[3.14] * 6' y='[3.14] * 2' z='[0.15] * 2'/>"
    "  <universe id='1' name='beer'>"
    "    <cell id='1' material='Beer' name='beer_cell' region='-1'/>"
    "  </universe>"
    "  <universe id='2' name='coffee'>"
    "    <cell id='1' material='Coffee' name='coffee_cell' region='-1'/>"
    "  </universe>"
    "  <universe id='3' name='milk'>"
    "    <cell id='1' material='Milk' name='milk_cell' region='-1'/>"
    "  </universe>"
    "</lattice>";

  _doc.Parse(test.c_str());

  _geo_input->buildGlobalPrimitives(_doc.FirstChildElement("global"));
  auto lattice = buildLattice(_doc.FirstChildElement("lattice"));

  // Test each attribute of the Lattice
  IntPyVec u_ids = {1, 2, 3};

  // Test lattice layout
  IntPyVec3D layout = {
    {
     {3, 3, 3, 3, 2, 2},
     {2, 2, 2, 2, 3, 1},
    },
    {
     {1, 1, 1, 2, 2, 2},
     {1, 1, 1, 2, 3, 3},
    },
  };

  Widths3D widths_oracles = {
    WidthVec{3.14} * 6,
    WidthVec{3.14} * 2,
    WidthVec{0.15} * 2,
  };

  testRecLattice(lattice, 101, "RecLatticeWithRows", &layout, &widths_oracles, &u_ids);

  delete lattice;
}

}/* namespace */
