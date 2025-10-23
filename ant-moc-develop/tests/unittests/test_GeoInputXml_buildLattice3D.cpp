/**
 * @file test_GeoInputXml_buildLattice3D.cpp
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
class test_GeoInputXml_buildLattice3D : public GeometryTestHarness {

protected:

  void SetUp() {
    // Add materials
    StringVec names = {
      "0", "1", "2", "Beer", "Coffee", "Milk"
    };

    for (auto &name : names)
      _mat_input->setMaterial(name, new Material(0, name.c_str()));
  }

};


TEST_F(test_GeoInputXml_buildLattice3D, recLatticeLayout2x2x3) {
  ROOT_ONLY();

  // Test input
  std::stringstream ss;
  ss << "<lattice id='100' name='RecLattice' type='Rectangle'>"
     << "  <row>              1 1 </row>"
     << "  <row>              1 2 </row>"
     << "  <row repeat='2'> 1 1 </row>"
     << "  <row repeat='2'> 2 2 </row>"
     << "  <widths x='[3.50364] * 2' y='[3.50364] * 2' z='[1.0] * 3'/>"
     << "  <format>"
     << "    <surface id='1' type='x-plane' coeffs='1.75182'/>"
     << "    <surface id='2' type='x-plane' coeffs='-1.75182'/>"
     << "    <surface id='3' type='y-plane' coeffs='1.75182'/>"
     << "    <surface id='4' type='y-plane' coeffs='-1.75182'/>"
     << "    <surface id='5' type='z-plane' coeffs='0.5'/>"
     << "    <surface id='6' type='z-plane' coeffs='-0.5'/>"
     << "  </format>"
     << "  <universe id='1' name='pin-cells'>"
     << "    <cell id='1' material='1' name='large' region='-7'>"
     << "      <surface id='7' type='z-cylinder' coeffs='0 0 0.3320'/>"
     << "    </cell>"
     << "  </universe>"
     << "  <universe id='2' name='pin-cells'>"
     << "    <cell id='2' material='1' name='small' region='+8'>"
     << "      <surface id='8' type='z-cylinder' coeffs='0 0 0.2867'/>"
     << "    </cell>"
     << "  </universe>"
     << "</lattice>";

  const std::string &xml = ss.str();
  _doc.Parse(xml.c_str());

  auto lattice = buildLattice(_doc.FirstChildElement("lattice"));

  // Test each attribute of the Lattice
  IntPyVec u_ids = {1, 2};

  IntPyVec3D layout = {
    {{1, 1},
     {1, 2}},
    {{1, 1},
     {1, 1}},
    {{2, 2},
     {2, 2}}
  };

  Widths3D widths_oracles = {
    {3.50364, 3.50364},
    {3.50364, 3.50364},
    {1.0, 1.0, 1.0}
  };

  testRecLattice(lattice, 100, "RecLattice", &layout, &widths_oracles, &u_ids);

  // Additional tests
  // 1. test removeUniverses, which removes all references to a
  //    universe from the layout
  u_ids = {2};
  layout = {
    {{-1, -1},
     {-1,  2}},
    {{-1, -1},
     {-1, -1}},
    {{ 2,  2},
     { 2,  2}}
  };
  auto u = lattice->getUniverse(0, 0, 2); // remove universe 1
  lattice->removeUniverse(u);
  testRecLattice(lattice, 100, "RecLattice", &layout, &widths_oracles, &u_ids);

  // 2. test updateUniverses, which update a universe at specified
  //    position
  u_ids = {1, 2};
  layout = {
    {{ 1,  1},
     { 1,  2}},
    {{-1, -1},
     {-1, -1}},
    {{ 2,  2},
     { 2,  2}}
  };
  lattice->updateUniverse(0, 0, 2, u);
  lattice->updateUniverse(0, 1, 2, u);
  lattice->updateUniverse(1, 1, 2, u);
  testRecLattice(lattice, 100, "RecLattice", &layout, &widths_oracles, &u_ids);

  delete lattice;
}


TEST_F(test_GeoInputXml_buildLattice3D, recLatticeExertedRegion) {
  ROOT_ONLY();

  const std::string test =
    "<cell id='1' fill='101' name='outer cell' region='-1 +2 -3 +4 -5 +6'>"
    "  <surface id='1' type='x-plane' coeffs='3.50364'/>"
    "  <surface id='2' type='x-plane' coeffs='-3.50364'/>"
    "  <surface id='3' type='y-plane' coeffs='3.50364'/>"
    "  <surface id='4' type='y-plane' coeffs='-3.50364'/>"
    "  <surface id='5' type='z-plane' coeffs='1.5'/>"
    "  <surface id='6' type='z-plane' coeffs='-1.5'/>"
    "  <lattice id='101' name='RecLattice' type='Rectangle'>"
    "    <row repeat='40'> 1 1 </row>"
    "    <widths x='[3.50364] * 2' y='[3.50364] * 2' z='[0.15] * 20'/>"
    "    <format region='-7 +8 -9 +10'>"
    "      <surface id='7'  type='x-plane' coeffs='1.75182'/>"
    "      <surface id='8'  type='x-plane' coeffs='-1.75182'/>"
    "      <surface id='9'  type='y-plane' coeffs='1.75182'/>"
    "      <surface id='10' type='y-plane' coeffs='-1.75182'/>"
    "    </format>"
    "    <universe id='1' name='pin-cell'>"
    "      <cell id='2' material='2' name='inner cell' region='+11'>"
    "        <surface id='11' type='z-cylinder' coeffs='0 0 0.3320'/>"
    "      </cell>"
    "    </universe>"
    "  </lattice>"
    "</cell>";

  _doc.Parse(test.c_str());

  std::string region = "-1 +2 -3 +4 -5 +6";
  IntPyVec sf_ids = {1, 2, 3, 4, 5, 6};

  auto outer_cell = buildCell(_doc.FirstChildElement("cell"));
  testCell(outer_cell, 1, "outer cell", false, "RecLattice", &region, &sf_ids);

  auto lattice = dynamic_cast<RecLattice*>(outer_cell->getFillUniverse());

  // Test each attribute of the Lattice
  IntPyVec u_ids = {1};

  IntPyVec3D layout = {
    {
     {1, 1},
     {1, 1}
    }
  };

  layout *= 20; // z-sections

  Widths3D widths_oracles = {
    {3.50364, 3.50364},
    {3.50364, 3.50364},
    WidthVec{0.15} * 20
  };

  testRecLattice(lattice, 101, "RecLattice", &layout, &widths_oracles, &u_ids);

  /* Test inner cells */
  auto univs = lattice->getUniqueUniverses();
  for(auto &univ_pair : univs) {
    Universe *u = univ_pair.second;
    for (auto &cell_pair : u->getCells()) {
      Cell *c = cell_pair.second;
      EXPECT_EQ   (2,            c->getInputId());
      //EXPECT_EQ   (2,           c->getFillMaterial()->getId());
      EXPECT_STREQ("inner cell", c->getName());
      ASSERT_EQ   (5,            c->getNumSurfaces());
    }
  }

  delete lattice;
}


TEST_F(test_GeoInputXml_buildLattice3D, refinedRecLatticeExertedWidths) {
  ROOT_ONLY();

  const std::string test =
    "<cell id='1' fill='101' name='outer cell' region='-1 +2 -3 +4 -5 +6'>"
    "  <surface id='1' type='x-plane' coeffs='3.50364'/>"
    "  <surface id='2' type='x-plane' coeffs='-3.50364'/>"
    "  <surface id='3' type='y-plane' coeffs='3.50364'/>"
    "  <surface id='4' type='y-plane' coeffs='-3.50364'/>"
    "  <surface id='5' type='z-plane' coeffs='1.5'/>"
    "  <surface id='6' type='z-plane' coeffs='-1.5'/>"
    "  <lattice id='101' name='RecLattice' type='Rectangle'>"
    "    <row repeat='4'> 1 1 </row>"
    "    <widths x='[3.50364] * 2' y='[3.50364] * 2' z='[0.15] * 2'/>"
    "    <refines x='5' y='5' z='5'/>"
    "    <universe id='1' name='pin-cell'>"
    "      <cell id='2' material='2' name='inner cell' region='+11'>"
    "        <surface id='11' type='z-cylinder' coeffs='0 0 0.3320'/>"
    "      </cell>"
    "    </universe>"
    "  </lattice>"
    "</cell>";
  _doc.Parse(test.c_str());

  constexpr int refines = 5;

  std::string region = "-1 +2 -3 +4 -5 +6";
  IntPyVec sf_ids = {1, 2, 3, 4, 5, 6};

  auto outer_cell = buildCell(_doc.FirstChildElement("cell"));
  testCell(outer_cell, 1, "outer cell", false, "RecLattice", &region, &sf_ids);

  auto lattice = dynamic_cast<RecLattice*>(outer_cell->getFillUniverse());

  // Test each attribute of the Lattice
  IntPyVec u_ids = {1};

  IntPyVec row = IntPyVec{1} * 2 * refines;
  IntPyVec2D plane = IntPyVec2D{row} * 2 * refines;
  IntPyVec3D layout = IntPyVec3D{plane} * 2 * refines;

  Widths3D widths_oracles = {
    WidthVec{3.50364 / refines} * 2 * refines,
    WidthVec{3.50364 / refines} * 2 * refines,
    WidthVec{0.15 / refines} * 2 * refines,
  };

  testRecLattice(lattice, 101, "RecLattice", &layout, &widths_oracles, &u_ids);

  /* Test inner cells */
  auto univs = lattice->getUniqueUniverses();
  for(auto &univ_pair : univs) {
    Universe *u = univ_pair.second;
    for (auto &cell_pair : u->getCells()) {
      Cell *c = cell_pair.second;
      EXPECT_EQ   (2,           c->getInputId());
      //EXPECT_EQ   (2,           c->getFillMaterial()->getId());
      EXPECT_STREQ("inner cell",    c->getName());
      ASSERT_EQ   (1,           c->getNumSurfaces());
    }
  }

  delete lattice;
}


TEST_F(test_GeoInputXml_buildLattice3D, axialRefinedRecLatticeSubdividedCells) {
  // Test input
  std::stringstream ss;
  ss << "<cell id='1' fill='101' name='outer cell'>"
     << "  <lattice id='101' name='RecLattice' type='Rectangle'>"
     << "    <row repeat='4'> 1 1 </row>"
     << "    <widths x='[3.50364] * 2' y='[3.50364] * 2' z='[0.15] * 2'/>"
     << "    <refines z='5'/>"
     << "    <universe id='1' name='pin-cell'>"
     << "      <cell id='2' material='2' name='inner cell' region='-1' sectors='4' rings='2'>"
     << "        <surface id='1' type='z-cylinder' coeffs='0 0 0.3320'/>"
     << "      </cell>"
     << "    </universe>"
     << "  </lattice>"
     << "</cell>";
  constexpr int refines = 5;

  const std::string &xml = ss.str();
  _doc.Parse(xml.c_str());

  auto outer_cell = buildCell(_doc.FirstChildElement("cell"));

  ASSERT_TRUE(outer_cell->getFillUniverse());
  auto lattice = dynamic_cast<RecLattice*>(outer_cell->getFillUniverse());

  EXPECT_STREQ("outer cell", outer_cell->getName());

  // Test each attribute of the Lattice
  IntPyVec u_ids = {1};

  IntPyVec2D plane = {
    {1, 1},
    {1, 1}
  };

  IntPyVec3D layout = IntPyVec3D{plane} * 2 * refines;

  Widths3D widths_oracles = {
    WidthVec{3.50364} * 2,
    WidthVec{3.50364} * 2,
    WidthVec{0.15 / refines} * 2 * refines,
  };

  testRecLattice(lattice, 101, "RecLattice", &layout, &widths_oracles, &u_ids);

  /* Test inner cells */
  auto univs = lattice->getUniqueUniverses();
  double radius = 3.50364 / std::sqrt(2);
  for(auto &univ_pair : univs) {
    Universe *u = univ_pair.second;
    u->subdivideCells(radius);
    for (auto &cell_pair : u->getCells()) {
      Cell *c = cell_pair.second;
      ASSERT_EQ   (2,           c->getInputId());
      ASSERT_STREQ("inner cell",c->getName());
      ASSERT_EQ   (1,           c->getNumSurfaces());
      // Loop over subdivided cells
      Universe *u2 = c->getFillUniverse();
      ASSERT_TRUE(u2);
      ASSERT_EQ(8, u2->getCells().size());
    }
  }

  delete lattice;
}


TEST_F(test_GeoInputXml_buildLattice3D, xyzRefinedRecLattice) {
  ROOT_ONLY();

  // Test input
  std::stringstream ss;
  ss << "<global>"
     << "  <surface id='1' type='z-cylinder' coeffs='0 0 0.3320'/>"
     << "</global>"
     << "<lattice id='101' name='RecLattice' type='Rectangle'>"
     << "  <row> 3 3 </row>"
     << "  <row> 3 2 </row>"

     << "  <row> 1 1 </row>"
     << "  <row> 1 2 </row>"
     << "  <widths x='[3.50364] * 2' y='[3.50364] * 2' z='[0.15] * 2'/>"
     << "  <refines x='2' y='2' z='2'/>"
     << "  <universe id='1' name='beer'>"
     << "    <cell id='1' material='Beer' name='beer_cell' region='-1'/>"
     << "  </universe>"
     << "  <universe id='2' name='coffee'>"
     << "    <cell id='1' material='Coffee' name='coffee_cell' region='-1'/>"
     << "  </universe>"
     << "  <universe id='3' name='milk'>"
     << "    <cell id='1' material='Milk' name='milk_cell' region='-1'/>"
     << "  </universe>"
     << "</lattice>";
  constexpr int refines = 2;

  const std::string &xml = ss.str();
  _doc.Parse(xml.c_str());

  _geo_input->buildGlobalPrimitives(_doc.FirstChildElement("global"));
  auto lattice = buildLattice(_doc.FirstChildElement("lattice"));

  // Test each attribute of the Lattice
  IntPyVec u_ids = {1, 2, 3};

  // Test lattice layout
  IntPyVec3D layout = {
    {
     {3, 3, 3, 3},
     {3, 3, 3, 3},
     {3, 3, 2, 2},
     {3, 3, 2, 2}
    },
    {
     {3, 3, 3, 3},
     {3, 3, 3, 3},
     {3, 3, 2, 2},
     {3, 3, 2, 2}
    },
    {
     {1, 1, 1, 1},
     {1, 1, 1, 1},
     {1, 1, 2, 2},
     {1, 1, 2, 2},
    },
    {
     {1, 1, 1, 1},
     {1, 1, 1, 1},
     {1, 1, 2, 2},
     {1, 1, 2, 2},
    }
  };

  Widths3D widths_oracles = {
    WidthVec{3.50364/refines} * 2 * refines,
    WidthVec{3.50364/refines} * 2 * refines,
    WidthVec{0.15   /refines} * 2 * refines,
  };

  testRecLattice(lattice, 101, "RecLattice", &layout, &widths_oracles, &u_ids);

  delete lattice;
}


TEST_F(test_GeoInputXml_buildLattice3D, xyzGlobalRefinedRecLattice) {
  ROOT_ONLY();

  // Test input
  std::stringstream ss;
  ss << "<lattice id='100' name='RecLattice' type='Rectangle'>"
     << "  <row> 1 </row>"
     << "  <widths x='3.50364' y='3.50364' z='0.15'/>"
     << "  <refines x='1' y='1' z='1'/>"
     << "  <universe id='1' name='beer'>"
     << "    <cell id='1' material='Beer' name='beer_cell'/>"
     << "  </universe>"
     << "</lattice>";
  const std::string &xml = ss.str();
  _doc.Parse(xml.c_str());

  _geo_input->setGlobalRefines({2, 2, 1});
  auto lattice = buildLattice(_doc.FirstChildElement("lattice"));

  IntPyVec u_ids = {1};

  IntPyVec3D layout = {
    {
     {1, 1},
     {1, 1},
    },
  };

  Widths3D widths_oracles = {
    WidthVec{3.50364/2} * 2,
    WidthVec{3.50364/2} * 2,
    {0.15}
  };

  testRecLattice(lattice, 100, "RecLattice", &layout, &widths_oracles, &u_ids);

  delete lattice;
}


TEST_F(test_GeoInputXml_buildLattice3D, partialGlobalRefinedRecLattice) {
  ROOT_ONLY();

  // Test input
  std::stringstream ss;
  ss << "<lattice id='100' name='RecLattice' type='Rectangle'>"
     << "  <row> 1 </row>"
     << "  <widths x='3.50364' y='3.50364' z='0.15'/>"
     << "  <refines x='1' y='1' z='1'/>"
     << "  <universe id='1' name='beer'>"
     << "    <cell id='1' material='Beer' name='beer_cell'/>"
     << "  </universe>"
     << "</lattice>";
  const std::string &xml = ss.str();
  _doc.Parse(xml.c_str());

  _geo_input->setGlobalRefines({0, 2, 0});
  auto lattice = buildLattice(_doc.FirstChildElement("lattice"));

  IntPyVec u_ids = {1};

  IntPyVec3D layout = {
    {
     {1},
     {1},
    },
  };

  Widths3D widths_oracles = {
    {3.50364},
    {1.75182, 1.75182},
    {0.15}
  };

  testRecLattice(lattice, 100, "RecLattice", &layout, &widths_oracles, &u_ids);

  delete lattice;
}



}/* namespace */
