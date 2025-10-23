/**
 * @file test_GeoInputXml_buildLattice2D.cpp
 * @brief Test Geometry reading and building process
 * @date April 20, 2019
 * @author An Wang, USTB (wangan@xs.ustb.edu.cn)
 */

#include "testing/geometry_testharness.h"
#include "antmoc/enum_types.h"
#include "antmoc/Surface.h"
#include "antmoc/Cell.h"
#include "antmoc/Universe.h"
#include "antmoc/Lattice.h"

#include <vector>

using namespace antmoc;

namespace {

// Test fixture
class test_GeoInputXml_buildLattice2D : public GeometryTestHarness {

protected:

  void SetUp() {
    // Add materials
    StringVec names = {
      "0", "1", "2"
    };

    for (auto &name : names)
      _mat_input->setMaterial(name, new Material(0, name.c_str()));
  }

};

#ifndef ENABLE_MPI_

TEST_F(test_GeoInputXml_buildLattice2D, buildLatticeRobust) {
  StringVec tests = {
    "<lattice name='assembly'/>",
    "<lattice name='assembly' type='sexangular'/>",
    "<lattice name='assembly' type='Mobius'/>",

    "<lattice name='assembly' type='rectangle'>"
    "  <name> assembly </name>" // duplicated
    "  <widths x='1.0' y='1.0' z='1.0'/>"
    "  <row> 1 </row>"
    "  <universe id='1'>"
    "    <cell id='1' material='0' name='cell 1'/>"
    "  </universe>"
    "</lattice>",

    "<lattice name='assembly' type='rectangle'>"
    "  <widths x='1.0' y='1.0' z='1.0'>"
    "    <x> 1.0 </x>" // duplicated
    "  </widths>"
    "  <row> 1 </row>"
    "  <universe id='1'>"
    "    <cell id='1' material='0' name='cell 1'/>"
    "  </universe>"
    "</lattice>",
  };


  for (auto &s : tests) {
    _doc.Parse(s.c_str());
    EXPECT_ANY_THROW(buildLattice(_doc.FirstChildElement("lattice")));
  }

  // These will pass
  StringVec pass = {
    "<lattice name='assembly' type='rectangle'>"
    "  <widths x='1.0' y='1.0' z='1.0'/>"
    "  <row> 1 </row>"
    "  <universe id='1'>"
    "    <cell id='1' material='0' name='cell 1'/>"
    "  </universe>"
    "</lattice>",
  };

  for (auto &s : pass) {
    _doc.Parse(s.c_str());
    EXPECT_NO_THROW(buildLattice(_doc.FirstChildElement("lattice")));
  }
}
#endif


TEST_F(test_GeoInputXml_buildLattice2D, extrudedRecLatticeLayout2x2x3) {
  ROOT_ONLY();

  // Test input
  std::stringstream ss;
  ss << "<lattice name='RecLattice' type='Rectangle' axial='extruded'>"
     << "  <id>  100 </id>"
     << "  <row id='1'> 1 1 </row>"
     << "  <row id='0'> 2 2 </row>"
     << "  <widths>"
     << "    <x> [3.50364] * 2 </x>"
     << "    <y> [3.50364] * 2 </y>"
     << "    <z> [1.0] * 3     </z>"
     << "  </widths>"
     << "  <format region='-1 +2 -3 +4 -5 +6'>"
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

  /* Test the layout of the Lattice
   * Note that the input is actually
   * 1 1
   * 2 2   Section 2
   * 1 1
   * 2 2   Section 1
   * 1 1
   * 2 2   Section 0 */

  IntPyVec2D plane = {
    {1, 1},
    {2, 2}
  };

  IntPyVec3D layout = IntPyVec3D{plane} * 3;

  Widths3D widths_oracles = {
    {3.50364, 3.50364},
    {3.50364, 3.50364},
    {1.0, 1.0, 1.0}
  };

  testRecLattice(lattice, 100, "RecLattice", &layout, &widths_oracles, &u_ids);

  delete lattice;
}


TEST_F(test_GeoInputXml_buildLattice2D, extrudedRecLatticeExertedRegion) {
  // Test input
  std::stringstream ss;
  ss << "<cell id='1' fill='101' name='outer cell' region='-1 +2 -3 +4 -5 +6'>"
     << "  <surface id='1' type='x-plane' coeffs='3.50364'/>"
     << "  <surface id='2' type='x-plane' coeffs='-3.50364'/>"
     << "  <surface id='3' type='y-plane' coeffs='3.50364'/>"
     << "  <surface id='4' type='y-plane' coeffs='-3.50364'/>"
     << "  <surface id='5' type='z-plane' coeffs='1.5'/>"
     << "  <surface id='6' type='z-plane' coeffs='-1.5'/>"
     << "  <lattice id='101' name='RecLattice' type='Rectangle' axial='extruded'>"
     << "    <row repeat='2'> 1 1 </row>"
     << "    <widths x='[3.50364] * 2' y='[3.50364] * 2' z='[0.15] * 20'/>"
     << "    <format region='-7 +8 -9 +10'>"
     << "      <surface id='7'  type='x-plane' coeffs='1.75182'/>"
     << "      <surface id='8'  type='x-plane' coeffs='-1.75182'/>"
     << "      <surface id='9'  type='y-plane' coeffs='1.75182'/>"
     << "      <surface id='10' type='y-plane' coeffs='-1.75182'/>"
     << "    </format>"
     << "    <universe id='1' name='pin-cell'>"
     << "      <cell id='2' material='2' name='inner cell' region='-11'>"
     << "        <surface id='11' type='z-cylinder' coeffs='0 0 0.3320'/>"
     << "      </cell>"
     << "    </universe>"
     << "  </lattice>"
     << "</cell>";

  const std::string &xml = ss.str();
  _doc.Parse(xml.c_str());

  auto outer_cell = buildCell(_doc.FirstChildElement("cell"));

  ASSERT_TRUE(outer_cell->getFillUniverse());
  auto lattice = dynamic_cast<RecLattice*>(outer_cell->getFillUniverse());

  EXPECT_STREQ("outer cell", outer_cell->getName());
  EXPECT_EQ("-1 +2 -3 +4 -5 +6", outer_cell->getRegion());
  EXPECT_EQ(6, outer_cell->getNumSurfaces());

  // Test each attribute of the Lattice
  IntPyVec u_ids = {1};

  IntPyVec2D plane = {
    {1, 1},
    {1, 1}
  };

  IntPyVec3D layout = IntPyVec3D{plane} * 20;

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
      EXPECT_EQ   (2,           c->getInputId());
      //EXPECT_EQ   (2,           c->getFillMaterial()->getId());
      EXPECT_STREQ("inner cell",    c->getName());
      ASSERT_EQ   (5,           c->getNumSurfaces());
    }
  }

  delete lattice;
}


TEST_F(test_GeoInputXml_buildLattice2D, extrudedRecLatticeExertedWidths) {
  ROOT_ONLY();

  // Test input
  std::stringstream ss;
  ss << "<cell id='1' fill='101' name='outer cell' region='-1 +2 -3 +4 -5 +6'>"
     << "  <surface id='1' type='x-plane' coeffs='3.50364'/>"
     << "  <surface id='2' type='x-plane' coeffs='-3.50364'/>"
     << "  <surface id='3' type='y-plane' coeffs='3.50364'/>"
     << "  <surface id='4' type='y-plane' coeffs='-3.50364'/>"
     << "  <surface id='5' type='z-plane' coeffs='1.5'/>"
     << "  <surface id='6' type='z-plane' coeffs='-1.5'/>"
     << "  <lattice id='101' name='RecLattice' type='Rectangle' axial='extruded'>"
     << "    <row> 1 1 </row>"
     << "    <row> 1 1 </row>"
     << "    <widths x='[3.50364] * 2' y='[3.50364] * 2' z='[0.15] * 2'/>"
     << "    <refines x='1' y='1' z='10'/>"
     << "    <universe id='1' name='pin-cell'>"
     << "      <cell id='2' material='2' name='inner cell' region='7'>"
     << "        <surface id='7' type='z-cylinder' coeffs='0 0 0.3320'/>"
     << "      </cell>"
     << "    </universe>"
     << "  </lattice>"
     << "</cell>";
  constexpr int refines = 10;

  const std::string &xml = ss.str();
  _doc.Parse(xml.c_str());

  auto outer_cell = buildCell(_doc.FirstChildElement("cell"));

  ASSERT_TRUE(outer_cell->getFillUniverse());
  auto lattice = dynamic_cast<RecLattice*>(outer_cell->getFillUniverse());

  EXPECT_STREQ("outer cell", outer_cell->getName());
  EXPECT_EQ("-1 +2 -3 +4 -5 +6", outer_cell->getRegion());
  EXPECT_EQ(6, outer_cell->getNumSurfaces());

  // Test each attribute of the Lattice
  IntPyVec u_ids = {1};

  IntPyVec2D plane = {
    {1, 1},
    {1 ,1}
  };

  IntPyVec3D layout = IntPyVec3D{plane} * 2 * refines;

  Widths3D widths_oracles = {
    {3.50364, 3.50364},
    {3.50364, 3.50364},
    WidthVec{0.15 / refines} * 2 * refines
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

}/* namespace */
