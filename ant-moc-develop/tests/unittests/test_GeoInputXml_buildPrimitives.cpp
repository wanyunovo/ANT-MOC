/**
 * @file test_GeoInputXml_buildPrimitives.cpp
 * @brief Test Geometry reading and building process
 * @date April 20, 2019
 * @author An Wang, USTB (wangan@xs.ustb.edu.cn)
 * @author Ya Fang, USTB (fangya201388@163.com)
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
class test_GeoInputXml_buildPrimitives : public GeometryTestHarness {

protected:

  void SetUp() {
    // Add materials
    StringVec names = {
      "1", "2", "3", "4", "m", "UO2"
    };

    for (auto &name : names)
      _mat_input->setMaterial(name, new Material(0, name.c_str()));
  }

};


#ifndef ENABLE_MPI_
TEST_F(test_GeoInputXml_buildPrimitives, buildUniverseRobust) {
  StringVec tests = {
    "<universe id='1'/>", // empty

    "<universe id='1'>"
    "  <id> 1 </id>"  // duplicated
    "  <cell id='1' material='m' name='c'/>"
    "</universe>",

    "<universe>"
    "  <id> 1 </id>"
    "  <id> 1 </id>"  // duplicated
    "  <cell id='1' material='m' name='c'/>"
    "</universe>",

    "<universe name='test'>"
    "  <name> test </name>"  // duplicated
    "  <cell id='1' material='m' name='c'/>"
    "</universe>",
  };

  for (auto &s : tests) {
    _doc.Parse(s.c_str());
    EXPECT_ANY_THROW(buildUniverse(_doc.FirstChildElement("universe")));
  }

  // These will pass
  StringVec pass = {
    "<universe id='1' name='test'>"
    "  <cell id='1' material='m' name='c'/>"
    "</universe>",
  };

  for (auto &s : pass) {
    _doc.Parse(s.c_str());
    EXPECT_NO_THROW(buildUniverse(_doc.FirstChildElement("universe")));
  }
}
#endif


TEST_F(test_GeoInputXml_buildPrimitives, buildUniverseWithManyCells) {
  const int id = 3;
  const char *name = "组件3";
  const int num_cells = 4;
  int c_id[] = {1, 4, 32, 512};
  const char *mat_name[] = {"1", "2", "3", "4"};
  float temp[] = {273.5, 500.0, 800., 10000.1};
  const char *c_name[] = {"C1", "C4", "C32", "C512"};
  std::stringstream ss;
  ss << "<universe id=\"" << id << "\" name=\"" << name << "\">";

  for (int i = 0; i < num_cells; i++) {
    ss << "  <cell id=\"" << c_id[i] << "\" material=\"" << mat_name[i] << "\""
       << "   temperature=\"" << temp[i] << "\" name=\"" << c_name[i] << "\">"
       << "  </cell>";
  }
  ss << "</universe>";

  const std::string &xml = ss.str();
  _doc.Parse(xml.c_str());

  Universe*  universe = buildUniverse(_doc.FirstChildElement("universe"));

  ASSERT_EQ(id, universe->getInputId());
  ASSERT_STREQ(name, universe->getName());

  std::map<int, Cell*> cells = universe->getCells();
  ASSERT_EQ(num_cells, cells.size());
  int i = 0;
  for (auto &c : cells){
     Cell *tem_cell = c.second;
     EXPECT_EQ(c_id[i], tem_cell->getInputId());
     EXPECT_STREQ(mat_name[i],tem_cell->getFillMaterial()->getName());
     EXPECT_EQ(temp[i], tem_cell->getTemperature());
     EXPECT_STREQ(c_name[i], tem_cell->getName());
     i++;
  }

}


#ifndef ENABLE_MPI_
TEST_F(test_GeoInputXml_buildPrimitives, buildCellRobust) {
  StringVec tests = {
    "<cell id='1'/>",
    "<cell id='1' fill='1' material='1'/>",
    "<cell id='1' fill='1'/>",
    //"<cell id='1' material='1'/>", FIXME, only gives a warning
    "<cell id='1' material='1' universes='1'>"
    "  <universes> 2 </universes>"  // duplicated
    "</cell>",

    "<cell id='1' material='1' region='-1'>"
    "  <region> -1 </region>"  // duplicated
    "  <surface id='1' name='circle' type='z-cylinder' coeffs='0 0 0.1'/>"
    "</cell>",

    "<cell id='1' material='1'>"
    "  <material> 2 </material>"  // duplicated
    "</cell>",

    "<cell id='1' material='1' temperature='273'>"
    "  <temperature> 273 </temperature>"  // duplicated
    "</cell>",

    "<cell id='1' material='1' sectors='3'>"
    "  <sectors> 5 </sectors>"  // duplicated
    "</cell>",
  };


  for (auto &s : tests) {
    _doc.Parse(s.c_str());
    EXPECT_ANY_THROW(buildCell(_doc.FirstChildElement("cell")));
  }
}
#endif


TEST_F(test_GeoInputXml_buildPrimitives, buildCellRegionByGlobal) {
  ROOT_ONLY();

  const std::string test =
    "<global>"
    "  <surface id='1' name='circle' type='z-cylinder' coeffs='0 0 0.3320'/>"
    "</global>"
    "<cell id='1' material='1' region='-1'/>";

  _doc.Parse(test.c_str());

  _geo_input->buildGlobalPrimitives(_doc.FirstChildElement("global"));
  auto cell = buildCell(_doc.FirstChildElement("cell"));

  ASSERT_EQ(1, cell->getNumSurfaces());

  // Check surface name
  auto surfaces = cell->getSurfaces();
  EXPECT_STREQ("circle", surfaces[1]->getName());

  Point p1(0.3, 0, 0);
  Point p2(0.4, 0, 0);
  EXPECT_TRUE(cell->containsPoint(&p1));
  EXPECT_FALSE(cell->containsPoint(&p2));
}


TEST_F(test_GeoInputXml_buildPrimitives, buildCellRegionByLocal) {
  ROOT_ONLY();

  const std::string test =
    "<cell id='1' material='1'>"
    "  <region> -1 +2 </region>"
    "  <surface id='1' name='outer circle' type='z-cylinder' coeffs='0 0 0.33'/>"
    "  <surface id='2' name='inner circle' type='z-cylinder' coeffs='0 0 0.11'/>"
    "</cell>";

  _doc.Parse(test.c_str());

  auto cell = buildCell(_doc.FirstChildElement("cell"));

  // Check global primitive maps
  EXPECT_EQ(0, _geo_input->getGlobalSurfaces().size());
  EXPECT_EQ(0, _geo_input->getGlobalCells().size());
  EXPECT_EQ(0, _geo_input->getGlobalUniverses().size());

  ASSERT_EQ(2, cell->getNumSurfaces());

  // Check surface name
  auto surfaces = cell->getSurfaces();
  EXPECT_STREQ("outer circle", surfaces[1]->getName());
  EXPECT_STREQ("inner circle", surfaces[2]->getName());

  Point p1(0.1, 0, 0), p2(0.4, 0, 0);
  Point p3(0.2, 0, 0);
  EXPECT_FALSE(cell->containsPoint(&p1));
  EXPECT_FALSE(cell->containsPoint(&p2));
  EXPECT_TRUE(cell->containsPoint(&p3));
}


TEST_F(test_GeoInputXml_buildPrimitives, buildCellRegionConflict) {
  ROOT_ONLY();

  const std::string test =
    "<global>"
    "  <surface id='1' name='circle' type='z-cylinder' coeffs='0 0 1.3320'/>"
    "</global>"
    "<cell id='1' material='1' region='-1'>"
    "  <surface id='1' name='circle' type='z-cylinder' coeffs='0 0 0.3320'/>"
    "</cell>";

  _doc.Parse(test.c_str());

  _geo_input->buildGlobalPrimitives(_doc.FirstChildElement("global"));
  auto cell = buildCell(_doc.FirstChildElement("cell"));

  ASSERT_EQ(1, cell->getNumSurfaces());

  // Check surface name
  auto surfaces = cell->getSurfaces();
  EXPECT_STREQ("circle", surfaces[1]->getName());

  Point p1(0.3, 0, 0);
  Point p2(0.4, 0, 0);
  EXPECT_TRUE(cell->containsPoint(&p1));
  EXPECT_FALSE(cell->containsPoint(&p2));
}


TEST_F(test_GeoInputXml_buildPrimitives, buildCellFilledWithMaterial) {
  ROOT_ONLY();

  const std::string test =
    "<cell id='1' name='pin' material='UO2'/>";

  _doc.Parse(test.c_str());
  auto cell = _geo_input->buildCell(_doc.FirstChildElement("cell"));

  ASSERT_TRUE(cell);
  EXPECT_EQ(1, cell->getInputId());
  EXPECT_STREQ("pin", cell->getName());
  EXPECT_STREQ("UO2", cell->getFillMaterial()->getName());
}


} // namespace
