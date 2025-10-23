/**
 * @file test_GeoInputXml_buildGlobal.cpp
 * @brief Test global primitives
 * @date Sep 4, 2019
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
class test_GeoInputXml_buildGlobal : public GeometryTestHarness {

protected:

  void SetUp() {
    // Add materials
    StringVec names = {
      "Beer", "Milk", "Cola", "Juice", "Bottle"
    };

    for (auto &name : names)
      _mat_input->setMaterial(name, new Material(0, name.c_str()));
  }

};


#ifndef ENABLE_MPI_

TEST_F(test_GeoInputXml_buildGlobal, robust) {
  ROOT_ONLY();
  StringVec tests = {
    "<global> <cell id='1' material='Beer' region='1'/> </global>",
    "<global> <cell id='1' fill='2'/> </global>",

    "<global>"
    "  <surface id='1' type='z-cylinder' coeffs='0 0 0.5'/>"
    "  <surface id='1' type='z-cylinder' coeffs='0 0 0.5'/>"
    "</global>",

    "<global>"
    "  <cell id='1' universes='1'/>"
    "  <cell id='1' universes='2'/>"
    "</global>",

    "<global>"
    "  <universe id='1'/>"
    "  <universe id='1'/>"
    "</global>",

  };

  for (auto &s : tests) {
    _doc.Parse(s.c_str());
    auto global = _doc.FirstChildElement("global");
    _geo_input->clear();
    EXPECT_ANY_THROW(_geo_input->buildGlobalPrimitives(global))
      << " Failed case:\n" << s << '\n';
  }
}

#endif  // ENABLE_MPI_


TEST_F(test_GeoInputXml_buildGlobal, recursively) {
  ROOT_ONLY();

  const std::string test =
    "<global>"
    "  <surface id='1' type='z-cylinder' coeffs='0 0 0.5'/>"
    "  <surface id='3' type='z-cylinder' coeffs='0 0 0.6'/>"

    "  <cell id='1' material='Beer'   name='beer_cell'   region='-1'/>"
    "  <cell id='2' material='Bottle' name='bottle_cell' region='+1 +3'/>"

    "  <universe id='1' name='beer'>"
    "    <cells> 1 2 </cells>"
    "  </universe>"
    "  <universe id='2' name='milk'>"
    "    <cells> 1 2 </cells>"
    "  </universe>"

    "  <lattice id='3' name='beers' type='rectangle'>"
    "    <widths x='[2.0]*2' y='2.0' z='5.0'/>"
    "    <row> 1 1 </row>"
    "  </lattice>"

    "  <universe id='4' name='milk box'>"
    "  <cell id='3' fill='5' name='milk_box'>"
    "    <lattice id='5' name='milk' type='rectangle'>"
    "      <widths x='[2.0]*2' y='2.0' z='5.0'/>"
    "      <row> 2 2 </row>"
    "    </lattice>"
    "  </cell>"
    "  </universe>"
    "</global>";

  _doc.Parse(test.c_str());

  _geo_input->buildGlobalPrimitives(_doc.FirstChildElement("global"));
  auto &surfaces   = _geo_input->getGlobalSurfaces();
  auto &cells      = _geo_input->getGlobalCells();
  auto &universes  = _geo_input->getGlobalUniverses();

  // Check surfaces
  std::vector<int>         ids     = {1, 3};
  std::vector<surfaceType> sftypes = {ZCYLINDER};
  std::vector<double>      radii   = {0.5, 0.6};

  ASSERT_EQ(ids.size(), surfaces.size());
  for (size_t i = 0; i < ids.size(); ++i) {
    const int id = ids[i];
    ASSERT_TRUE(surfaces.find(id) != surfaces.end());
    const auto sf = dynamic_cast<ZCylinder*>(surfaces[ids[i]]);
    ASSERT_TRUE(sf);
    EXPECT_EQ(sftypes[0], sf->getSurfaceType());
    EXPECT_EQ(radii[i],   sf->getRadius());
  }

  // Check cells
  ids = {1, 2, 3};
  StringVec names = {"beer_cell", "bottle_cell", "milk_box"};
  std::vector<bool> is_material = {true, true, false};
  StringVec fill_names = {"Beer", "Bottle", "milk"};
  StringVec regions = {
    "-1",
    "+1 +3",
    "",
  };
  IntPyVec2D sf_ids = {
    {1},
    {1, 3},
    {},
  };

  ASSERT_EQ(ids.size(), cells.size());
  for (size_t i = 0; i < ids.size(); ++i) {
    const int id = ids[i];
    ASSERT_TRUE(cells.find(id) != cells.end());
    const auto c = cells[ids[i]];
    testCell(c, id, names[i], is_material[i], fill_names[i], &regions[i], &sf_ids[i]);
  }

  // Separate lattices from universes
  RecLatticeMap lattices;
  auto it = universes.begin();
  while (it != universes.end()) {

    auto l = dynamic_cast<RecLattice*>(it->second);
    if (l) {
      lattices.insert({it->first, l});
      it = universes.erase(it);
    } else {
      ++it;
    }
  }

  // Check universes
  ids = {1, 2, 4};
  names = {"beer", "milk", "milk box"};

  IntPyVec2D c_ids = {
    {1, 2},
    {1, 2},
    {3},
  };

  ASSERT_EQ(ids.size(), universes.size());
  for (size_t i = 0; i < ids.size(); ++i) {
    const int id = ids[i];
    ASSERT_TRUE(universes.find(id) != universes.end());
    const auto u = universes[id];
    testUniverse(u, id, names[i], &c_ids[i]);
  }

  // Check lattices
  ids = {3, 5};
  names = {"beers", "milk"};

  std::vector<IntPyVec3D> layouts = {
    {{
     {1, 1},
    }},

    {{
     {2, 2},
    }},
  };

  std::vector<Widths3D> widths = {
    {
      {2.0, 2.0},
      {2.0},
      {5.0}
    },
    {
      {2.0, 2.0},
      {2.0},
      {5.0}
    }
  };

  IntPyVec2D u_ids = {
    {1},
    {2},
  };

  ASSERT_EQ(ids.size(), lattices.size());
  for (size_t i = 0; i < ids.size(); ++i) {
    const int id = ids[i];
    ASSERT_TRUE(lattices.find(id) != lattices.end());
    const auto l = lattices[id];
    testRecLattice(l, id, names[i], &layouts[i], &widths[i], &u_ids[i]);
  }
}


TEST_F(test_GeoInputXml_buildGlobal, lazyComputedUniverses) {
  ROOT_ONLY();

  const std::string test =
    "<global>"
    "  <surface id='1' type='z-cylinder' coeffs='0 0 0.5'/>"
    "  <surface id='3' type='z-cylinder' coeffs='0 0 0.6'/>"

    "  <cell id='1'  universes='1'   material='Beer'   name='beer'   region='-1'>"
    "    <sectors> 8 </sectors>"
    "    <rings>   4 </rings>"
    "  </cell>"

    "  <cell id='2'  universes='2'   material='Milk'   name='milk'   region='-1' rings='2'>"
    "    <sectors> 5 </sectors>"
    "  </cell>"

    "  <cell id='10'                 material='Bottle' name='bottle' region='+1 +3'>"
    "               <universes> 1 2 </universes>" // test 'universes' sub-element
    "  </cell>"

    "  <cell id='3'  universes='3'   fill='5'          name='beers_box' rings='2'>"
    "    <sectors> 5 </sectors>" // rings and sectors won't work for filled cells
    "  </cell>"

    "  <cell id='4' universes='4' fill='6' name='milk_box'>"
    "    <lattice id='6' name='milk' type='rectangle'>"
    "      <widths x='[2.0]*2' y='2.0' z='5.0'/>"
    "      <row> 2 2 </row>"
    "    </lattice>"
    "  </cell>"

    "  <lattice id='5' name='beers' type='rectangle'>"
    "    <widths x='[2.0]*2' y='2.0' z='5.0'/>"
    "    <row> 1 1 </row>"
    "  </lattice>"

    "  <lattice id='10' name='water' type='Hexagon' n_rings='2'>"
    "    <widths r='3.14' z='[1.0]*2'/>"
    "    <universes>"
    "          1     "
    "       1     1  "
    "          21    "
    "       1     1  "
    "          1     "
    "          2     "
    "       2     2  "
    "          1     "
    "       2     3  "
    "          3     "
    "    </universes>"
    "    <universe id='21' cells='10' name='global'/>"
    "    <universe         cells='1'  name='local1'/>"
    "    <universe         cells='1'  name='local2'/>"
    "  </lattice>"

    " <universe cells='10' name='anonymous'/>"
    "</global>";

  _doc.Parse(test.c_str());

  _geo_input->buildGlobalPrimitives(_doc.FirstChildElement("global"));
  auto &surfaces = _geo_input->getGlobalSurfaces();
  auto &cells      = _geo_input->getGlobalCells();
  auto &alluniverses  = _geo_input->getGlobalUniverses();

  // Check global surfaces
  std::vector<int>         ids     = {1, 3};
  std::vector<surfaceType> hstypes = {ZCYLINDER};
  std::vector<double>      radii   = {0.5, 0.6};

  ASSERT_EQ(ids.size(), surfaces.size());
  for (size_t i = 0; i < ids.size(); ++i) {
    const int id = ids[i];
    ASSERT_TRUE(surfaces.find(id) != surfaces.end());
    const auto sf = dynamic_cast<ZCylinder*>(surfaces[ids[i]]);
    ASSERT_TRUE(sf);
    EXPECT_EQ(hstypes[0], sf->getSurfaceType());
    EXPECT_EQ(radii[i],   sf->getRadius());
  }

  // Check global cells
  ids = {1, 2, 10, 3, 4};
  StringVec names = {"beer", "milk", "bottle", "beers_box", "milk_box"};
  std::vector<bool> is_material = {true, true, true, false, false};
  StringVec fill_names = {"Beer", "Milk", "Bottle", "beers", "milk"};
  StringVec regions = {
    "-1",
    "-1",
    "+1 +3",
    "",
    "",
  };
  IntPyVec2D sf_ids = {
    {1},
    {1},
    {1, 3},
    {},
    {},
  };

  IntPyVec sectors = { 8, 5, 0, 0, 0};
  IntPyVec rings   = { 4, 2, 0, 0, 0};

  ASSERT_EQ(ids.size(), cells.size());
  for (size_t i = 0; i < ids.size(); ++i) {
    const int id = ids[i];
    ASSERT_TRUE(cells.find(id) != cells.end());
    const auto c = cells[ids[i]];
    testCell(c, id, names[i], is_material[i], fill_names[i], &regions[i], &sf_ids[i], sectors[i], rings[i]);
  }

  // Separate lattices from universes
  UniverseMap universes;
  RecLatticeMap rec_lattices;
  HexLatticeMap hex_lattices;
  for (auto it = alluniverses.begin(); it != alluniverses.end(); ++it) {
    auto rec = dynamic_cast<RecLattice*>(it->second);
    auto hex = dynamic_cast<HexLattice*>(it->second);
    if (rec)
      rec_lattices.insert({it->first, rec});
    else if (hex)
      hex_lattices.insert({it->first, hex});
    else
      universes.insert(*it);
  }

  // Check global universes
  ids = {1, 2, 3, 4, 21};
  names = {"beer", "milk", "beers_box", "milk_box", "global"};

  IntPyVec2D c_ids = {
    {1, 10},
    {2, 10},
    {3},
    {4},
    {10},
  };

  ASSERT_EQ(ids.size(), universes.size());
  for (size_t i = 0; i < ids.size(); ++i) {
    const int id = ids[i];
    ASSERT_TRUE(universes.find(id) != universes.end());
    const auto u = universes[id];
    testUniverse(u, id, names[i], &c_ids[i]);
  }

  // Check rectangle lattices
  ids = {5, 6};
  names = {"beers", "milk"};

  std::vector<IntPyVec3D> layouts = {
    {{
     {1, 1},
    }},

    {{
     {2, 2},
    }},
  };

  std::vector<Widths3D> widths = {
    {
      {2.0, 2.0},
      {2.0},
      {5.0}
    },
    {
      {2.0, 2.0},
      {2.0},
      {5.0}
    }
  };

  IntPyVec2D u_ids = {
    {1},
    {2},
  };

  ASSERT_EQ(ids.size(), rec_lattices.size());
  for (size_t i = 0; i < ids.size(); ++i) {
    const int id = ids[i];
    ASSERT_TRUE(rec_lattices.find(id) != rec_lattices.end());
    const auto l = rec_lattices[id];
    testRecLattice(l, id, names[i], &layouts[i], &widths[i], &u_ids[i]);
  }

  // Check hexagon lattices
  ids = {10};
  names = {"water"};

  std::vector<IntPyVec> hex_layouts = {
    {3, 3, 2, 1, 2, 2, 2,
    1, 1, 1, 21, 1, 1, 1},
  };

  std::vector<double> hex_widths_r = {
    3.14,
  };

  std::vector<WidthVec> hex_widths_z = {
    {1.0, 1.0},
  };

  std::vector<int> num_r = {
    2
  };

  u_ids = {
    {1, 2, 3, 21},
  };

  ASSERT_EQ(ids.size(), hex_lattices.size());
  for (size_t i = 0; i < ids.size(); ++i) {
    const int id = ids[i];
    ASSERT_TRUE(hex_lattices.find(id) != hex_lattices.end());
    const auto l = hex_lattices[id];
    testHexLattice(l, id, names[i], &hex_layouts[i], &hex_widths_r[i],
                   &hex_widths_z[i], &num_r[i], &u_ids[i]);
  }

}


TEST_F(test_GeoInputXml_buildGlobal, useGlobalPrimitivesInGeometry) {
  ROOT_ONLY();

  const std::string test =
    "<global>"
    "  <surface id='1' type='z-cylinder' coeffs='0 0 0.5'/>"
    "  <surface id='3' type='z-cylinder' coeffs='0 0 0.6'/>"

    "  <cell id='1'  universes='1'   material='Beer'   name='beer'   region='-1'/>"
    "  <cell id='2'  universes='2'   material='Milk'   name='milk'   region='-1'/>"
    "  <cell id='3'  universes='1 2' material='Bottle' name='bottle' region='+1 +3'/>"
    "  <cell id='10'                 material='Cola'   name='cola'   region='-1'/>"
    "  <cell id='11'                 material='Juice'  name='juice'  region='-1'/>"
    "</global>";

  const std::string test_geo =
    "<cell id='1' fill='1' name='box'>"
    "  <lattice id='1' name='box_lattice' type='rectangle'>"
    "    <widths x='[2.0]*3' y='[2.0]*3' z='5.0'/>"
    "    <row> 1 2 1 </row>"
    "    <row> 2 3 2 </row>"
    "    <row> 4 2 4 </row>"

    "    <universe id='3' name='cola'>"
    "      <cells> 3 10 </cells>" // as sub-element
    "    </universe>"

    "    <universe id='4' name='juice' cells='3 11'/>"  // as attribute

    "  </lattice>"
    "</cell>";

  _doc.Parse(test.c_str());
  _geo_input->buildGlobalPrimitives(_doc.FirstChildElement("global"));
  _doc.Parse(test_geo.c_str());
  auto cell = _geo_input->buildCell(_doc.FirstChildElement("cell"));

  // Test the bounding cell
  ASSERT_TRUE(cell);
  testCell(cell, 1, "box", false, "box_lattice", nullptr);

  // Test the lattice
  auto lattice = dynamic_cast<RecLattice*>(cell->getFillUniverse());
  ASSERT_TRUE(lattice);

  IntPyVec u_ids = {1, 2, 3, 4};

  // Test lattice layout
  IntPyVec3D layout = {
    {
      {1, 2, 1},
      {2, 3, 2},
      {4, 2, 4},
    }
  };

  Widths3D widths = {
    {2.0, 2.0, 2.0},
    {2.0, 2.0, 2.0},
    {5.0},
  };

  testRecLattice(lattice, 1, "box_lattice", &layout, &widths, &u_ids);

  // Test the inner universe
  Universe *u1 = nullptr, *u2 = nullptr;
  const auto &universes = lattice->getUniqueUniverses();
  for (auto &u : universes) {
    if (3 == u.second->getInputId() && !u1)
      u1 = u.second;
    if (4 == u.second->getInputId() && !u2)
      u2 = u.second;
  }

  ASSERT_TRUE(u1);
  IntPyVec c_ids = {3, 10};
  testUniverse(u1, 3, "cola", &c_ids);

  ASSERT_TRUE(u2);
  c_ids = {3, 11};
  testUniverse(u2, 4, "juice", &c_ids);
}

}/* namespace */
