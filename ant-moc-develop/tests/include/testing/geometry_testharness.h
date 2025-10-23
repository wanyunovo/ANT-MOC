/// \file geometry_testharness.h
/// \brief A test harness for ANT-MOC
/// \date April 7, 2019
/// \author An Wang, USTB (wangan@xs.ustb.edu.cn)

#ifndef GEOMETRY_TESTHARNESS_H_
#define GEOMETRY_TESTHARNESS_H_

#include <sstream>

#include "testing/MockGeometry.h"
#include "testing/geometry_testutils.h"
#include "testing/test_utils.h"
#include "antmoc/Geometry.h"
#include "antmoc/GeoInputXml.h"
#include "antmoc/MaterialHandler.h"
#include "antmoc/Point.h"
#include "antmoc/string_utils.h"
#include "antmoc/xml_utils.h"

using namespace antmoc;


namespace {


///---------------------------------------------------------------------
/// \brief A simple track for testing
///---------------------------------------------------------------------
struct SimpleTrack {
  double azim;
  double polar;
  Point point;

  std::string toString() {
    std::stringstream ss;

    ss << "Starting point = " << point.toString()
       << ", azim = " << azim
       << ", polar= " << polar;

    return ss.str();
  }
};


///---------------------------------------------------------------------
/// \brief A test harness for geometry testing.
/// \details In googletest, this is actually a test fixture.
///---------------------------------------------------------------------
class GeometryTestHarness : public testing::Test {

public:

  /// \brief Initialize the test fixture.
  GeometryTestHarness():
    _mat_input(std::make_shared<MaterialHandler>()),
    _geo_input(std::make_shared<GeoInputXml>(&_geo_mock)) {

    antmoc::log::set_level("test");

    _geo_input->setMaterialHandler(_mat_input);
  }

  /// \brief Cleanup objects
  /// \details All of the unused objects will be deleted. So be careful with it
  ///          because it is tempting to delete allocated objects in individual
  ///          tests. For now, unused objects includes global primitives and
  ///          materials. That is, objects allocated by methods like newCell,
  ///          newSurface, etc. will not be handled by this method.
  void TearDown() {

    // Cleanup unused global primitives
    _geo_input->eraseUnusedPrimitives();

    // Cleanup geometry objects
    for (auto &s : _geo_mock.getAllSurfaces())
      delete s.second;
    for (auto &c : _geo_mock.getAllCells())
      delete c.second;
    for (auto &u : _geo_mock.getAllUniverses())
      delete u.second;

    // Cleanup all materials
    auto &materials = _mat_input->getAllMaterials();
    for (auto &e : materials)
      delete e.second;

    materials.clear();

  }


  /// \brief Build a surface and return the pointer to it.
  Surface* buildSurface(XMLElement *e) {
    return _geo_input->newSurface(e);
  }

  /// \brief Build a cell and return the pointer to it.
  Cell* buildCell(XMLElement *e,
                  std::string exerted_region = "",
                  SurfaceMap *exerted_surfaces = nullptr) {

    return _geo_input->buildCell(e, exerted_region, exerted_surfaces);
  }

  /// \brief Build a universe and return the pointer to it.
  Universe* buildUniverse(XMLElement *e,
                          std::string exerted_region = "",
                          SurfaceMap *exerted_surfaces = nullptr) {

    return _geo_input->buildUniverse(e, exerted_region, exerted_surfaces);
  }

  /// \brief Build a lattice and return the pointer to it.
  Lattice* buildLattice(XMLElement *e) {
    return _geo_input->buildLattice(e);
  }

  /// \brief Build global primitives, which can be automatically freed
  void buildGlobalPrimitives(XMLElement *e) {
    _geo_input->buildGlobalPrimitives(e);
  }


protected:

  tinyxml2::XMLDocument _doc;   ///< For reading and parsing xml files
  MockGeometry _geo_mock;       ///< Geometry
  std::shared_ptr<MaterialHandler> _mat_input;
  std::shared_ptr<GeoInputXml> _geo_input;

};


} // anonymous namespace

#endif
