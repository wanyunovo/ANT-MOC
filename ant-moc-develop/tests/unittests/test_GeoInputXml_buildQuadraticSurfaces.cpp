/// \file test_GeoInputXml_buildQuadraticSurfaces.cpp
/// \brief Test Geometry reading and building process
/// \date April 20, 2019
/// \author An Wang, USTB (wangan@xs.ustb.edu.cn)
/// \author Ya Fang, USTB (fangya201388@163.com)

#include "testing/geometry_testharness.h"
#include "antmoc/Surface.h"

using namespace antmoc;

namespace {

// Test fixture
class test_GeoInputXml_buildQuadraticSurfaces : public GeometryTestHarness{
  // do nothing
};


#ifndef ENABLE_MPI_
TEST_F(test_GeoInputXml_buildQuadraticSurfaces, robust) {
  StringVec tests = {
    "<surface id='1' coeffs='0 0 0.1'/>",
    "<surface id='1' type='z-cylinder'/>",
    "<surface id='1' type='undefined' coeffs='0 0 0.1'/>",

    "<surface id='2' coeffs='0 0 0.1'>"
    "  <type> </type>"
    "</surface>",

    "<surface id='3' type='z-cylinder' coeffs='0 0 0.1'>"
    "  <type> z-cylinder </type>" // duplicated
    "</surface>",
  };

  for (auto &s : tests) {
    _doc.Parse(s.c_str());
    EXPECT_ANY_THROW(buildSurface(_doc.FirstChildElement("surface")));
  }

  // These will pass
  StringVec pass = {
    "<surface id='1' type='z-cylinder' coeffs='0 0 0.1'/>",

    "<surface id='2' name='u2' coeffs='0 0 0.1'>"
    "  <type>      z-cylinder </type>"
    "  <boundary>  vacuum     </boundary>"
    "</surface>",
  };

  for (auto &s : pass) {
    _doc.Parse(s.c_str());
    try {
      buildSurface(_doc.FirstChildElement("surface"));
    }
    catch (const std::logic_error &e)
      { FAIL() << e.what(); }
  }
}
#endif

TEST_F(test_GeoInputXml_buildQuadraticSurfaces, createPlane) {
  ROOT_ONLY();

  const char *type = "plane";
  const char *name = "plane";
  const char *coeff = "1 1.2 3.45678 -9";
  double A, B, C, D;

  std::istringstream coeff_s(coeff);
  coeff_s >> A >> B >> C >> D;

  Surface *surface = _geo_input->newSurfaceQuadratic(type, coeff, 0, name);

  // Test if the surface is of type Plane
  Plane *plane = dynamic_cast<Plane *> (surface);
  ASSERT_TRUE(plane);
  EXPECT_STREQ(name, plane->getName());
  EXPECT_EQ(A, plane->getA());
  EXPECT_EQ(B, plane->getB());
  EXPECT_EQ(C, plane->getC());
  EXPECT_EQ(D, plane->getD());
}

TEST_F(test_GeoInputXml_buildQuadraticSurfaces, createXPlane) {
  ROOT_ONLY();

  const char *type = "x-plane";
  const char *name = "xplane";
  const char *coeff = "1.2";
  double X;

  std::istringstream coeff_s(coeff);
  coeff_s >> X;

  Surface *surface = _geo_input->newSurfaceQuadratic(type, coeff, 0, name);

  // Test if the surface is of type Plane
  XPlane *xplane = dynamic_cast<XPlane *> (surface);
  ASSERT_TRUE(xplane);
  EXPECT_STREQ(name, xplane->getName());
  EXPECT_EQ(X, xplane->getX());
}

TEST_F(test_GeoInputXml_buildQuadraticSurfaces, createYPlane) {
  ROOT_ONLY();

  const char *type = "y-plane";
  const char *name = "yplane";
  const char *coeff = "1.2";
  double Y;

  std::istringstream coeff_s(coeff);
  coeff_s >> Y;

  Surface *surface = _geo_input->newSurfaceQuadratic(type, coeff, 0, name);

  // Test if the surface is of type Plane
  YPlane *yplane = dynamic_cast<YPlane *> (surface);
  ASSERT_TRUE(yplane);
  EXPECT_STREQ(name, yplane->getName());
  EXPECT_EQ(Y, yplane->getY());
}

TEST_F(test_GeoInputXml_buildQuadraticSurfaces, createZPlane) {
  ROOT_ONLY();

  const char *type = "z-plane";
  const char *name = "zplane";
  const char *coeff = "1.2";
  double Z;

  std::istringstream coeff_s(coeff);
  coeff_s >> Z;

  Surface *surface = _geo_input->newSurfaceQuadratic(type, coeff, 0, name);

  // Test if the surface is of type Plane
  ZPlane *zplane = dynamic_cast<ZPlane *> (surface);
  ASSERT_TRUE(zplane);
  EXPECT_STREQ(name, zplane->getName());
  EXPECT_EQ(Z, zplane->getZ());
}

TEST_F(test_GeoInputXml_buildQuadraticSurfaces, createZCylinder) {
  ROOT_ONLY();

  const char *type = "z-cylinder";
  const char *name = "zcylinder";
  const char *coeff = "1 1.2 0.5";
  double X0, Y0, R;

  std::istringstream coeff_s(coeff);
  coeff_s >> X0 >> Y0 >> R;

  Surface *surface = _geo_input->newSurfaceQuadratic(type, coeff, 0, name);

  // Test if the surface is of type Plane
  ZCylinder *zcylinder = dynamic_cast<ZCylinder *> (surface);
  ASSERT_TRUE(zcylinder);
  EXPECT_STREQ(name, zcylinder->getName());
  EXPECT_EQ(X0, zcylinder->getX0());
  EXPECT_EQ(Y0, zcylinder->getY0());
  EXPECT_EQ(R, zcylinder->getRadius());
}

} // namespace
