/// \file test_GeoInputXml_buildHexLatticePrism.cpp
/// \brief Test Geometry reading and building process
/// \date March 22, 2020
/// \author An Wang, USTB (wangan@xs.ustb.edu.cn)

#include "testing/geometry_testharness.h"
#include "antmoc/Surface.h"

using namespace antmoc;

namespace {

// Test fixture
class test_GeoInputXml_buildHexLatticePrism : public GeometryTestHarness {
  // do nothing
};


TEST_F(test_GeoInputXml_buildHexLatticePrism, createFromXML) {
  ROOT_ONLY();

  StringVec tests = {
    // test default orientation
    "<surface id='1' name='hex-lattice-prism-1' type='hex-lattice-prism'>"
    "  <center> 0. 0. </center>"
    "  <pitch> 0.7 </pitch>"
    "  <n_rings> 5 </n_rings>"
    "</surface>",

    // test default n_rings
    "<surface id='2' name='hex-lattice-prism-2' type='hex-lattice-prism'>"
    "  <center> 0. 0. </center>"
    "  <pitch> 0.7 </pitch>"
    "  <orientation> x </orientation>"
    "</surface>",

    // test hexagonal prism
    "<surface id='3' name='hex-lattice-prism-3' type='hex-lattice-prism'>"
    "  <center> 0. 0. </center>"
    "  <pitch> 0.7 </pitch>"
    "  <n_rings> 5 </n_rings>"
    "  <orientation> x </orientation>"
    "</surface>",

    // test large HexLatticePrism
    "<surface id='4' name='hex-lattice-prism-4' type='hex-lattice-prism'>"
    "  <center> 0. 0. </center>"
    "  <pitch> 0.7 </pitch>"
    "  <n_rings> 100 </n_rings>"
    "  <orientation> y </orientation>"
    "</surface>",
  };

  PyVector<HexLatticePrism> prisms = {
    HexLatticePrism(0., 0., 0.7,   5, "y", 1, "hex-lattice-prism-1"),
    HexLatticePrism(0., 0., 0.7,   1, "x", 2, "hex-lattice-prism-2"),
    HexLatticePrism(0., 0., 0.7,   5, "x", 3, "hex-lattice-prism-3"),
    HexLatticePrism(0., 0., 0.7, 100, "y", 4, "hex-lattice-prism-4"),
  };

  for (size_t i = 0; i < tests.size(); ++i) {
    _doc.Parse(tests[i].c_str());
    Surface *surface = buildSurface(_doc.FirstChildElement("surface"));

    // compare the parsed surface with the manually initialized one
    EXPECT_EQ(surface->toString(), prisms[i].toString());

    // clean up
    delete surface;
  }
}


} // namespace
