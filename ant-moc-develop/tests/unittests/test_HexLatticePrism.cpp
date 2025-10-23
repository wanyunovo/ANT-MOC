/// \file test_HexLatticePrism.cpp
/// \brief Test surface HexLatticePrism
/// \date March 22, 2020
/// \author An Wang, USTB (wangan@xs.ustb.edu.cn)

#include "testing/test_utils.h"
#include "testing/geometry_testharness.h"
#include "antmoc/PyVector.h"
#include "antmoc/Surface.h"

#include <string>

using namespace antmoc;

namespace {

// Test fixture
class test_HexLatticePrism : public GeometryTestHarness {
  // do nothing
};


/// \brief Test evaluate()
/// \details This function won't actually evaluate anything for
///          a HexLatticePrism, but determines if a point is
///          inside/outside/on the surface.
TEST_F(test_HexLatticePrism, evaluate) {

  ROOT_ONLY();

  // Intialize the surface
  PyVector<HexLatticePrism> prisms = {
    HexLatticePrism(0., 0., 1.3, 1, "y"),
    HexLatticePrism(0., 0., 1.3, 2, "y"),
    HexLatticePrism(0., 0., 2.7, 2, "x"),
  };

  double tan30 = std::tan(M_PI/6);
  PyVector<PyVector<Point>> points = {
    {
      Point{0., 0., 0.},        // inside
      Point{1., 1.3, 0.},       // outside
      Point{0., 0.65, 0.},      // on surface
      Point{1.3*tan30, 0., 0.}, // on corner
    },
    {
      Point{0., 0.65, 0.},      // inside
      Point{1., 2.6, 0.},       // outside
      Point{0., 1.95, 0.},      // on surface
      Point{2.6*tan30, 0., 0.}, // on corner
    },
    {
      Point{0., 0., 0.},        // inside
      Point{5.4, 0., 0.},       // outside
      Point{4.05, 0., 0.},      // on surface
      Point{0., 5.4*tan30, 0.}, // on corner
    },
  };

  PyVector<double> oracles = {
    -1.,
    1.,
    0.,
    0.,
  };

  for (size_t i = 0; i < prisms.size(); ++i) {
    HexLatticePrism &prism = prisms[i];

    for (size_t j = 0; j < points[i].size(); ++j) {
      Point &point = points[i][j];
      double oracle = oracles[j];
      EXPECT_EQ(prism.evaluate(&point), oracle)
        << "HexLatticePrism " << i << ":\n" << prism.toString() << '\n'
        << "  Point " << j << ": " << point.toString() << std::endl;
    }
  }

}


/// \brief Test intersection() for points inside/outside the surface
/// \details Test intersection on a 2D plane
TEST_F(test_HexLatticePrism, intersection2D) {

  ROOT_ONLY();

  // Intialize the surface
  PyVector<HexLatticePrism> prisms = {
    HexLatticePrism(0., 0., 1.3, 1, "y"),
    HexLatticePrism(0., 0., 1.3, 2, "y"),
    HexLatticePrism(0., 0., 2.7, 2, "x"),
  };

  double cos30 = std::cos(M_PI/6);
  double sin30 = 0.5;
  PyVector<PyVector<SimpleTrack>> tracks = {
    {
      // inside
      SimpleTrack{  M_PI/6, M_PI/2, Point{0., 0., 3.33}},
      // outside the bounding cylinder, no intersection
      SimpleTrack{  M_PI/6, M_PI/2, Point{10, 10, 0.}},
      // outside, 1 intersection
      SimpleTrack{M_PI*5/6, M_PI/2, Point{1.3*cos30, -1.3*sin30, 0.}},
    },
    {
      // inside
      SimpleTrack{  M_PI/6, M_PI/2, Point{0., 0., 3.33}},
      // outside the bounding cylinder, no intersection
      SimpleTrack{  M_PI/6, M_PI/2, Point{10, 10, 0.}},
      // outside, 1 intersection
      SimpleTrack{M_PI*5/6, M_PI/2, Point{2.6*cos30, -2.6*sin30, 0.}},
    },
    {
      // inside
      SimpleTrack{      0, M_PI/2, Point{0., 0., 3.33}},
      // outside the bounding cylinder, no intersection
      SimpleTrack{  M_PI/6, M_PI/2, Point{10, 10, 0.}},
      // outside, 1 intersection
      SimpleTrack{-M_PI/3, M_PI/2, Point{-5.4*sin30, 5.4*cos30, 0.}},
    },
  };

  PyVector<PyVector<int>> oracle_num_intersections = {
    {
      1,
      0,
      1
    },
    {
      1,
      0,
      1
    },
    {
      1,
      0,
      1
    },
  };

  PyVector<PyVector<Point>> oracle_intersections = {
    {
      Point{0.65*cos30, 0.65*sin30, 3.33},
      Point{},
      Point{0.65*cos30, -0.65*sin30, 0.}
    },
    {
      Point{1.95*cos30, 1.95*sin30, 3.33},
      Point{},
      Point{1.95*cos30, -1.95*sin30, 0.}
    },
    {
      Point{4.05, 0., 3.33},
      Point{},
      Point{-4.05*sin30, 4.05*cos30, 0.}
    },
  };

  for (size_t i = 0; i < prisms.size(); ++i) {
    Point intersections[2];
    HexLatticePrism &prism = prisms[i];

    for (size_t j = 0; j < tracks[i].size(); ++j) {
      SimpleTrack &track = tracks[i][j];

      int num_intersections = prism.intersection(&(track.point),
                                                 track.azim,
                                                 track.polar,
                                                 intersections);
      // number of intersections
      ASSERT_EQ(num_intersections, oracle_num_intersections[i][j])
        << "HexLatticePrism " << i << ":\n" << prism.toString() << '\n'
        << "  Track: " << track.toString() << std::endl;

      // check the intersection
      if (num_intersections > 0) {
        EXPECT_EQ(intersections[0], oracle_intersections[i][j])
          << "HexLatticePrism " << i << ":\n" << prism.toString() << '\n'
          << "  Track: " << track.toString() << '\n'
          << "  Intersection " << j << ": " << intersections[0].toString() << '\n'
          << "  Oracle " << j << ": " << oracle_intersections[i][j].toString()
          << std::endl;
      }
    }
  }

}

} // namespace
