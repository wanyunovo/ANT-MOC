/// \date Apri 6, 2020
/// \author An Wang, USTB (wangan@xs.ustb.edu.cn)

#include "testing/test_utils.h"
#include "antmoc/log.h"
#include "antmoc/Point.h"
#include "antmoc/PyVector.h"
#include "antmoc/Surface.h"

#include <sstream>
#include <string>

using namespace antmoc;

namespace {

// Test fixture
class test_ZCylinder : public testing::Test {

protected:

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

};


/// \brief Test evaluate()
TEST_F(test_ZCylinder, evaluate) {

  ROOT_ONLY();

  // Intialize the surface
  PyVector<ZCylinder> cylinders = {
    ZCylinder(0., 0., 1.3),
    ZCylinder(2., 2., 1.3),
  };

  PyVector<PyVector<Point>> points = {
    {
      Point{0., 0., 0.},        // inside
      Point{1., 1.3, 0.},       // outside
      Point{0., 1.3, 0.},       // on surface
    },
    {
      Point{1.9, 1.9, 0.},      // inside
      Point{0.,   0., 0.},      // outside
      Point{0.7, 2.0, 0.},      // on surface
    },
  };

  PyVector<PyVector<double>> oracles = {
    {
      -1.69,
      1.,
      0.
    },
    {
      -1.67,
      6.31,
      0.
    }
  };

  for (size_t i = 0; i < cylinders.size(); ++i) {
    ZCylinder &cylinder = cylinders[i];

    for (size_t j = 0; j < points[i].size(); ++j) {
      Point &point = points[i][j];
      double oracle = oracles[i][j];
      EXPECT_DOUBLE_EQ(cylinder.evaluate(&point), oracle)
        << "ZCylinder " << i << ":\n" << cylinder.toString() << '\n'
        << "  Point " << j << ": " << point.toString() << std::endl;
    }
  }

}


/// \brief Test intersection() for points inside/outside the surface
/// \details Test intersection on a 2D plane
TEST_F(test_ZCylinder, intersection2D) {

  ROOT_ONLY();

  // Intialize the surface
  ZCylinder cylinder(0., 0., 1.3);

  double cos45 = std::cos(M_PI/4);
  double sin45 = cos45;
  PyVector<SimpleTrack> tracks = {
    // inside, parallel to x-axis
    SimpleTrack{      0, M_PI/2, Point{0., 0., 0.5} },
    // inside, parallel to y-axis
    SimpleTrack{ M_PI/2, M_PI/2, Point{0., 0., 0.5} },
    // inside
    SimpleTrack{ M_PI/4, M_PI/2, Point{0., 0., 0.5} },

    // outside, parallel to x-axis, no intersection
    SimpleTrack{      0, M_PI/2, Point{2., 0., 0.5} },
    // outside, no intersection
    SimpleTrack{ M_PI/4, M_PI/2, Point{0., 10., 0.5} },

    // outside, parallel to x-axis, 2 intersections
    SimpleTrack{      0, M_PI/2, Point{-10., 0., 0.5} },
    // outside, parallel to y-axis, 2 intersections
    SimpleTrack{ M_PI/2, M_PI/2, Point{0., -10., 0.5} },
    // outside, 2 intersections
    SimpleTrack{ M_PI/4, M_PI/2, Point{-1.3,-1.3, 0.5} },

    // outside, parallel to x-axis, 1 intersections (on surface)
    SimpleTrack{      0, M_PI/2, Point{-1.3,-1.3, 0.5} },
    // outside, parallel to y-axis, 1 intersections (on surface)
    SimpleTrack{ M_PI/2, M_PI/2, Point{1.3,-1.3, 0.5} },
    // outside, 1 intersections (on surface)
    SimpleTrack{ M_PI/4, M_PI/2, Point{0.,-1.3/cos45, 0.5} },
  };

  PyVector<int> oracle_num_intersections = {
    // a ray from an inside point can only have 1 intersection
    1,
    1,
    1,
    // outside, no intersection
    0,
    0,
    // a ray from an outside point can have two candidate intersections
    2,
    2,
    2,
    // on surface
    1,
    1,
    1
  };

  PyVector<PyVector<Point>> oracle_intersections = {
    // inside, parallel to x-axis
    {
      Point{ 1.3, 0., 0.5},
    },
    // inside, parallel to y-axis
    {
      Point{0., 1.3, 0.5},
    },
    // inside
    {
      Point{ 1.3*cos45, 1.3*sin45, 0.5},
    },
    // outside, parallel to x-axis, no intersection
    {
      Point{},
    },
    // outside, no intersection
    {
      Point{},
    },
    // outside, parallel to x-axis, 2 intersections
    {
      Point{ 1.3, 0., 0.5},
      Point{-1.3, 0., 0.5},
    },
    // outside, parallel to y-axis, 2 intersections
    {
      Point{0., 1.3, 0.5},
      Point{0.,-1.3, 0.5},
    },
    // outside, 2 intersections
    {
      Point{ 1.3*cos45, 1.3*sin45, 0.5},
      Point{-1.3*cos45,-1.3*sin45, 0.5},
    },
    // outside, parallel to x-axis, 1 intersections (on surface)
    {
      Point{0.,-1.3, 0.5},
    },
    // outside, parallel to y-axis, 1 intersections (on surface)
    {
      Point{1.3, 0., 0.5},
    },
    // outside, 1 intersections (on surface)
    {
      Point{ 1.3*cos45,-1.3*sin45, 0.5},
    }
  };

  Point intersections[2];
  for (size_t j = 0; j < tracks.size(); ++j) {
    SimpleTrack &track = tracks[j];

    int num_intersections = cylinder.intersection(&(track.point),
                                                  track.azim,
                                                  track.polar,
                                                  intersections);
    // number of intersections
    ASSERT_EQ(num_intersections, oracle_num_intersections[j])
      << "ZCylinder " << ":\n" << cylinder.toString() << '\n'
      << "  Track: " << j << ": " << track.toString() << std::endl;

    // check the intersection
    for (int k = 0; k < num_intersections; ++k) {
      EXPECT_EQ(intersections[k], oracle_intersections[j][k])
        << "ZCylinder " << ":\n" << cylinder.toString() << '\n'
        << "  Track " << j << ": " << track.toString() << '\n'
        << "  Intersection " << k << ": " << intersections[k].toString() << '\n'
        << "  Oracle " << k << ": " << oracle_intersections[j][k].toString()
        << std::endl;
    }
  }

}


/// \brief Test intersection() for points inside/outside the surface
/// \details Test intersection with 3-D tracks. Each of the track has
///          polar angle = pi/6.
TEST_F(test_ZCylinder, intersection3D) {

  ROOT_ONLY();

  // Intialize the surface
  ZCylinder cylinder(0., 0., 1.3);

  double cos45 = std::cos(M_PI/4);
  double sin45 = cos45;
  double tan30 = std::tan(M_PI/6);

  PyVector<SimpleTrack> tracks = {
    // inside, parallel to x-axis
    SimpleTrack{      0, M_PI/3, Point{0., 0., 0.5} },
    // inside, parallel to y-axis
    SimpleTrack{ M_PI/2, M_PI/3, Point{0., 0., 0.5} },
    // inside, parallel to y-axis (reverse)
    SimpleTrack{3*M_PI/2,2*M_PI/3, Point{0., 0., 0.5} },
    // inside
    SimpleTrack{ M_PI/4, M_PI/3, Point{0., 0., 0.5} },

    // outside, parallel to x-axis, no intersection
    SimpleTrack{      0, M_PI/3, Point{2., 0., 0.5} },
    // outside, no intersection
    SimpleTrack{ M_PI/4, M_PI/3, Point{0., 10., 0.5} },

    // outside, parallel to x-axis, 2 intersections
    SimpleTrack{      0, M_PI/3, Point{-10., 0., 0.5} },
    // outside, parallel to y-axis, 2 intersections
    SimpleTrack{ M_PI/2, M_PI/3, Point{ 0.,-10., 0.5} },
    // outside, parallel to y-axis, 2 intersections (reverse)
    SimpleTrack{3*M_PI/2,2*M_PI/3, Point{ 0., 10., 0.5} },
    // outside, 2 intersections
    SimpleTrack{ M_PI/4, M_PI/3, Point{-1.3,-1.3, 0.5} },

    // outside, parallel to x-axis, 1 intersections (on surface)
    SimpleTrack{      0, M_PI/3, Point{-1.3,-1.3, 0.5} },
    // outside, parallel to y-axis, 1 intersections (on surface)
    SimpleTrack{ M_PI/2, M_PI/3, Point{ 1.3,-1.3, 0.5} },
    // outside, parallel to y-axis, 1 intersections (on surface, reverse)
    SimpleTrack{3*M_PI/2,2*M_PI/3, Point{ 1.3, 1.3, 0.5} },
    // outside, 1 intersections (on surface)
    SimpleTrack{ M_PI/4, M_PI/3, Point{0.,-1.3/cos45, 0.5} },

    // inside, parallel to z-axis, intersect at infinity
    SimpleTrack{      0,      0, Point{ 0., 0., 0.5} },
    // outside, parallel to z-axis, intersect at infinity
    SimpleTrack{   M_PI,      0, Point{ 2., 0., 0.5} },
    // on-surface, parallel to z-axis, no intersection
    SimpleTrack{      0,2*M_PI/3, Point{1.3, 0., 0.5} },
  };

  PyVector<int> oracle_num_intersections = {
    // a ray from an inside point can only have 1 intersection
    1,
    1,
    1,
    1,
    // outside, no intersection
    0,
    0,
    // a ray from an outside point can have two candidate intersections
    2,
    2,
    2,
    2,
    // on surface
    1,
    1,
    1,
    1,
    // parallel to z-axis
    0,
    0,
    0,
  };

  PyVector<PyVector<Point>> oracle_intersections = {
    // inside, parallel to x-axis
    {
      Point{1.3, 0., 0.5 + 1.3*tan30},
    },
    // inside, parallel to y-axis
    {
      Point{0., 1.3, 0.5 + 1.3*tan30},
    },
    // inside, parallel to y-axis (reverse)
    {
      Point{0.,-1.3, 0.5 - 1.3*tan30},
    },
    // inside
    {
      Point{1.3*cos45, 1.3*sin45, 0.5 + 1.3*tan30},
    },
    // outside, parallel to x-axis, no intersection
    {
      Point{},
    },
    // outside, no intersection
    {
      Point{},
    },
    // outside, parallel to x-axis, 2 intersections
    {
      Point{ 1.3, 0., 0.5 + 11.3*tan30},
      Point{-1.3, 0., 0.5 + 8.7*tan30},
    },
    // outside, parallel to y-axis, 2 intersections
    {
      Point{0., 1.3, 0.5 + 11.3*tan30},
      Point{0.,-1.3, 0.5 + 8.7*tan30},
    },
    // outside, parallel to y-axis, 2 intersections (reverse)
    {
      Point{0., 1.3, 0.5 - 8.7*tan30},
      Point{0.,-1.3, 0.5 - 11.3*tan30},
    },
    // outside, 2 intersections
    {
      Point{ 1.3*cos45, 1.3*sin45, 0.5 + (1.3/cos45 + 1.3)*tan30},
      Point{-1.3*cos45,-1.3*sin45, 0.5 + (1.3/cos45 - 1.3)*tan30},
    },
    // outside, parallel to x-axis, 1 intersections (on surface)
    {
      Point{0.,-1.3, 0.5 + 1.3*tan30},
    },
    // outside, parallel to y-axis, 1 intersections (on surface)
    {
      Point{1.3, 0., 0.5 + 1.3*tan30},
    },
    // outside, parallel to y-axis, 1 intersections (on surface, reverse)
    {
      Point{1.3, 0., 0.5 - 1.3*tan30},
    },
    // outside, 1 intersections (on surface)
    {
      Point{ 1.3*cos45,-1.3*sin45, 0.5 + 1.3*tan30},
    },
    // inside, parallel to z-axis, no intersection (intersect at infinity)
    {
      Point{},
    },
    // outside, parallel to z-axis, no intersection (intersect at infinity)
    {
      Point{},
    },
    // on-surface, parallel to z-axis, no intersection
    {
      Point{},
    },
  };

  Point intersections[2];
  for (size_t j = 0; j < tracks.size(); ++j) {
    SimpleTrack &track = tracks[j];

    int num_intersections = cylinder.intersection(&(track.point),
                                                  track.azim,
                                                  track.polar,
                                                  intersections);
    // number of intersections
    ASSERT_EQ(num_intersections, oracle_num_intersections[j])
      << "ZCylinder " << ":\n" << cylinder.toString() << '\n'
      << "  Track: " << j << ": " << track.toString() << '\n'
      << (num_intersections > 0 ? "  1: " + intersections[0].toString() : "")
      << (num_intersections > 1 ? "  2: " + intersections[1].toString() : "")
      << std::endl;

    // check the intersection
    for (int k = 0; k < num_intersections; ++k) {
      EXPECT_EQ(intersections[k], oracle_intersections[j][k])
        << "ZCylinder " << ":\n" << cylinder.toString() << '\n'
        << "  Track " << j << ": " << track.toString() << '\n'
        << "  Intersection " << k << ": " << intersections[k].toString() << '\n'
        << "  Oracle " << k << ": " << oracle_intersections[j][k].toString()
        << std::endl;
    }
  }

}

} // namespace
