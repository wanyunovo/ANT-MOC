/// \file test_Cell_distance.cpp
/// \brief Test minimum distance from points to cells
/// \author An Wang, USTB (wangan@xs.ustb.edu.cn)

#include "testing/geometry_testharness.h"
#include "antmoc/Cell.h"


using namespace antmoc;

namespace {


// Test fixture
class test_Cell_distance: public GeometryTestHarness {

protected:

  void SetUp() {

    // Add materials
    _mat_input->setMaterial("water", new Material(0, "water"));

  }

};


/// \brief Test convex polygon
TEST_F(test_Cell_distance, inConvexRegion) {
  ROOT_ONLY();

  std::string test =
    "<global>"
    "  <surface id='1' type='x-plane' coeffs='1'/>"
    "  <surface id='2' type='x-plane' coeffs='-1'/>"
    "  <surface id='3' type='y-plane' coeffs='1'/>"
    "  <surface id='4' type='y-plane' coeffs='-1'/>"
    "  <surface id='5' type='z-plane' coeffs='1'/>"
    "  <surface id='6' type='z-plane' coeffs='-1'/>"
    "  <cell id='1' material='water' name='convex' region='-1 2 -3 4 -5 6'/>"
    "</global>";

  _doc.Parse(test.c_str());
  buildGlobalPrimitives(_doc.FirstChildElement("global"));

  auto &cells = _geo_input->getGlobalCells();
  ASSERT_EQ(cells.size(), 1);

  // The cell to be tested
  auto cell = cells[1];

  std::vector<SimpleTrack> tracks = {
    // parallel to the x axis
    SimpleTrack{       0, M_PI/2, Point{-0.5, 0., 0.}},
    // parallel to the y axis
    SimpleTrack{3*M_PI/2, M_PI/2, Point{0., 0.5, 0.}},
    // parallel to the z axis
    SimpleTrack{       0,      0, Point{0., 0., 0.}},
    // to the upper right corner
    SimpleTrack{  M_PI/4, M_PI/2, Point{0., 0., 0.}},
  };

  const std::vector<double> oracles = {
    1.5,
    1.5,
    1.0,
    1.0*std::sqrt(2),
  };

  for (size_t i = 0; i < tracks.size(); ++i) {
    auto &track = tracks[i];
    auto distance = cell->minSurfaceDist(&(track.point),
                                         track.azim,
                                         track.polar);
    EXPECT_DOUBLE_EQ(distance, oracles[i])
      << "Cell:\n" << cell->toString() << '\n'
      << "  Track " << i << ":\n" << track.toString() << '\n'
      << std::endl;
  }
}


/// \brief Test non-convex region
/// \details Non-convex region is not fully supported yet. The result of
///          minSurfaceDist(...) could be surprising. It is harmless to
///          treat non-convex regions as convex regions in most of cases.
///          For example, we usually want to find a cell for a point and
///          then find the next cell where the point will move. If the cell
///          contains a non-convex region, say, the outer space of a rectangle
///              ~(-1 2 -3 4)
///          the moving point may be stopped by some 'non-existent' boundaries
///          which are extra parts of the underlying surfaces. The result is
///          that the point remains in the same cell but a new segment has
///          been created. It seems like the point is on a superposition plane.
TEST_F(test_Cell_distance, inNonConvexRegion) {
  ROOT_ONLY();

  std::string test =
    "<global>"
    "  <surface id='1' type='x-plane' coeffs='1'/>"
    "  <surface id='2' type='x-plane' coeffs='-1'/>"
    "  <surface id='3' type='y-plane' coeffs='1'/>"
    "  <surface id='4' type='y-plane' coeffs='-1'/>"
    "  <surface id='5' type='z-plane' coeffs='1'/>"
    "  <surface id='6' type='z-plane' coeffs='-1'/>"
    "  <cell id='1' material='water' name='convex' region='~(-1 2 -3 4 -5 6)'/>"
    "</global>";

  _doc.Parse(test.c_str());
  buildGlobalPrimitives(_doc.FirstChildElement("global"));

  auto &cells = _geo_input->getGlobalCells();
  ASSERT_EQ(cells.size(), 1);

  // The cell to be tested
  auto cell = cells[1];

  std::vector<SimpleTrack> tracks = {
    // parallel to the x axis
    SimpleTrack{       0, M_PI/2, Point{-1.5, 0., 0.}},
    // parallel to the y axis
    SimpleTrack{3*M_PI/2, M_PI/2, Point{0., 1.5, 0.}},
    // parallel to the z axis
    SimpleTrack{       0,      0, Point{0., 0., -1.5}},
    // to the left boundary (actually being truncated by the bottom plane)
    SimpleTrack{std::atan(3), M_PI/2, Point{-1.5,-1.5, 0.}},
  };

  const std::vector<double> oracles = {
    0.5,
    0.5,
    0.5,
    std::sqrt(2.5)/3
  };

  for (size_t i = 0; i < tracks.size(); ++i) {
    auto &track = tracks[i];
    auto distance = cell->minSurfaceDist(&(track.point),
                                         track.azim,
                                         track.polar);
    EXPECT_DOUBLE_EQ(distance, oracles[i])
      << "Cell:\n" << cell->toString() << '\n'
      << "  Track " << i << ":\n" << track.toString() << '\n'
      << std::endl;
  }
}


} // namespace
