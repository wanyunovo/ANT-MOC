/// \brief Test the point is in the Cell or not when add the Set operation
/// \author An Wang, USTB (wangan@xs.ustb.edu.cn)

#include "testing/geometry_testharness.h"
#include "antmoc/Cell.h"


using namespace antmoc;

namespace {


// Test fixture
class test_Cell_containPoint: public GeometryTestHarness {

protected:

  void SetUp() {

    // Add materials
    _mat_input->setMaterial("1", new Material(0, "1"));

  }


  /// \brief A helper to test containsPoint
  void testCellContainsPoints(Cell *cell, std::vector<Point> in_points,
      std::vector<Point> out_points) {

    ASSERT_TRUE(cell);

    for (auto &p : in_points)
      EXPECT_TRUE(cell->containsPoint(&p))
        << " Cell id = " << cell->getInputId() << '\n'
        << " Expected to be inside the cell:\n" << p.toString();

    for (auto &p : out_points)
      EXPECT_FALSE(cell->containsPoint(&p))
        << " Cell id = " << cell->getInputId() << '\n'
        << " Expected to be outside the cell:\n" << p.toString();
  }

};


TEST_F(test_Cell_containPoint, singleHalfspaceSimple) {
  ROOT_ONLY();

  std::string test =
    " <cell id='1' material='1' name='large' region='-1'>"
    "   <surface id='1' type='z-cylinder' coeffs='0 0 1'/>" // the direction attribute can be added
    " </cell>";

  _doc.Parse(test.c_str());
  auto cell = buildCell(_doc.FirstChildElement("cell"));  // parse the cell

  std::vector<Point> in_points = {
   {0, 0, 0},
   {0.99, 0, 0},
   {-0.99, 0, 0},
   {0, 0.99, 0},
   {0, -0.99, 0},
   {0.5*std::sqrt(3), 0.49, 0},
  };

  std::vector<Point> out_points = {
   {1, 1, 0}
  };

  testCellContainsPoints(cell, in_points, out_points);

  delete cell;
}


// test the operator: complement(~)
TEST_F(test_Cell_containPoint, singleSurfaceComplement){
  ROOT_ONLY();

  std::string test =
    " <cell id='2' material='1' name='large' region='~(-2)'>" // actual the region = '+2'
    "   <surface id='2' type='z-cylinder' coeffs='0 0 1'/>"
    " </cell>";

  _doc.Parse(test.c_str());
  auto cell = buildCell(_doc.FirstChildElement("cell"));  // parse the cell

  std::vector<Point> in_points = {
   {1, 1, 0} // in
  };

  std::vector<Point> out_points = {
   {0, 0, 0}, // out
   {0.99, 0, 0},
   {-0.99, 0, 0},
   {0, 0.99, 0},
   {0, -0.99, 0},
   {0.5*std::sqrt(3), 0.49, 0},
  };

  testCellContainsPoints(cell, in_points, out_points);

  delete cell;
}


TEST_F(test_Cell_containPoint, localHalfspacesIntersection){
  ROOT_ONLY();

  const std::string test =
    " <cell id = '1' material='1' name='test1' region='8 -9'> "
    "   <surface id='8' type='z-cylinder' coeffs='0 0 1'/>"
    "   <surface id='9' type='z-cylinder' coeffs='0 0 2'/>"
    " </cell>";

  _doc.Parse(test.c_str());
  auto cell = buildCell(_doc.FirstChildElement("cell"));  // parse the cell

  std::vector<Point> in_points = {
    {1.999, 0, 0},
    {-1.999, 0, 0},
    {0, 1.999, 0},
    {0, -1.999, 0},
    {sqrt(3), 0.999, 0},
  };

  std::vector<Point> out_points = {
    {0, 0, 0},
    {0.999, 0, 0},
    {2.001, 0, 0},
    {sqrt(3), 1.001, 0}
  };

  testCellContainsPoints(cell, in_points, out_points);

  delete cell;
}


TEST_F(test_Cell_containPoint, globalHalfspacesSimple){
  ROOT_ONLY();

  const std::string test =
    " <global>"
    "   <surface id='1' type='z-cylinder' coeffs='0 0 1'/>"
    "   <cell id = '1' material='1' name='test1' region='+1 -2'>"
    "      <surface id='2' type='z-cylinder' coeffs='0 0 2'/>"
    "   </cell>"
    " </global>";

  _doc.Parse(test.c_str());
  _geo_input->buildGlobalPrimitives(_doc.FirstChildElement("global"));
  auto &cells = _geo_input->getGlobalCells();
  auto cell = cells[1];

  std::vector<Point> in_points = {
    {1.999, 0, 0},
    {-1.999, 0, 0},
    {0, 1.999, 0},
    {0, -1.999, 0},
    {sqrt(3), 0.999, 0},
  };

  std::vector<Point> out_points = {
    {0, 0, 0},
    {0.999, 0, 0},
    {2.001, 0, 0},
    {sqrt(3), 1.001, 0}
  };

  testCellContainsPoints(cell, in_points, out_points);
}


TEST_F(test_Cell_containPoint, globalHalfspacesComplex){
  ROOT_ONLY();

  const std::string test =
    " <global>"
    "   <surface id='3' type='x-plane' coeffs='-2'/>"
    "   <surface id='4' type='x-plane' coeffs='2'/>"
    "   <surface id='5' type='y-plane' coeffs='-2'/>"
    "   <surface id='6' type='y-plane' coeffs='2'/>"
    "   <surface id='7' type='z-cylinder' coeffs='0 0 1'/>"
    "   <cell id = '1' material='1' name='inner' region='-7'/>"
    "   <cell id = '2' material='1' name='middle' region='+3 -4 +5 -6 +7'/>"
    "   <cell id = '3' material='1' name='outer' region='~(+3 -4 +5 -6)'/>"
    " </global>";

  _doc.Parse(test.c_str());
  _geo_input->buildGlobalPrimitives(_doc.FirstChildElement("global"));
  auto &cells = _geo_input->getGlobalCells();

  // parse the global: complicated region
  std::vector<size_t> c_ids = {1, 2, 3};

  std::vector<Point> in_inner = {
    {0, 0, 0}, // in inner
    {0.999 ,0 ,0},
    {0.5*sqrt(3), 0.499, 0},
  };

  std::vector<Point> in_middle = {
    {1.001, 0, 0}, // in middle
    {0.1, -1.001, 0},
    {1.999, 1.999, 0},
  };

  std::vector<Point> in_outer = {
    {2.001, 0, 0}, // in outer
    {0, 2.001, 0},
    {3, 3, 0}
  };

  std::vector<std::vector<Point>> points =
    {in_inner, in_middle, in_outer};

  for (size_t i = 0; i < c_ids.size(); ++i) {
    auto id = c_ids[i];
    Cell *cell = cells[id];

    testCellContainsPoints(cell, points[i], std::vector<Point>{});
  }
}

} /** namespace */
