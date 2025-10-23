/// \date March 10, 2020
/// \author An Wang, USTB (wangan@xs.ustb.edu.cn)

#include "testing/test_utils.h"
#include "antmoc/log.h"
#include "antmoc/Point.h"

using namespace antmoc;

namespace {

// Test fixture
class test_Point_arithmetic : public testing::Test { };


TEST_F(test_Point_arithmetic, computedAssignment) {
  ROOT_ONLY();

  Point pt {0., 0., 0.}, pt1 {1.0001, 2.000002, 3.00000003};

  pt += pt1;

  EXPECT_EQ(pt.getX(), 1.0001);
  EXPECT_EQ(pt.getY(), 2.000002);
  EXPECT_EQ(pt.getZ(), 3.00000003);

  pt -= pt1;
  EXPECT_EQ(pt.getX(), 0.);
  EXPECT_EQ(pt.getY(), 0.);
  EXPECT_EQ(pt.getZ(), 0.);
}


TEST_F(test_Point_arithmetic, binaryOperation) {
  ROOT_ONLY();

  Point pt1, pt2;
  Point pt3 {0., 1., 2.}, pt4 {2., 1., 0.};

  pt1 = pt3 + pt4;
  pt2 = pt4 + pt3;

  EXPECT_EQ(pt1.getX(), 2.);
  EXPECT_EQ(pt1.getY(), 2.);
  EXPECT_EQ(pt1.getZ(), 2.);

  EXPECT_EQ(pt2.getX(), pt1.getX());
  EXPECT_EQ(pt2.getY(), pt1.getY());
  EXPECT_EQ(pt2.getZ(), pt1.getZ());

}


/// \brief Test equality operator
TEST_F(test_Point_arithmetic, equality) {
  ROOT_ONLY();

  Point pt1 {1., 2., 3.};
  Point pt2 {1., 1.9999999 + 0.0000001, 1. + 2.};
  Point pt3 {1., 2., 1. + 2.00000001};

  EXPECT_EQ(pt1, pt2);
  EXPECT_NE(pt1, pt3);
}

} // namespace
