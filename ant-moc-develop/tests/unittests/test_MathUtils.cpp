/// \file test_MathUtils.cpp
/// \brief Test math utility
/// \date Dec 17, 2019
/// \author An Wang, USTB (wangan@xs.ustb.edu.cn)

#include "testing/test_utils.h"
#include "antmoc/math_utils.h"

using namespace antmoc;

namespace {

/// Testing fixture
class test_MathUtils : public testing::Test {
};

TEST_F(test_MathUtils, greatestCommonDivisor) {
  ROOT_ONLY();

  EXPECT_EQ(antmoc::gcd(1, 2), 1);
  EXPECT_EQ(antmoc::gcd(2, 1), 1);
  EXPECT_EQ(antmoc::gcd(2, 2), 2);
  EXPECT_EQ(antmoc::gcd(6, 12), 6);
  EXPECT_EQ(antmoc::gcd(6, 15), 3);
}

#ifndef ENABLE_MPI_
/// Non-positive numbers are not supported
TEST_F(test_MathUtils, greatestCommonDivisorRobust) {

  EXPECT_ANY_THROW(antmoc::gcd(1, -1));
  EXPECT_ANY_THROW(antmoc::gcd(-1, 1));
  EXPECT_ANY_THROW(antmoc::gcd(-1, -1));
  EXPECT_ANY_THROW(antmoc::gcd(0, 1));

}
#endif  // ENABLE_MPI_


TEST_F(test_MathUtils, definitelyEqual) {

  double infinity  = std::numeric_limits<double>::infinity();
  float f_infinity = std::numeric_limits<float>::infinity();

  EXPECT_TRUE(definitelyEqual(1.0, 1.0 + FLT_EPSILON / 10));
  EXPECT_TRUE(definitelyEqual(1.0 + FLT_EPSILON / 10, 1.0));
  EXPECT_TRUE(definitelyEqual(infinity, infinity));
  EXPECT_TRUE(definitelyEqual(-infinity, -infinity));
  EXPECT_TRUE(definitelyEqual(INFINITY, INFINITY));
  EXPECT_TRUE(definitelyEqual(-INFINITY, -INFINITY));
  EXPECT_TRUE(definitelyEqual(INFINITY, f_infinity));
  EXPECT_TRUE(definitelyEqual(-INFINITY, -f_infinity));

  EXPECT_FALSE(definitelyEqual(1.0, 1.0 + FLT_EPSILON * 10));
  EXPECT_FALSE(definitelyEqual(1.0 + FLT_EPSILON * 10, 1.0));
  EXPECT_FALSE(definitelyEqual(infinity, -infinity));
  EXPECT_FALSE(definitelyEqual(-infinity, infinity));
  EXPECT_FALSE(definitelyEqual(INFINITY, -INFINITY));
  EXPECT_FALSE(definitelyEqual(-INFINITY, INFINITY));
  EXPECT_FALSE(definitelyEqual(INFINITY, -f_infinity));
  EXPECT_FALSE(definitelyEqual(-INFINITY, f_infinity));
}

} // namespace antmoc
