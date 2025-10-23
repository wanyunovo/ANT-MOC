/**
 * @date Aug 8, 2019
 * @author An Wang, USTB (wangan@xs.ustb.edu.cn)
 */

#include "testing/test_utils.h"
#include "antmoc/log.h"
#include "antmoc/PyVector.h"

using namespace antmoc;

namespace {

// Test fixture
class test_PyVector : public testing::Test {
  // do nothing
};


TEST_F(test_PyVector, computedAssignment) {
  ROOT_ONLY();

  PyVector<double> pv, pv0;
  PyVector<double> pv1{1.};

  pv += pv0;
  EXPECT_EQ(0, pv0.size());

  pv += pv1;
  EXPECT_EQ(1, pv.size());
  EXPECT_EQ(1.0, pv[0]);

  pv.clear();
  pv *= 2;
  EXPECT_EQ(0, pv.size());

  pv = {1., 2.};
  pv *= 0;
  EXPECT_EQ(0, pv.size());

  pv = {1., 2.};
  pv *= 2;
  EXPECT_EQ(4, pv.size());
  EXPECT_EQ(1.0, pv[0]);
  EXPECT_EQ(2.0, pv[1]);
  EXPECT_EQ(1.0, pv[2]);
  EXPECT_EQ(2.0, pv[3]);
}


TEST_F(test_PyVector, binaryOperation) {
  ROOT_ONLY();

  PyVector<double> pv, pv0;
  PyVector<double> pv1{1.}, pv2{2., 3.};

  EXPECT_EQ(0, (pv + pv0).size());

  auto v = pv + pv1;
  EXPECT_EQ(1, v.size());
  EXPECT_EQ(1.0, v[0]);

  v = pv1 + pv2;
  EXPECT_EQ(3, v.size());
  EXPECT_EQ(1.0, v[0]);
  EXPECT_EQ(2.0, v[1]);
  EXPECT_EQ(3.0, v[2]);

  EXPECT_EQ(0, (pv * 2).size());
  EXPECT_EQ(0, (pv1 * 0).size());

  pv = {1., 2.};
  v = pv * 2;
  EXPECT_EQ(4, v.size());
  EXPECT_EQ(1.0, v[0]);
  EXPECT_EQ(2.0, v[1]);
  EXPECT_EQ(1.0, v[2]);
  EXPECT_EQ(2.0, v[3]);
}

TEST_F(test_PyVector, commutativity) {
  ROOT_ONLY();

  PyVector<double> pv1{1.}, pv2{2.};

  auto v = pv1 + pv2;
  EXPECT_EQ(1.0, v[0]);
  EXPECT_EQ(2.0, v[1]);

  v = pv2 + pv1;
  EXPECT_EQ(2.0, v[0]);
  EXPECT_EQ(1.0, v[1]);

  v = pv1 + 2.0;
  EXPECT_EQ(1.0, v[0]);
  EXPECT_EQ(2.0, v[1]);

  v = 2.0 + pv1;
  EXPECT_EQ(2.0, v[0]);
  EXPECT_EQ(1.0, v[1]);
}

TEST_F(test_PyVector, implicitConversion) {
  ROOT_ONLY();

  PyVector<double> pv1{1.};

  auto v = pv1 + 2;
  EXPECT_EQ(1.0, v[0]);
  EXPECT_EQ(2.0, v[1]);

  v = pv1 + 2.0 + 3 + 4.0;
  EXPECT_EQ(1.0, v[0]);
  EXPECT_EQ(2.0, v[1]);
  EXPECT_EQ(3.0, v[2]);
  EXPECT_EQ(4.0, v[3]);

  std::vector<double> vec{1., 2.};
  PyVector<double> pv2(vec);

  v = pv2 + 3.;
  EXPECT_EQ(1.0, v[0]);
  EXPECT_EQ(2.0, v[1]);
  EXPECT_EQ(3.0, v[2]);
}

}/* namespace */
