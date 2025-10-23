/**
 * @brief Test Geometry reading and building process
 * @date Aug 7, 2019
 * @author An Wang, USTB (wangan@xs.ustb.edu.cn)
 */

#include "testing/test_utils.h"
#include "testing/MockGeometry.h"
#include "antmoc/Lattice.h"

#include <functional>
#include <vector>

using namespace antmoc;

namespace {

using DoubleVec = std::vector<double>;

// Test fixture
class test_Lattice_variableWidths : public testing::Test {
  // do nothing
};


TEST_F(test_Lattice_variableWidths, uniformRecLattice2x2x2) {
  ROOT_ONLY();
  RecLattice lattice;
  
  /* Construct a 1.2*1.2*1.2 lattice */
  lattice.setWidth(0.6, 0.6, 0.6);
  lattice.setNumX(2);
  lattice.setNumY(2);
  lattice.setNumZ(2);
  lattice.computeSizes();
  /* Boundaries */
  EXPECT_EQ(-0.6, lattice.getMinX());
  EXPECT_EQ(0.6, lattice.getMaxX());
  EXPECT_EQ(-0.6, lattice.getMinY());
  EXPECT_EQ(0.6, lattice.getMaxY());
  EXPECT_EQ(-0.6, lattice.getMinZ());
  EXPECT_EQ(0.6, lattice.getMaxZ());
  /* Points */
  Point pc(0., 0., 0.);
  Point px1(-0.6, 0., 0.), px2(-0.6, -0.01, 0.);
  Point py1(0., -0.6, 0.), py2(0., -0.6, -0.01);
  Point pz1(0., 0., -0.6), pz2(-0.01, 0., -0.6);
  Point pout(1.2, 1.2, 1.2);
  EXPECT_TRUE(lattice.containsPoint(&px1));
  EXPECT_TRUE(lattice.containsPoint(&px2));
  EXPECT_TRUE(lattice.containsPoint(&py1));
  EXPECT_TRUE(lattice.containsPoint(&py2));
  EXPECT_TRUE(lattice.containsPoint(&pz1));
  EXPECT_TRUE(lattice.containsPoint(&pz2));
  EXPECT_FALSE(lattice.containsPoint(&pout));
  /* Lattice cells */
  EXPECT_EQ(1, lattice.getLatX(&pc));
  EXPECT_EQ(1, lattice.getLatY(&pc));
  EXPECT_EQ(1, lattice.getLatZ(&pc));

  EXPECT_EQ(0, lattice.getLatX(&px1));
  EXPECT_EQ(1, lattice.getLatY(&px1));
  EXPECT_EQ(1, lattice.getLatZ(&px1));

  EXPECT_EQ(0, lattice.getLatX(&px2));
  EXPECT_EQ(0, lattice.getLatY(&px2));
  EXPECT_EQ(1, lattice.getLatZ(&px2));

  EXPECT_EQ(1, lattice.getLatX(&py1));
  EXPECT_EQ(0, lattice.getLatY(&py1));
  EXPECT_EQ(1, lattice.getLatZ(&py1));

  EXPECT_EQ(1, lattice.getLatX(&py2));
  EXPECT_EQ(0, lattice.getLatY(&py2));
  EXPECT_EQ(0, lattice.getLatZ(&py2));

  EXPECT_EQ(1, lattice.getLatX(&pz1));
  EXPECT_EQ(1, lattice.getLatY(&pz1));
  EXPECT_EQ(0, lattice.getLatZ(&pz1));

  EXPECT_EQ(0, lattice.getLatX(&pz2));
  EXPECT_EQ(1, lattice.getLatY(&pz2));
  EXPECT_EQ(0, lattice.getLatZ(&pz2));
}


TEST_F(test_Lattice_variableWidths, nonUniformRecLattice2x2x2) {
  ROOT_ONLY();
  RecLattice lattice;
  
  /* Construct a 1.2*1.2*1.2 lattice */
  DoubleVec widths = {0.3, 0.9};
  lattice.setWidths(widths, widths, widths);
  lattice.setNumX(2);
  lattice.setNumY(2);
  lattice.setNumZ(2);
  lattice.computeSizes();
  /* Boundaries */
  EXPECT_EQ(-0.6, lattice.getMinX());
  EXPECT_EQ(0.6, lattice.getMaxX());
  EXPECT_EQ(-0.6, lattice.getMinY());
  EXPECT_EQ(0.6, lattice.getMaxY());
  EXPECT_EQ(-0.6, lattice.getMinZ());
  EXPECT_EQ(0.6, lattice.getMaxZ());
  /* Points */
  Point pc(-0.3, -0.3, -0.3);
  Point px1(-0.6, -0.3, -0.3), px2(-0.6, -0.31, -0.3);
  Point py1(-0.3, -0.6, -0.3), py2(-0.3, -0.6, -0.31);
  Point pz1(-0.3, -0.3, -0.6), pz2(-0.31, -0.3, -0.6);
  Point pout(1.2, 1.2, 1.2);
  EXPECT_TRUE(lattice.containsPoint(&px1));
  EXPECT_TRUE(lattice.containsPoint(&px2));
  EXPECT_TRUE(lattice.containsPoint(&py1));
  EXPECT_TRUE(lattice.containsPoint(&py2));
  EXPECT_TRUE(lattice.containsPoint(&pz1));
  EXPECT_TRUE(lattice.containsPoint(&pz2));
  EXPECT_FALSE(lattice.containsPoint(&pout));
  /* Lattice cells */
  EXPECT_EQ(1, lattice.getLatX(&pc));
  EXPECT_EQ(1, lattice.getLatY(&pc));
  EXPECT_EQ(1, lattice.getLatZ(&pc));

  EXPECT_EQ(0, lattice.getLatX(&px1));
  EXPECT_EQ(1, lattice.getLatY(&px1));
  EXPECT_EQ(1, lattice.getLatZ(&px1));

  EXPECT_EQ(0, lattice.getLatX(&px2));
  EXPECT_EQ(0, lattice.getLatY(&px2));
  EXPECT_EQ(1, lattice.getLatZ(&px2));

  EXPECT_EQ(1, lattice.getLatX(&py1));
  EXPECT_EQ(0, lattice.getLatY(&py1));
  EXPECT_EQ(1, lattice.getLatZ(&py1));

  EXPECT_EQ(1, lattice.getLatX(&py2));
  EXPECT_EQ(0, lattice.getLatY(&py2));
  EXPECT_EQ(0, lattice.getLatZ(&py2));

  EXPECT_EQ(1, lattice.getLatX(&pz1));
  EXPECT_EQ(1, lattice.getLatY(&pz1));
  EXPECT_EQ(0, lattice.getLatZ(&pz1));

  EXPECT_EQ(0, lattice.getLatX(&pz2));
  EXPECT_EQ(1, lattice.getLatY(&pz2));
  EXPECT_EQ(0, lattice.getLatZ(&pz2));
}

}/* namespace */
