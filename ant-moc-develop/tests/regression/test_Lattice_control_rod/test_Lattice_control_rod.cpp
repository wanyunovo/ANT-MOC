/// \file test_Lattice_control_rod.cpp
/// \date 2020/03
/// \author An Wang, USTB (wangan@xs.ustb.edu.cn)
/// \purpose 测试只有部分控制棒的几何
///          轨迹密度的设置保证追踪出全部432个FSR。
///          在合适的轨迹分布下，该测试应该产生具有对称性的结果。


#ifdef USTB_

#include "testing/test_harness.h"

using namespace antmoc;


namespace {

class test_Lattice_control_rod : public TestHarness {

 public:

  test_Lattice_control_rod() {

    setTestDir("test_Lattice_control_rod");

    setArguments({
      "--primitives",   getPublicGeometryDir() + "/c5g7/simple-lattice.xml",
      // Dumping arguments
      "--mesh",         "[1.26]*3, [1.26]*3, [3.78]*4",
      "--dump-rx",      "F,T,PHI",
      "--dump-tracks",  "3D",
      // Ray tracing arguments
      "-f", "OTF Stacks",
      "-z", "-7.56, -3.78, 0, 3.78, 7.56",
      "-a", "16",
      "-s", "0.25",
      "-p", "2",
      "-l", "3.0",
      // Solver arguments
      "-q", "Equal Weight",
      "-i", "200",
    });

  }

};


/// \brief The basic test
/// \details Quadrature: Equal Weight
TEST_F(test_Lattice_control_rod, unrodded) {

  modifyArguments("-g", getTestDir() + "/geometry_unrodded.xml");

  std::string suffix = "unrodded";
  runTest(suffix, suffix);

}


TEST_F(test_Lattice_control_rod, roddedQuarter) {

  modifyArguments("-g", getTestDir() + "/geometry_quarter.xml");

  std::string suffix = "quarter";
  runTest(suffix, suffix);

}


TEST_F(test_Lattice_control_rod, roddedTwoQuarters) {

  modifyArguments("-g", getTestDir() + "/geometry_two_quarters.xml");

  std::string suffix = "two_quarters";
  runTest(suffix, suffix);

}


TEST_F(test_Lattice_control_rod, roddedThreeQuarters) {

  modifyArguments("-g", getTestDir() + "/geometry_three_quarters.xml");

  std::string suffix = "three_quarters";
  runTest(suffix, suffix);

}

} // anonymous namespace

#endif  // USTB_
