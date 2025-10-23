/// \file mpiTrack_test_boundary_reflective_3x3.cpp
/// \date 2020/03
/// \author An Wang, USTB (wangan@xs.ustb.edu.cn)
/// \purpose 回归测试。测试所有边界均为全反射边界的3x3棒束。
///          轨迹密度的设置保证追踪出全部216个FSR。


#ifdef USTB_

#include "testing/test_harness.h"

using namespace antmoc;


namespace {

class test_boundary_reflective_3x3 : public TestHarness {

 public:

  test_boundary_reflective_3x3() {

    setTestDir("test_boundary_reflective_3x3");

    setArguments({
      "--log-level",    "debug",
      "--num-domains",  "1,1,1,0",
      "--primitives",   getPublicGeometryDir() + "/c5g7/simple-lattice.xml",
      // Dumping arguments
      "--mesh",         "[1.26]*3, [1.26]*3, [3.78]*2",
      "--dump-rx",      "F,T,PHI",
      "--dump-tracks",  "3D",
      // Ray tracing arguments
      "-f", "OTF Stacks",
      "-z", "-3.78, 0, 3.78",
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
TEST_F(test_boundary_reflective_3x3, mpiTrack_EqualWeight) {

  setTestResultTypes({
    TestResultType::num_iterations,
    TestResultType::num_fsrs,
    TestResultType::keff
  });

  // This test is not allowed to overwrite the results
  setTestOverwrite(false);

  std::string suffix = "mpiTrack";
  runTest(suffix);  // compared with the basic results

}


} // anonymous namespace

#endif  // USTB_
