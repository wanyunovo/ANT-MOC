/// \file test_Lattice_tally_mesh.cpp
/// \date 2020/03
/// \author An Wang, USTB (wangan@xs.ustb.edu.cn)
/// \purpose 测试只有部分控制棒的几何，输出tallymesh。
///          轨迹密度的设置保证追踪出全部864个FSR。


#ifdef USTB_

#include "testing/test_harness.h"

using namespace antmoc;


namespace {

class test_Lattice_tally_mesh : public TestHarness {

 public:

  test_Lattice_tally_mesh() {

    setTestDir("test_Lattice_tally_mesh");

    setArguments({
      "--primitives",   getPublicGeometryDir() + "/c5g7/simple-lattice.xml",
      // Dumping arguments
      "--mesh",         "[1.26]*6, [1.26]*6, [3.78]*2",
      "--dump-rx",      "PHI",
      "--dump-tracks",  "3D",
      "--tally-groups", "0",
      // Ray tracing arguments
      "-f", "OTF Stacks",
      "-z", "-3.78, 0, 3.78",
      "-a", "8",
      "-s", "0.75",
      "-p", "2",
      "-l", "3.0",
      // Solver arguments
      "-q", "Equal Weight",
      "-i", "200",
    });

    setTestResultTypes({
      TestResultType::num_iterations,
      TestResultType::num_fsrs,
      TestResultType::num_tracks,
      TestResultType::num_chains,
      TestResultType::num_segments,
      TestResultType::keff,
    //  TestResultType::fluxes,
      TestResultType::mesh_data
    });
  }

};


/// \brief The basic test
TEST_F(test_Lattice_tally_mesh, tallyFluxes) {

  runTest("fluxes");

}


/// \brief The basic test
TEST_F(test_Lattice_tally_mesh, tallyAllRates) {

  modifyArguments("--dump-rx", "ALL");
  runTest("all_rates");

}


/// \brief The basic test
TEST_F(test_Lattice_tally_mesh, tallySelectedGroups) {

  modifyArguments("--tally-groups", "1,2:3");
  runTest("selected_groups");

}


} // anonymous namespace

#endif  // USTB_
