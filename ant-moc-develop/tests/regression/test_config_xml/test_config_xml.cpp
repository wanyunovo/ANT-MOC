/// \file test_config_xml.cpp
/// \date 2020/03
/// \author An Wang, USTB (wangan@xs.ustb.edu.cn)
/// \purpose 测试XML格式输入卡。
///          XML格式输入卡的优先级低于命令行输入参数，因此会被
///          命令行的参数覆盖。该测试要检查这种优先级。


#ifdef USTB_

#include "testing/test_harness.h"

using namespace antmoc;


namespace {

class test_config_xml : public TestHarness {

 public:

  test_config_xml() {

    setTestDir("test_config_xml");

    resetArguments({
      "-c",             getTestDir() + "/settings.xml",
      "--geometry",     getTestDir() + "/geometry.xml",
      "--primitives",   getPublicGeometryDir() + "/c5g7/simple-lattice.xml",
      "--materials",    getPublicMaterialDir() + "c5g7/mgxs.h5",
      // Dumping arguments
      "--dump-rx",      "PHI",
      "--dump-tracks",  "3D",
      "--tally-groups", "0",
      // Ray tracing arguments
      "-f", "OTF Stacks",
      "-z", "-3.78, 0, 3.78",
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
      TestResultType::fluxes,
    //  TestResultType::mesh_data
    });
  }

};


TEST_F(test_config_xml, overwriteArguments) {

  runTest();

}


} // anonymous namespace

#endif  // USTB_
