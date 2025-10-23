/// \file test_boundary_vacuum_3x3.cpp
/// \date 2020/03
/// \author An Wang, USTB (wangan@xs.ustb.edu.cn)
/// \purpose 回归测试。测试全真空边界的3x3棒束。
///          轨迹密度的设置保证追踪出全部216个FSR。


#ifdef USTB_

#include "testing/test_harness.h"

using namespace antmoc;


namespace {

class test_boundary_vacuum_3x3 : public TestHarness {

 public:

  test_boundary_vacuum_3x3() {

    setTestDir("test_boundary_vacuum_3x3");

    setArguments({
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
///          Tracing: OTF Stacks
TEST_F(test_boundary_vacuum_3x3, EqualWeight) {

  runTest();

}


/// \brief Testing tracing method OTF_Tracks
/// \details The result should be the same as the basic test
TEST_F(test_boundary_vacuum_3x3, OTFTracks) {

  // This test is not allowed to overwrite the results
  setTestOverwrite(false);

  modifyArguments("-f", "OTF Tracks");

  std::string suffix = "otf_tracks";
  runTest(suffix);  // compared with the basic results

}


/// \details Whe modules are set, the layout of tracks will be
///          changed so that k-eff and fluxes may not be the
///          same as the basic test.
TEST_F(test_boundary_vacuum_3x3, modules3x3) {

  modifyArguments("--modules", "3,3,1");

  std::string suffix = "modules_3x3";
  runTest(suffix, suffix);

}

/// \brief Testing computation method of k-eff
TEST_F(test_boundary_vacuum_3x3, KeffFromNeutronBalance) {

  // This test is not allowed to overwrite the results
  setTestOverwrite(false);

  setTestResultTypes({
    TestResultType::num_fsrs,
    TestResultType::num_tracks,
    TestResultType::num_chains,
    TestResultType::num_segments,
    TestResultType::keff,
    TestResultType::fluxes
  });

  modifyArguments({"--keff-neutron-balance"});

  std::string suffix = "keff_neutron_balance";
  runTest(suffix);  // compared with the basic results

}


TEST_F(test_boundary_vacuum_3x3, EqualAngle) {

  // This test is not allowed to overwrite the results
  setTestOverwrite(false);

  modifyArguments("-q", "Equal Angle");

  std::string suffix = "equal_angle";
  runTest(suffix);  // compared with the basic results

}


TEST_F(test_boundary_vacuum_3x3, GaussLegendre) {

  modifyArguments("-q", "Gauss Legendre");

  std::string suffix = "gauss_legendre";
  runTest(suffix, suffix);

}


TEST_F(test_boundary_vacuum_3x3, TabuchiYamamoto) {

  modifyArguments("-q", "Tabuchi Yamamoto");

  std::string suffix = "tabuchi_yamamoto";
  runTest(suffix, suffix);
}


TEST_F(test_boundary_vacuum_3x3, Leonard) {

  modifyArguments("-p", "4");
  modifyArguments("-q", "Leonard");

  std::string suffix = "leonard";
  runTest(suffix, suffix);

}


/// \brief Testing 'Compressed' XS file layout.
TEST_F(test_boundary_vacuum_3x3, compressedXSLayout) {

  // This test is not allowed to overwrite the results
  setTestOverwrite(false);

  modifyArguments("-m", getPublicMaterialDir() + "/c5g7/mgxs-compressed.h5");
  modifyArguments("--xs-layout", "Compressed");

  std::string suffix = "compressed_xs_layout";
  runTest(suffix);  // compared with the basic results

}


} // anonymous namespace

#endif  // USTB_
