/// \file test_reflector.cpp
/// \date 2020/03
/// \author An Wang, USTB (wangan@xs.ustb.edu.cn)
/// \purpose 回归测试。测试所有边界均为全反射边界的带反射层几何。
///          轨迹密度的设置保证追踪出全部FSR。


#ifdef USTB_

#include "testing/test_harness.h"

using namespace antmoc;


namespace {

class test_reflector : public TestHarness {

 public:

  test_reflector() {

    setTestDir("test_reflector");

    setArguments({
      "--primitives",   getPublicGeometryDir() + "/c5g7/simple-lattice.xml",
      // Dumping arguments
      "--mesh",         "[1.26]*6, [1.26]*6, [3.78]*2",
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
TEST_F(test_reflector, EqualWeight) {

  runTest();

}


/// \brief Testing tracing method OTF_Tracks
/// \details The result should be the same as the basic test
TEST_F(test_reflector, OTFTracks) {

  // This test is not allowed to overwrite the results
  setTestOverwrite(false);

  modifyArguments("-f", "OTF Tracks");

  std::string suffix = "otf_tracks";
  runTest(suffix);  // compared with the basic results

}


/// \details Whe modules are set, the layout of tracks will be
///          changed so that k-eff and fluxes may not be the
///          same as the basic test.
TEST_F(test_reflector, modules3x3) {

  modifyArguments("--modules", "3,3,1");

  std::string suffix = "modules_3x3";
  runTest(suffix, suffix);

}


/// \brief Testing computation method of k-eff
TEST_F(test_reflector, KeffFromNeutronBalance) {

  modifyArguments({"--keff-neutron-balance"});

  std::string suffix = "keff_neutron_balance";
  runTest(suffix, suffix);

}


TEST_F(test_reflector, EqualAngle) {

  // This test is not allowed to overwrite the results
  setTestOverwrite(false);

  modifyArguments("-q", "Equal Angle");

  std::string suffix = "equal_angle";
  runTest(suffix);  // compared with the basic results

}


TEST_F(test_reflector, GaussLegendre) {

  modifyArguments("-q", "Gauss Legendre");

  std::string suffix = "gauss_legendre";
  runTest(suffix, suffix);

}


TEST_F(test_reflector, TabuchiYamamoto) {

  modifyArguments("-q", "Tabuchi Yamamoto");

  std::string suffix = "tabuchi_yamamoto";
  runTest(suffix, suffix);

}


TEST_F(test_reflector, Leonard) {

  modifyArguments("-p", "4");
  modifyArguments("-q", "Leonard");

  std::string suffix = "leonard";
  runTest(suffix, suffix);

}


//----------------------------------------------------------------------
// Rodded reflectors
//----------------------------------------------------------------------
TEST_F(test_reflector, halfRodded) {

  modifyArguments("-g", getTestDir() + "/geometry_half_rodded.xml");

  std::string suffix = "half_rodded";
  runTest(suffix, suffix);

}


TEST_F(test_reflector, rodded) {

  modifyArguments("-g", getTestDir() + "/geometry_rodded.xml");

  std::string suffix = "rodded";
  runTest(suffix, suffix);

}


//----------------------------------------------------------------------
// Diffenrent boundary conditions
//----------------------------------------------------------------------
TEST_F(test_reflector, mixedBoundaryConditions) {

  modifyArguments("-g", getTestDir() + "/geometry_mixed_bc.xml");

  std::string suffix = "mixed_bc";
  runTest(suffix, suffix);

}


TEST_F(test_reflector, vacuumBoundaryConditions) {

  modifyArguments("-g", getTestDir() + "/geometry_vacuum_bc.xml");

  std::string suffix = "vacuum_bc";
  runTest(suffix, suffix);

}


} // anonymous namespace

#endif  // USTB_
