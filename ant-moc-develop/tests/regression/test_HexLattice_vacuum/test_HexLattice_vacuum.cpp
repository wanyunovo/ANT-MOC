/// \file test_HexLattice_vacuum.cpp
/// \date 2020/03
/// \author An Wang, USTB (wangan@xs.ustb.edu.cn)
/// \purpose 回归测试。测试简单六边形网格的反应率分布。
///          轨迹密度的设置保证追踪出全部64个FSR，其中棒束部分
///          56个FSR，棒束外8个FSR。
///          在轨迹密度合适的情况下，计算得到的反应率应该具有
///          对称性，周围6个栅元的反应率应该基本相同。反应率
///          分布对求积组类型非常敏感，有的求积组会导致计算结
///          果不对称。


#ifdef USTB_

#include "testing/test_harness.h"

using namespace antmoc;


namespace {

class test_HexLattice_vacuum : public TestHarness {

 public:

  test_HexLattice_vacuum() {

    setTestDir("test_HexLattice_vacuum");

    setArguments({
      // Dumping arguments
      "--mesh-type",    "Hexagon",
      "--mesh",         "2, 0.695, [1.89]",
      "--orientation",  "y",
      "--dump-rx",      "F,T,PHI",
      "--dump-tracks",  "3D",
      // Ray tracing arguments
      "-f", "OTF Stacks",
      "-z", "-0.945, 0.945",
      "-a", "16",
      "-s", "0.05",
      "-p", "2",
      "-l", "0.5",
      // Solver arguments
      "-q", "Equal Weight",
      "-i", "200",
    });

  }

};


/// \brief The basic test
/// \details Quadrature: Equal Weight
///          Tracing: OTF Stacks
TEST_F(test_HexLattice_vacuum, EqualWeight) {

  runTest();

}


/// \brief Testing tracing method OTF_Tracks
/// \details The result should be the same as the basic test
TEST_F(test_HexLattice_vacuum, OTFTracks) {

  // This test is not allowed to overwrite the results
  setTestOverwrite(false);

  modifyArguments("-f", "OTF Tracks");

  std::string suffix = "otf_tracks";
  runTest(suffix);  // compared with the basic results
}


TEST_F(test_HexLattice_vacuum, EqualAngle) {

  // This test is not allowed to overwrite the results
  setTestOverwrite(false);

  modifyArguments("-q", "Equal Angle");

  std::string suffix = "equal_angle";
  runTest(suffix);  // compared with the basic results

}


TEST_F(test_HexLattice_vacuum, GaussLegendre) {

  modifyArguments("-q", "Gauss Legendre");

  std::string suffix = "gauss_legendre";
  runTest(suffix, suffix);

}


TEST_F(test_HexLattice_vacuum, TabuchiYamamoto) {

  modifyArguments("-q", "Tabuchi Yamamoto");

  std::string suffix = "tabuchi_yamamoto";
  runTest(suffix, suffix);

}


TEST_F(test_HexLattice_vacuum, Leonard) {

  modifyArguments("-p", "4");
  modifyArguments("-q", "Leonard");

  std::string suffix = "leonard";
  runTest(suffix, suffix);

}


} // anonymous namespace

#endif  // USTB_
