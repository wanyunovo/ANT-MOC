/// \file test_fluxes_symmetry_3x3.cpp
/// \date 2019/12/10
/// \author An Wang, USTB (wangan@xs.ustb.edu.cn)
/// \purpose 回归测试。测试3x3棒束计算中的轨迹、线段数量和最终输出的标通量。
///
/// 问题描述：
/// 测试BEAVRS全堆时，裂变反应率与张庚提供的参考结果相差很大，分布规律不同。
/// 为测试代码中存在的问题，单独以91号组件为例，对比多个程序的计算结果：
///   1、ANT-MOC
///   2、OpenMOC
///   3、杨睿编写的MOC程序
/// ANT-MOC与OpenMOC对比C5G7单组件计算的keff和通量分布，发现结果相差很大。
/// 该测试的描述如下。
/// 该问题中长、宽、高均为3.78，轴向分为完全相等的2层。问题包含两种栅元：
/// 燃料栅元，控制棒栅元。所有材料均取自C5G7。该测试在程序错误未修复时可以
/// 准确复现以下问题：
///   1、线段数量增加，存在一些长度为0的线段；
///   2、标通量不对称，尤其在YMIN面靠近其与XMAX面交线的位置偏低。
///
/// 修复后程序与杨睿编写的程序对比，通量相对误差大约在千分之五以内(不同参数)。


#ifdef USTB_

#include "testing/test_harness.h"

using namespace antmoc;


namespace {

class test_fluxes_symmetry_3x3 : public TestHarness {

 public:

  test_fluxes_symmetry_3x3() {

    setTestDir("test_fluxes_symmetry_3x3");

    setArguments({
      "--log-level",    "debug",
      "--primitives",   getPublicGeometryDir() + "/c5g7/simple-lattice.xml",
      // Dumping arguments
      "--mesh",         "[1.26]*3, [1.26]*3, [3.78]",
      "--dump-rx",      "F,T,PHI",
      "--dump-tracks",  "3D",
      // Ray tracing arguments
      "-f", "OTF Stacks",
      "-z", "-1.89, 1.89",
      "-a", "4",
      "-s", "3.0",
      "-p", "2",
      "-l", "4.0",
      // Solver arguments
      "-q", "Equal Weight",
      "-i", "200",
    });

  }

};


/// \brief The basic test
/// \details Quadrature: Equal Weight
///          Tracing: OTF Stacks
TEST_F(test_fluxes_symmetry_3x3, EqualWeight) {

  runTest();

}


} // anonymous namespace

#endif  // USTB_
