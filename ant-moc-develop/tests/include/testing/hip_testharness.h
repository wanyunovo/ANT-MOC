/// \file   hip_testharness.h
/// \brief  A test harness for ANT-MOC
/// \date   December 31, 2020
/// \author An Wang, USTB (wangan@xs.ustb.edu.cn)

#ifndef HIP_TESTHARNESS_H_
#define HIP_TESTHARNESS_H_

#include "testing/test_utils.h"
#include "antmoc/hip/hip_utils.h"
#include "antmoc/log.h"

using namespace antmoc;

namespace {

class HipTestHarness : public testing::Test {

public:

  HipTestHarness() {
    log::set_level("verbose");
  }

  void SetUp() override {
    // Set device id by the MPI rank.
    hip::setDevice(hip::getGPUId());
  }

  void TearDown() override {
    HIP_ASSERT( hipDeviceSynchronize() );
  }

};  // class HipTestHarness

} // namespace

#endif  // HIP_TESTHARNESS_H_
