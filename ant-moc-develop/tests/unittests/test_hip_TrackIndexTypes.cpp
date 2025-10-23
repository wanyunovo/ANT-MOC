/// \file   test_hip_TrackIndexTypes.cpp
/// \brief  Test track indexing types.
/// \date   January 2, 2020
/// \author An Wang, USTB (wangan@xs.ustb.edu.cn)

#include "testing/hip_testharness.h"
#include "antmoc/TrackIndexTypes.h"

using namespace antmoc;

namespace {

/// Testing fixture
class test_hip_TrackIndexTypes : public HipTestHarness {};


TEST_F(test_hip_TrackIndexTypes, chainIndicesCopyHtoD) {
  // Prepare a host object
  TrackChainIndexes tci {1, 73, 2, 91, 11};

  // Allocate memory for the device object
  TrackChainIndexes *tci_d;
  HIP_ASSERT( hipMalloc(&tci_d, sizeof(tci)) );

  // Copy the host object to device
  HIP_ASSERT( hipMemcpyHtoD(tci_d, &tci, sizeof(tci)) );

  // Free device memory
  HIP_ASSERT( hipFree(tci_d) );

  SUCCEED();
}


TEST_F(test_hip_TrackIndexTypes, stackIndicesCopyHtoD) {
  // Prepare a host object
  TrackStackIndexes tsi {1, 133, 2, 91};

  // Allocate memory for the device object
  TrackStackIndexes *tsi_d;
  HIP_ASSERT( hipMalloc(&tsi_d, sizeof(tsi)) );

  // Copy the host object to device
  HIP_ASSERT( hipMemcpyHtoD(tsi_d, &tsi, sizeof(tsi)) );

  // Free device memory
  HIP_ASSERT( hipFree(tsi_d) );

  SUCCEED();
}

} // namespace
