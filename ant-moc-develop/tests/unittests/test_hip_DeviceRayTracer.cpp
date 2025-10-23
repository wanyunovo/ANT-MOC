/// \file   test_hip_DeviceRayTracer.cpp
/// \brief  Test DeviceRayTracer
/// \date   January 3, 2020
/// \author An Wang, USTB (wangan@xs.ustb.edu.cn)

#include "testing/hip_testharness.h"
#include "antmoc/hip/DeviceRayTracer.h"

using namespace antmoc;

namespace {

/// Testing fixture
class test_hip_DeviceRayTracer : public HipTestHarness {};


TEST_F(test_hip_DeviceRayTracer, allocation) {
  // Allocate memory for the device quadrature
  hip::DeviceRayTracer *ray_tracer;
  HIP_ASSERT( hipMalloc(&ray_tracer, sizeof(hip::DeviceRayTracer)) );

  // Copy memory from host to device.
  constexpr size_t num_fsrs = 4;
  FP_PRECISION fsr_volumes[num_fsrs] = {1.1, 1.2, 1.3, 1.4};
  hip::arrayCopyHtoD(ray_tracer->_FSR_volumes, fsr_volumes, num_fsrs);

  // Free device memory
  HIP_ASSERT( hipFree(ray_tracer) );

  SUCCEED();
}

} // namespace
