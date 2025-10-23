/// \file   test_hip_DevicePoint.cpp
/// \brief  Test DevicePoint
/// \date   December 31, 2020
/// \author An Wang, USTB (wangan@xs.ustb.edu.cn)

#include "testing/hip_testharness.h"
#include "antmoc/hip/DevicePoint.h"
#include "antmoc/Point.h"

using namespace antmoc;

namespace {

/// Testing fixture
class test_hip_DevicePoint : public HipTestHarness {

public:

  void SetUp() override {
    point_h = new Point(1.1, 2.2, 3.3);
  }

  void TearDown() override {
    delete point_h;
    point_h = nullptr;
  }

  Point *point_h = nullptr;
};


TEST_F(test_hip_DevicePoint, copyHtoD) {
  // Allocate memory for the device point
  hip::DevicePoint *point_d;
  HIP_ASSERT( hipMalloc(&point_d, sizeof(hip::DevicePoint)) );

  // Copy host point to device
  hip::pointCopyHtoD(point_d, point_h);

  // Print device point
  hipLaunchKernelGGL(
      hip::printObject, dim3(1), dim3(1), 0, 0,
      point_d);

  HIP_ASSERT( hipFree(point_d) );

  SUCCEED();
}

} // namespace
