/// \file   test_hip_DeviceQuadrature.cpp
/// \brief  Test DeviceQuadrature
/// \date   January 1, 2020
/// \author An Wang, USTB (wangan@xs.ustb.edu.cn)

#include "testing/hip_testharness.h"
#include "antmoc/hip/DeviceQuadrature.h"
#include "antmoc/Quadrature.h"

using namespace antmoc;

namespace {

/// Testing fixture
class test_hip_DeviceQuadrature : public HipTestHarness {

public:

  void SetUp() override {

    constexpr size_t num_azims  = 64;
    constexpr size_t num_polars = 12;
    constexpr double unified_xy_spacing = 0.005;
    constexpr double unified_z_spacing  = 0.075;

    // Create a host quadrature
    quadrature_h = new EqualWeightPolarQuad();
    quadrature_h->setNumAzimAngles(num_azims);
    quadrature_h->setNumPolarAngles(num_polars);

    for (size_t a = 0; a < num_azims/4; ++a) {
      quadrature_h->setAzimSpacing(unified_xy_spacing, a);

      for (size_t p = 0; p < num_polars/2; ++p)
        quadrature_h->setPolarSpacing(unified_z_spacing, a, p);
    }

    // Initialize angle values
    quadrature_h->initialize();

    // Compute weights
    quadrature_h->precomputeWeights(true);
  }

  void TearDown() override {
    delete quadrature_h;
  }

  Quadrature *quadrature_h = nullptr;
};


TEST_F(test_hip_DeviceQuadrature, copyHtoD) {
  // Allocate memory for the device quadrature
  hip::DeviceQuadrature *quadrature_d;
  HIP_ASSERT( hipMalloc(&quadrature_d, sizeof(hip::DeviceQuadrature)) );

  // Copy host quadrature to device
  hip::quadratureCopyHtoD(quadrature_d, quadrature_h);

  // Print device quadrature
  //hipLaunchKernelGGL(
  //    hip::printDeviceQuadrature, dim3(1), dim3(1), 0, 0,
  //    quadrature_d);

  HIP_ASSERT( hipFree(quadrature_d) );

  SUCCEED();
}

} // namespace
