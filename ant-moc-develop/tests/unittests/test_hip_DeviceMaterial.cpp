/// \file   test_hip_DeviceMaterial.cpp
/// \brief  Test DeviceMaterial
/// \date   December 31, 2020
/// \author An Wang, USTB (wangan@xs.ustb.edu.cn)

#include "testing/hip_testharness.h"
#include "antmoc/hip/DeviceMaterial.h"
#include "antmoc/Material.h"

using namespace antmoc;

namespace {

/// Testing fixture
class test_hip_DeviceMaterial : public HipTestHarness {

public:

  void SetUp() override {
    constexpr int num_groups = 4;
    constexpr int num_groups_squared = num_groups * num_groups;
    FP_PRECISION sigma[num_groups] = {1e-1, 1e-1, 1e-1, 1e-1};
    FP_PRECISION sigma_matrix[num_groups_squared] = {
      1e-1, 1e-1, 1e-1, 1e-1,
      1e-2, 1e-2, 1e-2, 1e-2,
      1e-3, 1e-3, 1e-3, 1e-3,
      1e-4, 1e-4, 1e-4, 1e-4,
    };

    // Create a host material
    material_h = new Material(0, "host material");
    material_h->setNumEnergyGroups(num_groups);
    material_h->setSigmaT(sigma, num_groups);
    material_h->setSigmaF(sigma, num_groups);
    material_h->setNuSigmaF(sigma, num_groups);
    material_h->setChi(sigma, num_groups);
    material_h->setSigmaS(sigma_matrix, num_groups_squared);
    material_h->buildFissionMatrix();
  }

  void TearDown() override {
    delete material_h;
    material_h = nullptr;
  }

  Material *material_h = nullptr;
};


TEST_F(test_hip_DeviceMaterial, copyHtoD) {
  // Allocate memory for the device material
  hip::DeviceMaterial *material_d;
  HIP_ASSERT( hipMalloc(&material_d, sizeof(hip::DeviceMaterial)) );

  // Copy host material to device
  hip::materialCopyHtoD(material_d, material_h);

  // Print device material
  hipLaunchKernelGGL(
      hip::printObject, dim3(1), dim3(1), 0, 0,
      material_d);

  HIP_ASSERT( hipFree(material_d) );

  SUCCEED();
}

} // namespace
