/// \file   test_hip_DeviceFSRDataTypes.cpp
/// \brief  Test device FSR data types.
/// \date   December 31, 2020
/// \author An Wang, USTB (wangan@xs.ustb.edu.cn)

#include "testing/hip_testharness.h"
#include "antmoc/hip/DeviceFSRDataTypes.h"
#include "antmoc/FSRDataTypes.h"
#include "antmoc/Point.h"

using namespace antmoc;

namespace {

/// Testing fixture
class test_hip_DeviceFSRDataTypes : public HipTestHarness {

public:

  void SetUp() override {
    // Initialize a FSRData for use.
    fsr_h = new FSRData();
    auto point    = new Point(1e-2, 1e-3, 1e-4);
    auto centroid = new Point(0.3, 0.3, 0.3);

    fsr_h->_fsr_id = 103;
    // We will set the material index for device FSRData manually.
    //fsr_h->_mat_id = 1;
    fsr_h->_point = point;
    fsr_h->_centroid = centroid;

    // Initialize an ExtrudedFSR for use.
    extfsr_h = new ExtrudedFSR();

    constexpr size_t num_fsrs = 4;
    auto fsr_ids    = new long[num_fsrs];
    auto axial_mesh = new double[num_fsrs];

    for (size_t i = 0; i < num_fsrs; ++i) {
      fsr_ids[i] = 1000 + i;
      axial_mesh[i] = 1e-2 * i;
    }

    extfsr_h->_fsr_id = 11;
    extfsr_h->_num_fsrs = num_fsrs;
    extfsr_h->_fsr_ids = fsr_ids;
    extfsr_h->_mesh = axial_mesh;
    // We will set material indices manually.
    //fsr_h->_materials = ;
  }

  void TearDown() override {
    delete fsr_h;
    delete extfsr_h;
  }

  FSRData *fsr_h = nullptr;
  ExtrudedFSR *extfsr_h = nullptr;
};


TEST_F(test_hip_DeviceFSRDataTypes, fsrDataCopyHtoD) {
  // Allocate memory for the device FSRData
  hip::DeviceFSRData *fsr_d;
  HIP_ASSERT( hipMalloc(&fsr_d, sizeof(hip::DeviceFSRData)) );

  // Copy host FSRData to device
  int material_index = 3;
  hip::fsrDataCopyHtoD(fsr_d, fsr_h, material_index);

  // Print device FSRData
  hipLaunchKernelGGL(
      hip::printObject, dim3(1), dim3(1), 0, 0,
      fsr_d);

  HIP_ASSERT( hipFree(fsr_d) );

  SUCCEED();
}


TEST_F(test_hip_DeviceFSRDataTypes, extrudedFSRCopyHtoD) {
  // Allocate memory for the device ExtrudedFSR
  hip::DeviceExtrudedFSR *extfsr_d;
  HIP_ASSERT( hipMalloc(&extfsr_d, sizeof(hip::DeviceExtrudedFSR)) );

  // Copy host ExtrudedFSR to device
  std::vector<int> material_indices;
  for (int i = 0; i < extfsr_h->_num_fsrs; ++i) {
    material_indices.push_back(2000 + i);
  }
  hip::extrudedFSRCopyHtoD(extfsr_d, extfsr_h, material_indices);

  // Print device ExtrudedFSR
  hipLaunchKernelGGL(
      hip::printObject, dim3(1), dim3(1), 0, 0,
      extfsr_d);

  HIP_ASSERT( hipFree(extfsr_d) );

  SUCCEED();
}

} // namespace
