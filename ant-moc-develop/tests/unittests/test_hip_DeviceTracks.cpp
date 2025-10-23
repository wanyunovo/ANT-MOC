/// \file   test_hip_DeviceTracks.cpp
/// \brief  Test device segment and tracks.
/// \date   January 1, 2020
/// \author An Wang, USTB (wangan@xs.ustb.edu.cn)

#include "testing/hip_testharness.h"
#include "antmoc/hip/DeviceTrack.h"
#include "antmoc/hip/DeviceTrack3D.h"
#include "antmoc/math_utils.h"
#include "antmoc/Material.h"
#include "antmoc/Track.h"

#include <map>
#include <vector>

using namespace antmoc;

namespace {

/// Testing fixture
class test_hip_DeviceTracks : public HipTestHarness {

public:

  void SetUp() override {
    // Create materials with automatic ids
    materials.push_back(new Material(0, "101"));
    materials.push_back(new Material(0, "102"));
    materials.push_back(new Material(0, "103"));

    for (size_t i = 0; i < materials.size(); ++i) {
      auto id = materials[i]->getId();
      material_IDs_to_indices[id] = i;
    }
  }

  void TearDown() override {
    for (auto &m : materials)
      delete m;
    materials.clear();
    material_IDs_to_indices.clear();
  }

protected:

  std::vector<Material*> materials;
  std::map<int, int> material_IDs_to_indices;
};


TEST_F(test_hip_DeviceTracks, trackCopyHtoD) {
  // Create a host Track
  Track *track_h = new Track;
  track_h->setValues(3., 0., 0., 0., M_PI);
  track_h->setUid(1);

  // Add segments to the Track
  for (size_t i = 0; i < materials.size(); ++i) {
    segment s;
    s._length = 1.0;
    s._material = materials[i];
    s._region_id = 3 - i;
    s._track_idx = 1;
    s._starting_position[0] = 3. - i;
    s._starting_position[1] = 0;

    track_h->addSegment(s);
    track_h->setNumSegments(i+1);
  }

  // Allocate memory for the device Track
  hip::DeviceTrack *track_d;
  HIP_ASSERT( hipMalloc(&track_d, sizeof(hip::DeviceTrack)) );

  // Copy host Track to device
  hip::trackCopyHtoD(track_d, track_h, material_IDs_to_indices);

  // Print device Track
  hipLaunchKernelGGL(
      hip::printObject, dim3(1), dim3(1), 0, 0,
      track_d);

  HIP_ASSERT( hipFree(track_d) );
  delete track_h;

  SUCCEED();
}


TEST_F(test_hip_DeviceTracks, createDeviceTrack3D) {

  // Allocate memory for the device Track3D
  hip::DeviceTrack3D *track3D_d;
  HIP_ASSERT( hipMalloc(&track3D_d, sizeof(hip::DeviceTrack3D)) );

  // Print device Track3D
  hipLaunchKernelGGL(
      hip::printObject, dim3(1), dim3(1), 0, 0,
      track3D_d);

  HIP_ASSERT( hipFree(track3D_d) );

  SUCCEED();
}


} // namespace
