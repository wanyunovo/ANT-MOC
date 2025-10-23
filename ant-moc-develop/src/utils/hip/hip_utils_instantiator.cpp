/// \file    hip_utils_instantiator.cpp
/// \details This is a temporary workaround for some templates. To make them
///          available in tests without any link error, we instantiate these
///          templates here.

#include "antmoc/hip/hip_utils.h"
#include "antmoc/hip/DeviceFSRDataTypes.h"
#include "antmoc/hip/DeviceMaterial.h"
#include "antmoc/hip/DevicePoint.h"
#include "antmoc/hip/DeviceQuadrature.h"
#include "antmoc/hip/DeviceTrack.h"
#include "antmoc/hip/DeviceTrack3D.h"

namespace antmoc {
namespace hip {


// Implementation
template <typename T>
__global__
void printObject(T *ptr) {
  if (ptr) {
    ptr->printString();
  }
  else {
    printf("Could not print empty pointer");
    abort();
  }
};

// Instantiation
template __global__ void printObject(DeviceFSRData *);
template __global__ void printObject(DeviceExtrudedFSR *);
template __global__ void printObject(DeviceMaterial *);
template __global__ void printObject(DevicePoint *);
template __global__ void printObject(DeviceTrack *);
template __global__ void printObject(DeviceTrack3D *);


} // namespace hip
} // namespace antmoc
