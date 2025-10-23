#include "antmoc/hip/hip_info.h"
#include "antmoc/hip/hip_utils.h"
#include "antmoc/log.h"
#include "antmoc/mpi_utils.h"

#include "hip/hip_runtime.h"

namespace antmoc {

namespace hip {

void printDeviceProperties() {

  auto gpu_id = getGPUId();

  hipDeviceProp_t prop;
  HIP_ASSERT( hipGetDeviceProperties(&prop, gpu_id) );

  log::verbose_once(
      "Properties of device {}:\n"
      "  name               = {}\n"
      "  GCN architecture   = {}\n"
      "  CUs                = {}\n"
      "  compute capability = {}.{}\n"
      "  warp size          = {}\n"
      "  global mem         = {:.2f} MiB\n"
      "  shared mem/block   = {:.2f} KiB\n"
      "  registers/block    = {}\n"
      "  concurrent kernels = {}\n"
      "  is multi-gpu board = {}\n"
      , gpu_id
      , prop.name
      , prop.gcnArch
      , prop.multiProcessorCount
      , prop.major, prop.minor
      , prop.warpSize
      , prop.totalGlobalMem / double(1 << 20)
      , prop.sharedMemPerBlock / double(1 << 10)
      , prop.regsPerBlock
      , prop.concurrentKernels
      , prop.isMultiGpuBoard
  );
}


bool containsGPU() {
  return getNumGPUs() > 0;
}


int getNumCUs() {
  hipDeviceProp_t dev_prop;
  HIP_ASSERT( hipGetDeviceProperties(&dev_prop, getGPUId()) );

  return dev_prop.multiProcessorCount;
}


}   // namespace

}   // namespace antmoc
