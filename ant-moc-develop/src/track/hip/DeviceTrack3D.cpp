#include "antmoc/hip/DeviceTrack3D.h"


namespace antmoc {
namespace hip {

__host__ __device__
DeviceTrack3D::DeviceTrack3D():

  DeviceTrack(),

  _theta(0.),
  _polar_index(-1),
  _z_index(-1),
  _lz_index(-1),
  _cycle_index(-1),
  _cycle_track_index(-1),
  _train_index(-1)
  { }


__device__ void DeviceTrack3D::printString() {

  printf("Track3D:\n"
         "theta = %f\n"
         "polar_index = %d, z_index = %d, lz_index = %d\n",
          _theta, _polar_index, _z_index, _lz_index);
}


} // namespace hip
} // namespace antmoc
