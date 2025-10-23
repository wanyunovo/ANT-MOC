#include "antmoc/hip/DeviceRayTracer.h"
#include "antmoc/hip/DeviceFSRDataTypes.h"
#include "antmoc/hip/DeviceMaterial.h"
#include "antmoc/hip/DevicePoint.h"
#include "antmoc/hip/DeviceQuadrature.h"
#include "antmoc/hip/DeviceTrack3D.h"

namespace antmoc {
namespace hip {

__host__ __device__
DeviceRayTracer::DeviceRayTracer() :
  _num_groups(0),
  _num_2D_tracks(0),
  _quadrature(nullptr),
  _tracks_2D_array(nullptr),
  _FSR_volumes(nullptr),
  _fsr_centroids(nullptr),
  _extruded_fsrs(nullptr),
  _global_z_mesh(nullptr)
{}


__host__ __device__
DeviceRayTracer::~DeviceRayTracer() {
  if (_quadrature)
    free(_quadrature);

  if (_tracks_2D_array)
    free(_tracks_2D_array);

  if (_FSR_volumes)
    free(_FSR_volumes);

  if (_fsr_centroids)
    free(_fsr_centroids);

  if (_extruded_fsrs)
    free(_extruded_fsrs);

  if (_global_z_mesh)
    free(_global_z_mesh);
}


__device__
int DeviceRayTracer::getNumEnergyGroups() const {
  return _num_groups;
}


__device__
long DeviceRayTracer::getNum2DTracks() const {
  return _num_2D_tracks;
}


__device__
DeviceTrack* DeviceRayTracer::get2DTrack(long uid) {
  return _tracks_2D_array + uid;
}


__device__
DeviceQuadrature* DeviceRayTracer::getQuadrature() {
  return _quadrature;
}


__device__
FP_PRECISION* DeviceRayTracer::getFSRVolumesBuffer() {
  return _FSR_volumes;
}

} // namespace hip
} // namespace antmoc
