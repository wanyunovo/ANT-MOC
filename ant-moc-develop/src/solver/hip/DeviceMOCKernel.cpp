#include "antmoc/hip/DeviceMOCKernel.h"
#include "antmoc/hip/DeviceMaterial.h"
#include "antmoc/hip/DeviceQuadrature.h"
#include "antmoc/hip/DeviceRayTracer.h"
#include "antmoc/hip/DeviceTrack3D.h"

namespace antmoc {
namespace hip {


__device__
DeviceMOCKernel::DeviceMOCKernel(DeviceRayTracer *ray_tracer) :
  _count(0),
  _num_groups(ray_tracer->getNumEnergyGroups())
{}


__device__
int DeviceMOCKernel::getCount() {
  return _count;
}


//------------------------------------------------------------------------------
// DeviceVolumeKernel
//------------------------------------------------------------------------------
__device__
DeviceVolumeKernel::DeviceVolumeKernel(DeviceRayTracer *ray_tracer) :
  DeviceMOCKernel(ray_tracer) {

  _FSR_volumes = ray_tracer->getFSRVolumesBuffer();
  _quadrature = ray_tracer->getQuadrature();
  _weight = 0;
}


__device__
void DeviceVolumeKernel::newTrack(DeviceTrack *track) {

  /* Compute the Track azimuthal weight */
  auto azim_index = track->_azim_index;
  _weight = _quadrature->getAzimSpacing(azim_index)
            * _quadrature->getAzimWeight(azim_index);

  /* Multiply by polar weight if 3D */
  auto track_3D = dynamic_cast<DeviceTrack3D*>(track);
  if (track_3D) {
    auto polar_index = track_3D->_polar_index;
    _weight *= _quadrature->getPolarSpacing(azim_index, polar_index)
               * _quadrature->getPolarWeight(azim_index, polar_index);
  }
}


__device__
void DeviceVolumeKernel::execute(
        FP_PRECISION length, DeviceMaterial* mat, long fsr_id, int track_idx,
        int cmfd_surface_fwd, int cmfd_surface_bwd,
        FP_PRECISION x_start, FP_PRECISION y_start, FP_PRECISION z_start,
        FP_PRECISION phi, FP_PRECISION theta) {

  // Add value to buffer
  atomicAdd(_FSR_volumes + fsr_id, _weight * length);

  // Increment the count of the number of times the kernel is executed.
  _count++;
}


} // namespace hip
} // namespace antmoc
