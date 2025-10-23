/// \file   DeviceMOCKernel.h
/// \brief  MOCKernel objects on device.
/// \date   January 3, 2020
/// \author An Wang, USTB (wangan.cs@gmail.com)

#ifndef DEVICE_MOCKERNEL_H_
#define DEVICE_MOCKERNEL_H_

#include "antmoc/hip/hip_utils.h"

namespace antmoc {
namespace hip {

// Forward declarations
struct DeviceMaterial;
struct DeviceTrack;
struct DeviceQuadrature;
class DeviceRayTracer;


/// \class   DeviceMOCKernel MOCKernel.h "src/MOCKernel.h"
/// \brief   An MOCKernel object specifies a functionality to apply to MOC
///          segments.
/// \details See MOCKernel.h
class DeviceMOCKernel {

protected:

  ///< Count referring to the segment number.
  int _count;

  ///< Number of energy groups in the current problem
  int _num_groups;

public:

  /// \brief Constructor for the DeviceMOCKernel assigns default values.
  /// \param ray_tracer The DeviceRayTracer used to pull relevant tracking
  ///        data from.
  __device__
  DeviceMOCKernel(DeviceRayTracer *ray_tracer);

  /// \brief Destructor
  __device__
  virtual ~DeviceMOCKernel() = default;

  /// \brief   Returns current segment count
  /// \details MOC kernels count how many times they are accessed. This value
  ///          returns the value of the counter (number of execute accesses)
  ///          since kernel creation or last call to newTrack.
  /// \return  _count the counter value
  __device__
  int getCount();

  /// \brief Prepares DeviceMOCKernel for handling a new track
  __device__
  virtual void newTrack(DeviceTrack *track);

  /// \brief Executing function describes kernel behavior
  __device__
  virtual void execute(FP_PRECISION length, DeviceMaterial* mat, long fsr_id,
                       int track_idx, int cmfd_surface_fwd, int cmfd_surface_bwd,
                       FP_PRECISION x_start, FP_PRECISION y_start, FP_PRECISION z_start,
                       FP_PRECISION phi, FP_PRECISION theta) = 0;

};


/// \class   DeviceVolumeKernel
/// \brief   Calculates the volume in FSRs by adding weighted segment lengths.
/// \details See class VolumeKernel for more details.
class DeviceVolumeKernel: public DeviceMOCKernel {

private:

  ///< Pointer to array of FSR volumes.
  FP_PRECISION *_FSR_volumes;

  ///< The associated quadrature from which weights are derived.
  DeviceQuadrature *_quadrature;

  ///< The Track's volume weight.
  FP_PRECISION _weight;

public:

  /// \brief Constructor for the DeviceVolumeKernel assigns default values,
  ///        calls the MOCKernel constructor, and pulls references FSR volumes
  ///        from the provided DeviceRayTracer.
  /// \param ray_tracer The DeviceRayTracer used to pull relevant tracking
  ///        data from.
  __device__
  DeviceVolumeKernel(DeviceRayTracer *ray_tracer);

  /// \brief   Prepares a DeviceVolumeKernel for a new Track.
  /// \details Resets the segment count and updates the weight for the new Track.
  /// \param   track The new Track the DeviceMOCKernel prepares to handle.
  __device__
  void newTrack(DeviceTrack *track) final;

  /// \brief   Adds segment contribution to the FSR volume.
  /// \details The DeviceVolumeKernel execute function adds the product of the
  ///          track length and track weight to the buffer array at index id,
  ///          referring to the array of FSR volumes.
  /// \param   length Segment length.
  /// \param   fsr_id The FSR ID of the FSR associated with the segment.
  __device__
  void execute(FP_PRECISION length, DeviceMaterial* mat, long fsr_id,
               int track_idx, int cmfd_surface_fwd, int cmfd_surface_bwd,
               FP_PRECISION x_start, FP_PRECISION y_start, FP_PRECISION z_start,
               FP_PRECISION phi, FP_PRECISION theta) final;
};


} // namespace hip
} // namespace antmoc

#endif /* DEVICE_MOCKERNEL_H_ */
