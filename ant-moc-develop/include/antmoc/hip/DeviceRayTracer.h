/// \file DeviceRayTracer.h
/// \brief The DeviceRayTracer class.
/// \date December 24, 2020
/// \author An Wang, USTB (wangan.cs@gmail.com)

#ifndef DEVICE_OTF_TRACER_H_
#define DEVICE_OTF_TRACER_H_

#include "antmoc/hip/hip_utils.h"

namespace antmoc {
namespace hip {

// Forward declarations
struct DevicePoint;
struct DeviceTrack;
struct DeviceExtrudedFSR;
struct DeviceQuadrature;

/// \class   DeviceRayTracer
/// \brief   Data and helper methods for ray tracing on device.
/// \details This class is initialized by GPUTrackGenerator.
class DeviceRayTracer {

// Temporarily make this block public.
public:

  ///< Number of energy groups.
  int _num_groups;

  ///< Number of 2D tracks for all azimuthal angles.
  long _num_2D_tracks;

  ///< Quadrature for MOC.
  DeviceQuadrature *_quadrature;

  ///< Array of 2D device tracks.
  DeviceTrack *_tracks_2D_array;

  ///< A buffer holding the computed FSR volumes.
  FP_PRECISION  *_FSR_volumes;

  DevicePoint *_fsr_centroids;

  DeviceExtrudedFSR *_extruded_fsrs;

  double *_global_z_mesh;


public:

  /// \brief Default constructor.
  __host__ __device__
  DeviceRayTracer();

  /// \brief Free either host or device memory.
  __host__ __device__
  virtual ~DeviceRayTracer();

  /// \brief Returns the number of energy groups.
  __device__
  int getNumEnergyGroups() const;

  /// \brief Returns the number of 2D trakcs.
  __device__
  long getNum2DTracks() const;

  /// \brief Returns the pointer to a 2D track.
  __device__
  DeviceTrack* get2DTrack(long uid);

  /// \brief Returns an initialized DeviceQuadrature.
  __device__
  DeviceQuadrature* getQuadrature();

  /// \brief Returns the FSR volumes array.
  __device__
  FP_PRECISION* getFSRVolumesBuffer();

  __device__
  bool containsGlobalZMesh() const;

  __device__
  bool containsFSRCentroids() const;

};  // class DeviceRayTracer

} // namespace hip
} // namespace antmoc

#endif  // DEVICE_OTF_TRACER_H_
