/// \file   DeviceTrack.h
/// \brief  Track, segment and point structures on a GPU.
/// \date   December 19, 2020
/// \author An Wang, USTB (wangan.cs@gmail.com)
/// \author William Boyd, MIT, Course 22 (wboyd@mit.edu, 2012)

#ifndef DEVICE_TRACK_H_
#define DEVICE_TRACK_H_

#include "antmoc/hip/DevicePoint.h"
#include "antmoc/hip/hip_utils.h"
#include "antmoc/enum_types.h"
#include <map>

namespace antmoc {

// Forward declarations
class Track;

namespace hip {

/// \struct DeviceSegment
/// \brief A DeviceSegment represents a line segment within a single flat source
///        region along a track.
struct DeviceSegment {

  /** The length of the segment (cm) */
  double _length;

  /** An index into the materials array that contains Material pointers */
  int _material_index;

  /** The ID for flat source region in which this segment resides */
  long _region_id;

  /** The ID of the track */
  int _track_idx;

  /** The starting point of the segment relative to the FSR centroid */
  FP_PRECISION _start[3];

  /// \brief Print myself.
  __device__ void printString();
};


/// \struct DeviceTrack
/// \brief A Track represents a characteristic line across the geometry.
/// \details A Track has particular starting and ending points on the
///          boundaries of the geometry and an azimuthal and polar angle.
struct DeviceTrack {

  /** A monotonically increasing unique ID for each Track created */
  long _uid;

  /** The Track's start point */
  DevicePoint _start;

  /** The azimuthal angle for the Track */
  double _phi;

  /** A dynamically sized vector of segments making up this Track */
  DeviceSegment *_segments;

  /** Number of segments recorded during volume calculation */
  int _num_segments;

  /* Indices that are used to locate the track in the various track arrays */
  int _azim_index;
  long _xy_index;
  long _link_index;
  long _chain_index;

  /** Pointers to next, reflective, and periodic Tracks in the forward and
   *  reverse directions */
  long _track_next_fwd;
  long _track_next_bwd;
  long _track_prdc_fwd;
  long _track_prdc_bwd;
  long _track_refl_fwd;
  long _track_refl_bwd;

  /** Booleans to indicate wheter the reflective Tracks in the forward and
   *  and backward direction enter into Tracks pointed in the forward
   *  direction. */
  bool _next_fwd_is_fwd;
  bool _next_bwd_is_fwd;

  /** An enum to indicate whether the outgoing angular flux along this
   *  Track's "forward" direction should be zeroed out for vacuum boundary
   *  conditions or sent to a periodic or reflective track. */
  boundaryType _bc_fwd;

  /** An enum to indicate whether the outgoing angular flux along this
   *  Track's "reverse" direction should be zeroed out for vacuum boundary
   *  conditions or sent to a periodic or reflective track. */
  boundaryType _bc_bwd;

  /// \brief Constructor
  __host__ __device__
  DeviceTrack();

  /// \brief Destructor
  __host__ __device__
  virtual ~DeviceTrack();

  /// \brief Print myself
  __device__ void printString();
};


/// \brief   Given a pointer to a Track on the host, a DeviceTrack on
///          the GPU, and the map of material IDs to indices in the
///          _materials array, copy all of the class attributes and
///          segments from the Track object on the host to the GPU.
/// \details This routine is called by the GPUSolver::initializeTracks()
///          private class method and is not intended to be called
///          directly.
/// \param track_d Pointer to a DeviceTrack on the GPU
/// \param track_h Pointer to a Track on the host
/// \param material_IDs_to_indices Map of material IDs to indices
///        in the _materials array.
void trackCopyHtoD(DeviceTrack *track_d, Track *track_h,
                   std::map<int, int> &material_IDs_to_indices);

} // namespace hip
} // namespace antmoc

#endif /* DEVICE_TRACK_H_ */
