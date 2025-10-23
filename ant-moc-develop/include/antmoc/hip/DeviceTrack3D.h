/// \file   DeviceTrack3D.h
/// \brief  3D Track structure on a GPU.
/// \date   December 19, 2020
/// \author An Wang, USTB (wangan.cs@gmail.com)
//
#ifndef DEVICE_TRACK3D_H_
#define DEVICE_TRACK3D_H_

#include "antmoc/hip/DeviceTrack.h"

namespace antmoc {
namespace hip {

/// \struct  DeviceTrack3D
/// \brief   A 3D Track represents a characteristic line across the geometry.
/// \details A 3D Track has particular starting and ending points on the
///          boundaries of the geometry and an azimuthal and polar angle.
struct DeviceTrack3D : public DeviceTrack {

  /** The polar angle for the Track */
  double _theta;

  /* Indices that are used to locate the track in the various track arrays */
  int _polar_index;
  int _z_index;
  int _lz_index;
  int _cycle_index;
  int _cycle_track_index;
  int _train_index;

  /// \brief Default constructor
  __host__ __device__
  DeviceTrack3D();

  /// \brief Print myself
  __device__ void printString();

};  // class DeviceTrack3D

} // namespace hip
} // namespace antmoc

#endif  // DEVICE_TRACK3D_H_
