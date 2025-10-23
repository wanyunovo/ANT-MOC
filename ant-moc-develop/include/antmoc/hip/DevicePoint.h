/// \file   DevicePoint.h
/// \brief  Point on device.
/// \date   December 19, 2020
/// \author An Wang, USTB (wangan.cs@gmail.com)

#ifndef DEVICE_POINT_H_
#define DEVICE_POINT_H_

#include "antmoc/hip/hip_utils.h"

namespace antmoc {

// Forward declaration
class Point;

namespace hip {

/// \struct DevicePoint
/// \brief A 3-D point on GPU.
struct DevicePoint {

  FP_PRECISION _x;
  FP_PRECISION _y;
  FP_PRECISION _z;

  __host__ __device__
  DevicePoint(FP_PRECISION x = .0, FP_PRECISION y = .0, FP_PRECISION z = .0):
    _x(x), _y(y), _z(z) { }

  /// \brief Compute the distance between two points.
  __device__
  FP_PRECISION distance(const FP_PRECISION x, const FP_PRECISION y, const FP_PRECISION z) const;

  /// \brief Compute the distance between two points.
  __device__
  FP_PRECISION distanceToPoint(const DevicePoint *point);

  /// \brief Print myself.
  __device__ void printString();

};  // class DevicePoint


/// \brief Copy a Point object from host to device.
/// \param point_d Pointer to a DevicePoint object.
/// \param point_h Pointer to a Point object.
void pointCopyHtoD(DevicePoint *point_d, Point *point_h);

} // namespace hip
} // namespace antmoc

#endif  // DEVICE_POINT_H_
