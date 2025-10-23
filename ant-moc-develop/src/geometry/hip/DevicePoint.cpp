#include "antmoc/hip/DevicePoint.h"
#include "antmoc/Point.h"

namespace antmoc{
namespace hip{

__device__
FP_PRECISION DevicePoint::
distance(const FP_PRECISION x, const FP_PRECISION y, const FP_PRECISION z) const {
  auto deltax = _x - x;
  auto deltay = _y - y;
  auto deltaz = _z - z;
  // sqrt() accepts only type double, so a promotion may occur.
  return sqrt(deltax*deltax + deltay*deltay + deltaz*deltaz);
}


__device__
FP_PRECISION DevicePoint::
distanceToPoint(const DevicePoint *point) {
  return distance(point->_x, point->_y, point->_z);
}


__device__
void DevicePoint::printString() {
  printf("Point: x = %f, y = %f, z = %f\n", _x, _y, _z);
}


//------------------------------------------------------------------------------
// Helper functions
//------------------------------------------------------------------------------
void pointCopyHtoD(DevicePoint *point_d, Point *point_h) {
  FP_PRECISION x = point_h->getX();
  FP_PRECISION y = point_h->getY();
  FP_PRECISION z = point_h->getZ();

  HIP_ASSERT( hipMemcpyHtoD(&point_d->_x, &x, sizeof(FP_PRECISION)) );
  HIP_ASSERT( hipMemcpyHtoD(&point_d->_y, &y, sizeof(FP_PRECISION)) );
  HIP_ASSERT( hipMemcpyHtoD(&point_d->_z, &z, sizeof(FP_PRECISION)) );
}

} // namespace hip
} // namespace antmoc
