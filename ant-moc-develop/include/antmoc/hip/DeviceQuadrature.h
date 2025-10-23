/// \file DeviceQuadrature.h
/// \brief Quadrature on device. See Quadrature.h.
/// \date December 31, 2020
/// \author An Wang, USTB (wangan.cs@gmail.com)

#ifndef DEVICE_QUADRATURE_H_
#define DEVICE_QUADRATURE_H_

#include "hip/hip_runtime.h"


namespace antmoc {

// Forward declaration
class Quadrature;

namespace hip {

/// \struct DeviceQuadrature
/// \brief  Quadrature data.
struct DeviceQuadrature {

  ///< The number of azimuthal angles in (0, 2*PI)
  size_t _num_azim;

  ///< The number of polar angles in (0, PI)
  size_t _num_polar;

  ///< An array of the sines of quadrature polar angles indexed by (azim, polar)
  double* _sin_thetas;

  ///< An array of the quadrature polar angles indexed by (azim, polar)
  double* _thetas;

  ///< An array of the quadrature azimuthal angles
  double* _phis;

  ///< The actual track azimuthal spacing (cm) by azimuthal angle
  double* _azim_spacings;

  ///< An array of the quadrature azimuthal weights
  double* _azim_weights;

  ///< The actual track polar spacing (cm) by (azim, polar)
  double* _polar_spacings;

  ///< An array of the quadrature polar weights indexed by (azim, polar)
  double* _polar_weights;

  ///< An array of the total weights indexed by (azim, polar)
  double* _total_weights;

  __host__ __device__
  DeviceQuadrature();

  __host__ __device__
  ~DeviceQuadrature();

  __device__ size_t getNumPolarAngles() const;
  __device__ size_t getNumAzimAngles() const;
  __device__ double getSinTheta(size_t azim, size_t polar) const;
  __device__ double getTheta(size_t azim, size_t polar) const;
  __device__ double getPhi(size_t azim) const;
  __device__ double getAzimWeight(size_t azim) const;
  __device__ double getPolarWeight(size_t azim, size_t polar) const;
  __device__ double getWeight(size_t azim, size_t polar) const;
  __device__ const double* getSinThetas() const;
  __device__ const double* getThetas() const;
  __device__ const double* getPhis() const;
  __device__ const double* getAzimWeights() const;
  __device__ const double* getPolarWeights() const;
  __device__ const double* getAzimSpacings() const;
  __device__ double getAzimSpacing(size_t azim) const;
  __device__ const double* getPolarSpacings() const;
  __device__ double getPolarSpacing(size_t azim, size_t polar) const;
};


/// \brief Copy quadrature data from host to device.
/// \param quad_d Pointer to device quadrature.
/// \param quad_h Pointer to host quadrature.
void quadratureCopyHtoD(DeviceQuadrature *quad_d, Quadrature *quad_h);

} // namespace hip
} // namespace antmoc

#endif /* DEVICE_QUADRATURE_H_ */
