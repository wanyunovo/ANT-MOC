#include "antmoc/hip/DeviceQuadrature.h"
#include "antmoc/hip/hip_utils.h"
#include "antmoc/Quadrature.h"

namespace antmoc {
namespace hip {

__host__ __device__
DeviceQuadrature::DeviceQuadrature() :
  _num_azim(0),
  _num_polar(0),
  _sin_thetas(nullptr),
  _thetas(nullptr),
  _phis(nullptr),
  _azim_spacings(nullptr),
  _azim_weights(nullptr),
  _polar_spacings(nullptr),
  _polar_weights(nullptr),
  _total_weights(nullptr)
{}


__host__ __device__
DeviceQuadrature::~DeviceQuadrature() {
  if (_sin_thetas)
    free(_sin_thetas);

  if (_thetas)
    free(_thetas);

  if (_phis)
    free(_phis);

  if (_azim_spacings)
    free(_azim_spacings);

  if (_azim_weights)
    free(_azim_weights);

  if (_polar_spacings)
    free(_polar_spacings);

  if (_polar_weights)
    free(_polar_weights);

  if (_total_weights)
    free(_total_weights);
}


__device__
size_t DeviceQuadrature::getNumPolarAngles() const {
  return _num_polar;
}


__device__
size_t DeviceQuadrature::getNumAzimAngles() const {
  return _num_azim;
}


__device__
double DeviceQuadrature::getSinTheta(size_t azim, size_t polar) const {

  if (azim >= _num_azim/2)
    azim = _num_azim - azim - 1;

  return _sin_thetas[_num_azim/2 * azim + polar];
}


__device__
double DeviceQuadrature::getTheta(size_t azim, size_t polar) const {

  if (azim >= _num_azim/2)
    azim = _num_azim - azim - 1;

  return _thetas[_num_azim/2 * azim + polar];
}


__device__
double DeviceQuadrature::getPhi(size_t azim) const {

  if (azim >= _num_azim/2)
    azim = _num_azim - azim - 1;

  return _phis[azim];
}


__device__
double DeviceQuadrature::getAzimWeight(size_t azim) const {

  if (azim >= _num_azim/2)
    azim = _num_azim - azim - 1;

  return _azim_weights[azim];
}


__device__
double DeviceQuadrature::getPolarWeight(size_t azim, size_t polar) const {

  if (azim >= _num_azim/2)
    azim = _num_azim - azim - 1;

  return _polar_weights[_num_azim/2 * azim + polar];
}


__device__
double DeviceQuadrature::getWeight(size_t azim, size_t polar) const {

  if (azim >= _num_azim/2)
    azim = _num_azim - azim - 1;

  return _total_weights[_num_azim/2 * azim + polar];
}


__device__
const double* DeviceQuadrature::getSinThetas() const {

  return _sin_thetas;
}


__device__
const double* DeviceQuadrature::getThetas() const {

  return _thetas;
}


__device__
const double* DeviceQuadrature::getPhis() const {

  return _phis;
}


__device__
const double* DeviceQuadrature::getAzimWeights() const {

  return _azim_weights;
}


__device__
const double* DeviceQuadrature::getPolarWeights() const {

  return _polar_weights;
}


__device__
const double* DeviceQuadrature::getAzimSpacings() const {
  return _azim_spacings;
}


__device__
double DeviceQuadrature::getAzimSpacing(size_t azim) const {
  if (azim >= _num_azim/2)
    azim = _num_azim - azim - 1;
  return _azim_spacings[azim];
}


__device__
const double* DeviceQuadrature::getPolarSpacings() const {
  return _polar_spacings;
}


__device__
double DeviceQuadrature::getPolarSpacing(size_t azim, size_t polar) const {
  if (azim >= _num_azim/2)
    azim = _num_azim - azim - 1;
  return _polar_spacings[_num_azim/2 * azim + polar];
}


void quadratureCopyHtoD(DeviceQuadrature *quad_d, Quadrature *quad_h) {
  size_t num_azims    = quad_h->getNumAzimAngles();
  size_t num_polars   = quad_h->getNumPolarAngles();

  // Copy data from host to device
  HIP_ASSERT( hipMemcpyHtoD(&quad_d->_num_azim, &num_azims, sizeof(size_t)) );
  HIP_ASSERT( hipMemcpyHtoD(&quad_d->_num_polar, &num_polars, sizeof(size_t)) );
  vectorCopyHtoD(quad_d->_phis, quad_h->getPhis());
  vectorCopyHtoD(quad_d->_azim_spacings,   quad_h->getAzimSpacings());
  vectorCopyHtoD(quad_d->_azim_weights,    quad_h->getAzimWeights());
  matrixCopyHtoD(quad_d->_thetas,          quad_h->getThetas());
  matrixCopyHtoD(quad_d->_sin_thetas,      quad_h->getSinThetas());
  matrixCopyHtoD(quad_d->_polar_spacings,  quad_h->getPolarSpacings());
  matrixCopyHtoD(quad_d->_polar_weights,   quad_h->getPolarWeights());
  matrixCopyHtoD(quad_d->_total_weights,   quad_h->getTotalWeights());
}

} // namespace hip
} // namespace antmoc
