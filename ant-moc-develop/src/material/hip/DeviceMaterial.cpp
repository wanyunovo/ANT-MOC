#include "antmoc/hip/DeviceMaterial.h"
#include "antmoc/Material.h"

namespace antmoc {
namespace hip {


__host__ __device__
DeviceMaterial::DeviceMaterial():
  _id(-1),
  _num_groups(0),
  _sigma_t(nullptr),
  _sigma_s(nullptr),
  _sigma_a(nullptr),
  _sigma_f(nullptr),
  _nu_sigma_f(nullptr),
  _chi(nullptr),
  _fiss_matrix(nullptr)
  { }


__host__ __device__
DeviceMaterial::~DeviceMaterial() {
  if (_sigma_t)
    free(_sigma_t);

  if (_sigma_s)
    free(_sigma_s);

  if (_sigma_a)
    free(_sigma_a);

  if (_sigma_f)
    free(_sigma_f);

  if (_nu_sigma_f)
    free(_nu_sigma_f);

  if (_chi)
    free(_chi);

  if (_fiss_matrix)
    free(_fiss_matrix);
}


__device__
FP_PRECISION DeviceMaterial::getSigmaTByGroup(int group) {
  return _sigma_t[group-1];
}


__device__
FP_PRECISION DeviceMaterial::getSigmaSByGroup(int origin, int destination) {
  return _sigma_s[(destination-1)*_num_groups + (origin-1)];
}


__device__
FP_PRECISION DeviceMaterial::getSigmaAByGroup(int group) {

  /* If not initialized, compute _sigma_a the absorption cross section */
  if (_sigma_a == nullptr) {
    _sigma_a = new FP_PRECISION[_num_groups];
    for (int g=0; g < _num_groups; g++) {
      _sigma_a[g] = _sigma_t[g];
      for (int gp=0; gp < _num_groups; gp++)
        _sigma_a[g] -= _sigma_s[gp*_num_groups + g];
    }
  }

  return _sigma_a[group-1];
}


__device__
FP_PRECISION DeviceMaterial::getSigmaFByGroup(int group) {
  return _sigma_f[group-1];
}


__device__
FP_PRECISION DeviceMaterial::getNuSigmaFByGroup(int group) {
  return _nu_sigma_f[group-1];
}


__device__
FP_PRECISION DeviceMaterial::getChiByGroup(int group) {
  return _chi[group-1];
}


__device__
FP_PRECISION DeviceMaterial::getFissionMatrixByGroup(int origin, int destination) {
  return _fiss_matrix[(destination-1)*_num_groups + (origin-1)];
}


__device__ void DeviceMaterial::printString() {

  printf("Material id = %d, num_groups = %d\n", _id , _num_groups);

  printf("Sigma_t = \n");
  for (int e = 0; e < _num_groups; e++)
    printf("%f, ", _sigma_t[e]);

  printf("\nSigma_f = \n");
  for (int e = 0; e < _num_groups; e++)
    printf("%f, ", _sigma_f[e]);

  printf("\nnu_sigma_f = \n");
  for (int e = 0; e < _num_groups; e++)
    printf("%f, ", _nu_sigma_f[e]);

  printf("\nSigma_s = \n");
  for (int G = 0; G < _num_groups; G++) {
    for (int g = 0; g < _num_groups; g++)
      printf("%f\t\t", _sigma_s[G+g*_num_groups]);
    printf("\n");
  }

  printf("\nChi = \n");
  for (int e = 0; e < _num_groups; e++)
    printf("%f, ", _chi[e]);

  printf("\nFiss. Matrix = \n");
  for (int G = 0; G < _num_groups; G++) {
    for (int g = 0; g < _num_groups; g++)
      printf("%f\t\t", _fiss_matrix[G+g*_num_groups]);
    printf("\n");
  }
}


//------------------------------------------------------------------------------
// Helper functions
//------------------------------------------------------------------------------
void materialCopyHtoD(DeviceMaterial *material_d, Material *material_h) {
  // Copy simple values from host to device.
  auto id          = material_h->getId();
  auto max_sigma_t = material_h->getMaxSigmaT();
  auto num_groups  = material_h->getNumEnergyGroups();

  HIP_ASSERT( hipMemcpyHtoD(&material_d->_id, &id, sizeof(id)) );
  HIP_ASSERT( hipMemcpyHtoD(&material_d->_max_sigma_t, &max_sigma_t, sizeof(max_sigma_t)) );
  HIP_ASSERT( hipMemcpyHtoD(&material_d->_num_groups,  &num_groups,  sizeof(num_groups)) );

  // Data array sizes.
  auto vec_len = static_cast<size_t>(num_groups);
  auto mat_len = static_cast<size_t>(num_groups * num_groups);

  // Copy Material data from host to arrays on the device
  arrayCopyHtoD(material_d->_sigma_t,     material_h->getSigmaT(),    vec_len);
  arrayCopyHtoD(material_d->_sigma_s,     material_h->getSigmaS(),    mat_len);
  arrayCopyHtoD(material_d->_sigma_a,     material_h->getSigmaA(),    vec_len);
  arrayCopyHtoD(material_d->_sigma_f,     material_h->getSigmaF(),    vec_len);
  arrayCopyHtoD(material_d->_nu_sigma_f,  material_h->getNuSigmaF(),  vec_len);
  arrayCopyHtoD(material_d->_chi,         material_h->getChi(),       vec_len);
  arrayCopyHtoD(material_d->_fiss_matrix, material_h->getFissionMatrix(), mat_len);
}

void materialCopyDtoH(Material *material_h, Material *material_d) {
}

}   // namespace hip
}   // namespace antmoc
