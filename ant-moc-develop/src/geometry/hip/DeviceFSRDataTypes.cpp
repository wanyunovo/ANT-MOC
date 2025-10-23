#include "antmoc/hip/DeviceFSRDataTypes.h"
#include "antmoc/hip/DevicePoint.h"
#include "antmoc/FSRDataTypes.h"
#include "antmoc/Material.h"


namespace antmoc {
namespace hip {

__host__ __device__
DeviceFSRData::DeviceFSRData() :
  _id(-1),
  _material_index(-1),
  _point(nullptr),
  _centroid(nullptr)
{}


__host__ __device__
DeviceFSRData::~DeviceFSRData() {
  if (_point) {
    free(_point);
  }

  if (_centroid) {
    free(_centroid);
  }
}


__device__
void DeviceFSRData::printString() {

  printf("DeviceFSRData: id = %ld, material_index = %d\n",
          _id, _material_index);
  printf("characteristic point:\n");
  if (_point)
    _point->printString();

  printf("centroid:\n");
  if (_centroid)
    _centroid->printString();

  printf("\n");
}


__host__ __device__
DeviceExtrudedFSR::DeviceExtrudedFSR() :
  _id(-1),
  _num_fsrs(0),
  _material_indices(nullptr),
  _fsr_ids(nullptr),
  _mesh(nullptr)
{}


__host__ __device__
DeviceExtrudedFSR::~DeviceExtrudedFSR() {
  if (_material_indices)
    free(_material_indices);

  if (_fsr_ids)
    free(_fsr_ids);

  if (_mesh)
    free(_mesh);
}


__device__
void DeviceExtrudedFSR::printString() {

  printf("DeviceExtrudedFSR: id = %d, num_fsrs = %d\n", _id, _num_fsrs);

  printf("material indices:\n");
  for (int i = 0; i < _num_fsrs; i++) {
    printf("%d, ", _material_indices[i]);
  }

  printf("\nfsr ids:\n");
  for (int i = 0; i < _num_fsrs; i++) {
    printf("%ld, ", _fsr_ids[i]);
  }

  printf("\naxial mesh:\n");
  for (int i = 0; i < _num_fsrs; i++) {
    printf("%f, ", _mesh[i]);
  }

  printf("\n");
}


//------------------------------------------------------------------------------
// Helper functions
//------------------------------------------------------------------------------
void fsrDataCopyHtoD(DeviceFSRData *fsr_d, FSRData *fsr_h, int material_index) {

  long id = fsr_h->_fsr_id;

  HIP_ASSERT( hipMemcpyHtoD(&fsr_d->_id, &id, sizeof(id)) );
  HIP_ASSERT( hipMemcpyHtoD(&fsr_d->_material_index, &material_index, sizeof(material_index)) );

  // Allocate memory for device points
  DevicePoint *point, *centroid;
  HIP_ASSERT( hipMalloc(&point,    sizeof(DevicePoint)) );
  HIP_ASSERT( hipMalloc(&centroid, sizeof(DevicePoint)) );

  // Copy points from host to device
  pointCopyHtoD(point, fsr_h->_point);
  pointCopyHtoD(centroid, fsr_h->_centroid);
  HIP_ASSERT( hipMemcpyHtoD(&fsr_d->_point,    &point,    sizeof(DevicePoint*)) );
  HIP_ASSERT( hipMemcpyHtoD(&fsr_d->_centroid, &centroid, sizeof(DevicePoint*)) );
}


void fsrDataCopyHtoD(DeviceFSRData *fsr_d, FSRData *fsr_h,
                     std::map<int, int> &material_IDs_to_indices) {

  fsrDataCopyHtoD(fsr_d, fsr_h, material_IDs_to_indices[fsr_h->_mat_id]);
}


void extrudedFSRCopyHtoD(DeviceExtrudedFSR *fsr_d, ExtrudedFSR *fsr_h,
                         std::vector<int> &material_indices) {

  int id = fsr_h->_fsr_id;
  int num_fsrs = fsr_h->_num_fsrs;

  HIP_ASSERT( hipMemcpyHtoD(&fsr_d->_id, &id, sizeof(id)) );
  HIP_ASSERT( hipMemcpyHtoD(&fsr_d->_num_fsrs, &num_fsrs, sizeof(num_fsrs)) );

  // Copy arrays from host to device.
  vectorCopyHtoD(fsr_d->_material_indices, material_indices);
  arrayCopyHtoD(fsr_d->_fsr_ids, fsr_h->_fsr_ids, num_fsrs);

  if (fsr_h->_mesh) {
    arrayCopyHtoD(fsr_d->_mesh, fsr_h->_mesh, num_fsrs);
  }
}


void extrudedFSRCopyHtoD(DeviceExtrudedFSR *fsr_d, ExtrudedFSR *fsr_h,
                         std::map<int, int> &material_IDs_to_indices) {

  // Setup a buffer for material indices.
  std::vector<int> material_indices;
  material_indices.resize(fsr_h->_num_fsrs);

  // Fill the buffer.
  for (int i = 0; i < fsr_h->_num_fsrs; ++i) {
    auto material_id = fsr_h->_materials[i]->getId();
    material_indices.push_back(material_IDs_to_indices[material_id]);
  }

  extrudedFSRCopyHtoD(fsr_d, fsr_h, material_indices);
}


} // namespace hip
} // namespace antmoc
