/// \file   DeviceFSRDataTypes.h
/// \brief  FSR data structures on device.
/// \date   December 31, 2020
/// \author An Wang (wangan@xs.ustb.edu.cn)

#ifndef DEVICE_FSR_DATA_TYPES_H_
#define DEVICE_FSR_DATA_TYPES_H_

#include "antmoc/hip/hip_utils.h"
#include <cstddef>
#include <map>
#include <vector>

namespace antmoc {

// Forward declarations
struct ExtrudedFSR;
struct FSRData;
class Material;

namespace hip {

// Forward declarations
struct DevicePoint;


/// \struct DeviceFSRData
/// \brief  See FSRDataTypes.h
struct DeviceFSRData {

  ///< The FSR ID
  long _id;

  ///< The Material index in the contiguous material array
  int _material_index;

  ///< Characteristic point in Root Universe that lies in FSR
  DevicePoint *_point;

  ///< Global numerical centroid in Root Universe
  DevicePoint *_centroid;

  /// \brief Constructor
  __host__ __device__
  DeviceFSRData();

  /// \brief Destructor
  __host__ __device__
  ~DeviceFSRData();

  /// \brief Print myself.
  __device__ void printString();
};


/// \struct DeviceExtrudedFSR
/// \brief  See FSRDataTypes.h
struct DeviceExtrudedFSR {

  ///< Axial extruded FSR ID
  int _id;

  ///< Number of FSRs in the axially extruded FSR
  int _num_fsrs;

  ///< Array of material indices for each FSR
  int* _material_indices;

  ///< Array of 3D FSR IDs
  long* _fsr_ids;

  ///< Array defining the axial mesh
  double* _mesh;

  /// \brief Constructor
  __host__ __device__
  DeviceExtrudedFSR();

  /// \brief Destructor
  __host__ __device__
  ~DeviceExtrudedFSR();

  /// \brief Print myself.
  __device__ void printString();
};


/// \brief   Copy an FSRData from host to device.
/// \details This method provides a way to set the material index directly.
/// \param   fsr_d Pointer to a DeviceFSRData object.
/// \param   fsr_h Pointer to an FSRData object.
/// \param   material_index Index of the material the FSR contains.
void fsrDataCopyHtoD(DeviceFSRData *fsr_d, FSRData *fsr_h, int material_index);


/// \brief Copy an FSRData from host to device.
/// \param fsr_d Pointer to a DeviceFSRData object.
/// \param fsr_h Pointer to an FSRData object.
/// \param material_IDs_to_indices Map of material IDs to contiguous indices.
void fsrDataCopyHtoD(DeviceFSRData *fsr_d, FSRData *fsr_h,
                     std::map<int, int> &material_IDs_to_indices);


/// \brief   Copy an ExtrudedFSR from host to device.
/// \details This method provides a way to set material indices directly.
/// \param   fsr_d Pointer to a DeviceExtrudedFSR object.
/// \param   fsr_h Pointer to an ExtrudedFSR object.
/// \param   material_indices Contiguous indices of materials the FSR contains.
void extrudedFSRCopyHtoD(DeviceExtrudedFSR *fsr_d, ExtrudedFSR *fsr_h,
                         std::vector<int> &material_indices);


/// \brief Copy an ExtrudedFSR from host to device.
/// \param fsr_d Pointer to a DeviceExtrudedFSR object.
/// \param fsr_h Pointer to an ExtrudedFSR object.
/// \param material_IDs_to_indices Map of material IDs to contiguous indices.
void extrudedFSRCopyHtoD(DeviceExtrudedFSR *fsr_d, ExtrudedFSR *fsr_h,
                         std::map<int, int> &material_IDs_to_indices);


} // namespace hip
} // namespace antmoc

#endif /* DEVICE_FSR_DATA_TYPES_H_ */
