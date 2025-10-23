/// \file DeviceMaterial.h
/// \brief A material containing cross-sections.
/// \date December 19, 2020
/// \author An Wang, USTB (wangan.cs@gmail.com)

#ifndef DEVICE_MATERIAL_H_
#define DEVICE_MATERIAL_H_

#include "antmoc/hip/hip_utils.h"


namespace antmoc {

// Forward declarations
class Material;

namespace hip {


/// \class Material
/// \brief The Material class represents a unique material and its relevant
///        nuclear data (i.e., multigroup cross-sections) for neutron transport.
struct DeviceMaterial {

  /** A user-defined ID for each Material created */
  int _id;

  /** The number of energy groups */
  int _num_groups;

  /** Max total cross section */
  FP_PRECISION _max_sigma_t;

  /** An array of the total cross-sections for each energy group */
  FP_PRECISION *_sigma_t;

  /** A 2D array of the scattering cross-section matrix from/into each group */
  FP_PRECISION *_sigma_s;

  /** An array of the absorption cross-sections for each energy group */
  FP_PRECISION *_sigma_a;

  /** An array of the fission cross-sections for each energy group */
  FP_PRECISION *_sigma_f;

  /** An array of the fission cross-sections multiplied by nu \f$ \nu \f$
   *  for each energy group */
  FP_PRECISION *_nu_sigma_f;

  /** An array of the chi \f$ \chi \f$ values for each energy group */
  FP_PRECISION *_chi;

  /** A 2D array of the fission matrix from/into each group */
  FP_PRECISION *_fiss_matrix;

  /** A boolean representing whether or not this Material contains a non-zero
   *  fission cross-section and is fissionable */
  bool _fissionable;

  /// \brief Constructor for DeviceMaterial.
  __host__ __device__
  DeviceMaterial();

  // FIXME: not sure if this is the correct way for memory management on device.
  /// \brief Destructor for DeviceMaterial to free GPU memory.
  __host__ __device__
  ~DeviceMaterial();

  /// \brief Get the Material's total cross section for some energy group.
  /// \param group The energy group
  /// \return The total cross section
  __device__
  FP_PRECISION getSigmaTByGroup(int group);

  /// \brief Get the Material's scattering cross section for some energy group.
  /// \param origin      The incoming energy group
  /// \param destination The outgoing energy group
  /// \return The scattering cross section
  __device__
  FP_PRECISION getSigmaSByGroup(int origin, int destination);

  /// \brief Get the Material's absorption cross section for some energy group.
  /// \param group The energy group
  /// \return The absorption cross section
  __device__
  FP_PRECISION getSigmaAByGroup(int group);

  /// \brief Get the Material's fission cross section for some energy group.
  /// \param group The energy group
  /// \return The fission cross section
  __device__
  FP_PRECISION getSigmaFByGroup(int group);

  /// \brief Get the Material's nu-fission cross section for some energy group.
  /// \param group The energy group
  /// \return The nu-fission cross section
  __device__
  FP_PRECISION getNuSigmaFByGroup(int group);

  /// \brief  Get the Material's fission spectrum for some energy group.
  /// \param  group The energy group
  /// \return The fission spectrum
  __device__
  FP_PRECISION getChiByGroup(int group);

  /// \brief  Get the Material's fission matrix for some energy group.
  /// \param  origin      The incoming energy group \f$ E_{0} \f$
  /// \param  destination The outgoing energy group \f$ E_{1} \f$
  /// \return The fission matrix entry \f$ \nu\Sigma_{f}(E_{0}) * \chi(E_{1})\f$
  __device__
  FP_PRECISION getFissionMatrixByGroup(int origin, int destination);

  /// \brief Print myself.
  __device__ void printString();

};  // class DeviceMaterial


/// \brief Copies a material from host to device.
/// \param material_d  Pointer to device material.
/// \param material_h  Pointer to host material.
void materialCopyHtoD(DeviceMaterial *material_d, Material *material_h);

/// \brief Copies a material from device to host.
/// \param material_h  Pointer to host material.
/// \param material_d  Pointer to device material.
void materialCopyDtoH(Material *material_h, Material *material_d);

}   // namespace hip
}   // namespace antmoc

#endif /* DEVICE_MATERIAL_H_ */
