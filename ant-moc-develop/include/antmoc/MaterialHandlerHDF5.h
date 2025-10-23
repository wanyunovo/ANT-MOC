/// \file MaterialHandlerHDF5.h
/// \brief Handling materials data in HDF5 format.
/// \date April 22, 2019
/// \author Gen Wang, USTB
/// \author Ya Fang, USTB (fangya201388@163.com)
/// \author An Wang, USTB (wangan@xs.ustb.edu.cn)


#ifndef MATERIAL_HANDLER_HDF5_H_
#define MATERIAL_HANDLER_HDF5_H_


#include <map>
#include <string>
#include <vector>

#include "antmoc/HDF5Handler.h"
#include "antmoc/MaterialHandler.h"

namespace antmoc
{


///---------------------------------------------------------------------
/// \class MaterialHandlerHDF5
/// \brief Methods for writing and reading materials in HDF5 format
///---------------------------------------------------------------------
class MaterialHandlerHDF5 : public HDF5Handler, public MaterialHandler {

  ///< Data layout of the material file
  XSFileLayout _layout = XSFileLayout::NAMED;

  ///< Name of the top group
  std::string _domain_type = "material";

  ///< ID of the top group
  hid_t _domain_id;

public:

  /// \brief Create a material handler for HDF5 files.
  /// \param file The path to the materials file.
  /// \param mode File access mode.
  /// \param layout Data layout of the material file.
  MaterialHandlerHDF5(std::string file, HDF5Mode mode = HDF5Mode::ReadOnly,
                      XSFileLayout layout = XSFileLayout::NAMED);
  virtual ~MaterialHandlerHDF5();

  /// \brief Read all Multi-group XS of a material from an HDF5 file.
  void readAllMaterials();

  /// \brief Read Multi-group XS of a material from an HDF5 file.
  /// \param name Name of the material.
  void readMaterial(std::string name);

  /// \brief Check the material map
  bool hasMaterials();

  /// \brief Check the existence of a specific material
  /// \param name Name of the material.
  bool hasMaterial(std::string name);

  /// \brief Return the total number of materials in the file
  size_t getTotalNumMaterials();

private:

  /// \brief Read all materials which are stored in a compressed format.
  void readAllCompressedMaterials();

  /// \brief Read a material by name.
  void readCompressedMaterial(std::string mat_name);

  /// \brief Read all materials which are stored in an OpenMOC-compatible format.
  void readAllNamedMaterials();

  /// \brief Read a material by name;
  void readNamedMaterial(std::string mat_name);

  /// \brief A helper method for compressed material data.
  /// \param domain_mat ID of the top group.
  std::map<std::string, int8_t> readIdxMapAttributes(hid_t domain_id);

  /// \brief Check whether a specified XS dataset exists.
  /// \param mat_name Name of the material, which is also the name of an HDF5 group.
  /// \param xs_name Name of the XS, which is also the name of an HDF5 dataset.
  bool existsXS(const std::string &mat_name, const std::string &xs_name);

  /// \brief Retrieve specified XS data from a dataset.
  /// \details The number of dimensions of the dataset should not be less than ndim.
  ///          It is allowed to read part of the data. If that is so, there will be
  ///          a warning.
  /// \param mat_name Name of the material, which is also the name of an HDF5 group.
  /// \param xs_name Name of the XS, which is also the name of an HDF5 dataset.
  /// \param offsets The starting index of each dimensions of the dataset.
  /// \param shape The size of each dimensions of the dataset.
  /// \param buffer Where to store the XS data (1-D array).
  /// \return Whether the dataset exists.
  void readXS(std::string mat_name, std::string xs_name,
              std::vector<hsize_t> offsets,
              std::vector<hsize_t> shape, double *buffer);
};

} // namespace antmoc

#endif  // MATERIAL_HANDLER_HDF5_H_
