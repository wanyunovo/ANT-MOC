#include "antmoc/MaterialHandlerHDF5.h"
#include "antmoc/log.h"
#include "antmoc/Material.h"
#include "antmoc/string_utils.h"

#include <vector>

namespace antmoc
{

MaterialHandlerHDF5::
MaterialHandlerHDF5(std::string file, HDF5Mode mode, XSFileLayout layout)
  : HDF5Handler(file, mode, false, false),
    _layout(layout) {

  // Get the id of H5 file
  auto h5_file = getFileId();

  // Get the number of energy groups from the file
  readScalarAttribute(h5_file, "# groups", _num_groups_in_file);
  setNumEnergyGroups(getNumEnergyGroupsInFile());

  // Open the top group and get the number of materials
  _domain_id = openH5Group(h5_file, _domain_type);
}


MaterialHandlerHDF5::~MaterialHandlerHDF5() {

  H5Gclose(_domain_id);

}


/// \details The handler has materials if any material has been retrieved or
///          the file has defined any material name.
bool MaterialHandlerHDF5::hasMaterials() {
  if (MaterialHandler::hasMaterials())
    return true;
  else
    return getTotalNumMaterials() > 0;
}


/// \details The handler has a material if the material name can be found
///          in the file.
bool MaterialHandlerHDF5::hasMaterial(std::string name) {
  if (MaterialHandler::hasMaterial(name))
    return true;
  else
    return existsH5Group(_domain_id, name);
}


size_t MaterialHandlerHDF5::getTotalNumMaterials() {
  StringVec mat_names;
  discoverNames(_domain_id, mat_names);
  return mat_names.size();
}


/// \details XS data is organized in a specified way, which is specified
///          in the constructor.
void MaterialHandlerHDF5::readAllMaterials() {

  switch (_layout) {
    case XSFileLayout::NAMED :
      readAllNamedMaterials();
      break;
    case XSFileLayout::COMPRESSED :
      readAllCompressedMaterials();
      break;
  }

  log::info("Read {} materials during readAllMaterials", _materials.size());
}


/// \details XS data is organized in a specified way, which is specified
///          in the constructor.
void MaterialHandlerHDF5::readMaterial(std::string name) {

  switch (_layout) {
    case XSFileLayout::NAMED :
      readNamedMaterial(name);
      break;
    case XSFileLayout::COMPRESSED :
      readCompressedMaterial(name);
      break;
  }
}


/// \details Call readCompressedMaterial() repeatedly.
void MaterialHandlerHDF5::readAllCompressedMaterials() {

  // Discover material names
  StringVec mat_names;
  discoverNames(_domain_id, mat_names);

  // Loop over all materials
  for (auto &mat_name : mat_names)
    readCompressedMaterial(mat_name);
}


/// \details XS data is organized as follows
///          The top group is named 'material'. Cross-section data of each
///          material is stored in an subgroup of 'material' with the name
///          of the material as the group name. Different kinds of XS are
///          separated by datasets.
///          Groups:
///          /material/name1
///          /material/name2
///          ...
///
///          Datasets:
///          /material/name1/reactions
///          /material/name1/scattering
///
///          'reactions' includes absorption, fission, nu-fission,
///          transport and chi.
///          The index of each of them are stored in an attribute:
///          /material/idx absorption
///          /material/idx fission
///          /material/idx nu-fission
///          /material/idx transport
///          /material/idx chi
void MaterialHandlerHDF5::readCompressedMaterial(std::string mat_name) {

  // Do nothing if there is a material with the same name, or
  // the required material doesn't exist.
  if (_materials.count(mat_name) || !existsH5Group(_domain_id, mat_name))
    return;

  // Get the number of energy groups
  const hsize_t num_groups = getNumEnergyGroups();
  validateNumEnergyGroups();

  log::debug("Reading {}-group material '{}' in compressed format...", num_groups, mat_name);

  const std::string name_reactions = "reactions";
  const std::string name_scatter   = "scattering";

  // Check if there are enough datasets to read
  if (!existsXS(mat_name, name_reactions) ||
      !existsXS(mat_name, name_scatter)) {
    log::error("Material {} doesn't has dataset '{}' or '{}'", mat_name, name_reactions, name_scatter);
  }

  //------------------------------------------------------------------
  // Create a material object
  //------------------------------------------------------------------
  Material *material = new Material(0, mat_name.c_str());
  _materials.emplace(mat_name, material);

  material->setNumEnergyGroups(num_groups);

  //------------------------------------------------------------------
  // Reading dataset 'reactions' from this data group
  //------------------------------------------------------------------
  // Buffer
  std::string xs_name;
  auto sigma = new double[num_groups * num_groups]();

  // The indices must exist for locating XS
  auto idx_map = readIdxMapAttributes(_domain_id);

  // Read 'reactions' dataset
  std::vector<hsize_t> offsets_required = {0, 0};
  std::vector<hsize_t> shape_required   = {num_groups, 1};

  auto readReactionXS = [&](double *buffer, const std::string &xs_name) {
    log::debug("Reading '{}/{}/{}', idx = {}", mat_name, name_reactions, xs_name, idx_map[xs_name]);
    offsets_required[1] = idx_map[xs_name];
    readXS(mat_name, name_reactions, offsets_required, shape_required, buffer);
  };

  if (xs_name = "chi", idx_map.count(xs_name)) {
    readReactionXS(sigma, xs_name);
    material->setChi(sigma, num_groups);
  }
  else {
    log::warn("No 'chi' MGXS found for '{}/{}'", _domain_type, mat_name);
  }

  if (xs_name = "nu-fission", idx_map.count(xs_name)) {
    readReactionXS(sigma, xs_name);
    material->setNuSigmaF(sigma, num_groups);
  }
  else {
    log::warn("No 'nu-fission' MGXS found for '{}/{}'", _domain_type, mat_name);
  }

  if (xs_name = "total", idx_map.count(xs_name)) {
    readReactionXS(sigma, xs_name);
    material->setSigmaT(sigma, num_groups);
  }
  else if (xs_name = "transport", idx_map.count(xs_name)) {
    readReactionXS(sigma, xs_name);
    material->setSigmaT(sigma, num_groups);
  }
  else {
    log::warn("No 'total' or 'transport' MGXS found for '{}/{}'", _domain_type, mat_name);
  }

  if (xs_name = "fission", idx_map.count(xs_name)) {
    readReactionXS(sigma, xs_name);
    material->setSigmaF(sigma, num_groups);
  }
  else if (xs_name = "absorption", idx_map.count(xs_name)) {
    readReactionXS(sigma, xs_name);
    material->setSigmaA(sigma, num_groups);
  }
  else {
    log::warn("No 'fission' or 'absorption' MGXS found for '{}/{}'", _domain_type, mat_name);
  }

  //------------------------------------------------------------------
  // Reading dataset 'scattering', currently no high-order scattering rates
  //------------------------------------------------------------------
  offsets_required = {0, 0};
  shape_required   = {num_groups, num_groups};

  readXS(mat_name, name_scatter, offsets_required, shape_required, sigma);
  material->setSigmaS(sigma, num_groups *num_groups);

  // Delete the buffer
  delete [] sigma;
}


/// \details Cass readNamedMaterial() repeatedly.
void MaterialHandlerHDF5::readAllNamedMaterials() {

  // Discover material names
  StringVec mat_names;
  discoverNames(_domain_id, mat_names);

  // Loop over all materials
  for (auto &mat_name : mat_names)
    readNamedMaterial(mat_name);
}


/// \details XS data is organized as follows
///          The top group is named 'material'. Cross-section data of each
///          material is stored in an subgroup of 'material' with the name
///          of the material as the group name. Different kinds of XS are
///          separated by datasets.
///          Groups:
///          /material/name1
///          /material/name2
///          ...
///
///          Datasets:
///          /material/name1/absorption
///          /material/name1/fission
///          /material/name1/nu-fission
///          /material/name1/total
///              or 'transport'
///          /material/name1/chi
///          /material/name1/scatter matrix
///              or 'nu-scatter matrix'
///              or 'consistent scatter matrix'
void MaterialHandlerHDF5::readNamedMaterial(std::string mat_name) {

  // Do nothing if there is a material with the same name, or
  // the required material doesn't exist.
  if (_materials.count(mat_name) || !existsH5Group(_domain_id, mat_name))
    return;

  // Get the number of energy groups
  const hsize_t num_groups = getNumEnergyGroups();
  validateNumEnergyGroups();

  log::debug("Reading {}-group material '{}' in named format...", num_groups, mat_name);

  // Create a new material object
  Material *material = new Material(0, mat_name.c_str());
  _materials.emplace(mat_name, material);

  material->setNumEnergyGroups(num_groups);

  // Data buffers, no high-order scattering rates
  std::string xs_name;
  auto sigma = new double[num_groups * num_groups]();
  std::vector<hsize_t> offsets_required = {0};
  std::vector<hsize_t> shape_required   = {num_groups};

  auto readXSSimple = [&] (double *buffer, const std::string &xs_name) {
    readXS(mat_name, xs_name, offsets_required, shape_required, buffer);
  };

  // Search for the total/transport cross section
  if ( xs_name = "total", existsXS(mat_name, xs_name) ) {
    readXSSimple(sigma, xs_name);
    material->setSigmaT(sigma, num_groups);
  }
  else if ( xs_name = "transport", existsXS(mat_name, xs_name) ) {
    readXSSimple(sigma, xs_name);
    material->setSigmaT(sigma, num_groups);
  }
  else {
    log::warn("No 'total' or 'transport' MGXS found for '{}/{}'", _domain_type, mat_name);
  }

  // Search for the fission production cross section
  if ( xs_name = "nu-fission", existsXS(mat_name, xs_name) ) {
    readXSSimple(sigma, xs_name);
    material->setNuSigmaF(sigma, num_groups);
  }
  else {
    log::warn("No 'nu-fission' MGXS found for '{}/{}'", _domain_type, mat_name);
  }

  // Search for chi (fission spectrum)
  if ( xs_name = "chi", existsXS(mat_name, xs_name) ) {
    readXSSimple(sigma, xs_name);
    material->setChi(sigma, num_groups);
  }
  else {
    log::warn("No 'chi' MGXS found for '{}/{}'", _domain_type, mat_name);
  }

  // Search for optional cross sections
  if ( xs_name = "fission", existsXS(mat_name, xs_name) ) {
    readXSSimple(sigma, xs_name);
    material->setSigmaF(sigma, num_groups);
  }
  else if ( xs_name = "absorption", existsXS(mat_name, xs_name) ) {
    readXSSimple(sigma, xs_name);
    material->setSigmaA(sigma, num_groups);
  }
  else {
    log::warn("No 'fission' or 'absorption' MGXS found for '{}/{}'", _domain_type, mat_name);
  }

  // Search for the scattering matrix cross section
  shape_required = {num_groups * num_groups};
  if ( xs_name = "nu-scatter matrix", existsXS(mat_name, xs_name) ) {
    readXSSimple(sigma, xs_name);
    material->setSigmaS(sigma, num_groups * num_groups);
  }
  else if ( xs_name = "scatter matrix", existsXS(mat_name, xs_name) ) {
    readXSSimple(sigma, xs_name);
    material->setSigmaS(sigma, num_groups * num_groups);
  }
  else if ( xs_name = "consistent nu-scatter matrix", existsXS(mat_name, xs_name) ) {
    readXSSimple(sigma, xs_name);
    material->setSigmaS(sigma, num_groups * num_groups);
  }
  else {
    log::warn("No 'scatter matrix' found for '{}/{}'", _domain_type, mat_name);
  }

  delete[] sigma;
}


bool MaterialHandlerHDF5::existsXS(const std::string &mat_name,
                                   const std::string &xs_name) {
  return existsH5Dataset(_domain_id, fmt::format("{}/{}", mat_name, xs_name));
}


void MaterialHandlerHDF5::readXS(std::string mat_name, std::string xs_name,
                                 std::vector<hsize_t> offsets,
                                 std::vector<hsize_t> shape, double *buffer) {

  if (shape.empty())
    return;

  // Get the number of dimensions required
  int ndim_required = shape.size();
  offsets.resize(ndim_required, 0);

  log::debug("Reading XS dataset '{}/{}': require {} dim(s), shape = ({}), offsets = ({})",
              mat_name, xs_name, ndim_required,
              stringutils::join(shape, ","), stringutils::join(offsets, ","));

  // Open the required dataset
  auto dataset_id = openH5Dataset(_domain_id, fmt::format("{}/{}", mat_name, xs_name));

  // Read shape from the file and check the consistency
  hsize_t shape_in_file[3] = {1, 1, 1};
  const int ndim_in_file = getH5SpaceSimpleExtentDims(dataset_id, shape_in_file);

  if (ndim_required < ndim_in_file) {
    // Resize the shape and offset vectors when needed
    shape.resize(ndim_in_file, 1);
    offsets.resize(ndim_in_file, 0);
  }
  else if (ndim_required > ndim_in_file) {
    log::error("Ill-formed XS file: dataset '{}/{}' has only {} dimensions but {} are required.",
                mat_name, xs_name, ndim_in_file, ndim_required);
  }


  // Reading part of energy groups is allowed with a warning.
  if (shape[0] < shape_in_file[0]) {
    log::warn_once("Trying to truncate the XS dataset: take the first {} elements.", shape[0]);
  }
  else if (shape[0] > shape_in_file[0]) {
    log::error("When reading '{}/{}': not enough elements to be read in.", mat_name, xs_name);
  }

  // Construct data spaces for memory and dataset layout.
  hsize_t mspace_count = 1;
  for (int i = 0; i < ndim_required; ++i) {
    mspace_count *= shape[i];
  }

  auto mspace = H5Screate_simple(1, &mspace_count, NULL);
  auto dspace = H5Dget_space(dataset_id);

  // Select the data in the dataset
  H5Sselect_hyperslab(dspace, H5S_SELECT_SET, offsets.data(), NULL, shape.data(), NULL);

  // Read the dataset into the buffer
  auto mem_datatype = getH5MemoryDatatype<double>();
  readH5Dataset(dataset_id, mem_datatype, mspace, dspace, buffer);

  // Close H5 objects
  H5Sclose(mspace);
  H5Sclose(dspace);
  H5Dclose(dataset_id);
}


std::map<std::string, int8_t>
MaterialHandlerHDF5::readIdxMapAttributes(hid_t domain_id) {

  std::map<std::string, int8_t> idx_map;

  // Read indices from attributes
  int8_t idx;
  auto readAttr =
    [this, &domain_id, &idx]
    (std::string name) {
      if (existsH5Attribute(domain_id, name)) {
        readScalarAttribute(domain_id, name, idx);
        return true;
      }
      else {
        return false;
      }
    };

  std::string xs_name;

  // Search for the total/transport cross section
  if ( xs_name = "total", readAttr("idx " + xs_name) )
    idx_map.insert({xs_name, idx});
  else if ( xs_name = "transport", readAttr("idx " + xs_name) )
    idx_map.insert({xs_name, idx});
  else {
    log::warn("No 'idx total' or 'idx transport' found");
  }

  // Search for the fission production cross section
  if ( xs_name = "nu-fission", readAttr("idx " + xs_name) )
    idx_map.insert({xs_name, idx});
  else {
    log::warn("No 'idx nu-fission' found");
  }

  // Search for chi (fission spectrum)
  if ( xs_name = "chi", readAttr("idx " + xs_name) )
    idx_map.insert({xs_name, idx});
  else {
    log::warn("No 'idx chi' found");
  }

  // Search for optional cross sections
  if ( xs_name = "fission", readAttr("idx " + xs_name) )
    idx_map.insert({xs_name, idx});
  else if ( xs_name = "absorption", readAttr("idx " + xs_name) )
    idx_map.insert({xs_name, idx});
  else {
    log::warn("No 'idx fission' or 'idx absorption' found");
  }

  return idx_map;
}


} // namespace antmoc
