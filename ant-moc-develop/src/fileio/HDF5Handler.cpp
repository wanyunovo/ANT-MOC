#include "antmoc/HDF5Handler.h"
#include "antmoc/file_utils.h"
#include "antmoc/log.h"

#include <sstream>

namespace antmoc
{


/// \details The file id is managed by constructor/destructor
HDF5Handler::HDF5Handler(std::string file, HDF5Mode mode, bool is_parallel, bool for_fsrs)
  : _file(file),
    _is_parallel(is_parallel),
    _for_fsrs(for_fsrs) {

  // Create a new file
  _fapl_id = createFileAPL(is_parallel, for_fsrs);

  switch (mode) {
    case HDF5Mode::Truncate :
      _file_id = createH5File(file, H5F_ACC_TRUNC, _fapl_id);
      break;
    case HDF5Mode::ReadOnly :
      _file_id = openH5File(file, H5F_ACC_RDONLY, _fapl_id);
      break;
    case HDF5Mode::ReadWrite :
      _file_id = openH5File(file, H5F_ACC_RDWR, _fapl_id);
      break;
    default :
      log::error("Undefined HDF5 file mode");
  }

  // A flag
  _file_opened = true;

  log::debug("HDF5Handler accepts file: {}", _file);
}


/// \brief Close HDF5 objects
HDF5Handler::~HDF5Handler() {
  H5Fclose(_file_id);
  H5Pclose(_fapl_id);
}


//----------------------------------------------------------------------
// HDF5 utility
//----------------------------------------------------------------------

/// \brief Discover element names of a group
/// \param groupid id of the group
/// \param names buffer
void HDF5Handler::discoverNames(hid_t groupid, StringVec &names) {

  // Define the operation to be taken during the iteration
  auto op =
    [](hid_t group, const char *name,
       const H5L_info_t *info, void *context)
      {
        auto ctx = (StringVec *)context;
        ctx->push_back(name);
        return 0;
      };

  // Iterate through subelements of the group
  H5Literate(groupid, H5_INDEX_NAME, H5_ITER_NATIVE, NULL,
             op, (void *) &names);
}


/// \brief Create a file access property list
/// \details This method is intended to handle both serial and parallel I/O.
hid_t HDF5Handler::createFileAPL(bool is_parallel, bool for_fsrs) {
  hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);

#ifdef ENABLE_MPI_
  if (is_parallel) {
    log::profile("Creating file access property list for MPI I/O");

    MPI_Comm comm;
    if (for_fsrs)
      comm = getCommForFSRs();
    else
      comm = getCommForTracks();

    MPI_Info info = MPI_INFO_NULL;
    H5Pset_fapl_mpio(plist_id, comm, info);

    // Show the size of the communicator
    int size;
    MPI_Comm_size(comm, &size);
    log::verbose_once("Number of process participating in MPI I/O = {}", size);
  }
#endif

  return plist_id;
}


/// \brief Create a data transfer property list
/// \details This method is intended to handle both serial and parallel I/O.
hid_t HDF5Handler::createDatasetXPL(bool is_parallel) {
  hid_t plist_id = H5Pcreate(H5P_DATASET_XFER);

#ifdef ENABLE_MPI_
  if (is_parallel)
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
#endif

  return plist_id;
}


//--------------------------------------
// File interfaces
//--------------------------------------

/// \brief Create an H5 file and return its id
/// \param flags defaults to H5F_ACC_TRUNC
/// \param fapl_id id of the file access property list (defaults to
///                H5P_DEFAULT)
hid_t HDF5Handler::createH5File(std::string file, unsigned flags,
                                hid_t fapl_id) {
  log::verbose_once("Creating H5 file '{}'", file);
  auto file_id = H5Fcreate(file.c_str(),
                           flags,
                           H5P_DEFAULT,
                           fapl_id);
  return file_id;
}


/// \brief Open an H5 file and return its id
/// \param flags defaults to H5F_ACC_RDONLY
/// \param fapl_id id of the file access property list (defaults to
///                H5P_DEFAULT)
hid_t HDF5Handler::openH5File(std::string file, unsigned flags,
                              hid_t fapl_id) {
  if (!fileutils::existsFile(file))
    log::error("Cannot find HDF5 file: '{}'", file);

  auto file_id = H5Fopen(file.c_str(),
                         flags,
                         fapl_id);

  if (file_id < 0)
    log::error("Failed to open HDF5 file: '{}'", file);

  return file_id;
}


//--------------------------------------
// Group interfaces
//--------------------------------------

/// \brief Create an H5 group and return its id
hid_t HDF5Handler::createH5Group(hid_t loc_id, std::string path) {
  auto group_id = H5Gcreate(loc_id,
                            path.c_str(),
                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  return group_id;
}


/// \brief Open an H5 group and return its id
hid_t HDF5Handler::openH5Group(hid_t loc_id, std::string path) {
  auto group_id = H5Gopen(loc_id,
                          path.c_str(),
                          H5P_DEFAULT);
  return group_id;
}


/// \brief Determines whether an H5 group exists
/// \param loc_id a group
/// \param path relative path to the group
hid_t HDF5Handler::existsH5Group(hid_t loc_id, std::string path) {
  htri_t ret = H5Lexists(loc_id,
                         path.c_str(),
                         H5P_DEFAULT);
  return (ret > 0);
}


//--------------------------------------
// Dataset interfaces
//--------------------------------------

/// \brief Create an H5 dataset and return its id
hid_t HDF5Handler::createH5Dataset(hid_t loc_id, std::string path,
                                  hid_t type_id, hid_t space_id) {
  auto dataset_id = H5Dcreate(loc_id, path.c_str(),
                              type_id, space_id,
                              H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  return dataset_id;
}


/// \brief Open an H5 dataset and return its id
hid_t HDF5Handler::openH5Dataset(hid_t loc_id, std::string path) {
  auto dataset_id = H5Dopen(loc_id,
                            path.c_str(),
                            H5P_DEFAULT);
  return dataset_id;
}


/// \brief Determines whether an H5 dataset exists
/// \param loc_id a group
/// \param path relative path to the dataset
hid_t HDF5Handler::existsH5Dataset(hid_t loc_id, std::string path) {
  htri_t ret = H5Lexists(loc_id,
                         path.c_str(),
                         H5P_DEFAULT);
  return (ret > 0);
}


/// \details Transfer property list for this I/O operation is pre-defined.
///          If the I/O is parallel, the xpl will be parallel as well.
herr_t HDF5Handler::writeH5Dataset(hid_t dataset_id, hid_t mem_type_id,
                                   hid_t mem_space_id, hid_t file_space_id,
                                   const void *buf) {
  auto dataset_xpl= createDatasetXPL();
  auto ret = H5Dwrite(dataset_id,
                      mem_type_id,
                      mem_space_id, file_space_id,
                      dataset_xpl,
                      buf);
  H5Pclose(dataset_xpl);
  return ret;
}


/// \details See writeH5Dataset(...).
herr_t HDF5Handler::readH5Dataset(hid_t dataset_id, hid_t mem_type_id,
                                  hid_t mem_space_id, hid_t file_space_id,
                                  void *buf) {
  auto dataset_xpl= createDatasetXPL();
  auto ret = H5Dread(dataset_id,
                     mem_type_id,
                     mem_space_id, file_space_id,
                     dataset_xpl,
                     buf);
  H5Pclose(dataset_xpl);
  return ret;
}


//--------------------------------------
// Attribute interfaces
//--------------------------------------

/// \brief Create an H5 attribute and return its id
hid_t HDF5Handler::createH5Attribute(hid_t loc_id, std::string name,
                                     hid_t type_id, hid_t space_id) {
  auto attr_id = H5Acreate(loc_id,
                           name.c_str(),
                           type_id,
                           space_id,
                           H5P_DEFAULT, H5P_DEFAULT);
  return attr_id;
}


/// \brief Writes data to an attribute
herr_t HDF5Handler::writeH5Attribute(hid_t attr_id, hid_t mem_type_id,
                                     const void *buf) {
  return H5Awrite(attr_id, mem_type_id, buf);
}


/// \brief Reads data from an attribute
herr_t HDF5Handler::readH5Attribute(hid_t attr_id, hid_t mem_type_id,
                                    void *buf) {
  return H5Aread(attr_id, mem_type_id, buf);
}


/// \brief Open an H5 attribute and return its id
hid_t HDF5Handler::openH5Attribute(hid_t obj_id, std::string name) {
  auto attr_id = H5Aopen(obj_id,
                         name.c_str(),
                         H5P_DEFAULT);
  return attr_id;
}


/// \brief Determines whether an H5 attribute exists
bool HDF5Handler::existsH5Attribute(hid_t obj_id, std::string name) {
  // Note that if this function fails, ret < 0
  htri_t ret = H5Aexists(obj_id, name.c_str());
  return (ret > 0);
}


//--------------------------------------
// Dataspace interfaces
//--------------------------------------

int HDF5Handler::getH5SpaceSimpleExtentDims(hid_t dataset_id, hsize_t *dims) {
  auto space_id = H5Dget_space(dataset_id);
  auto ndims = H5Sget_simple_extent_dims(space_id, dims, NULL);

  H5Sclose(space_id);

  return ndims;
}


int HDF5Handler::getH5SpaceSimpleExtentDims(hid_t loc_id,
                                            const std::string &dataset_path,
                                            hsize_t *dims) {
  auto dataset_id = openH5Dataset(loc_id, dataset_path);
  auto space_id = H5Dget_space(dataset_id);
  auto ndims = H5Sget_simple_extent_dims(space_id, dims, NULL);

  H5Sclose(space_id);
  H5Dclose(dataset_id);

  return ndims;
}

//----------------------------------------------------------------------
// Parallel I/O utility
//----------------------------------------------------------------------

#ifdef ENABLE_MPI_

/// \brief Get a communicator within which processes have their own FSRs
/// \details For now, the spatial domain decomposition won't be enabled
///          at the same time the periodic track decomposition is enabled.
MPI_Comm HDF5Handler::getCommForFSRs() {
  return mpi::getCommUniqueDomains();
}


/// \brief Get a communicator within which processes have their own tracks
/// \details For now, the spatial domain decomposition won't be enabled
///          at the same time the periodic track decomposition is enabled.
MPI_Comm HDF5Handler::getCommForTracks() {
  return mpi::getMPIComm();
}


/// \brief Get the number of processes for writing FSR data
int HDF5Handler::getNumProcsForFSRs() {
  return mpi::getNumUniqueDomains();
}


/// \brief Get my rank for writing FSR data
int HDF5Handler::getMPIRankForFSRs() {
  return mpi::getRankUniqueDomains();
}

#endif

} // namespace antmoc
