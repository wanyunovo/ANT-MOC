/// \file HDF5Handler.h
/// \brief An I/O class for HDF5 files
/// \date March 6, 2020
/// \author An Wang, USTB (wangan@xs.ustb.edu.cn)

#ifndef HDF5_HANDLER_H_
#define HDF5_HANDLER_H_

#include <string>
#include <typeindex>
#include <typeinfo>
#include <unordered_map>

#include "hdf5.h"

#include "antmoc/log.h"
#include "antmoc/mpi_utils.h"
#include "antmoc/string_utils.h"


namespace antmoc
{


///---------------------------------------------------------------------
/// \enum HDF5 file mode
///---------------------------------------------------------------------
enum class HDF5Mode {
  ///< Create a file, and truncate the file if is exists
  Truncate,

  ///< Open a file in read-only mode
  ReadOnly,

  ///< Open a file in read-write mode
  ReadWrite
};


///---------------------------------------------------------------------
/// \class HDF5Handler
/// \brief Methods for writing HDF5 files
///---------------------------------------------------------------------
class HDF5Handler {

protected:

  ///< Path to the file
  std::string _file;

  ///< HDF5 file id
  hid_t _file_id;

  ///< File access property list
  hid_t _fapl_id;

  ///< A flag instructing whether MPI I/O is enabled
  bool _is_parallel;

  ///< A flag instructing wheter the file is created for FSRs
  bool _for_fsrs;

  ///< A flag to indicate file status
  bool _file_opened;

public:

  /// \brief Constructor of HDF5Handler
  /// \param mode a flag instructing whether to create or open a file
  /// \param is_parallel decide whether MPI I/O is enabled
  /// \param for_fsrs instruct whether the file is created or opened for FSRs
  HDF5Handler(std::string file, HDF5Mode mode = HDF5Mode::ReadOnly,
              bool is_parallel = false, bool for_fsrs = true);

  HDF5Handler(const HDF5Handler&) = delete;
  virtual ~HDF5Handler();

  //--------------------------------------------------------------------
  // HDF5 utility
  //--------------------------------------------------------------------
  static void discoverNames(hid_t groupid, StringVec &names);

  hid_t getFileId() { return _file_id; }
  bool isParallel() { return _is_parallel; }

  /// \brief Returns HDF5 pre-defined native datatype
  /// \param T An integral or floating-point type.
  template <typename T> hid_t getH5MemoryDatatype();

  /// \brief Returns HDF5 pre-defined file datatype
  /// \param T An integral or floating-point type.
  template <typename T> hid_t getH5FileDatatype();

  hid_t createFileAPL(bool is_parallel = false, bool for_fsrs = true);
  hid_t createDatasetXPL(bool is_parallel = false);

  //--------------------------------------
  // File interfaces
  //--------------------------------------
  hid_t createH5File(std::string file, unsigned flags = H5F_ACC_TRUNC,
                     hid_t fapl_id = H5P_DEFAULT);
  hid_t openH5File(std::string file, unsigned flags = H5F_ACC_RDONLY,
                   hid_t fapl_id = H5P_DEFAULT);

  //--------------------------------------
  // Group interfaces
  //--------------------------------------
  hid_t createH5Group(hid_t loc_id, std::string path);
  hid_t openH5Group(hid_t loc_id, std::string path);
  hid_t existsH5Group(hid_t loc_id, std::string path);

  //--------------------------------------
  // Dataset interfaces
  //--------------------------------------
  hid_t createH5Dataset(hid_t loc_id, std::string path,
                        hid_t type_id, hid_t space_id);
  hid_t openH5Dataset(hid_t loc_id, std::string path);
  hid_t existsH5Dataset(hid_t loc_id, std::string path);


  /// \brief Writes data to an H5 dataset.
  /// \param dataset_id ID of the H5 dataset.
  /// \param mem_type_id ID of the memory datatype.
  /// \param mem_space_id ID of the memory dataspace.
  /// \param file_space_id ID of the dataset's dataspace in the file.
  /// \param buf Buffer with data to be written to the file. 
  herr_t writeH5Dataset(hid_t dataset_id, hid_t mem_type_id,
                        hid_t mem_space_id, hid_t file_space_id,
                        const void *buf);

  /// \brief Reads data from an H5 dataset.
  herr_t readH5Dataset(hid_t dataset_id, hid_t mem_type_id,
                       hid_t mem_space_id, hid_t file_space_id,
                       void *buf);

  //--------------------------------------
  // Attribute interfaces
  //--------------------------------------
  hid_t createH5Attribute(hid_t loc_id, std::string name,
                          hid_t type_id, hid_t space_id);
  hid_t openH5Attribute(hid_t obj_id, std::string name);
  herr_t writeH5Attribute(hid_t attr_id, hid_t mem_type_id,
                          const void *buf);
  herr_t readH5Attribute(hid_t attr_id, hid_t mem_type_id,
                         void *buf);
  bool existsH5Attribute(hid_t obj_id, std::string name);

  //--------------------------------------
  // Dataspace interfaces
  //--------------------------------------
  /// \brief Retrieves dataspace dimension size of a dataset.
  /// \param dataset_id ID of the dataset.
  /// \param dims Pointer to array to store the size of each dimension.
  /// \return The number of dimensions.
  int getH5SpaceSimpleExtentDims(hid_t dataset_id, hsize_t *dims);

  /// \brief Retrieves dataspace dimension size of a dataset.
  /// \param loc_id ID of the H5 object to which the dataset is attached
  /// \param dataset_path The relative path of the dataset
  /// \param dims Pointer to array to store the size of each dimension.
  /// \return The number of dimensions.
  int getH5SpaceSimpleExtentDims(hid_t loc_id, const std::string &dataset_path,
                                 hsize_t *dims);

  //--------------------------------------------------------------------
  // ANTMOC utility
  //--------------------------------------------------------------------
  /// \brief Writes an attribute of type T into an H5 object. T must be an
  ///        integral type or an floating-point type.
  /// \param obj_id ID of the H5 object to which attribute is attached
  /// \param name Attribute name.
  /// \param value A value of type T.
  template <typename T>
  void writeScalarAttribute(hid_t obj_id, std::string name, T value);

  /// \brief Reads an attribute of type T from an H5 object. T must be an
  ///        integral type or an floating-point type.
  template <typename T>
  void readScalarAttribute(hid_t obj_id, std::string name, T &value);


  //--------------------------------------------------------------------
  // Parallel I/O utility
  //--------------------------------------------------------------------
#ifdef ENABLE_MPI_
  MPI_Comm getCommForFSRs();
  MPI_Comm getCommForTracks();
  int getNumProcsForFSRs();
  int getMPIRankForFSRs();

#endif

};


/// \details This method maps typeids of c++ integral types and floating-point
///          types to HDF5 pre-defined native datatypes.
template <typename T>
hid_t HDF5Handler::getH5MemoryDatatype() {
  static_assert(std::is_integral<T>::value || std::is_floating_point<T>::value,
                "Only integral datatypes and floating-point datatypes are "
                "pre-defined by HDF5");

  std::unordered_map<std::type_index, hid_t> type_map = {
    {typeid(int8_t),              H5T_NATIVE_INT8},
    {typeid(int),                 H5T_NATIVE_INT},
    {typeid(short),               H5T_NATIVE_SHORT},
    {typeid(long),                H5T_NATIVE_LONG},
    {typeid(long long),           H5T_NATIVE_LLONG},
    {typeid(uint8_t),             H5T_NATIVE_UINT8},
    {typeid(unsigned int),        H5T_NATIVE_UINT},
    {typeid(unsigned short),      H5T_NATIVE_USHORT},
    {typeid(unsigned long),       H5T_NATIVE_ULONG},
    {typeid(unsigned long long),  H5T_NATIVE_ULLONG},
    {typeid(float),               H5T_NATIVE_FLOAT},
    {typeid(double),              H5T_NATIVE_DOUBLE}
  };

  hid_t mem_type_id = H5T_NATIVE_INT;

  if ( type_map.count(std::type_index(typeid(T))) )
    mem_type_id = type_map[typeid(T)];
  else
    log::ferror("Undefined memory datatype '%s'", typeid(T).name());

  return mem_type_id;
}


/// \details This method maps typeids of c++ integral types and floating-point
///          types to HDF5 pre-defined standard types.
template <typename T>
hid_t HDF5Handler::getH5FileDatatype() {
  static_assert(std::is_integral<T>::value || std::is_floating_point<T>::value,
                "Only integral datatypes and floating-point datatypes are "
                "pre-defined by HDF5");

  std::unordered_map<std::type_index, hid_t> type_map = {
    {typeid(int8_t),              H5T_STD_I8LE},
    {typeid(int),                 H5T_STD_I32LE},
    {typeid(short),               H5T_STD_I16LE},
    {typeid(long),                H5T_STD_I64LE},
    {typeid(long long),           H5T_STD_I64LE},
    {typeid(uint8_t),             H5T_STD_U8LE},
    {typeid(unsigned int),        H5T_STD_U32LE},
    {typeid(unsigned short),      H5T_STD_U16LE},
    {typeid(unsigned long),       H5T_STD_U64LE},
    {typeid(unsigned long long),  H5T_STD_U64LE},
    {typeid(float),               H5T_IEEE_F32LE},
    {typeid(double),              H5T_IEEE_F64LE}
  };

  hid_t type_id = H5T_STD_I32LE;

  if ( type_map.count(std::type_index(typeid(T))) )
    type_id = type_map[typeid(T)];
  else
    log::ferror("Undefined memory datatype '%s'", typeid(T).name());

  return type_id;
}

//----------------------------------------------------------------------
// ANTMOC utility
//----------------------------------------------------------------------

/// \details The attribute and data space are managed by this method.
///          Integral types will be written as LONG, and floating-point
///          types will be written as DOUBLE.
template <typename T>
void HDF5Handler::writeScalarAttribute(hid_t obj_id,
                                       std::string name,
                                       T value) {

  static_assert(std::is_integral<T>::value || std::is_floating_point<T>::value,
                "Method writeScalarAttribute supports only integral or floating-"
                "point attribute");

  // Set datatypes
  hid_t mem_type_id = getH5MemoryDatatype<T>();
  hid_t type_id     = getH5FileDatatype<T>();

  // Create the data space for the file
  hsize_t dim = 1;
  auto dataspace_id = H5Screate_simple(1, &dim, NULL);

  // Create the attribute
  auto attr_id = createH5Attribute(obj_id, name, type_id, dataspace_id);

  // Write the attribute data
  writeH5Attribute(attr_id, mem_type_id, &value);

  // Close H5 objects
  H5Aclose(attr_id);
  H5Sclose(dataspace_id);
}


template <typename T>
void HDF5Handler::readScalarAttribute(hid_t obj_id,
                                      std::string name,
                                      T &value) {

  static_assert(std::is_integral<T>::value || std::is_floating_point<T>::value,
                "Method readScalarAttribute supports only integral or floating-"
                "point attribute");

  // Set datatypes
  hid_t mem_type_id = getH5MemoryDatatype<T>();

  // Open attribute
  auto attr_id = openH5Attribute(obj_id, name);

  // Read a scalar
  readH5Attribute(attr_id, mem_type_id, &value);

  // Close H5 objects
  H5Aclose(attr_id);
}


} // namespace antmoc

#endif  // HDF5_HANDLER_H_
