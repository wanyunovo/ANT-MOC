/// \file FSRDataHandlerHDF5.h
/// \brief An I/O class for FSR data
/// \date March 6, 2020
/// \author An Wang, USTB (wangan@xs.ustb.edu.cn)

#ifndef FSRDATA_HANDLER_HDF5_H_
#define FSRDATA_HANDLER_HDF5_H_


#include <memory>
#include <set>
#include <string>
#include <unordered_map>
#include <vector>

#include "antmoc/HDF5Handler.h"
#include "antmoc/FSRDataHandler.h"


namespace antmoc
{

///---------------------------------------------------------------------
/// \class FSRDataHandlerHDF5
/// \brief Methods for writing and reading FSR data in HDF5 format
///---------------------------------------------------------------------
class FSRDataHandlerHDF5 : public HDF5Handler, public FSRDataHandler {

public:

  /// \brief Construct an object of FSRDataHandlerHDF5
  FSRDataHandlerHDF5(std::string file, HDF5Mode mode, SolverPtr solver,
                     bool is_parallel = true);
  ~FSRDataHandlerHDF5() = default;

  //--------------------------------------------------------------------
  // ANTMOC utility
  //--------------------------------------------------------------------
  /// \brief Dump reaction rates, cross-sections or volumes. Data will
  ///        written into a group named '/FSR'.
  /// \param types Types of the tallied reaction rates.
  /// \param energy_groups Selected energy groups (defaults to all).
  void dumpFSRData(std::set<TallyType> types,
                   std::set<int> energy_groups = std::set<int>());

  /// \brief A helper method for dumpFSRData(...).
  void writeFieldArrays(hid_t domain_id,
                        std::set<TallyType> types,
                        std::set<int> energy_groups = std::set<int>());

  /// \brief Read all of the field arrays from file.
  /// \return A map from dataset names to data.
  std::unordered_map<std::string, std::vector<FP_PRECISION>>
  readFieldArrays(hid_t domain_id,
                  std::set<TallyType> types,
                  std::set<int> energy_groups);

  /// \brief Write the reaction rate of a specific type and a certain group.
  /// \param group_id ID of the H5 group.
  /// \param type Type of the tallied reaction rate
  /// \param group A certain energy group, and "0" means the total rates.
  void writeSingleFieldArray(hid_t group_id, TallyType type, const int group);

  /// \brief Read the reaction rate of a specific type and a certain group.
  std::vector<FP_PRECISION> readSingleFieldArray(hid_t group_id,
                                                 TallyType type,
                                                 const int group);

  /// \brief Write characteristic points or centroids of FSRs to the H5 file.
  /// \param loc_id Parent object of the point group.
  /// \param group_path Relative path to the point group.
  /// \param fsr_points A reference to the data.
  void writeFSRPoints(hid_t loc_id,
                      std::string dataset_path,
                      std::vector<std::vector<FP_PRECISION>> &fsr_points);

  void writeFSRVolumes(hid_t loc_id);
  void writeMetaData();

  /// \brief Write a simple dataset to a specific place.
  /// \param loc_id Parent object of the dataset.
  /// \param dataset_path Path to the H5 dataset.
  /// \param data A pointer to the data.
  /// \param count Number of elements.
  void writeDataSetForFSRs(hid_t loc_id,
                           std::string dataset_path,
                           FP_PRECISION *data,
                           hsize_t count);

  /// \brief Read a simple dataset from a file or group
  /// \param loc_id Parent object of the dataset.
  /// \param dataset_path Path to the H5 dataset.
  /// \param count Number of elements.
  /// \return Data vector.
  std::vector<FP_PRECISION> readDataSetForFSRs(hid_t loc_id,
                                               std::string dataset_path,
                                               hsize_t count);

  hsize_t getTotalCountForFSRs(hsize_t);
  hsize_t getFirstOffsetForFSRs(hsize_t);
  std::vector<FP_PRECISION> getFSRDataArray(TallyType, const int);

};


} // namespace antmoc

#endif  // FSRDATA_HANDLER_HDF5_H_
