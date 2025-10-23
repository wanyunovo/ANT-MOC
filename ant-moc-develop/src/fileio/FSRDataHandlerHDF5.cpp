#include "antmoc/FSRDataHandlerHDF5.h"
#include "antmoc/Geometry.h"
#include "antmoc/Material.h"
#include "antmoc/PyVector.h"
#include "antmoc/Solver.h"
#include "antmoc/TrackGenerator3D.h"

#include <sstream>

namespace antmoc
{


/// \details The file id is managed by constructor/destructor
FSRDataHandlerHDF5::FSRDataHandlerHDF5(std::string file, HDF5Mode mode,
                                       SolverPtr solver,
                                       bool is_parallel)
  : HDF5Handler(file, mode, is_parallel, true), FSRDataHandler(solver)
{ }


//----------------------------------------------------------------------
// ANTMOC utility
//----------------------------------------------------------------------

/// \brief Write metadata of the solver to the HDF5 file
void FSRDataHandlerHDF5::writeMetaData() {
  // Write the number of groups
  auto num_groups = getNumEnergyGroups();
  writeScalarAttribute(_file_id, "# groups", num_groups);

  // Write the total number of FSRs
  hsize_t num_fsrs = getTotalCountForFSRs(getNumFSRs());
  writeScalarAttribute(_file_id, "# fsrs", num_fsrs);
}


/// \details The data is assumed to be distributed over comm_unique_domains.
///          Datasets from different processes will be given different offsets.
void FSRDataHandlerHDF5::writeDataSetForFSRs(hid_t loc_id,
                                             std::string dataset_path,
                                             FP_PRECISION *data,
                                             hsize_t count) {

  if (!data)
    log::error("In writeDataSetForFSRs(...), null pointer is passed in");

  if (count < 1) {
    log::verbose("In writeDataSetForFSRs(...), nothing to write: count == 0");
    return;
  }

  // Create memory dataspace
  auto memspace_id = H5Screate_simple(1, &count, NULL);

  // Create file dataspace and map local elements into the shared file
  auto total_count  = getTotalCountForFSRs(count);
  auto filespace_id = H5Screate_simple(1, &total_count, NULL);
  auto offset       = getFirstOffsetForFSRs(count);
  H5Sselect_hyperslab(filespace_id, H5S_SELECT_SET, &offset,
                      NULL, &count, NULL);

  // Create the dataset
  auto file_datatype = getH5FileDatatype<FP_PRECISION>();
  auto dataset_id = createH5Dataset(loc_id, dataset_path,
                                    file_datatype, filespace_id);

  // Write the dataset
  auto mem_datatype = getH5MemoryDatatype<FP_PRECISION>();
  writeH5Dataset(dataset_id, mem_datatype,
                 memspace_id, filespace_id, data);

  // Close H5 objects
  H5Sclose(memspace_id);
  H5Sclose(filespace_id);
  H5Dclose(dataset_id);

  log::debug("Writing {} values at offset {} of dataset '{}'",
             count, offset, dataset_path);
}


/// \details The data is assumed to be distributed over comm_unique_domains.
///          Datasets from different processes will be given different offsets.
///          Note that the order of FSRs could change between different runs,
///          So that we cannot tell if the IDs of FSRs in memory are equivalent
///          to the IDs of FSRs in the file.
std::vector<FP_PRECISION>
FSRDataHandlerHDF5::
readDataSetForFSRs(hid_t loc_id, std::string dataset_path, hsize_t count) {

  // Create memory dataspace
  auto memspace_id = H5Screate_simple(1, &count, NULL);

  // Create file dataspace and map local elements into the shared file
  auto total_count  = getTotalCountForFSRs(count);
  auto filespace_id = H5Screate_simple(1, &total_count, NULL);
  auto offset       = getFirstOffsetForFSRs(count);
  H5Sselect_hyperslab(filespace_id, H5S_SELECT_SET, &offset,
                      NULL, &count, NULL);

  // Open the dataset
  auto dataset_id = openH5Dataset(loc_id, dataset_path);

  // Read the dataset to a vector
  auto mem_datatype = getH5MemoryDatatype<FP_PRECISION>();
  auto buf = new FP_PRECISION[count]();
  readH5Dataset(dataset_id, mem_datatype,
                memspace_id, filespace_id, (void*)buf);

  // Copy into a vector
  std::vector<FP_PRECISION> data(buf, buf + count);

  // Close H5 objects
  H5Sclose(memspace_id);
  H5Sclose(filespace_id);
  H5Dclose(dataset_id);

  return data;
}


/// \details The hierarchy of the data file is as following
///          group      /
///          group      /FSR
///          group      /FSR/Points
///          dataset    /FSR/Points/X
///          dataset    /FSR/Points/Y
///          dataset    /FSR/Points/Z
///          group      /FSR/Centroids
///          dataset    /FSR/Centroids/X
///          dataset    /FSR/Centroids/Y
///          dataset    /FSR/Centroids/Z
///          dataset    /FSR/Volumes
///          group      /FSR/Fission RX
///          dataset    /FSR/Fission RX/g1
///          dataset    /FSR/Fission RX/g2
///          dataset    /FSR/Fission RX/sum
///          group      /FSR/Scalar Flux
///          dataset    /FSR/Scalar Flux/g1
///          dataset    /FSR/Scalar Flux/g2
///          dataset    /FSR/Scalar Flux/sum
///          ...
void FSRDataHandlerHDF5::dumpFSRData(std::set<TallyType> types,
                                     std::set<int> energy_groups) {

#ifdef ENABLE_MPI_
  // Only processes in the main communicator can participate in
  // writing the FSR data file
  if (!mpi::isInMainCommUniqueDomains()) {
    return;
  }
#endif

  log::result("Dumping FSR reaction rates, xs, and fluxes...");

  // Default
  if (energy_groups.empty())
    energy_groups = getEnergyGroupSet();

  // Always print aggregated data
  energy_groups.insert(0);

  // Create the top-level group
  std::string domain_path = "/FSR";
  hid_t domain_id = createH5Group(_file_id, domain_path);

  // Write characteristic points and centroids of FSRs
  auto fsr_points = getFSRPoints();
  writeFSRPoints(domain_id, "Points", fsr_points);
  fsr_points = getFSRCentroids();
  writeFSRPoints(domain_id, "Centroids", fsr_points);
  fsr_points.clear();

  // Write FSR volumes
  writeFSRVolumes(domain_id);

  // Write all of the field arrays
  writeFieldArrays(domain_id, types, energy_groups);

  // Close H5 objects
  H5Gclose(domain_id);

  log::result("Finished writing file: {}", _file);
}


/// \details Selected groups and datasets will be created to store the data.
///          The hierarchy is shown as dumpFSRData(...).
void FSRDataHandlerHDF5::writeFieldArrays(hid_t domain_id,
                                          std::set<TallyType> types,
                                          std::set<int> energy_groups) {

  // Write all of the reaction rates
  auto num_groups = getNumEnergyGroups();
  PyVector<int> skipped_groups;
  for (auto type : types) {
    auto type_name = tallyutils::getTallyTypeName(type);

    // Create a group for the data
    auto group_id = createH5Group(domain_id, type_name);

    log::debug("Writing data into group '{}'", type_name);

    for (auto group : energy_groups) {
      if (group <= num_groups) {
        // Write a single field as a dataset
        writeSingleFieldArray(group_id, type, group);
      } else {
        skipped_groups.push_back(group);
      }
    }

    // Create the H5 group
    H5Gclose(group_id);
  }

  if (!skipped_groups.empty()) {
    log::verbose("The number of invalid energy groups is {}, which have "
                 "skipped during FSR data dumping: {}",
                 skipped_groups.size(),
                 skipped_groups.toString());
  }
}


std::unordered_map<std::string, std::vector<FP_PRECISION>>
FSRDataHandlerHDF5::readFieldArrays(hid_t domain_id,
                                     std::set<TallyType> types,
                                     std::set<int> energy_groups) {

  std::unordered_map<std::string, std::vector<FP_PRECISION>> fields;

  // Write all of the reaction rates
  auto num_groups = getNumEnergyGroups();
  PyVector<int> skipped_groups;
  for (auto type : types) {
    auto type_name = tallyutils::getTallyTypeName(type);

    // Create a group for the data
    auto group_id = openH5Group(domain_id, type_name);

    for (auto group : energy_groups) {
      if (group <= num_groups) {

        // Prepare dataset path
        std::string array_name = type_name + "/" +
                                 (group == 0 ? "sum" :
                                               "g" + std::to_string(group));

        // Read a single field and store it in the map
        fields.insert({array_name, readSingleFieldArray(group_id, type, group)});
      } else {
        skipped_groups.push_back(group);
      }
    }

    // Create the H5 group
    H5Gclose(group_id);
  }

  if (!skipped_groups.empty()) {
    log::verbose("The number of invalid energy groups is {}, which have "
                 "skipped during FSR data dumping: {}",
                 skipped_groups.size(),
                 skipped_groups.toString());
  }

  return fields;
}


/// \details A dataset will be created to store the array. Path to the
///          dataset is relative to the given group and is determined
///          by the method itself.
void FSRDataHandlerHDF5::writeSingleFieldArray(hid_t group_id,
                                               TallyType type,
                                               const int group) {
  auto num_groups = getNumEnergyGroups();
  if (group > num_groups)
    log::error("Group id {} is out of range", group);

  // Prepare dataset path
  std::ostringstream buf;
  if (group == 0)
    buf << "sum";
  else
    buf << "g" << group;

  std::string dataset_path = buf.str();

  // Prepare data and count
  auto fsr_data_array = getFSRDataArray(type, group);
  hsize_t num_fsrs = getNumFSRs();

  // Write data to the dataset
  writeDataSetForFSRs(group_id, dataset_path, fsr_data_array.data(), num_fsrs);
}


std::vector<FP_PRECISION>
FSRDataHandlerHDF5::
readSingleFieldArray(hid_t group_id, TallyType type, const int group) {

  auto num_groups = getNumEnergyGroups();
  if (group > num_groups)
    log::error("Group id {} is out of range", group);

  // Prepare dataset path
  std::ostringstream buf;
  if (group == 0)
    buf << "sum";
  else
    buf << "g" << group;

  std::string dataset_path = buf.str();

  // Prepare data and count
  hsize_t num_fsrs = getNumFSRs();

  // Write data to the dataset
  return readDataSetForFSRs(group_id, dataset_path, num_fsrs);
}


void FSRDataHandlerHDF5::writeFSRPoints(hid_t loc_id,
                                        std::string dataset_path,
                                        std::vector<std::vector<FP_PRECISION>> &fsr_points) {
  // Create a group for centroids
  hid_t group_id = createH5Group(loc_id, dataset_path);

  hsize_t num_fsrs = getNumFSRs();

  // Write x-coordinates
  writeDataSetForFSRs(group_id, "X",
                      fsr_points[0].data(), num_fsrs);

  // Write y-coordinates
  writeDataSetForFSRs(group_id, "Y",
                      fsr_points[1].data(), num_fsrs);

  // Write z-coordinates
  writeDataSetForFSRs(group_id, "Z",
                      fsr_points[2].data(), num_fsrs);

  H5Gclose(group_id);
}


/// \brief Write volumes of FSRs to the H5 file
/// \param loc_id parent object of the centroid dataset
void FSRDataHandlerHDF5::writeFSRVolumes(hid_t loc_id) {
  hsize_t num_fsrs = getNumFSRs();
  auto fsr_volumes = getFSRVolumes();
  writeDataSetForFSRs(loc_id, "Volumes", fsr_volumes, num_fsrs);
}


/// \brief Get the number of global elements
hsize_t FSRDataHandlerHDF5::getTotalCountForFSRs(hsize_t count) {
  hsize_t total_count = count;

#ifdef ENABLE_MPI_
  if (mpi::isSpatialDecomposed()) {
    MPI_Comm comm = getCommForFSRs();
    MPI_Allreduce(&count, &total_count, 1, mpi::getDatatype<long>(),
                  MPI_SUM, comm);
  }
#endif

  return total_count;
}


/// \brief Get the global offset of the first local elements
hsize_t FSRDataHandlerHDF5::getFirstOffsetForFSRs(hsize_t count) {
  hsize_t offset = 0;

#ifdef ENABLE_MPI_
  if (mpi::isSpatialDecomposed()) {
    int num_procs = getNumProcsForFSRs();
    MPI_Comm comm = getCommForFSRs();

    hsize_t *all_counts = new hsize_t[num_procs]();

    MPI_Allgather(&count, 1, mpi::getDatatype<long>(),
                  all_counts, 1, mpi::getDatatype<long>(),
                  comm);

    int my_rank = getMPIRankForFSRs();
    for (int r = 0; r < my_rank; ++r)
      offset += all_counts[r];
  }
#endif

  return offset;
}


/// \brief Extract the reaction rate of a certain type
/// \param type type of the tallied reaction rate
/// \param group a certain energy group, and "-1" means the total rates
std::vector<FP_PRECISION>
FSRDataHandlerHDF5::getFSRDataArray(TallyType type, const int group) {

  long num_fsrs = getNumFSRs();
  auto num_groups = getNumEnergyGroups();

  // Loop over all flat source regions
  std::vector<FP_PRECISION> fsr_data_array;
  fsr_data_array.reserve(num_fsrs);
  for (long r = 0; r < num_fsrs; ++r) {

    // Determine the volume and cross-sections of the FSR
    FP_PRECISION volume = getFSRVolume(r);

    // Define the lambda for computing values in different cases
    auto compute_value = [&](int g) {
      double xs = getFSRXS(r, g, type);
      double value = 0.;
      if (tallyutils::isRXTallyType(type)) {
        value = getFSRFlux(r, g) * xs;
      } else if (tallyutils::isXSTallyType(type)) {
        value = xs;
      } else if (type == +TallyType::Volume) {
        value = volume;
      }
      return value;
    };

    FP_PRECISION value = 0.;

    // Tally the mesh data summed across all groups
    if (group == 0) {
      for (int g = 1; g <= num_groups; g++)
        value += compute_value(g);
    }
    else { // A single group
      value = compute_value(group);
    }

    fsr_data_array.push_back(value);

  } // End loop for FSRs

  return fsr_data_array;
}


} // namespace antmoc
