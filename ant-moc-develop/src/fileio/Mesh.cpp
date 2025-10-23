#include "antmoc/Mesh.h"
#include "antmoc/FSRDataTypes.h"
#include "antmoc/Geometry.h"
#include "antmoc/Lattice.h"
#include "antmoc/Material.h"
#include "antmoc/PyVector.h"
#include "antmoc/Solver.h"
#include "antmoc/TrackGenerator3D.h"
#include "antmoc/Universe.h"

#include <iomanip>
#include <ios>
#include <sstream>
#include <cstdlib>


namespace antmoc
{

/// \brief The Mesh constructor.
/// \details If no lattice is given, a default lattice can be constructed with
///          the Mesh::createLattice function.
/// \param solver The solver from which scalar fluxes and cross-sections are
///        extracted
/// \param lattice An optional parameter for the lattice across which reaction
///        rates are tallied
Mesh::Mesh(SolverPtr solver, Lattice* lattice):
  FSRDataHandler(solver),
  _lattice(lattice),
  _lattice_allocated(false),
  _num_outsides(0)
  { }


/// \brief The Mesh destructor deletes its lattice if the lattice was allocated
///        internally.
Mesh::~Mesh() {
  delete _lattice;
  _lattice = nullptr;

  _cur_elements.clear();

  for (auto e : _docs)
    delete e.second;
}


/// \brief Creates an internal lattice over which to tally reaction rates with
///        the user-input dimensions.
/// \param num_x the number of mesh cells in the x-direction
/// \param num_y the number of mesh cells in the y-direction
/// \param num_z the number of mesh cells in the z-direction
void Mesh::createLattice(int num_x, int num_y, int num_z) {

  /* Delete the current lattice if currently allocated */
  if (_lattice_allocated)
    delete _lattice;

  /* Allocate a new lattice */
  auto lattice = new RecLattice();
  lattice->setNumX(num_x);
  lattice->setNumY(num_y);
  lattice->setNumZ(num_z);

  /* Get the root universe */
  Geometry* geometry = _solver->getGeometry();
  Universe* root_universe = geometry->getRootUniverse();

  /* Determine the geometry widths in each direction */
  double width_x = (root_universe->getMaxX() - root_universe->getMinX())/num_x;
  double width_y = (root_universe->getMaxY() - root_universe->getMinY())/num_y;
  double width_z = (root_universe->getMaxZ() - root_universe->getMinZ())/num_z;

  /* Determine the center-point of the geometry */
  double offset_x = (root_universe->getMinX() + root_universe->getMaxX()) / 2;
  double offset_y = (root_universe->getMinY() + root_universe->getMaxY()) / 2;
  double offset_z = (root_universe->getMinZ() + root_universe->getMaxZ()) / 2;

  /* Create the Mesh lattice */
  lattice->setWidth(width_x, width_y, width_z);
  lattice->setOffset(offset_x, offset_y, offset_z);
  lattice->computeSizes();

  _lattice = lattice;
  _lattice_allocated = true;
}


/// \brief Creates an internal lattice over which to tally reaction rates with
///        the user-input dimensions.
/// \param widths_x the widths of mesh cells in the x-direction
/// \param widths_y the widths of mesh cells in the y-direction
/// \param widths_z the widths of mesh cells in the z-direction
void Mesh::createLattice(const PyVector<double> &widths_x,
                         const PyVector<double> &widths_y,
                         const PyVector<double> &widths_z) {

  /* Delete the current lattice if currently allocated */
  if (_lattice_allocated)
    delete _lattice;

  /* Allocate a new lattice */
  auto lattice = new RecLattice();
  lattice->setNumX(widths_x.size());
  lattice->setNumY(widths_y.size());
  lattice->setNumZ(widths_z.size());
  lattice->setWidths(widths_x.get(), widths_y.get(), widths_z.get());

  /* Get the root universe */
  Geometry* geometry = _solver->getGeometry();
  Universe* root_universe = geometry->getRootUniverse();

  /* Determine the center-point of the geometry */
  double offset_x = (root_universe->getMinX() + root_universe->getMaxX()) / 2;
  double offset_y = (root_universe->getMinY() + root_universe->getMaxY()) / 2;
  double offset_z = (root_universe->getMinZ() + root_universe->getMaxZ()) / 2;

  /* Create the Mesh lattice */
  lattice->setOffset(offset_x, offset_y, offset_z);
  lattice->computeSizes();

  _lattice = lattice;
  _lattice_allocated = true;
}


/// \brief Creates an internal hexlattice over which to tally reaction rates
///         with the user-input dimensions.
/// \param num_r the number of radial lattice cells
/// \param width_r the radial pitch
/// \param widths_z the axial pitches
/// \param orientation hexlattice orientation (x or y)
/// \param offset the offset with respect to the center of the geometry
void Mesh::createLattice(int num_r, double width_r,
                         const PyVector<double> &widths_z,
                         std::string orientation,
                         const Point &offset) {

  // Delete the current lattice if currently allocated
  if (_lattice_allocated)
    delete _lattice;

  // Allocate a new lattice
  auto lattice = new HexLattice();
  lattice->setNumR(num_r);
  lattice->setNumZ(widths_z.size());
  lattice->setWidths(width_r, widths_z.get());
  lattice->setOrientation(orientation);

  // Get the root universe
  Geometry* geometry = _solver->getGeometry();
  Universe* root_universe = geometry->getRootUniverse();

  // Determine the center-point of the geometry
  double offset_x = (root_universe->getMinX() + root_universe->getMaxX()) / 2;
  double offset_y = (root_universe->getMinY() + root_universe->getMaxY()) / 2;
  double offset_z = (root_universe->getMinZ() + root_universe->getMaxZ()) / 2;

  lattice->setOffset(offset_x + offset.getX(),
                     offset_y + offset.getY(),
                     offset_z + offset.getZ());

  lattice->computeSizes();

  _lattice = lattice;
  _lattice_allocated = true;
}


/// \brief Set the _lattice of a mesh to be an existing one, for which the user
///        inputs the dimensions.
/// \param lattice the existing lattice to be set to the Mesh
void Mesh::setLattice(Lattice* lattice) {
  /* Delete the current lattice if currently allocated */
  if (_lattice_allocated)
    delete _lattice;

  _lattice = lattice;
  _lattice_allocated = false;
}


/// \brief Tallies data of the given type over the Mesh lattice
/// \details Mesh data could be reaction rates, cross-sections, volumes, etc.
///          The returned vector contains all of the lattice cells, including
///          the invalid ones in HexLattices.
/// \param tally The type of mesh data to tally
/// \param volume_average whether the data should be volume averaged
/// \return The mesh data in a 1D vector indexed by the lattice cell IDs
std::vector<FP_PRECISION> Mesh::getMeshData(TallyType tally,
                                            bool volume_average,
                                            const int group) {

  /* Check that the Mesh contains a lattice */
  if (_lattice == NULL)
    log::error("A Lattice must be set or created to get mesh data "
               "form a Mesh object");

  /* Create a 1D array of mesh data with the appropriate size */
  std::vector<FP_PRECISION> mesh_data;
  std::vector<FP_PRECISION> volumes_lattice;
  int size = _lattice->getMaxNumUniverses();  // includes invalid ones
  mesh_data.resize(size, 0.);
  if (volume_average)
    volumes_lattice.resize(size, 0.);

  /* Extract the number of groups */
  int num_groups = getNumEnergyGroups();
  if (group > num_groups)
    log::error("Group id {} is out of range", group);

  // Reset the number of outside points
  resetNumOutsides(0);

  /* Loop over all flat source regions */
  long num_fsrs = getNumFSRs();
  for (long r = 0; r < num_fsrs; ++r) {

    /* Determine the FSR material and which Mesh cell contains it */
    Point* pt = getFSRPoint(r);

    // If the point is outside the lattice, skips it and warns the user.
    // It is necessary to skip these points because we need only the
    // points inside the lattice.
    if (!_lattice->containsPoint(pt)) {
      increaseNumOutsides();
      continue;
    }

    int lat_cell = _lattice->getLatticeCell(pt);

    /* Determine the volume and cross-sections of the FSR */
    FP_PRECISION volume = getFSRVolume(r);

    // Compute values in different cases
    auto compute_value = [&](int g) {
      double xs = getFSRXS(r, g, tally);
      double value = 0.;
      if (tallyutils::isRXTallyType(tally)) {
        value = getFSRFlux(r, g) * volume * xs;
      } else if (tallyutils::isXSTallyType(tally)) {
        value = volume * xs;
      } else if (tally == +TallyType::Volume) {
        value = volume;
      }
      return value;
    };

    /* Tally the mesh data summed across all groups */
    double fsr_mesh_data = 0.;
    if (group == 0) {
      for (int g = 1; g <= num_groups; g++) {
        // Compute the value of current mesh cell
        auto value = compute_value(g);
        mesh_data.at(lat_cell) += value;
        fsr_mesh_data += value;
      }
    }
    else { // A single group
      // Compute the value of current mesh cell
      auto value = compute_value(group);
      mesh_data.at(lat_cell) += value;
      fsr_mesh_data += value;
    }

    /* Tally fsr volume to cell volume, only for a non-zero reaction rate */
    if (std::abs(fsr_mesh_data) > FLT_EPSILON && volume_average)
      volumes_lattice.at(lat_cell) += volume;
  }

#ifdef ENABLE_MPI_
  if (mpi::isSpatialDecomposed()) {

    // If domain decomposed, do a reduction for mesh data
    MPI_Comm comm = mpi::getCommUniqueDomains();
    FP_PRECISION* mesh_data_array = mesh_data.data();
    FP_PRECISION* mesh_data_send = new FP_PRECISION[size]();
    memcpy(mesh_data_send, mesh_data_array,
           size * sizeof(FP_PRECISION));

    MPI_Allreduce(mesh_data_send,
                  mesh_data_array,
                  size,
                  mpi::getDatatype<FP_PRECISION>(),
                  MPI_SUM,
                  comm);
    delete [] mesh_data_send;

    log::debug("Mesh data array reduced");

    // If domain decomposed, do a reduction for mesh volumes
    if (volume_average) {
      FP_PRECISION* volume_data_array = volumes_lattice.data();
      FP_PRECISION* volume_data_send = new FP_PRECISION[size]();
      memcpy(volume_data_send, volume_data_array,
             size * sizeof(FP_PRECISION));

      MPI_Allreduce(volume_data_send,
                    volume_data_array,
                    size,
                    mpi::getDatatype<FP_PRECISION>(),
                    MPI_SUM,
                    comm);
      delete [] volume_data_send;

      log::debug("FSR volumes array reduced");
    }
  }
#endif

  /* If volume average requested, divide by volume */
  if (volume_average)
    for (size_t i = 0; i < mesh_data.size(); ++i)
      if (volumes_lattice.at(i) > FLT_EPSILON)
        mesh_data.at(i) /= volumes_lattice.at(i);

  return mesh_data;
}


/// \brief Tallies data of the given type over the Mesh lattice
/// \param type The type of reaction to tally
/// \param volume_average whether the mesh data should be volume averaged
/// \return The mesh data in a 3D vector indexed by the lattice cell
///         x, y, and z indexes
Vector3D Mesh::getFormattedMeshData(TallyType type, bool volume_average,
                                    const int group) {

  /* Extract mesh data*/
  Vector3D mesh_data;
  auto mesh_data_array = getMeshData(type, volume_average, group);

  /* Format mesh data into a 3D array */
  // FIXME, only for RecLattice
  auto lattice = dynamic_cast<RecLattice*>(_lattice);
  int num_x = lattice->getNumX();
  for (int i=0; i < num_x; i++) {
    int num_y = lattice->getNumY();
    std::vector<std::vector<FP_PRECISION> > vector_2D;
    for (int j=0; j < num_y; j++) {
      int num_z = lattice->getNumZ();
      std::vector<FP_PRECISION> vector_1D;
      for (int k=0; k < num_z; k++) {
        int idx = k * num_x * num_y + j * num_x + i;
        vector_1D.push_back(mesh_data_array[idx]);
      }
      vector_2D.push_back(vector_1D);
    }
    mesh_data.push_back(vector_2D);
  }
  return mesh_data;
}


//----------------------------------------------------------------------
// Utility
//----------------------------------------------------------------------
/// \brief Write a named XML document to a file
/// \details If the file path is not specified, data will be print to
///          standard output. When the dumping is finished, the XML
///          document will be deleted.
void Mesh::writeXMLDocToFile(std::string doc_name, std::string file) {
  if (_docs.find(doc_name) == _docs.end()) {
    log::error("Trying to dump an nonexistent document: {}", doc_name);
  }

  auto doc = _docs.at(doc_name);
  if (file.empty()) {
    doc->Print();
  }
  else {
    FILE *fp = fopen(file.c_str(), "w");
    XMLPrinter printer(fp);
    doc->Print(&printer);

    // Inform the user
    log::debug("Partially writing file: {}", file);
    log::result("Finished writing file: {}", file);
  }
  doc->Clear();
  _docs.erase(doc_name);
}


/// \details This method print a document to memory, which introduces extra
///          memory cost.
std::string Mesh::printDocToString(std::string doc_name) {

  if (_docs.find(doc_name) == _docs.end()) {
    log::error("Document {} doesn't exist. printMeshDataToString() "
               "must be called after initializing the document", doc_name);
  }

  auto &doc = _docs.at(doc_name);
  XMLPrinter printer;
  doc->Print(&printer);

  std::string str(printer.CStr());
  printer.ClearBuffer();

  return str;
}


//----------------------------------------------------------------------
// Dumping mesh data
//----------------------------------------------------------------------
/// \details This method first writes all of the data into an XML document
///          in memory, and then writes it into a file.
void Mesh::dumpMeshDataToFile(std::string file, std::set<TallyType> types,
                              std::set<int> energy_groups) {

  log::result("Dumping mesh reaction rates, xs, and fluxes...");

  std::string doc_name = "mesh";
  printMeshDataToXML(doc_name, types, energy_groups);

  // The root will write the file since the data has been reduced
  if (mpi::isMPIRoot())
    writeXMLDocToFile(doc_name, file);
}


/// \details The unstructured grid is used to dump mesh data such as
///          reaction, cross-sections, and volumes for both RecLattices
///          and HexLattices, cross-sections and volumes.
///          Once the initialization successes, a pointer of XMLElement
///          will point to the current data element.
void Mesh::initializeUnstructuredDoc(std::string doc_name) {

  // Creates XML elements
  if (_docs.find(doc_name) != _docs.end())
    _docs.erase(doc_name);
  _docs.insert({doc_name, new XMLDocument});

  auto doc = _docs.at(doc_name);
  auto e_root     = doc->NewElement("VTKFile");
  auto e_grid     = doc->NewElement("UnstructuredGrid");
  auto e_piece    = doc->NewElement("Piece");
  auto e_points   = doc->NewElement("Points");
  auto e_cells    = doc->NewElement("Cells");
  auto e_celldata = doc->NewElement("CellData");
  auto e_da_pt    = doc->NewElement("DataArray"); // Data array of points
  auto e_da_con   = doc->NewElement("DataArray"); // Connectivity for Cells
  auto e_da_off   = doc->NewElement("DataArray"); // Offsets for Cells
  auto e_da_type  = doc->NewElement("DataArray"); // Offsets for Cells
  auto e_da_indx  = doc->NewElement("DataArray"); // Lattice cell indices

  // Creates the hierarchy
  // root
  //   --UnstructuredGrid
  //     --Piece
  //       --Points
  //         --data
  //       --Cells
  //         --connectivity
  //         --offsets
  //         --types
  //       --CellData
  //         --data
  doc->InsertEndChild(e_root);
  e_root->InsertEndChild(e_grid);
  e_grid->InsertEndChild(e_piece);
  e_piece->InsertEndChild(e_points);
  e_piece->InsertEndChild(e_cells);
  e_piece->InsertEndChild(e_celldata);
  e_points->InsertEndChild(e_da_pt);
  e_cells->InsertEndChild(e_da_con);
  e_cells->InsertEndChild(e_da_off);
  e_cells->InsertEndChild(e_da_type);
  e_celldata->InsertEndChild(e_da_indx);

  // Initializes the root element
  e_root->SetAttribute("type", "UnstructuredGrid");
  e_root->SetAttribute("version", "0.1");

  //--------------------------------------------------------------------
  // Sets vertices and indices
  //--------------------------------------------------------------------
  // Computes the number of cells and vertices
  // (only for 3-D geometry
  int max_num_cells = _lattice->getMaxNumUniverses();
  int num_cells = _lattice->getNumLatticeCells();
  int vertices_per_cell;
  int cell_type;

  if (_lattice->getLatticeType() == latticeType::Rectangle) {
    vertices_per_cell = 8;
    cell_type = 12;
  } else {
    vertices_per_cell = 12;
    cell_type = 16;
  }

  // Total number of vertices
  int num_points = num_cells * vertices_per_cell;
  // Number of cells on a z-section
  int num_cells_xy = num_cells / _lattice->getNumZ();

  e_piece->SetAttribute("NumberOfPoints", std::to_string(num_points).c_str());
  e_piece->SetAttribute("NumberOfCells", std::to_string(num_cells).c_str());
  e_piece->SetAttribute("NumberOfCellsXY", std::to_string(num_cells_xy).c_str());

  // Records information about the tally mesh shape
  if (_lattice->getLatticeType() == latticeType::Rectangle)
    e_piece->SetAttribute("type", "Rectangular Mesh");
  else
    e_piece->SetAttribute("type", "Hexagonal Mesh");

  std::ostringstream extent;
  // count invalid cells if its a hexagonal mesh
  extent << _lattice->getNumX() << ' '
         << _lattice->getNumY() << ' '
         << _lattice->getNumZ();
  e_piece->SetAttribute("Extent", extent.str().c_str());

  // String buffers
  std::ostringstream ss_pt, ss_indx;
  ss_pt << '\n' << std::fixed << std::setprecision(6);
  ss_indx << '\n';

  // Loop over lattice cells
  // Note that getMaxNumUniverses() returns the number of all universes,
  // including invalid ones.
  for (int idx = 0; idx < max_num_cells; ++idx) {
    // Skips invalid indices
    if (_lattice->isValidIndex(idx)) {
      ss_indx << ' ' << idx;  // records valid indices
      auto vertices = _lattice->getLatticeCellVertices(idx);
      for (auto &pt : vertices)
        ss_pt << ' ' << pt.getX() << ' ' << pt.getY() << ' ' << pt.getZ();
    }
  }

  // Writes points to the doc
  ss_pt << '\n';
  e_da_pt->SetAttribute("Name", "Vertices");
  e_da_pt->SetAttribute("type", "Float64");
  e_da_pt->SetAttribute("NumberOfComponents", "3");
  e_da_pt->SetAttribute("format", "ascii");
  e_da_pt->SetText(ss_pt.str().c_str());

  // Writes valid indices to the doc
  ss_indx << '\n';
  e_da_indx->SetAttribute("Name", "Valid Indices");
  e_da_indx->SetAttribute("type", "Int32");
  e_da_indx->SetAttribute("format", "ascii");
  e_da_indx->SetText(ss_indx.str().c_str());


  //--------------------------------------------------------------------
  // Sets data related to cells:
  //   connectivity, offsets, and cell types
  //--------------------------------------------------------------------
  // Writes connectivity and offsets
  // Connectivity shows how you connect vertices to form a cell, e.g.
  //   0 1 2 3 4 5 6 7 8  ---- a hexahedron
  //   0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16  ---- a hexagonal prism
  // Offsets record the past-the-end point of cells, e.g.
  //   8 16 24   ---- three hexahedrons
  //   16 32 48  ---- three hexagonal prisms
  // Cell types are integers. Only two of them are used in the following
  // code.
  //   12  ---- hexahedrons
  //   16  ---- hexagonal prism

  std::ostringstream ss_con, ss_off; // Connectivity and offsets
  std::ostringstream ss_type;
  ss_con << '\n';
  ss_off << '\n';
  ss_type << '\n';

  for (int i = 0; i < num_points; ++i) {
    ss_con << ' ' << i;
    if (i % vertices_per_cell == 0) {
      ss_off << ' ' << i + vertices_per_cell;
      ss_type << ' ' << cell_type;
    }
  }

  // Initializes the DataArray of connectivity
  ss_con << '\n';
  e_da_con->SetAttribute("Name", "connectivity");
  e_da_con->SetAttribute("type", "Int32");
  e_da_con->SetAttribute("format", "ascii");
  e_da_con->SetText(ss_con.str().c_str());

  // Initializes the DataArray of offsets
  ss_off << '\n';
  e_da_off->SetAttribute("Name", "offsets");
  e_da_off->SetAttribute("type", "Int32");
  e_da_off->SetAttribute("format", "ascii");
  e_da_off->SetText(ss_off.str().c_str());

  // Initializes the DataArray of cell types
  ss_type << '\n';
  e_da_type->SetAttribute("Name", "types");
  e_da_type->SetAttribute("type", "UInt8");
  e_da_type->SetAttribute("format", "ascii");
  e_da_type->SetText(ss_type.str().c_str());

  // Set current element to <CellData>
  _cur_elements[doc_name] = e_celldata;
}


/// \details This method should be called after _cur_elements is initialized.
void Mesh::printSingleMeshDataToXML(std::string doc_name, std::string data_name,
                              TallyType type,
                              bool volume_average,
                              const int group,
                              unsigned precision) {

  // Gets the data array
  auto mesh_data = getMeshData(type, volume_average, group);

  size_t max_num_cells = _lattice->getMaxNumUniverses();
  if (mesh_data.size() != max_num_cells) {
    log::error("The number of mesh data points ({}) is not equal "
               "to the maximum number of lattice cells ({})",
               mesh_data.size(), max_num_cells);
  }

  // Initializes element CellData
  auto cur_element = _cur_elements.at(doc_name);

  // Creates a data array
  auto doc = _docs.at(doc_name);
  auto e_da = doc->NewElement("DataArray");
  cur_element->InsertEndChild(e_da);

  std::ostringstream buf;
  buf << data_name;
  if (group > 0 &&
      (tallyutils::isRXTallyType(type) ||
       tallyutils::isXSTallyType(type)))
    buf << " g" << group;

  e_da->SetAttribute("Name", buf.str().c_str());
  e_da->SetAttribute("type", "Float64");
  e_da->SetAttribute("format", "ascii");

  // Data string buffer
  std::ostringstream ss;
  ss << '\n';

  // Set the precision of mesh data to
  ss << std::fixed << std::setprecision(precision);

  // VTK files follows the order of x-y-z. This is the very order
  // we stored the data in underlying lattice layout, which means
  // we should start from the lower left corner at z=0.
  int num_x = _lattice->getNumX();
  int num_y = _lattice->getNumY();
  int num_z = _lattice->getNumZ();

  for (int z = 0; z < num_z; ++z) {
    for (int y = 0; y < num_y; ++y) {
      for (int x = 0; x < num_x; ++x) {
        // Skips invalid indices
        if (_lattice->areValidIndices(x, y, z)) {
          auto indx = _lattice->getLatticeCell(x, y, z);
          ss << mesh_data[indx] << ' ';
        }
      }
      ss << '\n'; // line breaks
    }
  }

  //std::cout << ss.str();
  e_da->SetText(ss.str().c_str());
}


/// \details This method initializes a new document object and print mesh
///          data to it. If the document has be created before, it will be
///          erased. If this method succeeds, a new document filled with
///          data will be stored in the document map.
///          Only data tallied on rectangular or hexagonal lattices could
///          be dumped properly by this method. Energy groups must be
///          specified to dump data with selected groups.
void Mesh::printMeshDataToXML(std::string doc_name, std::set<TallyType> types,
                              std::set<int> energy_groups) {
  // Default
  if (energy_groups.empty())
    energy_groups = getEnergyGroupSet();

  // Always print aggregated data
  energy_groups.insert(0);

  // Initializes the XML document
  initializeUnstructuredDoc(doc_name);

  // Print volumes to XML doc
  printSingleMeshDataToXML(doc_name, "Volume", TallyType::Volume, false, 1);

  auto num_groups = getNumEnergyGroups();
  PyVector<int> skipped_groups;
  for (auto type : types) {
    auto type_name = tallyutils::getTallyTypeName(type);

    for (auto group : energy_groups) {
      if (group <= num_groups) {
        printSingleMeshDataToXML(doc_name, type_name, type, false, group);
        // Volume-averaged
        printSingleMeshDataToXML(doc_name, "Avg " + type_name, type, true, group);
      } else {
        skipped_groups.push_back(group);
      }
    }
  }

  // Reports the number of points outside the mesh
  if (getNumOutsides() > 0) {
    log::verbose_once("There are {} points outside the tally mesh, which have been "
                      "skipped during tallying reaction rates. You are safe to "
                      "ignore this message if this is the expected behaviour.",
                       getNumOutsides());
  }

  // Reports skipped energy groups
  if (!skipped_groups.empty()) {
    log::verbose("The number of invalid energy groups is {}, which have "
                 "skipped during mesh data dumping: {}",
                 skipped_groups.size(), skipped_groups.toString());
  }
}


//----------------------------------------------------------------------
// Dumping tracks
//----------------------------------------------------------------------
/// \brief Dump track mesh to files
/// \param tg TrackGenerator which holds tracks
/// \param types types of tallied mesh data (All, None, Tracks_2D or Tracks_3D)
/// \param dir directory of output files
/// \param prefix prefix of output files
void Mesh::dumpTrackMeshToDir(TrackGeneratorPtr tg, std::set<TallyType> types,
                              std::string dir) {

#ifdef ENABLE_MPI_
  mpi::mpiBarrier();
#endif

  if (dir.empty()) {
    log::error("Dir cannot be empty when dumping tracks");
  }

  // Expand TallyType::None and TallyType::All
  tallyutils::expandTallyTypeForTrack(types);
  if (types.empty()) {
    return;
  }

  log::result("Dumping tracks...");

  auto r = mpi::getMPIRank();          // MPI rank

  for (auto t : types) {
    if (!tallyutils::isTracksTallyType(t)) {
      log::error("Only tracks should be dumped to track dumping files");
    }

    auto type_name = tallyutils::getTallyTypeName(t); // Get the name of TallyType for identification
    std::string file = dir + "/" + type_name;  // Common part of the path

    // Write a main file for poly data
    if (mpi::isMPIRoot()) {
      std::string pdoc_name = "PPoly " + type_name; // Initializes the XML document
      initializePPolyDataDoc(pdoc_name, type_name); // An XML document holding PPolyData
      writeXMLDocToFile(pdoc_name, file + ".pvtp");
    }
    // Write pieces
    std::string doc_name = "Piece " + type_name;    // Initializes the XML document
    printTracksMeshToXML(tg, doc_name, type_name, t);
    writeXMLDocToFile(doc_name, getPieceName(file, r));
  }

#ifdef ENABLE_MPI_
  mpi::mpiBarrier();
#endif
}


/// \brief Creates an XML document for the ppoly data
/// \details This is a helpler method for dumping track mesh.
///          Output files are named automatically according to a given file
///          name and ids of processes.
/// \param doc_name the name of the XML document
/// \param prefix the prefix of dumping files
void Mesh::initializePPolyDataDoc(std::string doc_name, std::string prefix) {

  // Creates XML elements
  if (_docs.find(doc_name) != _docs.end())
    _docs.erase(doc_name);
  _docs.insert({doc_name, new XMLDocument});

  auto doc = _docs.at(doc_name);
  auto e_root = doc->NewElement("VTKFile");
  auto e_ppoly = doc->NewElement("PPolyData");
  auto e_ppoints = doc->NewElement("PPoints");
  auto e_pda_ppt = doc->NewElement("PDataArray");
  auto e_ids = doc->NewElement("PCellData");
  auto e_pda_ids = doc->NewElement("PDataArray");

  // Creates the hierarchy
  doc->InsertEndChild(e_root);
  e_root->InsertEndChild(e_ppoly);
  e_ppoly->InsertEndChild(e_ppoints);
  e_ppoints->InsertEndChild(e_pda_ppt);
  e_ppoly->InsertEndChild(e_ids);
  e_ids->InsertEndChild(e_pda_ids);

  // Initializes the root element
  e_root->SetAttribute("type", "PPolyData");
  e_root->SetAttribute("version", "0.1");
  e_root->SetAttribute("byte_order", "LittleEndian");

  // Initializes the PPolyData element
  e_ppoly->SetAttribute("GhostLevel", "0");

  // Initializes the DataArray of PPoints element
  e_pda_ppt->SetAttribute("type", "Float64");
  e_pda_ppt->SetAttribute("NumberOfComponents", "3");

  // Initializes elements for IDs
  e_ids->SetAttribute("Scalars", "track_ids");
  e_pda_ids->SetAttribute("type", "Int64");
  e_pda_ids->SetAttribute("Name", "track_ids");

  // Initializes the Piece elements
  auto num_pieces = mpi::getNumProcs();
  for (int r = 0; r < num_pieces; ++r) {
    auto e_piece = doc->NewElement("Piece");
    e_piece->SetAttribute("Source", getPieceName(prefix, r).c_str());
    e_ppoly->InsertEndChild(e_piece);
  }
}


/// \brief Creates and dump an XML document to files
/// \details Points and connectivities are dumped to files under a given
///          directory. Each file is a serial VTK polydata file containing
///          a piece of tracks.
void Mesh::printTracksMeshToXML(TrackGeneratorPtr tg, std::string doc_name,
                                std::string data_name, TallyType t) {

  // Get array data from track generator
  std::vector<TrackMeshData> track_coords;
  switch (t) {
    case TallyType::Tracks_2D:
      track_coords = tg->TrackGenerator::getFormattedTracks();
      break;
    case TallyType::Tracks_3D:
      if (std::dynamic_pointer_cast<TrackGenerator3D>(tg)) {
        track_coords = tg->getFormattedTracks();
      }
      break;

    default:
      log::error("Unrecognized type for dumping tracks");
  }

  // Creates XML elements
  if (_docs.find(doc_name) == _docs.end())
    _docs.insert({doc_name, new XMLDocument});

  auto doc = _docs.at(doc_name);
  auto e_root = doc->NewElement("VTKFile");
  auto e_poly = doc->NewElement("PolyData");
  auto e_piece = doc->NewElement("Piece");
  auto e_points = doc->NewElement("Points");
  auto e_da_p = doc->NewElement("DataArray");  // Point data
  auto e_lines = doc->NewElement("Lines");
  auto e_da_con = doc->NewElement("DataArray");  // Connectivity
  auto e_da_off = doc->NewElement("DataArray");  // Offsets
  auto e_ids = doc->NewElement("CellData");
  auto e_da_ids = doc->NewElement("DataArray");  // Track IDs

  // Creates the hierarchy
  // root
  //   --polydata
  //     --piece
  //       --points
  //         --data
  //       --lines
  //         --connectivity
  //         --offsets
  //       --ids(celldata)
  //         --data
  doc->InsertEndChild(e_root);
  e_root->InsertEndChild(e_poly);
  e_poly->InsertEndChild(e_piece);
  e_piece->InsertEndChild(e_points);
  e_piece->InsertEndChild(e_lines);
  e_piece->InsertEndChild(e_ids);
  e_points->InsertEndChild(e_da_p);
  e_lines->InsertEndChild(e_da_con);
  e_lines->InsertEndChild(e_da_off);
  e_ids->InsertEndChild(e_da_ids);

  // Initializes the root element
  e_root->SetAttribute("type", "PolyData");
  e_root->SetAttribute("version", "0.1");
  e_root->SetAttribute("byte_order", "LittleEndian");

  //--------------------------------------------------------------------
  // Initializes the Piece element, fills in data
  //--------------------------------------------------------------------

  // Initializes the DataArray of Points element and IDs element
  std::ostringstream ss_pt, ss_id; // Points
  ss_pt << '\n';
  ss_id << '\n';

  // Set the precision of point coordinates to 6
  ss_pt << std::fixed << std::setprecision(6);

  for (auto &line : track_coords) {
    // String of point coordinates
    auto &b = line._start;
    auto &e = line._end;
    ss_pt << ' ' << b.getX() << ' ' << b.getY() << ' ' << b.getZ()
          << ' ' << e.getX() << ' ' << e.getY() << ' ' << e.getZ();

    // String of track IDs
    ss_id << ' ' << line._uid;
  }

  ss_pt << '\n';
  e_da_p->SetAttribute("type", "Float64");
  e_da_p->SetAttribute("Name", data_name.c_str());
  e_da_p->SetAttribute("NumberOfComponents", "3");
  e_da_p->SetAttribute("format", "ascii");
  e_da_p->SetText(ss_pt.str().c_str());

  ss_id << '\n';
  e_da_ids->SetAttribute("type", "Int64");
  e_da_ids->SetAttribute("Name", "track_ids");
  e_da_ids->SetAttribute("format", "ascii");
  e_da_ids->SetText(ss_id.str().c_str());

  // Initializes connectivity and offsets
  size_t num_points = 0, num_lines = 0;
  std::ostringstream ss_con, ss_off; // Connectivity and offsets
  ss_con << '\n';
  ss_off << '\n';

  num_lines = track_coords.size();
  for (size_t t = 0; t < num_lines; ++t) {
    // Strings of connectivity and offsets
    size_t indx = num_points;
    ss_off << ' ' <<  (num_points += 2);  // Number of points
    for (; indx < num_points; ++indx)
      ss_con << ' ' << indx;
  }

  ss_con << '\n';
  ss_off << '\n';

  e_piece->SetAttribute("NumberOfPoints", std::to_string(num_points).c_str());
  e_piece->SetAttribute("NumberOfLines", std::to_string(num_lines).c_str());

  // Initializes the DataArray of connectivity
  e_da_con->SetAttribute("type", "Float64");
  e_da_con->SetAttribute("Name", "connectivity");
  e_da_con->SetAttribute("format", "ascii");
  e_da_con->SetText(ss_con.str().c_str());

  // Initializes the DataArray of offsets
  e_da_off->SetAttribute("type", "Float64");
  e_da_off->SetAttribute("Name", "offsets");
  e_da_off->SetAttribute("format", "ascii");
  e_da_off->SetText(ss_off.str().c_str());

}


} /* namespace antmoc */

