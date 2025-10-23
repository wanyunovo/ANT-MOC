#include "antmoc/Cmfd.h"
#include "antmoc/ConfigInputFile.h"
#include "antmoc/debug_utils.h"
#include "antmoc/enum_types.h"
#include "antmoc/Factory.h"
#include "antmoc/file_utils.h"
#include "antmoc/FSRDataHandlerHDF5.h"
#include "antmoc/GeoInputXml.h"
#include "antmoc/MaterialHandlerHDF5.h"
#include "antmoc/Mesh.h"
#include "antmoc/mpi_utils.h"
#include "antmoc/log.h"
#include "antmoc/Timer.h"
#include "antmoc/TrackGenerator3D.h"

#ifdef ENABLE_HIP_
#include "antmoc/hip/hip_info.h"
#endif

using namespace antmoc;

void printTimerReport(Timer &);
std::shared_ptr<Mesh> createSimpleMesh(Factory::ConfInputPtr, Factory::SolverPtr);
void dumpMeshData(Factory::ConfInputPtr, std::shared_ptr<Mesh>);
void dumpDataToH5File(Factory::ConfInputPtr, Factory::SolverPtr);


int main(int argc, char* argv[]){

#ifdef ENABLE_MPI_
  int provided;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_SERIALIZED, &provided);
  // Initialize MPI static variables to default values
  mpi::defaultInitialize();
#endif
  // Initialize the logger
  log::initialize();

  Timer timer;
  timer.startTimer();

  //------------------------------------------------------------
  // Read simulation arguments
  //------------------------------------------------------------
  timer.startTimer();
  auto conf = Factory::getConfInput<ConfigInputFile>(argc, argv);

  // Print the help message when needed
  if (argc < 2) {
    conf->showHelp();
  }
  else {
    conf->showHelpAsNeeded();
  }

  // Set log level
  conf->initializeLogger();

  // Print the title
  log::title("ANT-MOC: an advanced neutron transport code from the CVR suite");

  // Report input arguments
  conf->printArgumentsReport();

  // Validate input parameters
  conf->validateArguments();

  timer.stopTimer("Read Settings");

  // Print pid and hostname with log::profile
#ifdef ENABLE_MPI_
  mpi::mpiBarrier();
#endif
  printHostnames();

#ifdef ENABLE_HIP_
  hip::printDeviceProperties();
#endif

#ifdef ENABLE_MPI_
  // Report MPI datatypes
  mpi::showMPIDatatypes();

  //------------------------------------------------------------
  // Domain decomposition
  //------------------------------------------------------------
  timer.startTimer();
  auto domains = conf->getNumDomains();
  mpi::setDomainDecomposition(MPI_COMM_WORLD,
                              domains[0],
                              domains[1],
                              domains[2],
                              domains[3]);
  timer.stopTimer("Domain Decomposition");
#endif

  //------------------------------------------------------------
  // Initialize the MaterialHandler for later use
  //------------------------------------------------------------
  log::info("Creating geometry and materials...");

  timer.startTimer();

  Factory::MaterialHandlerPtr mat_input;
  mat_input = Factory::getMaterialHandler<MaterialHandlerHDF5>(conf);

  // Set the number of energy group for debugging
  {
    auto n = conf->getRequiredEnergyGroups();
    if (n > 0)
      mat_input->setNumEnergyGroups(n);
  }
  //mat_input->readAllMaterials();

  /* Create Non-uniform CMFD mesh */
  /*
  log::info("Creating CMFD mesh...");
  int axial_refines = 15;
  Cmfd* cmfd = new Cmfd();
  cmfd->setSORRelaxationFactor(1.5);
  cmfd->setCMFDRelaxationFactor(0.7);
  std::vector<std::vector<double>> widths = {
    {1.26,1.26,1.26},
    {1.26,1.26,1.26},
    {10.71,10.71,10.71,10.71}
  };
  cmfd->setWidths(widths);
  */
  
  /* Create CMFD mesh 创建CMFD网格 */
  log::fdebug("Creating CMFD mesh...");
  Cmfd* cmfd = new Cmfd();
  cmfd->setSORRelaxationFactor(1.5);
  //cmfd->setCMFDRelaxationFactor(0.7);
  cmfd->setLatticeStructure(3, 3, 30);//设置CMFD网格的逻辑结构  每一个cmfd网格都是一个lattice 重要 

  if(!conf->isHexTallyMesh()){  //设置MOC能群和CMFD能群的对应关系 四边形  能群压缩
    std::vector<std::vector<int> > cmfd_group_structure;
    cmfd_group_structure.resize(2);
    for (int g=0; g<3; g++)
      cmfd_group_structure.at(0).push_back(g+1);//这里g+1是为了便于统计总共有多少个MOC能群数，以及便于后续判断CMFD对应的MOC能群数是否是线性递增的
    for (int g=3; g<7; g++)
      cmfd_group_structure.at(1).push_back(g+1);
    cmfd->setGroupStructure(cmfd_group_structure);
  }
  cmfd->setKNearest(3);//设置更新MOC通量所关联的径向CMFD网格数量

  // set for HexLattice 000000000
  cmfd->setHexLatticeEnable(conf->isHexTallyMesh());
  // CMFD Lattice aligns with mesh
  if(conf->isHexTallyMesh()){
    log::fdebug("Create CMFD Lattice with Hex ...");
    cmfd->setHexGroups(10);
    auto widths = conf->getTallyMesh();
    cmfd->setOrientation(conf->getTallyMeshOrientation());
    cmfd->setNumR(widths[0][0]);
    cmfd->setWidthR(widths[1][0]);
    cmfd->setNumZ(10);
    //cmfd->setWidthsZ(widths[2]);
  }


  //------------------------------------------------------------
  // Create the geometry
  //------------------------------------------------------------
  auto nmods = conf->getNumModules();
  Geometry geometry;
  geometry.setNumDomainModules(nmods[0], nmods[1], nmods[2]);

  auto geo_input =
    Factory::getGeoInput<GeoInputXml> (&geometry, mat_input, conf);

  // Try to read primitives
  geo_input->readGlobalPrimitives(conf->getGlobalPrimitivesPath());
  // Read the geometry
  geo_input->readGeometryFromFile(conf->getGeoInputPath());

  // Dump geometry settings as needed
  if (conf->doesDumpSettings()) {
    geo_input->dumpSettings(conf->getOutputDirectory());
  }

  timer.stopTimer("Read Geometry");

  // Remove unused CSG objects and delete temporary containers
  timer.startTimer();
  if (conf->doesCleanupPrimitives())
    geo_input->clear();
  timer.stopTimer("Clean Up");

#ifdef ENABLE_MPI_
  timer.startTimer();
  if (mpi::isSpatialDecomposed()) {
    geometry.setDomainDecomposition();
  }
  timer.stopTimer("Domain Decomposition");
#endif
  
  //set CMFD for
  geometry.setCmfd(cmfd);//将cmfd指针传递给geometry  关联几何

  // Initialize FSR
  timer.startTimer();
  geometry.initializeFlatSourceRegions();//涉及CMFD的初始化
  timer.stopTimer("Meshing");

  // Remove unused materials
  timer.startTimer();
  if (conf->doesCleanupMaterials())
    geo_input->eraseUnusedMaterials();
  timer.stopTimer("Clean Up");

  // Print the number of materials and delete the pointer
  geo_input->printReport();
  geo_input.reset();

  // Print the memory usage of materials
  geometry.printMemUsageMaterials();

  // FIXME: HDF5 interfaces must be called before mpi finalizing
  timer.startTimer();
  mat_input.reset();
  timer.stopTimer("Clean Up");

  timer.stopTimer("Input");
  printTimerReport(timer);

  //------------------------------------------------------------
  // Generate tracks
  //------------------------------------------------------------
  log::info("Initializing the track generator...");
  auto quad = Factory::getQuadrature(conf);
  auto tg   = Factory::getTrackGenerator(&geometry, conf);
  tg->setQuadrature(quad);
  tg->generateTracks();

  //------------------------------------------------------------
  // Run simulation
  //------------------------------------------------------------
  log::info("Running simulation...");
  SolverPtr solver = Factory::getSolver(tg, conf);

  solver->computeEigenvalue(conf->getMaxIterations());
  solver->printTimerReport();

  //------------------------------------------------------------
  // Output reaction rates, fluxes and volumes
  //------------------------------------------------------------
  if (conf->doesDumpVisualizationData()) {
    auto mesh_simple = createSimpleMesh(conf, solver);
    dumpMeshData(conf, mesh_simple);
  }

  // Dump output data to HDF5 files
  if (conf->doesDumpFSRData() && mpi::isInMainCommUniqueDomains())
    dumpDataToH5File(conf, solver);

  log::header("Finished");

#ifdef ENABLE_MPI_
  MPI_Finalize();
#endif
  return 0;
}


//----------------------------------------------------------------------
// Helpers
//----------------------------------------------------------------------
/// \brief Print a timing report for pre-processing
void printTimerReport(Timer &timer) {

#ifdef ENABLE_MPI_
  timer.reduceTimer(mpi::getMPIComm());
#endif

  log::header("Timing Report - Input (average)");

  const std::string tot_string = "Input";
  timer.printSplit(tot_string, "Total Input Time");

  timer.printSplit("Read Settings", "Time to read settings",
                   1, tot_string);

  timer.printSplit("Read Geometry", "Time to read geometry and materials",
                   1, tot_string);

#ifdef ENABLE_MPI_
  timer.printSplit("Domain Decomposition", "Time of domain decomposition",
                   1, tot_string);
#endif

  timer.printSplit("Meshing", "Time of meshing and FSRs initialization",
                   1, tot_string);

  timer.printSplit("Clean Up", "Time to remove unused objects",
                   1, tot_string);

  log::separator("-");
}


/// \brief Create a Mesh object with the shape of a simple lattice.
/// \param conf A pointer to ConfInput object.
/// \param solver A pointer to Solver object.
/// \return A pointer to the new Mesh object.
std::shared_ptr<Mesh> createSimpleMesh(Factory::ConfInputPtr conf, Factory::SolverPtr solver) {

  log::info("Creating simple tally mesh...");

  // Define the mesh object
  std::shared_ptr<Mesh> mesh = nullptr;

  // Retrieve the shape of the lattice
  auto widths = conf->getTallyMesh();

  // Initialize the mesh object
  mesh = std::make_shared<Mesh>(solver);

  if (conf->isHexTallyMesh()) { // Hexagonal
    int num_r = (int) widths[0][0];
    double width_r = widths[1][0];
    mesh->createLattice(num_r, width_r, widths[2],
                        conf->getTallyMeshOrientation(),
                        conf->getTallyMeshOffset());
  } else {  // Rectangular
    mesh->createLattice(widths[0], widths[1], widths[2]);
  }

  return mesh;
}


/// \brief Dump visualization data to files.
/// \details This method dumps rates, cross-sections, volumes and
///          tracks data on demand.
/// \param conf A pointer to ConfInput object.
/// \param mesh A mesh object to dump visualization data.
void dumpMeshData(Factory::ConfInputPtr conf, std::shared_ptr<Mesh> mesh) {

  // Create the directory if it doesn't exist
  std::string dump_dir = conf->getOutputDirectory();
  fileutils::createDirectory(dump_dir);

  // Cross-sections are ignored because I don't not how to sum them up across FSRs
  auto field_types = conf->getDumpRXTypes();
  tallyutils::expandTallyTypeForRX(field_types);

  // Selected energy groups
  auto energy_groups = conf->getTallyEnergyGroups();
  mesh->dumpMeshDataToFile(dump_dir + "/reaction_rates.vtu", field_types, energy_groups);

  // Dump 2D or 3D tracks
  // FIXME: need to be refactored to take advantage of PHDF5
  //mesh.dumpTrackMeshToDir(tg, conf->getDumpTracksTypes(), dump_dir);
}


/// \brief Dump visualization data to HDF5 files.
/// \details This method dumps rates, cross-sections, volumes and
///          tracks data on demand.
/// \param conf a pointer to ConfInput object
/// \param solver a pointer to Solver object
void dumpDataToH5File(Factory::ConfInputPtr conf, Factory::SolverPtr solver) {

  // Create the directory if it doesn't exist
  std::string dump_dir = conf->getOutputDirectory();
  std::string file = dump_dir + "/fsr_data.h5";
  fileutils::createDirectory(dump_dir);

  // Selected field types
  auto field_types = tallyutils::combineRXAndXS(conf->getDumpRXTypes(),
                                                conf->getDumpXSTypes());

  // Create a file and write data to it
  FSRDataHandlerHDF5 h5_writer(file, HDF5Mode::Truncate, solver);

  // Selected energy groups
  auto energy_groups = conf->getTallyEnergyGroups();

  h5_writer.writeMetaData();
  h5_writer.dumpFSRData(field_types, energy_groups);
}
