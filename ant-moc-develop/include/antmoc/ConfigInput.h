/// \file include/ConfigInput.h
/// \brief Handling configuration input.
/// \date June 13, 2019
/// \author An Wang, USTB (wangan@xs.ustb.edu.cn)

#ifndef CONFIGINPUT_H_
#define CONFIGINPUT_H_

#include <set>
#include <string>
#include <vector>

#include "antmoc/enum_types.h"
#include "antmoc/lattice_utils.h"
#include "antmoc/log.h"
#include "antmoc/mpi_utils.h"
#include "antmoc/Option.h"
#include "antmoc/Point.h"
#include "antmoc/tally_utils.h"

namespace antmoc {

// Forwarding declarations
class Geometry;
class Quadrature;
class Solver;
class TrackGenerator;

///---------------------------------------------------------------------
/// \class ConfigInput
/// \brief A reader for configuration file
/// \details A configuration reader must provide several basic interfaces
///          such as:
///            get arguments (objects)
///            validate arguments
///            print help messages
///
///          Most of the arguments are supposed to be checked when they
///          are set because they may not be used until a certain call
///          to some functional procedure. For example, if it doesn't
///          validate the options of tally mesh, the program won't be
///          terminated until data output. This is usually harmful since
///          we would expect the program to dump data if no exception
///          was throwed during a simulation.
///---------------------------------------------------------------------
class ConfigInput : public Option {

public:

  ConfigInput() = default;
  ConfigInput(int &argc, char **argv): Option(argc, argv) { }
  ConfigInput(const StringVec &argv): Option(argv) { }
  ConfigInput(std::initializer_list<std::string> ilist): Option(ilist) { }
  virtual ~ConfigInput() = default;

  //--------------------------------------
  // Interfaces
  //--------------------------------------
  virtual void showHelp() = 0;
  virtual void showHelpAsNeeded() = 0;
  virtual void showExtraHelpAsNeeded() = 0;
  virtual void validateArguments();

  //--------------------------------------
  // Logging and reporting
  //--------------------------------------
  void initializeLogger();
  virtual void printArgumentsReport();

  log::level getLogLevel();
  size_t getLogLineLength();
  std::string getLogFile();
  std::string getOutputDirectory();
  std::set<int> getLoggingRanks();

  //--------------------------------------
  // Parallel computing options
  //--------------------------------------
  int getNumThreads();
  std::vector<int> getNumDomains();
  std::string getDefaultNumDomains();
  trackMappingType getTrackMappingType();

  //-------------------------------------------
  // Configuration, geometry, materials options
  //-------------------------------------------
  std::string getConfInputPath();
  std::string getGeoInputPath();
  std::string getGlobalPrimitivesPath();
  std::string getMatInputPath();
  XSFileLayout getXSFileLayout();
  std::vector<int> getLatticeRefines();
  int getCellSectors();
  int getCellRings();
  bool doesCleanupPrimitives();
  bool doesCleanupMaterials();

  //--------------------------------------
  // Dumping options
  //--------------------------------------
  bool doesDumpSettings();
  bool doesDumpVisualizationData();
  bool doesDumpFSRData();

  std::string getTallyMeshStr();
  std::string getTallyMeshOffsetStr();
  std::string getTallyEnergyGroupsStr();

  tallyMeshType getTallyMeshType();
  std::vector<WidthVec> getTallyMesh();
  Point getTallyMeshOffset();
  std::string getTallyMeshOrientation();
  std::set<TallyType> getDumpRXTypes();
  std::set<TallyType> getDumpXSTypes();
  std::set<TallyType> getDumpTracksTypes();
  std::set<int> getTallyEnergyGroups();

  std::string getTimeStamp();
  bool isHexTallyMesh();

  //--------------------------------------
  // Ray tracing options
  //--------------------------------------
  std::string getZonesStr();
  std::string getZMeshTypeStr();

  size_t getNumAzim();
  size_t getNumPolar();
  double getAzimSpacing();
  double getPolarSpacing();
  segmentationType getSegmentationType();
  std::vector<FP_PRECISION> getZones();
  std::vector<int> getNumModules();
  bool isGlobalZMesh();

  //--------------------------------------
  // Solver options
  //--------------------------------------
  std::string getKeffTypeStr();

  solverType getSolverType();
  quadratureType getQuadratureType();
  bool getKeffFromNeutronBalance();
  int getMaxIterations();
  double getTolerance();
  stabilizationType getStabilizationType();
  double getStabilizationFactor();

  /// \brief Infos the solver to check MGXS.
  bool doesCheckXS();

  /// \brief Infos the solver to check and correct MGXS.
  bool doesCorrectXS();

  /// \brief Log level of XS checking messages. (defaults to error)
  log::level getCheckXSLogLevel();

  //-------------------------------------
  // Experimental options
  //-------------------------------------
  std::string getRequiredEnergyGroupsStr();

  int getRequiredEnergyGroups();


protected:

  const std::string _l_loglevel = "log-level";
  const std::string _l_logline  = "log-line";
  const std::string _l_logfile  = "log-file";
  const std::string _a_logfile  = "o,log-file";
  const std::string _l_logranks  = "logging-ranks";

  //-------------------------------------
  // Parallel computing
  //-------------------------------------
  const std::string _l_threads = "num-omp-threads";
  const std::string _l_domains = "num-domains";
  const std::string _l_track_mp = "track-mapping";

  //-------------------------------------
  // I/O options
  //-------------------------------------
  const std::string _l_conf     = "config";
  const std::string _a_conf     = "c,config";
  const std::string _l_geo      = "geometry";
  const std::string _a_geo      = "g,geometry";
  const std::string _l_prim     = "primitives";
  const std::string _l_mat      = "materials";
  const std::string _a_mat      = "m,materials";
  const std::string _l_xsl      = "xs-layout";

  const std::string _l_lat_refines  = "lattice-refines";
  const std::string _l_cell_sectors = "cell-sectors";
  const std::string _l_cell_rings   = "cell-rings";

  const std::string _l_prim_cl  = "cleanup-primitives";
  const std::string _l_mat_cl   = "cleanup-materials";

  const std::string _l_outputdir = "output-directory";
  const std::string _a_outputdir = "d,output-directory";
  const std::string _l_dump_settings = "dump-settings";
  const std::string _l_dump_fsrs     = "dump-fsrs";
  const std::string _l_dumptracks    = "dump-tracks";

  const std::string _l_meshtype     = "mesh-type";
  const std::string _l_mesh         = "mesh";
  const std::string _l_meshoffset   = "mesh-offset";
  const std::string _l_orientation  = "orientation";
  const std::string _l_dumprx       = "dump-rx";
  const std::string _l_dumpxs       = "dump-xs";
  const std::string _l_tallygroups  = "tally-groups";

  //-------------------------------------
  // Quadrature
  //-------------------------------------
  const std::string _l_quad     = "quadrature";
  const std::string _a_quad     = "q,quadrature";
  const std::string _l_nazims   = "num-azims";
  const std::string _a_nazims   = "a,num-azims";
  const std::string _l_npolars  = "num-polars";
  const std::string _a_npolars  = "p,num-polars";
  const std::string _l_sazim    = "azim-spacing";
  const std::string _a_sazim    = "s,azim-spacing";
  const std::string _l_spolar   = "z-spacing";
  const std::string _a_spolar   = "l,z-spacing";

  //-------------------------------------
  // Ray tracing options
  //-------------------------------------
  const std::string _l_segform    = "formation";
  const std::string _a_segform    = "f,formation";
  const std::string _l_segzones   = "zones";
  const std::string _a_segzones   = "z,zones";
  const std::string _l_modules    = "modules";
  const std::string _l_zmesh      = "z-mesh";

  //-------------------------------------
  // Solver options
  //-------------------------------------
  const std::string _l_solver   = "solver";
  const std::string _l_niters   = "max-iters";
  const std::string _a_niters   = "i,max-iters";
  const std::string _l_keff     = "keff-neutron-balance";
  const std::string _l_tol      = "tolerance";
  const std::string _a_tol      = "t,tolerance";
  const std::string _l_stab     = "stabilization";
  const std::string _l_stab_factor = "stabilization-factor";
  const std::string _l_check_xs   = "check-xs";
  const std::string _l_correct_xs = "correct-xs";
  const std::string _l_check_xs_log_level = "check-xs-log-level";

  //-------------------------------------
  // Experimental options
  //-------------------------------------
  const std::string _l_requiredgroups = "required-groups";
};


} // namespace antmoc

#endif  // CONFIGINPUT_H_
