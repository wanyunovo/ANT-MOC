/// \file src/ConfigInput.cpp

#include "antmoc/ConfigInput.h"
#include "antmoc/VecParser.h"

#include <cmath>
#include <ctime>
#include <locale>

namespace antmoc {

/// \brief Check user-defined arguments
/// \details It is usually unnecessary to check all arguments.
///          Besides, most of the arguments have been checked when they
///          are set.
void ConfigInput::validateArguments() {
  // Check file paths
  if (stringutils::isSpaces( getGeoInputPath() ))
      log::error("Path to the geometry file is empty!");
  if (stringutils::isSpaces( getMatInputPath() ))
      log::error("Path to the material file is empty!");

  // Check tally mesh
  getTallyMesh();
}


//----------------------------------------------------------------------
// Logging and reporting
//----------------------------------------------------------------------

/// \brief Initialize the logger
void ConfigInput::initializeLogger() {

  log::set_logging_ranks(getLoggingRanks());
  log::set_line_length(getLogLineLength());
  log::set_level(getLogLevel());

  auto path = getLogFile();

  // Create the log file only as needed
  if (!path.empty()) {
    log::set_path(path);
  }

}


/// \brief Returns the log level
log::level ConfigInput::getLogLevel() {

  log::level e = log::level::info;
  auto s = getOptionValue(_l_loglevel);

  // Keep the default if there is an empty string
  if (!stringutils::isSpaces(s)) {
    try {
      e = log::level::_from_string_nocase(s.c_str());
    } catch (...) {
      log::error("Invalid log level: {}", s);
    }
  }

  return e;

}


/// \brief Returns the log line length
/// \return The length of a log line (between 64 and 256)
size_t ConfigInput::getLogLineLength() {

  size_t length = getOptionValueSizet(_l_logline);

  if (length < 64 || length > 256) {
    log::error("Log line length must be a number within [64, 256]: {}" , length);
  }

  return length;

}


/// \brief Returns the log file path
std::string ConfigInput::getLogFile() {

  auto path = stringutils::trim( getOptionValue(_l_logfile) );

  if (stringutils::isSpaces(path)) {
    path = "";
  }
  else {

    auto pos = path.find_last_of("/");

    if (pos == path.size() - 1) {
      // No file name specified, take the default
      path = path + "antmoc.0.out";
    }
  }
  return path;
}


/// \brief Returns the directory for dumping data
std::string ConfigInput::getOutputDirectory() {
  return stringutils::trim( getOptionValue(_l_outputdir) );
}


/// \brief Returns the list of logging ranks
std::set<int> ConfigInput::getLoggingRanks() {

  std::set<int> groups;
  auto s = stringutils::trim(getOptionValue(_l_logranks));

  if (!stringutils::isSpaces(s)) {
    auto strvec = stringutils::splitString(s, ",");
    for (auto &w : strvec) {
      stringutils::trim(w);
      auto int_set = stringutils::toIntegerSet(w, ":");
      groups.insert(int_set.begin(), int_set.end());
    }
  }

  return groups;
}


/// \brief Report input arguments
void ConfigInput::printArgumentsReport() {

  const int name_width = 24;

  log::verbose_once("Parsing arguments: {}", *this);

  auto num_procs = mpi::getNumProcs();
  auto domains = getNumDomains();
  int num_procs_shared_domain = num_procs / domains[0]
                                / domains[1]
                                / domains[2];

  log::info("\nReporting runtime parameters...");
  log::info("{:-<{}}", "", name_width);
  log::info("{:<{}} = {}", "Log level",         name_width, stringutils::toUpper(enumToString(getLogLevel())));
  log::info("{:<{}} = {}", "Log line length",   name_width, getLogLineLength());
  log::info("{:<{}} = {}", "Log file",          name_width, getLogFile());

  // Report MPI and OpenMP status
  log::info("{:-<{}}", "", name_width);
  log::info("{:<{}} = {}", "Ranks",             name_width, num_procs);
  log::info("{:<{}} = {}", "Threads per rank",  name_width, getNumThreads());
  log::info("{:<{}} = {}", "Spatial Domains",   name_width, stringutils::join(getNumDomains(), " x "));
  log::info("{:<{}} = {}", "Sub track domains", name_width, num_procs_shared_domain);
  log::info("{:<{}} = {}", "Track mapping",     name_width, enumToString(getTrackMappingType()));
  log::info("{:-<{}}", "", name_width);

  // Report input parameters
  log::info("{:<{}} = {}", "Config file",     name_width, getConfInputPath());
  log::info("{:<{}} = {}", "Geometry file",   name_width, getGeoInputPath());
  log::info("{:<{}} = {}", "Primitives file", name_width, getGlobalPrimitivesPath());
  log::info("{:<{}} = {}", "Materials file",  name_width, getMatInputPath());
  log::info("{:<{}} = {}", "XS file layout",  name_width, enumToString(getXSFileLayout()));
  log::info("{:<{}} = {}", "Global lattice refines",  name_width, stringutils::join(getLatticeRefines(), ", "));

  auto sectors = getCellSectors();
  auto rings = getCellRings();
  log::info("{:<{}} = {}", "Global cell sectors",     name_width, sectors < 0 ? "" : std::to_string(sectors));
  log::info("{:<{}} = {}", "Global cell rings",       name_width, rings < 0 ? "" : std::to_string(rings));

  // Report the quadrature
  log::info("{:-<{}}", "", name_width);
  log::info("{:<{}} = {}", "Quadrature",        name_width, enumToString(getQuadratureType()));
  log::info("{:<{}} = {}", "Azimuthal angles",  name_width, getNumAzim());
  log::info("{:<{}} = {}", "Azimuthal spacing", name_width, getAzimSpacing());
  log::info("{:<{}} = {}", "Polar angles",      name_width, getNumPolar());
  log::info("{:<{}} = {}", "Polar spacing",     name_width, getPolarSpacing());
  log::info("{:<{}} = {}", "Segment formation", name_width, enumToString(getSegmentationType()));
  log::info("{:<{}} = {}", "Modules",           name_width, stringutils::join(getNumModules(), ", "));
  log::info("{:<{}} = {}", "Axial zones",       name_width, getZonesStr());
  log::info("{:<{}} = {}", "Z mesh type",       name_width, getZMeshTypeStr());

  log::info("{:-<{}}", "", name_width);
  log::info("{:<{}} = {}", "Solver",                name_width, enumToString(getSolverType()));
  log::info("{:<{}} = {}", "Max iterations",        name_width, getMaxIterations());
  log::info("{:<{}} = {:.4E}", "Tolerance",         name_width, getTolerance());
  log::info("{:<{}} = {}", "K-eff",                 name_width, getKeffTypeStr());
  log::info("{:<{}} = {}", "Stabilization",         name_width, enumToString(getStabilizationType()));
  if (getStabilizationType() != +stabilizationType::NONE)
    log::info("{:<{}} = {}", "Stabilization factor",  name_width, getStabilizationFactor());

  if (doesCheckXS()) {
    log::info("{:<{}} = {}", "MGXS checks",     name_width, doesCheckXS());
    log::info("{:<{}} = {}", "MGXS correction", name_width, doesCorrectXS());
    log::info("{:<{}} = {}", "MGXS log level",  name_width, getCheckXSLogLevel());
  }

  log::info("{:-<{}}", "", name_width);


  log::info("\nReporting output configuration...");
  log::info("{:-<{}}", "", name_width);

  log::info("{:<{}} = {}", "Output directory",    name_width, getOutputDirectory());
  log::info("{:<{}} = {}", "Output settings ",    name_width, doesDumpSettings());
  log::info("{:<{}} = {}", "Output tally mesh",   name_width, doesDumpVisualizationData());
  log::info("{:<{}} = {}", "Output FSR data",     name_width, doesDumpFSRData());
  log::info("{:<{}} = {}", "Output track types",  name_width, tallyutils::join(getDumpTracksTypes(), ", "));

  if (doesDumpVisualizationData() || doesDumpFSRData()) {
    log::info("{:<{}} = {}", "Tally RX types",    name_width, tallyutils::join(getDumpRXTypes(), ", "));
    log::info("{:<{}} = {}", "Tally XS types",    name_width, tallyutils::join(getDumpXSTypes(), ", "));
    log::info("{:<{}} = {}", "Selected groups",   name_width, getTallyEnergyGroupsStr());

    if (doesDumpVisualizationData()) {
      log::info("{:<{}} = {}", "Tally mesh type",   name_width, enumToString(getTallyMeshType()));
      log::info("{:<{}} = {}", "Tally mesh shape",  name_width, getTallyMeshStr());
      log::info("{:<{}} = {}", "Tally mesh offset", name_width, getTallyMeshOffsetStr());

      if (isHexTallyMesh())
        log::info("{:<{}} = {}", "Mesh orientation", name_width, getTallyMeshOrientation());
    }
  }
  else {
    log::info("No output specified.");
  }

  log::info("{:-<{}}", "", name_width);


  log::info("\nReporting other options...");
  log::info("{:-<{}}", "", name_width);
  log::info("{:<{}} = {}", "Required energy group", name_width, getRequiredEnergyGroupsStr());
  log::info("{:-<{}}", "", name_width);
  log::info("");
}


//----------------------------------------------------------------------
// Parallel computing options
//----------------------------------------------------------------------

/// \brief Returns the number of omp threads
int ConfigInput::getNumThreads() {
  return getOptionValueInt(_l_threads);
}


/// \brief Returns the number of domains in each direction
std::vector<int> ConfigInput::getNumDomains() {

  auto s = getOptionValue(_l_domains);
  if (stringutils::isSpaces(s)) {
#ifdef MIT_
    s = getDefaultNumDomains();
#else
    s = "1, 1, 1";
#endif
  }

  std::vector<int> domains;
  auto strvec = stringutils::splitString(s, ",");

  for (auto &w : strvec) {
    stringutils::trim(w);
    domains.push_back(std::stoi(w));
  }

  if (domains.size() < 3 || domains.size() > 4)
    log::ferror("Number of domains must be specified in "
                "3 or 4 directions: %s", s);

  if (domains[0] < 1 || domains[1] < 1 || domains[2] < 1) {
    log::ferror("Number of domains must be no less than 1"
                " : %s", s);
  }

  // Default number of track domains
  auto num_procs = mpi::getNumProcs();
  int ns = num_procs / domains[0] / domains[1] / domains[2];
  if (domains.size() == 3)
    domains.push_back(ns);
  else if (domains[3] < 1)
    domains[3] = ns;

  return domains;
}


/// \brief Get the default number of domains if purely domain
///        decomposition is performed
std::string ConfigInput::getDefaultNumDomains() {
  // Handle the number of domains
  auto num_procs = mpi::getNumProcs();
  int exp = std::log2(num_procs);
  int nxyz = std::pow(2, exp/3);

  std::vector<int> domains = {nxyz, nxyz, nxyz};
  switch (exp % 3) {
    case 1 : domains[0] <<= 1;
             break;
    case 2 : domains[1] = (domains[0] <<= 1);
             break;
    default: break;
  }

  // FIXME
  std::ostringstream ss;
  ss << domains[0] << ',' << domains[1] << ',' << domains[2];

  return ss.str();
}


/// \brief Set the track mapping type
trackMappingType ConfigInput::getTrackMappingType() {

  auto s = getOptionValue(_l_track_mp);

  // Keep the default
  if (stringutils::isSpaces(s))
    s = "Auto";

  auto type = stringutils::spaceToUnderscore(s);
  return trackMappingType::_from_string_nocase(type.c_str());

}


//----------------------------------------------------------------------
// Configuration, geometry, materials options
//----------------------------------------------------------------------

/// \brief Returns the path of settings file
std::string ConfigInput::getConfInputPath() {
  return getOptionValue(_l_conf);
}

/// \brief Returns the path of geometry file
std::string ConfigInput::getGeoInputPath() {
  return getOptionValue(_l_geo);
}

/// \brief Returns the path of primitives file
std::string ConfigInput::getGlobalPrimitivesPath() {
  return getOptionValue(_l_prim);
}

/// \brief Returns the path of materials file
std::string ConfigInput::getMatInputPath() {
  return getOptionValue(_l_mat);
}


/// \brief Returns the type of xs file layout
XSFileLayout ConfigInput::getXSFileLayout() {

  auto s = getOptionValue(_l_xsl);

   // Keep the default
  if (stringutils::isSpaces(s))
    s = "NAMED";

  auto type = stringutils::spaceToUnderscore(s);
  return XSFileLayout::_from_string_nocase(type.c_str());

}


/// \brief Indicate whether unused primitives should be deleted
bool ConfigInput::doesCleanupPrimitives() {
  return getOptionValueBool(_l_prim_cl);
}


/// \brief Indicate whether unused materials should be deleted
bool ConfigInput::doesCleanupMaterials() {
  return getOptionValueBool(_l_mat_cl);
}


/// \brief Returns the global lattice refines
std::vector<int> ConfigInput::getLatticeRefines() {

  auto s = getOptionValue(_l_lat_refines);
  if (stringutils::isSpaces(s))
    s = "0, 0, 0";

  std::vector<int> refines;
  auto strvec = stringutils::splitString(s, ",");

  for (auto &w : strvec) {
    stringutils::trim(w);

    try {
      refines.push_back(std::stoi(w));
    } catch (...) {
      log::error("Invalid global refines: {}", s);
    }
  }

  if (refines.size() != 3)
    log::error("Number of global refines must be specified in 3 directions: {}", s);

  if (refines[0] < 0 || refines[1] < 0 || refines[2] < 0) {
    log::error("Number of global refines must be positive: {}", s);
  }

  return refines;
}


/// \brief Returns the global number of cell sectors.
/// \return The number of sectors (non-negative) or -1 (do nothing).
int ConfigInput::getCellSectors() {

  auto s = getOptionValue(_l_cell_sectors);
  int sectors = -1;

  if (!stringutils::isSpaces(s)) {
    try {
      sectors = std::stoi(s);
    } catch (...) {
      log::error("Invalid number of cell sectors: {}", s);
    }

    if (sectors < 0)
      log::error("Unable to set a negative number of sectors: {}", sectors);
  }

  return sectors;
}


/// \brief Returns the global number of cell rings.
/// \return The number of rings (non-negative) or -1 (do nothing).
int ConfigInput::getCellRings() {

  auto s = getOptionValue(_l_cell_rings);
  int rings = -1;

  if (!stringutils::isSpaces(s)) {
    try {
      rings = std::stoi(s);
    } catch (...) {
      log::error("Invalid number of cell rings: {}", s);
    }

    if (rings < 0)
      log::error("Unable to set a negative number of rings: {}", rings);
  }

  return rings;
}


//----------------------------------------------------------------------
// Dumping options
//----------------------------------------------------------------------
/// \brief Returns a boolean indicating dumping settings or not.
/// \details A directory is required.
bool ConfigInput::doesDumpSettings() {
  return getOptionValueBool(_l_dump_settings)
         && !getOutputDirectory().empty();
}


/// \brief Returns a boolean indicating dumping visualization data or not.
/// \details A directory, at least one tally type and a tally mesh are required.
bool ConfigInput::doesDumpVisualizationData() {
  auto field_types = getDumpRXTypes();
  tallyutils::expandTallyTypeForRX(field_types);

  return !(getOutputDirectory().empty() ||
           getTallyMesh().empty() ||
           field_types.empty());
}


/// \brief Returns a boolean indicating dumping FSR data or not.
/// \details A directory and at least one tally type are required.
bool ConfigInput::doesDumpFSRData() {
  auto field_types = tallyutils::combineRXAndXS(getDumpRXTypes(),
                                                getDumpXSTypes());

  return getOptionValueBool(_l_dump_fsrs)
         && !getOutputDirectory().empty()
         && !field_types.empty();
}


/// \brief Returns the type of tally mesh
tallyMeshType ConfigInput::getTallyMeshType() {

  auto s = getOptionValue(_l_meshtype);

  // Keep the default
  if (stringutils::isSpaces(s))
    s = "RECTANGLE";

  auto type = stringutils::spaceToUnderscore(s);
  return tallyMeshType::_from_string_nocase(type.c_str());

}

/// \brief Indicate whether we have a hexagonal tally mesh 这段代码的主要功能是判断当前的网格类型是否为六边形网格（HEXAGON），并返回一个布尔值。
bool ConfigInput::isHexTallyMesh() {
  return getTallyMeshType() == +tallyMeshType::HEXAGON;
}


/// \brief Returns the tally mesh
std::vector<WidthVec> ConfigInput::getTallyMesh() {

  std::vector<WidthVec> tmesh;
  auto s = getOptionValue(_l_mesh);

  if (!stringutils::isSpaces(s)) {

    auto strvec = stringutils::splitString(s, ",");
    for (auto &w : strvec)
      stringutils::trim(w);

    if (strvec.size() != 3) {
      log::ferror("Tally mesh must be defined in 3 directions "
                  "at the same time: %s", s);
    }

    VecParser p;
    for (auto &w : strvec)
      tmesh.push_back(p.parse(w));

    if (isHexTallyMesh()) {
      if (tmesh[0].size() != 1 || tmesh[1].size() != 1)
        log::ferror("Bad format of hexgonal tally mesh. It must "
                    "be set to 'num_r, width_r, widths_z'");
    }
  }

  return tmesh;
}


/// \brief Returns the tally mesh in a string
std::string ConfigInput::getTallyMeshStr() {
  return getOptionValue(_l_mesh);
}


/// \brief Returns the tally mesh offset
Point ConfigInput::getTallyMeshOffset() {

  auto s = getOptionValue(_l_meshoffset);
  if (stringutils::isSpaces(s))
    s = "0, 0, 0";

  std::vector<double> offset;
  auto strvec = stringutils::splitString(s, ",");

  for (auto &w : strvec) {
    stringutils::trim(w);
    offset.push_back(std::stod(w));
  }

  if (offset.size() != 3) {
    log::ferror("Tally mesh offset must have 3 components: %s", s);
  }

  return Point{offset[0], offset[1], offset[2]};
}


/// \brief Returns the tally mesh offset in a string
std::string ConfigInput::getTallyMeshOffsetStr() {
  auto p = getTallyMeshOffset();
  std::stringstream ss;
  ss << '(' << p.getX() << ", " << p.getY() << ", " << p.getZ() << ')';
  return ss.str();
}


/// \brief Returns the orientation of tally mesh (only for hexagon mesh)
std::string ConfigInput::getTallyMeshOrientation() {

  auto type = getOptionValue(_l_orientation);

  if (stringutils::isSpaces(type))
    type = "Y";

  stringutils::toUpper(stringutils::trim(type));

  if (type != "Y" && type != "X") {
    if (!isHexTallyMesh())  // the type is invalid and will be reset
      type = "Y";
    else
      log::ferror("Undefined tally mesh orientation: %s", type);
  }

  return type;
}


/// \brief Returns tally types of reaction rates
std::set<TallyType> ConfigInput::getDumpRXTypes() {

  auto s = getOptionValue(_l_dumprx);
  if (stringutils::isSpaces(s))
    s = "PHI";

  std::set<TallyType> types;
  auto strvec = stringutils::splitString(s, ",");

  for (auto &w : strvec) {
    stringutils::trim(w);

    // Get an object of TallyType
    auto t = tallyutils::codeToTallyType(w, "RX");
    if (!tallyutils::isRXTallyType(t)) {
      log::ferror("Cannot dump other mesh data to RX data files: %s", w);
    }

    types.insert(t);
  }

  return types;
}


/// \brief Returns tally types
std::set<TallyType> ConfigInput::getDumpXSTypes() {

  auto s = getOptionValue(_l_dumpxs);
  if (stringutils::isSpaces(s))
    s = "None";

  std::set<TallyType> types;
  auto strvec = stringutils::splitString(s, ",");

  for (auto &w : strvec) {
    stringutils::trim(w);

    // Get an object of TallyType
    auto t = tallyutils::codeToTallyType(w, "XS");
    if (!tallyutils::isXSTallyType(t)) {
      log::ferror("Cannot dump other mesh data to XS data files: %s", w);
    }

    types.insert(t);
  }

  return types;
}


/// \brief Parses and returns the string of tallied data types
std::set<TallyType> ConfigInput::getDumpTracksTypes() {

  auto s = getOptionValue(_l_dumptracks);
  if (stringutils::isSpaces(s))
    s = "None";

  std::set<TallyType> types;
  auto strvec = stringutils::splitString(s, ",");

  for (auto &w : strvec) {
    stringutils::trim(w);

    // Get an object of TallyType
    auto t = tallyutils::codeToTallyType(w, "TRACKS");
    if (!tallyutils::isTracksTallyType(t)) {
      log::ferror("Cannot dump other mesh data to tracks data files: %s", w);
    }

    types.insert(t);
  }

  return types;
}


/// \brief Returns the tally mesh offset
std::set<int> ConfigInput::getTallyEnergyGroups() {

  auto s = stringutils::trim(getOptionValue(_l_tallygroups));

  if (stringutils::isSpaces(s))
    s = "ALL";

  std::set<int> groups;

  // Special case
  if (s.size() == 3) {
    std::string str(s);
    stringutils::toUpper(str);
    if (str == "ALL")
      return groups;
  }

  auto strvec = stringutils::splitString(s, ",");
  for (auto &w : strvec) {
    stringutils::trim(w);
    auto int_set = stringutils::toIntegerSet(w, ":");
    groups.insert(int_set.begin(), int_set.end());
  }

  return groups;
}


std::string ConfigInput::getTallyEnergyGroupsStr() {
  auto s = getTallyEnergyGroups();
  if (s.empty())
    return "ALL";
  else
    return stringutils::join(s, ",");
}


/// \brief Get the timestamp
std::string ConfigInput::getTimeStamp() {
  // Time buffer
  char buf[100];

  std::time_t t = std::time(nullptr);
  std::strftime(buf, sizeof(buf), "%Y%m%d-%H%M%S", std::localtime(&t));
#ifdef ENABLE_MPI_
  MPI_Bcast(buf, sizeof(buf), MPI_CHAR, 0, mpi::getMPIComm());
#endif

  return std::string(buf);
}


//----------------------------------------------------------------------
// Ray tracing options
//----------------------------------------------------------------------

/// \brief Set the segmentation type
segmentationType ConfigInput::getSegmentationType() {

  auto s = getOptionValue(_l_segform);

  // Keep the default
  if (stringutils::isSpaces(s))  //检查变量 s 是否为空白（即仅包含空格字符）
    s = "OTF Stacks";

  auto type = stringutils::spaceToUnderscore(s);//将字符串 s 中的空格转换为下划线
  return segmentationType::_from_string_nocase(type.c_str());//将字符串 type 转换为对应的 segmentationType 枚举值，并返回该枚举值

}


/// \brief Returns the segmentation zones
std::vector<FP_PRECISION> ConfigInput::getZones() {

  std::vector<FP_PRECISION> zones;
  const auto s = getZonesStr();

  if (stringutils::toUpper(s) != "AUTO") {

    auto strvec = stringutils::splitString(s, ",");

    for (auto &w : strvec) {
      stringutils::trim(w);
      zones.push_back(std::stod(w));
    }

    if (zones.size() < 2) {
      log::ferror("Segmentation zones must have at least 2 numbers: '%s'", s);
    }
  }

  return zones;
}


/// \brief Returns the segmentation zones in a string
std::string ConfigInput::getZonesStr() {

  const auto s = stringutils::trim(getOptionValue(_l_segzones));

  if (s.empty() || stringutils::toUpper(s) == "AUTO")
    return "AUTO";
  else
    return s;
}


/// \brief Returns the type of axial mesh
std::string ConfigInput::getZMeshTypeStr() {

  if (isGlobalZMesh())
    return "GLOBAL";
  else
    return "LOCAL";
}


/// \brief Determine if the Z Mesh is global
bool ConfigInput::isGlobalZMesh() {

  auto s = getOptionValue(_l_zmesh);
  stringutils::toUpper(stringutils::trim(s));

  if (stringutils::isSpaces(s))
    s = "LOCAL";

  bool result = false;
  if (s == "GLOBAL")
    result = true;
  else if (s == "LOCAL")
    result = false;
  else
    log::ferror("Undefined z-mesh type: %s", s);

  return result;
}


/// \brief Returns the number of modules in each direction
std::vector<int> ConfigInput::getNumModules() {

  auto s = getOptionValue(_l_modules);
  if (stringutils::isSpaces(s))
    s = "1, 1, 1";

  std::vector<int> modules;
  auto strvec = stringutils::splitString(s, ",");

  for (auto &w : strvec) {
    stringutils::trim(w);
    modules.push_back(std::stoi(w));
  }

  if (modules.size() != 3) {
    log::ferror("Number of modules must be specified in "
                "3 directions: %s", s);
  }

  if (modules[0] < 1 || modules[1] < 1 || modules[2] < 1) {
    log::ferror("Number of modules must be no less than 1: %s", s);
  }

  return modules;
}


/// \brief Returns the number of azimuthal angles
size_t ConfigInput::getNumAzim() {

  int num_azim = getOptionValueInt(_l_nazims);

  if (num_azim < 0) {
    log::ferror("Unable to set a negative number of azimuthal "
                "angles %d", num_azim);
  }

  if (num_azim % 4 != 0) {
    log::ferror("Unable to set # azimuthal angles to %d "
                "since it is not a multiple of 4", num_azim);
  }

  return static_cast<size_t>(num_azim);
}


/// \brief Returns the number of polar angles
size_t ConfigInput::getNumPolar() {

  int num_polar = getOptionValueInt(_l_npolars);

  if (num_polar < 0) {
    log::ferror("Unable to set a negative number of polar angles %d",
                num_polar);
  }

  if (num_polar % 2 != 0) {
    log::ferror("Unable to set # polar angles to %d "
                "since it is not a multiple of 2", num_polar);
  }

  return static_cast<size_t>(num_polar);
}


/// \brief Returns the spacing of azimuthal tracks
double ConfigInput::getAzimSpacing() {

  double azim_spacing = getOptionValueDouble(_l_sazim);

  if (azim_spacing < 0) {
    log::ferror("Unable to set a negative spacing '%f'"
                " for azimuthal tracks", azim_spacing);
  }

  return azim_spacing;
}


/// \brief Returns the spacing of axial tracks
double ConfigInput::getPolarSpacing() {

  double polar_spacing = getOptionValueDouble(_l_spolar);

  if (polar_spacing < 0) {
    log::ferror("Unable to set a negative spacing '%f'"
                " for axial tracks", polar_spacing);
  }

  return polar_spacing;
}


//----------------------------------------------------------------------
// Solver options
//----------------------------------------------------------------------

/// \brief Returns the solver type
solverType ConfigInput::getSolverType() {

  auto s = getOptionValue(_l_solver);

  // Keep the default
  if (stringutils::isSpaces(s))
    s = "CPU Solver";

  auto type = stringutils::spaceToUnderscore(s);
  return solverType::_from_string_nocase(type.c_str());

}


/// \brief Returns the quadrature type
quadratureType ConfigInput::getQuadratureType() {

  auto s = getOptionValue(_l_quad);

  // Keep the default
  if (stringutils::isSpaces(s))
    s = "Equal Angle";

  auto type = stringutils::spaceToUnderscore(s);
  return quadratureType::_from_string_nocase(type.c_str());

}


/// \brief Indicate the type of keff formulation
bool ConfigInput::getKeffFromNeutronBalance() {
  return getOptionValueBool(_l_keff);
}


/// \brief Returns the type of keff formulation
std::string ConfigInput::getKeffTypeStr() {
  if (getKeffFromNeutronBalance())
    return "Neutron Balance";
  else
    return "Successive Fission Rates";
}


/// \brief Returns the maximum number of iterations
int ConfigInput::getMaxIterations() {

  int max_iters = getOptionValueInt(_l_niters);

  if (max_iters < 0)
    log::ferror("Unable to set a negative number of maximum "
               "iterations %d", max_iters);

  return max_iters;
}


/// \brief Returns the tolerance of source iterations
double ConfigInput::getTolerance() {

  double tolerance = getOptionValueDouble(_l_tol);

  if (tolerance <= 0)
    log::ferror("Unable to set a non-positive tolerance: %e", tolerance);

  return tolerance;
}


/// \brief Returns the stabilization type
stabilizationType ConfigInput::getStabilizationType() {

  auto s = getOptionValue(_l_stab);

  // Keep the default
  if (stringutils::isSpaces(s))
    s = "NONE";

  return stabilizationType::_from_string_nocase(s.c_str());

}


bool ConfigInput::doesCheckXS() {
  return getOptionValueBool(_l_check_xs) || doesCorrectXS();
}


bool ConfigInput::doesCorrectXS() {
  return getOptionValueBool(_l_correct_xs);
}


log::level ConfigInput::getCheckXSLogLevel() {

  log::level e = log::level::error;
  auto s = getOptionValue(_l_check_xs_log_level);

  // Keep the default if there is an empty string
  if (!stringutils::isSpaces(s)) {
    try {
      e = log::level::_from_string_nocase(s.c_str());
    } catch (...) {
      log::error("Invalid log level for checking MGXS: {}", s);
    }
  }

  return e;
}


/// \brief Returns the stabilization factor
double ConfigInput::getStabilizationFactor() {

  double factor = getOptionValueDouble(_l_stab_factor);

  if (factor < 0)
    log::error("Unable to set a negative stabilization factor : {}", factor);

  return factor;
}


//----------------------------------------------------------------------
// Experimental options
//----------------------------------------------------------------------
/// \brief Returns the number of required energy groups
/// \details Empty string and positive numbers are valid. If the user
///          leave the choice to antmoc, 0 will be returned.
int ConfigInput::getRequiredEnergyGroups() {
  auto s = stringutils::trim(getOptionValue(_l_requiredgroups));

  if (stringutils::isSpaces(s))
    s = "ALL";

  // Take 0 as 'All'
  std::string str(s);
  stringutils::toUpper(str);
  if (str == "ALL")
    return 0;

  int num_groups = 0;
  try {
    num_groups = std::stoi(s);
    if (num_groups < 0)
      throw std::logic_error("the number of energy group must be non-negative");
  }
  catch (std::logic_error &e) {
    log::error("Invalid number of energy group '{}': {}", s, e.what());
  }

  return num_groups;
}


/// \brief Returns the number of required energy groups in string
std::string ConfigInput::getRequiredEnergyGroupsStr() {
  auto num_groups = getRequiredEnergyGroups();
  if (num_groups == 0)
    return "ALL";
  else
    return std::to_string(num_groups);
}


} // namespace antmoc

