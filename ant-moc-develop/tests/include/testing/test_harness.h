/// \file test_harness.h
/// \brief A test harness for ANT-MOC
/// \date Dec 10, 2019
/// \author An Wang, USTB (wangan@xs.ustb.edu.cn)

#ifndef TEST_HARNESS_H_
#define TEST_HARNESS_H_

#define STRINGIFY(s)           #s
#define STRINGIFY_MACRO(macro) STRINGIFY(macro)

#ifndef MATERIAL_DIR
#  define MATERIAL_DIR .
#endif
#ifndef GEOMETRY_DIR
#  define GEOMETRY_DIR .
#endif

#define TEST_WORKING_DIR STRINGIFY_MACRO(WORKING_DIR)
#define DEFAULT_MATERIAL_DIR STRINGIFY_MACRO(MATERIAL_DIR)
#define DEFAULT_GEOMETRY_DIR STRINGIFY_MACRO(GEOMETRY_DIR)

#include "sys/stat.h"
#include "unistd.h"

#include <algorithm>
#include <iterator>
#include <fstream>
#include <sstream>

#include "testing/test_utils.h"
#include "antmoc/container_utils.h"
#include "antmoc/enum_types.h"
#include "antmoc/file_utils.h"
#include "antmoc/log.h"
#include "antmoc/string_utils.h"
#include "antmoc/Factory.h"
#include "antmoc/FSRDataHandlerHDF5.h"
#include "antmoc/math_utils.h"
#include "antmoc/Mesh.h"

using namespace antmoc;


namespace {

///---------------------------------------------------------------------
/// \brief Result types
///---------------------------------------------------------------------
BETTER_ENUM(TestResultType, char,
  num_iterations,
  num_fsrs,
  num_chains,
  num_tracks,
  num_segments,
  keff,
  fluxes,
  mesh_data
)


///---------------------------------------------------------------------
/// \brief A test harness for regression testing.
/// \details In googletest, this is actually a test fixture.
///---------------------------------------------------------------------
class TestHarness : public testing::Test {

 public:

  /// \brief Initialize the test fixture.
  /// \details The constructor is where we initialize the options
  ///          handler, or we have to initialize them manually.
  TestHarness() {

    antmoc::log::set_level("info");

    // Set the directory of results
    // setTestDir("./");

    // Set acceptance criteria
    setTestCriteriaKeff(10);
    setTestCriteriaFlux(5e-3);

    // Overwrite results or not
    setTestOverwrite(true);

    // Default test result types
    setTestResultTypes({
      TestResultType::num_iterations,
      TestResultType::num_fsrs,
      TestResultType::num_tracks,
      TestResultType::num_chains,
      TestResultType::num_segments,
      TestResultType::keff,
      TestResultType::fluxes,
    // TestResultType::mesh_data
    });

  }


 protected:

  /// \brief Create materials, geometry, tracks and solver.
  virtual void SetUp() {
    // Reset auto ids to ensure the same behaviour of hashmaps
    antmoc::reset_auto_ids();
  }

  /// \brief Clean output files.
  virtual void TearDown() {

    ROOT_ONLY();

    // Delete results_test.dat if it exists
    std::string results = getDataFilePath("test");

    if (fileutils::existsFile(results)) {
      remove(results.c_str());
    }

  }

  /// \brief Run the test and compare results.
  /// \details This method is supposed to be called by the developer.
  /// \param test_suffix suffix of the filename for test results
  /// \param oracle_suffix suffix of the filename for oracle results
  virtual void runTest(std::string test_suffix = "",
                       std::string oracle_suffix = "") {

#ifdef ENABLE_MPI_
    mpi::mpiBarrier();
#endif

    auto &test_options = TestOptions::get();

    setTestSuffix(test_suffix);
    setOracleSuffix(oracle_suffix);

    modifyArguments("--log-level", test_options.getLogLevel());
    parseArguments();

#ifdef ENABLE_MPI_
    setDomainDecomposition();
#endif
    createMaterials();
    createGeometry();
    createTrackGenerator();
    generateTracks();
    createSolver();

    solve();

    if (mpi::isInMainCommUniqueDomains())
      dumpResults();

    if (test_options.visualize()) {
      dumpVTKData();
      dumpFSRData();
    }

    if (test_options.updateResults())
      overwriteResults();
    else
      compareResults();

#ifdef ENABLE_MPI_
    mpi::mpiBarrier();
#endif
  }


  //--------------------------------------------------------------------
  // Sets and reads TestHarness attributes
  //--------------------------------------------------------------------
  /// \brief Set the directory for results (relative to the working directory).
  /// \details In the current fixture, developers have to call this
  ///          method explicitly to start testing.
  void setTestDir(const std::string &dir) { _test_dir = dir; }
  std::string getTestDir()                { return getWorkingDir() + _test_dir + "/"; }

  /// \brief Return the working directory
  std::string getWorkingDir()        { return std::string(TEST_WORKING_DIR) + "/"; }
  /// \brief Return a direcory set by CMake
  std::string getPublicMaterialDir() { return std::string(DEFAULT_MATERIAL_DIR) + "/"; }
  /// \brief Return a direcory set by CMake
  std::string getPublicGeometryDir() { return std::string(DEFAULT_GEOMETRY_DIR) + "/"; }

  /// \brief Set the filename for test results.
  void setTestSuffix(const std::string &suffix) { _test_suffix = suffix; }
  std::string getTestSuffix()                   { return _test_suffix; }

  /// \brief Set the filename for oracle file.
  void setOracleSuffix(const std::string &suffix) { _oracle_suffix = suffix; }
  std::string getOracleSuffix()                   { return _oracle_suffix; }

  /// \brief Set the result types to be written into files
  void setTestResultTypes(std::set<TestResultType> types) {
    _test_result_types = types;
  }

  std::set<TestResultType> getTestResultTypes() { return _test_result_types; }

  /// \brief Set acceptance criteria for k-eff (+/- pcm)
  void setTestCriteriaKeff(unsigned criteria) { _test_criteria_keff = criteria; }
  unsigned getTestCriteriaKeff()              { return _test_criteria_keff; }

  /// \brief Set acceptance criteria for FSR flux (+/-)
  void setTestCriteriaFlux(FP_PRECISION criteria) { _test_criteria_flux = criteria; }
  FP_PRECISION getTestCriteriaFlux()              { return _test_criteria_flux; }

  /// \brief Set a flag for overwriting results
  void setTestOverwrite(bool overwrite) { _test_overwrite = overwrite; }
  bool getTestOverwrite()               { return _test_overwrite; }

  /// \brief Set default arguments for testing.
  /// \details This method initializes an object of ConfInput. Default
  ///          parameters are set and will be overwritten by user inputs.
  virtual void setDefaultArguments() {

    // Default arguments
    static const StringVec default_argv = {

      // By default, the number of threads is set to 1 to ensure
      // the order of FSRs.
      "--output-directory", "",
      //"--num-omp-threads",  "1",
      "--num-domains",      "1,1,1,0",

      "--materials",        getPublicMaterialDir() + "c5g7/mgxs.h5",
      "--geometry",         getTestDir() + "geometry.xml",
      "--dump-rx",          "F",
      "--tally-groups",     "0",

      "--num-azims",        "4",
      "--num-polars",       "2",
      "--azim-spacing",     "0.5",
      "--z-spacing",        "0.75",

      "--solver",           "CPU Solver",
      "--quadrature",       "Equal Weight",
      "--tolerance",        "1E-5",

    };

    // Initialize a runtime options handler.
    // All of the options are supposed to be hard-coded except geometry
    // files and the material file.
    if (_conf_in == nullptr) {
      _conf_in = Factory::getConfInput<ConfigInputFile>(default_argv);
    }

  }


  /// \brief Initialize default arguments and append new ones
  virtual void setArguments(const StringVec &argv) {
    setDefaultArguments();
    _conf_in->appendArguments(argv);
  }


  /// \brief Clear arguments and set them to new ones
  virtual void resetArguments(const StringVec &argv) {
    setDefaultArguments();
    _conf_in->setArguments(argv);
  }


  /// \brief Append a new argument or modify an existing argument
  virtual void modifyArguments(const std::string &name, const std::string &value) {
    setDefaultArguments();
    _conf_in->appendArguments({name, value});
  }

  /// \brief Append a list of arguments
  virtual void modifyArguments(const StringVec &argv) {
    setDefaultArguments();
    _conf_in->appendArguments(argv);
  }

  /// \brief Parse arguments
  virtual void parseArguments() {
    _conf_in->initializeLogger();

    // Reports
    log::test("Arguments passed in: {}", _conf_in->toString());
    _conf_in->printArgumentsReport();
  }


  //--------------------------------------------------------------------
  // Methods for driving tests
  //--------------------------------------------------------------------

  Factory::ConfInputPtr getConfigInput()            { return _conf_in; }
  Factory::MaterialHandlerPtr getMaterialHandler()  { return _mat_in; }
  Geometry *getGeometry()                           { return &_geometry; }
  Factory::TrackGeneratorPtr getTrackGenerator()    { return _track_generator; }
  Factory::SolverPtr getSolver()                    { return _solver; }


#ifdef ENABLE_MPI_
  /// \brief Set domain decomposition
  virtual void setDomainDecomposition() {
    auto comm = mpi::getMPIComm();
    auto domains = _conf_in->getNumDomains();
    mpi::setDomainDecomposition(comm, domains[0], domains[1], domains[2], domains[3]);
  }
#endif

  /// \brief Create materials
  virtual void createMaterials() {

    // Read materials from an HDF5 file and initialize a handler
    _mat_in = Factory::getMaterialHandler<MaterialHandlerHDF5>(_conf_in);
  }

  /// \brief Create the geometry
  /// \details This method initializes an object of Geometry and return
  ///          a smart pointer to it. The geometry is constructed from
  ///          an XML file rather than string. When the reading is done,
  ///          the geometry will be meshed.
  virtual void createGeometry() {

    // Initialize a handler to read the geometry from a file
    auto geo_in = Factory::getGeoInput<GeoInputXml>(&_geometry, _mat_in, _conf_in);

    // Set modules
    auto nmods = _conf_in->getNumModules();
    _geometry.setNumDomainModules(nmods[0], nmods[1], nmods[2]);

    // Read the geometry
    geo_in->readGlobalPrimitives(_conf_in->getGlobalPrimitivesPath());
    geo_in->readGeometryFromFile(_conf_in->getGeoInputPath());

    // Clean up unused primitives
    geo_in->eraseUnusedPrimitives();

    // Meshing
    _geometry.initializeFlatSourceRegions();

    // Clean up unused materials
    geo_in->eraseUnusedMaterials();

    // Print the number of materials
    geo_in->printReport();
  }

  /// \brief Create a quadrature and a track generator.
  /// \details This method initializes an object of Quadrature and an
  ///          object of TrackGenerator. The type of the TrackGenerator
  ///          object depends on user-defined inputs.
  virtual void createTrackGenerator() {

    // Initialize a quadrature
    auto quad = Factory::getQuadrature(_conf_in);

    // Initialize a track generator
    _track_generator = Factory::getTrackGenerator(&_geometry, _conf_in);
    _track_generator->setQuadrature(quad);

  }

  /// \brief Generate 2-D or 3-D tracks
  virtual void generateTracks() {
    _track_generator->generateTracks();
  }

  /// \brief Create a solver.
  /// \details This method initializes an object of Solver whose type
  ///          depends on user-defined inputs. An object of TrackGenerator
  ///          must be initialized before calling this method.
  virtual void createSolver() {
    _solver = Factory::getSolver(_track_generator, _conf_in);
  }

  /// \brief Call the solver
  virtual void solve() {
    _solver->computeEigenvalue(_conf_in->getMaxIterations());
  }


  //--------------------------------------------------------------------
  // Methods for manipulating results
  //--------------------------------------------------------------------

  /// \brief Returns mesh data in string
  std::string getMeshData();

  /// \brief Returns testing output filename
  /// \param root_name the root name of the file (test or oracle)
  virtual std::string getDataFileName(std::string root_name = "test") {
    std::string filename = "results_" + root_name;
    std::string suffix = (root_name == "oracle" ? _oracle_suffix : _test_suffix);
    if (stringutils::isSpaces(suffix))
      filename += ".dat";
    else
      filename += "_" + suffix + ".dat";

    return filename;
  }


  /// \brief Returns the path to the data file
  /// \param root_name the root name of the path (test or oracle)
  virtual std::string getDataFilePath(std::string root_name = "test") {
    return getTestDir() + "/" + getDataFileName(root_name);
  }

  /// \brief Dump the results in HDF5
  virtual void dumpResults();

  /// \brief Read and compare results in HDF5
  virtual void compareResults();

  /// \brief Compare a single array of reaction rates
  /// \details This method computes relative error on each FSR.
  /// \param name A descriptive name of the data array, usually the tally type.
  virtual void compareFSRRates(std::string name,
                               std::vector<FP_PRECISION> rates_test,
                               std::vector<FP_PRECISION> rates_oracle);

  /// \brief Compare statistices of reaction rates
  /// \details This method computes max, min, mean, and stddev.
  /// \param name A descriptive name of the data array, usually the tally type.
  virtual void compareFSRRatesStatistics(std::string name,
                                         std::vector<FP_PRECISION> rates_test,
                                         std::vector<FP_PRECISION> rates_oracle);

  /// \brief Overwrite results files
  virtual void overwriteResults() {

    ROOT_ONLY();

    if (_test_overwrite) {
      std::string file1(getDataFilePath("test"));
      std::string file2(getDataFilePath("oracle"));

      if (fileutils::existsFile(file1)) {
        rename(file1.c_str(), file2.c_str());
      }
    }
  }


  //--------------------------------------------------------------------
  // Utility
  //--------------------------------------------------------------------

  /// \brief Write results to XML files
  void dumpVTKData();

  /// \brief Write FSR data to H5 files
  void dumpFSRData();


private:

  std::string _test_dir;           ///< Directory of the current test
  std::string _test_suffix;        ///< User-defined suffix of output files
  std::string _oracle_suffix;      ///< User-defined suffix of output files
  std::set<TestResultType> _test_result_types;
  unsigned _test_criteria_keff;      ///< The acceptance criteria for k-eff
  FP_PRECISION _test_criteria_flux;  ///< The acceptance criteria for FSR flux

  bool _test_overwrite;            ///< A flag for overwriting results

  Factory::ConfInputPtr _conf_in;  ///< Arguments handler
  Factory::MaterialHandlerPtr _mat_in;    ///< Materials handler
  Geometry _geometry;              ///< Geometry
  Factory::TrackGeneratorPtr _track_generator; ///< Track generator
  Factory::SolverPtr _solver;      ///< Solver

};  // class TestHarness




//----------------------------------------------------------------------
// Inline functions
//----------------------------------------------------------------------
inline void TestHarness::dumpResults() {
  // Set booleans
  const auto &types = _test_result_types;
  bool keff         = types.find(TestResultType::keff) != types.end();
  bool num_iters    = types.find(TestResultType::num_iterations) != types.end();
  bool num_fsrs     = types.find(TestResultType::num_fsrs) != types.end();
  bool num_chains   = types.find(TestResultType::num_chains) != types.end();
  bool num_tracks   = types.find(TestResultType::num_tracks) != types.end();
  bool num_segments = types.find(TestResultType::num_segments) != types.end();
  bool fluxes       = types.find(TestResultType::fluxes) != types.end();
  bool mesh_data    = types.find(TestResultType::mesh_data) != types.end();

  // Determine the file path
  std::string file_path = getDataFilePath("test");
  FSRDataHandlerHDF5 h5_writer(file_path, HDF5Mode::Truncate, _solver, true);

  auto file_id = h5_writer.getFileId();

  // Write out the eigenvalue
  if (keff) {
    double k = _solver->getKeff();
    h5_writer.writeScalarAttribute(file_id, "keff", k);
  }

  // Write out the number of iterations
  if (num_iters) {
    auto n = (size_t)_solver->getNumIterations();
    h5_writer.writeScalarAttribute(file_id, "# Iterations", n);
  }

  // Write out the number of FSRs
  if (num_fsrs) {
    auto n = (size_t)_geometry.getNumFSRs();
    h5_writer.writeScalarAttribute(file_id, "# FSRs", n);
  }

  // Write out the number of 2D chains
  if (num_chains) {
    auto n = (size_t)_track_generator->getNum2DChains();
    h5_writer.writeScalarAttribute(file_id, "# 2D Chains", n);
  }

  // Write out the number of tracks
  if (num_tracks) {
    auto n = (size_t)_track_generator->getNumTracks();
    h5_writer.writeScalarAttribute(file_id, "# Tracks", n);
  }

  // Write out the number of segments
  if (num_segments) {
    auto n = (size_t)_track_generator->getNumSegments();
    h5_writer.writeScalarAttribute(file_id, "# Segments", n);
  }

  // Write out the scalar fluxes
  if (fluxes) {
    std::set<TallyType> field_types = {TallyType::Scalar_Flux};
    std::set<int> groups = _conf_in->getTallyEnergyGroups();
    h5_writer.writeFieldArrays(file_id, field_types, groups);
  }

  // Write out the tallied data
  if (mesh_data) {
  }
}


/// \details Entries are read and compared one by one. Some of them have
///          pre-defined criteria.
///          This method is supposed to be invoked by the root.
inline void TestHarness::compareResults() {
  ROOT_ONLY();

  auto file_path1 = getDataFilePath("test");
  auto file_path2 = getDataFilePath("oracle");

  if (!fileutils::existsFile(file_path1) || !fileutils::existsFile(file_path2)) {
    log::ferror("Results files may not exist: '%s' and '%s'",
                      file_path1, file_path2);
  }

  // Set booleans
  const auto &types = _test_result_types;
  bool keff         = types.find(TestResultType::keff) != types.end();
  bool num_iters    = types.find(TestResultType::num_iterations) != types.end();
  bool num_fsrs     = types.find(TestResultType::num_fsrs) != types.end();
  bool num_chains   = types.find(TestResultType::num_chains) != types.end();
  bool num_tracks   = types.find(TestResultType::num_tracks) != types.end();
  bool num_segments = types.find(TestResultType::num_segments) != types.end();
  bool fluxes       = types.find(TestResultType::fluxes) != types.end();
  bool mesh_data    = types.find(TestResultType::mesh_data) != types.end();

  // Initialize H5 readers
  FSRDataHandlerHDF5 h5_reader1(file_path1, HDF5Mode::ReadOnly, _solver, false);
  FSRDataHandlerHDF5 h5_reader2(file_path2, HDF5Mode::ReadOnly, _solver, false);

  auto file1_id = h5_reader1.getFileId();
  auto file2_id = h5_reader2.getFileId();

  // Compare the eigenvalue
  if (keff) {
    double keff_test, keff_oracle;
    h5_reader1.readScalarAttribute(file1_id, "keff", keff_test);
    h5_reader2.readScalarAttribute(file2_id, "keff", keff_oracle);

    unsigned delta_keff = static_cast<unsigned>(std::abs(keff_test - keff_oracle) * 1e5);

    EXPECT_LE(delta_keff, _test_criteria_keff)
      << "k-eff (test)   = " << keff_test << '\n'
      << "k-eff (oracle) = " << keff_oracle
      << std::endl;
  }

  // Compare the number of iterations
  if (num_iters) {
    unsigned int num_iters_test, num_iters_oracle;
    h5_reader1.readScalarAttribute(file1_id, "# Iterations", num_iters_test);
    h5_reader2.readScalarAttribute(file2_id, "# Iterations", num_iters_oracle);
    EXPECT_EQ(num_iters_test, num_iters_oracle);
  }

  // Compare the number of FSRs
  if (num_fsrs) {
    unsigned long num_fsrs_test, num_fsrs_oracle;
    h5_reader1.readScalarAttribute(file1_id, "# FSRs", num_fsrs_test);
    h5_reader2.readScalarAttribute(file2_id, "# FSRs", num_fsrs_oracle);
    EXPECT_EQ(num_fsrs_test, num_fsrs_oracle);
  }

  // Compare the number of 2D chains
  if (num_chains) {
    unsigned long num_chains_test, num_chains_oracle;
    h5_reader1.readScalarAttribute(file1_id, "# 2D Chains", num_chains_test);
    h5_reader2.readScalarAttribute(file2_id, "# 2D Chains", num_chains_oracle);
    EXPECT_EQ(num_chains_test, num_chains_oracle);
  }

  // Compare the number of tracks
  if (num_tracks) {
    unsigned long num_tracks_test, num_tracks_oracle;
    h5_reader1.readScalarAttribute(file1_id, "# Tracks", num_tracks_test);
    h5_reader2.readScalarAttribute(file2_id, "# Tracks", num_tracks_oracle);
    EXPECT_EQ(num_tracks_test, num_tracks_oracle);
  }

  // Compare the number of segments
  if (num_segments) {
    unsigned long num_segments_test, num_segments_oracle;
    h5_reader1.readScalarAttribute(file1_id, "# Segments", num_segments_test);
    h5_reader2.readScalarAttribute(file2_id, "# Segments", num_segments_oracle);
    EXPECT_EQ(num_segments_test, num_segments_oracle);
  }

  // Compare the scalar fluxes
  if (fluxes) {
    std::set<TallyType> field_types = {TallyType::Scalar_Flux};
    std::set<int> groups = _conf_in->getTallyEnergyGroups();
    auto fluxes_test   = h5_reader1.readFieldArrays(file1_id, field_types, groups);
    auto fluxes_oracle = h5_reader2.readFieldArrays(file2_id, field_types, groups);

    auto keys_test   = mapExtractKeys(fluxes_test);
    auto keys_oracle = mapExtractKeys(fluxes_oracle);
    std::sort(keys_test.begin(), keys_test.end());
    std::sort(keys_oracle.begin(), keys_oracle.end());

    // Keys of test results must be contained in keys of oracle results
    ASSERT_LE(keys_test, keys_oracle);

    for (const auto &key : keys_test) {
      ASSERT_TRUE(contains(keys_oracle, key));
      // Compare statistics
      compareFSRRatesStatistics(key, fluxes_test[key], fluxes_oracle[key]);
      // Compare reaction rates
      //compareFSRRates(key, fluxes_test[key], fluxes_oracle[key]);
    }
  }

  // Compare the tallied data
  if (mesh_data) {
  }
}


inline void TestHarness::compareFSRRates(std::string name,
                                         std::vector<FP_PRECISION> rates_test,
                                         std::vector<FP_PRECISION> rates_oracle) {
  ASSERT_EQ(rates_test.size(), rates_oracle.size());

  PyVector<FP_PRECISION> diff_vec;
  diff_vec.reserve(rates_test.size());

  FP_PRECISION max_rel_diff = 0;

  for (size_t r = 0; r < rates_test.size(); ++r) {

    const auto &test   = rates_test[r];
    const auto &oracle = rates_oracle[r];

    FP_PRECISION diff = test - oracle;

    // Spacial cases
    if (oracle == 0.) {
      if (test != 0.)
        diff = std::numeric_limits<FP_PRECISION>::infinity();
      else
        diff = 0.;
    } else {
      diff = std::abs(test/oracle - 1);
    }

    // Save the error for debugging
    diff_vec.push_back(diff);

    // Compute the maximum relative error
    max_rel_diff = std::max(max_rel_diff, diff);
  }

  EXPECT_LE(max_rel_diff, _test_criteria_flux)
    << "Comparing " << name << " by FSRs, relative errors:\n"
    << diff_vec.toString()
    << std::endl;
}


inline void TestHarness::compareFSRRatesStatistics(std::string name,
                                                   std::vector<FP_PRECISION> rates_test,
                                                   std::vector<FP_PRECISION> rates_oracle) {
  ASSERT_EQ(rates_test.size(), rates_oracle.size());

  // Statistics of the test
  auto rates_test_min    = mathutils::min(rates_test);
  auto rates_test_max    = mathutils::max(rates_test);
  auto rates_test_mean   = mathutils::mean(rates_test);
  auto rates_test_stddev = mathutils::stddev(rates_test);

  // Statistics of the oracle
  auto rates_oracle_min    = mathutils::min(rates_oracle);
  auto rates_oracle_max    = mathutils::max(rates_oracle);
  auto rates_oracle_mean   = mathutils::mean(rates_oracle);
  auto rates_oracle_stddev = mathutils::stddev(rates_oracle);

  // Relative error
  auto compute_relative_error = [](const FP_PRECISION &test, const FP_PRECISION &oracle) {
      FP_PRECISION diff = 0.;

      if (oracle != 0.) {
        diff = std::abs(test/oracle - 1);
      }
      else {
        // If both of the results are definitely 0, we let it pass.
        if (mathutils::definitelyEqual(test, oracle))
          diff = 0.;
        else
          diff = std::numeric_limits<FP_PRECISION>::infinity();
      }

      return diff;
    };

  auto rel_error_min    = compute_relative_error(rates_test_min, rates_oracle_min);
  auto rel_error_max    = compute_relative_error(rates_test_max, rates_oracle_max);
  auto rel_error_mean   = compute_relative_error(rates_test_mean, rates_oracle_mean);
  auto rel_error_stddev = compute_relative_error(rates_test_stddev, rates_oracle_stddev);

  EXPECT_LE(rel_error_min, _test_criteria_flux)
    << fmt::format("Comparing {}, relative error of min rates ({}, {}) is {}",
                    name, rates_test_min, rates_oracle_min, rel_error_min)
    << std::endl;

  EXPECT_LE(rel_error_max, _test_criteria_flux)
    << fmt::format("Comparing {}, relative error of max rates ({}, {}) is {}",
                    name, rates_test_max, rates_oracle_max, rel_error_max)
    << std::endl;

  EXPECT_LE(rel_error_mean, _test_criteria_flux)
    << fmt::format("Comparing {}, relative error of mean rates ({}, {}) is {}",
                    name, rates_test_mean, rates_oracle_mean, rel_error_mean)
    << std::endl;

  EXPECT_LE(rel_error_stddev, _test_criteria_flux)
    << fmt::format("Comparing {}, relative error of stddev rates ({}, {}) is {}",
                    name, rates_test_stddev, rates_oracle_stddev, rel_error_stddev)
    << std::endl;
}


// FIXME: use HDF5
inline std::string TestHarness::getMeshData() {
  Mesh mesh(_solver);

  auto widths = _conf_in->getTallyMesh();
  if (_conf_in->isHexTallyMesh()) {
    int num_r = (int) widths[0][0];
    double width_r = widths[1][0];
    mesh.createLattice(num_r, width_r, widths[2],
                       _conf_in->getTallyMeshOrientation(),
                       _conf_in->getTallyMeshOffset());
  } else {
    mesh.createLattice(widths[0], widths[1], widths[2]);
  }

  // Tallied types
  auto field_types = _conf_in->getDumpRXTypes();
  tallyutils::expandTallyTypeForRX(field_types);
  // Tallied groups
  std::set<int> groups = _conf_in->getTallyEnergyGroups();
  mesh.printMeshDataToXML("mesh", field_types, groups);

  return mesh.printDocToString("mesh");
}


inline void TestHarness::dumpVTKData() {

  Mesh mesh(_solver);

  auto widths = _conf_in->getTallyMesh();
  if (_conf_in->isHexTallyMesh()) {
    int num_r = (int) widths[0][0];
    double width_r = widths[1][0];
    mesh.createLattice(num_r, width_r, widths[2],
                       _conf_in->getTallyMeshOrientation(),
                       _conf_in->getTallyMeshOffset());
  } else {
    mesh.createLattice(widths[0], widths[1], widths[2]);
  }

  // Determine the file path
  std::string dump_dir = getTestDir();
  std::string file_path = dump_dir + "/reaction_rates" +
                          (stringutils::isSpaces(_test_suffix) ? "" : "_" + _test_suffix)
                          + ".vtu";

  // Dump mesh data
  auto field_types = _conf_in->getDumpRXTypes();
  tallyutils::expandTallyTypeForRX(field_types);
  mesh.dumpMeshDataToFile(file_path, field_types);

  // Dump 2D or 3D tracks
  auto track_types = _conf_in->getDumpTracksTypes();
  mesh.dumpTrackMeshToDir(_track_generator, track_types, dump_dir);
}


inline void TestHarness::dumpFSRData() {

  // Determine the file path
  std::string dump_dir = getTestDir();
  std::string file_path = dump_dir + "/fsr_data" +
                          (stringutils::isSpaces(_test_suffix) ? "" : "_" + _test_suffix)
                          + ".h5";

  // Reaction rates
  auto field_types = _conf_in->getDumpRXTypes();
  tallyutils::expandTallyTypeForRX(field_types);

  // Dump FSR data
  FSRDataHandlerHDF5 h5_writer(file_path, HDF5Mode::Truncate, _solver);
  h5_writer.dumpFSRData(field_types);

}

} // anonymous namespace

#endif
