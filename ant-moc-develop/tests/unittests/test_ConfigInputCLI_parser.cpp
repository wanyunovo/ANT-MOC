/**
 * @file test_ConfigInputCLI_parser
 * @brief Test ConfigInputCLI functions
 * @date June 8, 2019
 * @author An Wang, USTB (wangan@xs.ustb.edu.cn)
 */

#include "testing/test_utils.h"
#include "antmoc/log.h"
#include "antmoc/ConfigInputCLI.h"

#include <vector>
#include <string>


using namespace antmoc;

namespace {

// Test fixture
class test_ConfigInputCLI_parser : public testing::Test {
  // do nothing
};

TEST_F(test_ConfigInputCLI_parser, checkDefaultParameters) {
  ROOT_ONLY();

  ConfigInputCLI input;

  EXPECT_EQ("", input.getConfInputPath());
  EXPECT_EQ("", input.getGeoInputPath());
  EXPECT_EQ("", input.getGlobalPrimitivesPath());
  EXPECT_EQ("", input.getMatInputPath());

  EXPECT_EQ("ALL", input.getTallyEnergyGroupsStr());
  EXPECT_EQ("ALL", input.getRequiredEnergyGroupsStr());
}


TEST_F(test_ConfigInputCLI_parser, parseLongCommandLineOptions) {
  ROOT_ONLY();

  const std::vector<std::string> argv = {
    "--config", "config.xml",
    "--geometry", "geometry.xml",
    "--primitives", "geometry.xml",
    "--materials", "materials.h5",
  };

  ConfigInputCLI input(argv);

  EXPECT_EQ("config.xml", input.getConfInputPath());
  EXPECT_EQ("geometry.xml", input.getGeoInputPath());
  EXPECT_EQ("geometry.xml", input.getGlobalPrimitivesPath());
  EXPECT_EQ("materials.h5", input.getMatInputPath());
}


TEST_F(test_ConfigInputCLI_parser, parseShortCommandLineOptions) {
  ROOT_ONLY();

  const std::vector<std::string> argv = {
    "-c", "config.xml",
    "-g", "geometry.xml",
    "-m", "materials.h5",
    "-o", "output_dir",
  };

  ConfigInputCLI input(argv);

  EXPECT_EQ(input.getConfInputPath(), "config.xml");
  EXPECT_EQ(input.getGeoInputPath(), "geometry.xml");
  EXPECT_EQ(input.getMatInputPath(), "materials.h5");

}


/// \brief Testing log level
/// \details Log Level:
///             debug, profile, profile_once, verbose, verbose_once, info,
///             node, separator, header, title, warn, warn_once, critical,
///             result, test, error
TEST_F(test_ConfigInputCLI_parser, parseLogLevel) {
  ROOT_ONLY();

  // Object to be tested
  ConfigInputCLI input;

  // Normal functionality
  const std::vector<StringVec> tests = {
    {"--log-level", ""},
    {"--log-level=verbose"},
    {"--log-level", "debug"},
    {"--log-level", "profile"},
    {"--log-level", "profile_once"},
    {"--log-level", "verbose"},
    {"--log-level", "verbose_once"},
    {"--log-level", "info"},
    {"--log-level", "node"},
    {"--log-level", "warn"},
    {"--log-level", "warn_once"},
    {"--log-level", "critical"},
    {"--log-level", "result"},
    {"--log-level", "test"},
    {"--log-level", "error"},
  };

  const std::vector<log::level> oracles = {
    +log::level::info,
    +log::level::verbose,
    +log::level::debug,
    +log::level::profile,
    +log::level::profile_once,
    +log::level::verbose,
    +log::level::verbose_once,
    +log::level::info,
    +log::level::node,
    +log::level::warn,
    +log::level::warn_once,
    +log::level::critical,
    +log::level::result,
    +log::level::test,
    +log::level::error,
  };

  for (size_t i = 0; i < tests.size(); ++i) {
    input.setArguments(tests[i]);
    EXPECT_EQ(oracles[i], input.getLogLevel());
  }

#ifndef ENABLE_MPI_
  // Robust
  const std::vector<StringVec> fails = {
    {"--log-level", "Undefined"},
    {"--log-level", "1.0"},
  };

  for (size_t i = 0; i < fails.size(); ++i) {
    input.setArguments(fails[i]);
    EXPECT_ANY_THROW(input.getLogLevel())
      << "Case " << i << ": " << fails[i][0]
      << " " << fails[i][1] << std::endl;
  }
#endif
}


/// \brief Testing log file
/// \details Log file name: non-empty string
TEST_F(test_ConfigInputCLI_parser, parseLogFile) {
  ROOT_ONLY();

  // Object to be tested
  ConfigInputCLI input;

  // Normal functionality
  const std::vector<StringVec> tests = {
    {"-o", " "},  // whitespace will be ignored
    {"-o", "log"},
    {"--log-file", " "},  // whitespace will be ignored
    {"--log-file", "log/"},  // filename not specified
    {"--log-file", "log"},
    {"--log-file", "antmoc.log.out"},
    {"--log-file=antmoc.log.out"},
    {"--log-file=log/test/antmoc.log.out"},
    {"--log-file=/opt/log/antmoc.log.out"},
  };

  const StringVec oracles = {
    "",
    "log",
    "",
    "log/antmoc.0.out",
    "log",
    "antmoc.log.out",
    "antmoc.log.out",
    "log/test/antmoc.log.out",
    "/opt/log/antmoc.log.out",
  };

  for (size_t i = 0; i < tests.size(); ++i) {
    input.setArguments(tests[i]);
    EXPECT_EQ(oracles[i], input.getLogFile());
  }
}


/// \brief Testing output directory path
/// \details Output directory: non-empty string
TEST_F(test_ConfigInputCLI_parser, parseOutputDirectory) {
  ROOT_ONLY();

  // Object to be tested
  ConfigInputCLI input;

  // Normal functionality
  const std::vector<StringVec> tests = {
    {"-d", "dir"},
    {"-d", " "},  // whitespace will be ignored
    {"--output-directory", "output.d"},
    {"--output-directory", ""},
    {"--output-directory=output.d"},
  };

  const StringVec oracles = {
    "dir",
    "",
    "output.d",
    "",
    "output.d",
  };

  for (size_t i = 0; i < tests.size(); ++i) {
    input.setArguments(tests[i]);
    EXPECT_EQ(oracles[i], input.getOutputDirectory());
  }
}


/// \brief Testing the number of azimuthal angles
/// \details # of azims: integer, multiples of 4
TEST_F(test_ConfigInputCLI_parser, parseNumAzim) {
  ROOT_ONLY();

  // Object to be tested
  ConfigInputCLI input;

  // Normal functionality
  const std::vector<StringVec> tests = {
    {"-a", "0"},
    {"-a", "4"},
    {"-a", "12"},
    {"--num-azims", "64"},
    {"--num-azims=64"},
  };

  const std::vector<int> oracles = {
    0,
    4,
    12,
    64,
    64,
  };

  for (size_t i = 0; i < tests.size(); ++i) {
    input.setArguments(tests[i]);
    EXPECT_EQ(oracles[i], input.getNumAzim());
  }

#ifndef ENABLE_MPI_
  // Robust
  const std::vector<StringVec> fails = {
    {"-a", "-1"},
    {"-a", "7"},
    {"-a", " "},
    {"-a", ""},
  };

  for (size_t i = 0; i < fails.size(); ++i) {
    input.setArguments(fails[i]);
    EXPECT_ANY_THROW(input.getNumAzim())
      << "Case " << i << ": " << fails[i][0]
      << " " << fails[i][1] << std::endl;
  }
#endif
}


TEST_F(test_ConfigInputCLI_parser, parseNumPolar) {
  ROOT_ONLY();

  // Object to be tested
  ConfigInputCLI input;

  // Normal functionality
  const std::vector<StringVec> tests = {
    {"-p", "0"},
    {"-p", "2"},
    {"-p", "6"},
    {"--num-polars", "14"},
    {"--num-polars=14"},
  };

  const std::vector<int> oracles = {
    0,
    2,
    6,
    14,
    14,
  };

  for (size_t i = 0; i < tests.size(); ++i) {
    input.setArguments(tests[i]);
    EXPECT_EQ(oracles[i], input.getNumPolar());
  }

#ifndef ENABLE_MPI_
  // Robust
  const std::vector<StringVec> fails = {
    {"-p", "-1"},
    {"-p", "7"},
    {"-p", " "},
    {"-p", ""},
  };

  for (size_t i = 0; i < fails.size(); ++i) {
    input.setArguments(fails[i]);
    EXPECT_ANY_THROW(input.getNumPolar())
      << "Case " << i << ": " << fails[i][0]
      << " " << fails[i][1] << std::endl;
  }
#endif
}


TEST_F(test_ConfigInputCLI_parser, parseAzimSpacing) {
  ROOT_ONLY();

  // Object to be tested
  ConfigInputCLI input;

  // Normal functionality
  const std::vector<StringVec> tests = {
    {"-s", "0"},
    {"-s", "1."},
    {"-s", ".1"},
    {"-s", "1.1"},
    {"--azim-spacing", "0.5"},
    {"--azim-spacing=0.5"},
  };

  const std::vector<double> oracles = {
    0,
    1.,
    .1,
    1.1,
    0.5,
    0.5,
  };

  for (size_t i = 0; i < tests.size(); ++i) {
    input.setArguments(tests[i]);
    EXPECT_EQ(oracles[i], input.getAzimSpacing());
  }

#ifndef ENABLE_MPI_
  // Robust
  const std::vector<StringVec> fails = {
    {"-s", "-1"},
    {"-s", " "},
    {"-s", ""},
  };

  for (size_t i = 0; i < fails.size(); ++i) {
    input.setArguments(fails[i]);
    EXPECT_ANY_THROW(input.getAzimSpacing())
      << "Case " << i << ": " << fails[i][0]
      << " " << fails[i][1] << std::endl;
  }
#endif
}


TEST_F(test_ConfigInputCLI_parser, parsePolarSpacing) {
  ROOT_ONLY();

  // Object to be tested
  ConfigInputCLI input;

  // Normal functionality
  const std::vector<StringVec> tests = {
    {"-l", "0"},
    {"-l", "1."},
    {"-l", ".1"},
    {"-l", "1.1"},
    {"--z-spacing", "0.5"},
    {"--z-spacing=0.5"},
  };

  const std::vector<double> oracles = {
    0,
    1.,
    .1,
    1.1,
    0.5,
    0.5,
  };

  for (size_t i = 0; i < tests.size(); ++i) {
    input.setArguments(tests[i]);
    EXPECT_EQ(oracles[i], input.getPolarSpacing());
  }

#ifndef ENABLE_MPI_
  // Robust
  const std::vector<StringVec> fails = {
    {"-l", "-1"},
    {"-l", " "},
    {"-l", ""},
  };

  for (size_t i = 0; i < fails.size(); ++i) {
    input.setArguments(fails[i]);
    EXPECT_ANY_THROW(input.getPolarSpacing())
      << "Case " << i << ": " << fails[i][0]
      << " " << fails[i][1] << std::endl;
  }
#endif
}


TEST_F(test_ConfigInputCLI_parser, parseMaxIterations) {
  ROOT_ONLY();

  // Object to be tested
  ConfigInputCLI input;

  // Normal functionality
  const std::vector<StringVec> tests = {
    {"-i", "0"},
    {"-i", "1"},
    {"--max-iters", "1000"},
    {"--max-iters=1000"},
  };

  const std::vector<int> oracles = {
    0,
    1,
    1000,
    1000,
  };

  for (size_t i = 0; i < tests.size(); ++i) {
    input.setArguments(tests[i]);
    EXPECT_EQ(oracles[i], input.getMaxIterations());
  }


#ifndef ENABLE_MPI_
  // Robust
  const std::vector<StringVec> fails = {
    {"-i", "-1"},
    {"-i", " "},
    {"-i", ""},
  };

  for (size_t i = 0; i < fails.size(); ++i) {
    input.setArguments(fails[i]);
    EXPECT_ANY_THROW(input.getMaxIterations())
      << "Case " << i << ": " << fails[i][0]
      << " " << fails[i][1] << std::endl;
  }
#endif
}


TEST_F(test_ConfigInputCLI_parser, parseTolerance) {
  ROOT_ONLY();

  // Object to be tested
  ConfigInputCLI input;

  // Normal functionality
  const std::vector<StringVec> tests = {
    {"--tolerance", "1"},
    {"--tolerance", "1."},
    {"--tolerance", "1E-0"},
    {"--tolerance", "1E-1"},
    {"--tolerance=1e-5"},
  };

  const std::vector<double> oracles = {
    1,
    1.,
    1E-0,
    1E-1,
    1E-5,
  };

  for (size_t i = 0; i < tests.size(); ++i) {
    input.setArguments(tests[i]);
    EXPECT_EQ(oracles[i], input.getTolerance());
  }

#ifndef ENABLE_MPI_
  // Robust
  const std::vector<StringVec> fails = {
    {"--tolerance", "-1"},
    {"--tolerance", "0"},
    {"--tolerance", " "},
    {"--tolerance", ""},
  };

  for (size_t i = 0; i < fails.size(); ++i) {
    input.setArguments(fails[i]);
    EXPECT_ANY_THROW(input.getTolerance())
      << "Case " << i << ": " << fails[i][0]
      << " " << fails[i][1] << std::endl;
  }
#endif
}


TEST_F(test_ConfigInputCLI_parser, parseSegmentationType) {
  ROOT_ONLY();

  // Object to be tested
  ConfigInputCLI input;

  // Normal functionality
  const std::vector<StringVec> tests = {
    {"-f", ""}, // empty
    {"-f", "Explicit 2D"},
    {"-f", "Explicit 3D"},
    {"-f", "OTF tracks"},
    {"-f", "OTF stacks"},
    {"-f", "OTF_stacks"},
    {"--formation", "OTF tracks"},
    {"--formation=OTF tracks"},
  };

  const std::vector<segmentationType> oracles = {
    +segmentationType::OTF_STACKS,
    +segmentationType::EXPLICIT_2D,
    +segmentationType::EXPLICIT_3D,
    +segmentationType::OTF_TRACKS,
    +segmentationType::OTF_STACKS,
    +segmentationType::OTF_STACKS,
    +segmentationType::OTF_TRACKS,
    +segmentationType::OTF_TRACKS,
  };

  for (size_t i = 0; i < tests.size(); ++i) {
    input.setArguments(tests[i]);
    EXPECT_EQ(oracles[i], input.getSegmentationType());
  }

#ifndef ENABLE_MPI_
  // Robust
  const std::vector<StringVec> fails = {
    {"-f", "Undefined"},
    {"-f", "1.0"},
  };

  for (size_t i = 0; i < fails.size(); ++i) {
    input.setArguments(fails[i]);
    EXPECT_ANY_THROW(input.getSegmentationType())
      << "Case " << i << ": " << fails[i][0]
      << " " << fails[i][1] << std::endl;
  }
#endif
}


TEST_F(test_ConfigInputCLI_parser, parseQuadratureType) {
  ROOT_ONLY();

  // Object to be tested
  ConfigInputCLI input;

  // Normal functionality
  const std::vector<StringVec> tests = {
    {"-q", ""},
    {"-q", "Tabuchi Yamamoto"},
    {"-q", "tabuchi yamamoto"},
    {"-q", "leonard"},
    {"-q", "gauss legendre"},
    {"-q", "equal weight"},
    {"-q", "equal angle"},
    {"-q", "equal_angle"},
    {"--quadrature", "Equal Weight"},
    {"--quadrature=Equal Weight"},
  };

  const std::vector<quadratureType> oracles = {
    +quadratureType::EQUAL_ANGLE,
    +quadratureType::TABUCHI_YAMAMOTO,
    +quadratureType::TABUCHI_YAMAMOTO,
    +quadratureType::LEONARD,
    +quadratureType::GAUSS_LEGENDRE,
    +quadratureType::EQUAL_WEIGHT,
    +quadratureType::EQUAL_ANGLE,
    +quadratureType::EQUAL_ANGLE,
    +quadratureType::EQUAL_WEIGHT,
    +quadratureType::EQUAL_WEIGHT,
  };

  for (size_t i = 0; i < tests.size(); ++i) {
    input.setArguments(tests[i]);
    EXPECT_EQ(oracles[i], input.getQuadratureType());
  }

#ifndef ENABLE_MPI_
  // Robust
  const std::vector<StringVec> fails = {
    {"-q", "Undefined"},
    {"-q", "1.0"},
  };

  for (size_t i = 0; i < fails.size(); ++i) {
    input.setArguments(fails[i]);
    EXPECT_ANY_THROW(input.getQuadratureType())
      << "Case " << i << ": " << fails[i][0]
      << " " << fails[i][1] << std::endl;
  }
#endif
}


TEST_F(test_ConfigInputCLI_parser, parseSolverType) {
  ROOT_ONLY();

  // Object to be tested
  ConfigInputCLI input;

  // Normal functionality
  const std::vector<StringVec> tests = {
    {"--solver", ""},
    {"--solver", "CPU Solver"},
    {"--solver", "cpu solver"},
    {"--solver", "CPU LS Solver"},
    {"--solver", "GPU Solver"},
    {"--solver", "GPU_Solver"},
    {"--solver=CPU Solver"},
  };

  const std::vector<solverType> oracles = {
    +solverType::CPU_SOLVER,
    +solverType::CPU_SOLVER,
    +solverType::CPU_SOLVER,
    +solverType::CPU_LS_SOLVER,
    +solverType::GPU_SOLVER,
    +solverType::GPU_SOLVER,
    +solverType::CPU_SOLVER,
  };

  for (size_t i = 0; i < tests.size(); ++i) {
    input.setArguments(tests[i]);
    EXPECT_EQ(oracles[i], input.getSolverType());
  }

#ifndef ENABLE_MPI_
  // Robust
  const std::vector<StringVec> fails = {
    {"--solver", "Undefined"},
    {"--solver", "1.0"},
  };

  for (size_t i = 0; i < fails.size(); ++i) {
    input.setArguments(fails[i]);
    EXPECT_ANY_THROW(input.getSolverType())
      << "Case " << i << ": " << fails[i][0]
      << " " << fails[i][1] << std::endl;
  }
#endif
}


TEST_F(test_ConfigInputCLI_parser, parseXSFileLayout) {
  ROOT_ONLY();

  // Object to be tested
  ConfigInputCLI input;

  // Normal functionality
  const std::vector<StringVec> tests = {
    {"--xs-layout", ""},
    {"--xs-layout", "named"},
    {"--xs-layout=named"},
    {"--xs-layout", "compressed"},
    {"--xs-layout=compressed"},
  };

  const std::vector<XSFileLayout> oracles = {
    +XSFileLayout::NAMED,
    +XSFileLayout::NAMED,
    +XSFileLayout::NAMED,
    +XSFileLayout::COMPRESSED,
    +XSFileLayout::COMPRESSED,
  };

  for (size_t i = 0; i < tests.size(); ++i) {
    input.setArguments(tests[i]);
    EXPECT_EQ(oracles[i], input.getXSFileLayout());
  }

#ifndef ENABLE_MPI_
  // Robust
  const std::vector<StringVec> fails = {
    {"--xs-layout", "Undefined"},
    {"--xs-layout", "1.0"},
  };

  for (size_t i = 0; i < fails.size(); ++i) {
    input.setArguments(fails[i]);
    EXPECT_ANY_THROW(input.getXSFileLayout())
      << "Case " << i << ": " << fails[i][0]
      << " " << fails[i][1] << std::endl;
  }
#endif
}


TEST_F(test_ConfigInputCLI_parser, parseNumModules) {
  ROOT_ONLY();

  ConfigInputCLI input;

  const std::vector<StringVec> tests = {
    {"--modules", ""},
    {"--modules", "1,2,3"},
    {"--modules=3,3,3"},
  };

  const std::vector<std::vector<int>> oracles = {
    {1, 1, 1},
    {1, 2, 3},
    {3, 3, 3},
  };

  for (size_t i = 0; i < tests.size(); ++i) {
    input.setArguments(tests[i]);
    EXPECT_EQ(oracles[i], input.getNumModules())
      << "Case " << i << std::endl;
  }

#ifndef ENABLE_MPI_
  // Robust
  const std::vector<StringVec> fails = {
    {"--modules", "Undefined"},
    {"--modules", "0, 0, 0"},
    {"--modules", "0, 1, 1"},
    {"--modules", "1, 0, 0"},
    {"--modules", "3"},
    {"--modules", "3,  2"},
    {"--modules", "3, -2, 1"},
    {"--modules", "3,  2, a"},
    {"--modules", "3,  2, 1, 1"},
  };

  for (size_t i = 0; i < fails.size(); ++i) {
    input.setArguments(fails[i]);
    EXPECT_ANY_THROW(input.getNumModules())
      << "Case " << i << ": " << fails[i][0]
      << " " << fails[i][1] << std::endl;
  }
#endif
}


TEST_F(test_ConfigInputCLI_parser, parseLatticeRefines) {
  ROOT_ONLY();

  ConfigInputCLI input;

  const std::vector<StringVec> tests = {
    {"--lattice-refines", ""},
    {"--lattice-refines", "1,2,3"},
    {"--lattice-refines", "0, 0, 0"},
    {"--lattice-refines", "0, 1, 1"},
    {"--lattice-refines", "1, 0, 0"},
    {"--lattice-refines", "3.4, 2, 1"},// will be truncated
    {"--lattice-refines=5,5,5"},
  };

  const std::vector<std::vector<int>> oracles = {
    {0, 0, 0},
    {1, 2, 3},
    {0, 0, 0},
    {0, 1, 1},
    {1, 0, 0},
    {3, 2, 1},
    {5, 5, 5},
  };

  for (size_t i = 0; i < tests.size(); ++i) {
    input.setArguments(tests[i]);
    EXPECT_EQ(oracles[i], input.getLatticeRefines())
      << "Case " << i << std::endl;
  }

#ifndef ENABLE_MPI_
  // Robust
  const std::vector<StringVec> fails = {
    {"--lattice-refines", "Undefined"},
    {"--lattice-refines", "3"},
    {"--lattice-refines", "3, 2"},
    {"--lattice-refines", "3, -2, 1"},
    {"--lattice-refines", "3,  2, a"},
    {"--lattice-refines", "3,  2, 1, 1"},
  };

  for (size_t i = 0; i < fails.size(); ++i) {
    input.setArguments(fails[i]);
    EXPECT_ANY_THROW(input.getLatticeRefines())
      << "Case " << i << ": " << fails[i][0]
      << " " << fails[i][1] << std::endl;
  }
#endif
}


TEST_F(test_ConfigInputCLI_parser, parseCellSectors) {
  ROOT_ONLY();

  ConfigInputCLI input;

  const std::vector<StringVec> tests = {
    {"--cell-sectors", ""},
    {"--cell-sectors", "0"},
    {"--cell-sectors", "1"},
    {"--cell-sectors", "10"},
    {"--cell-sectors", "10.7"}, // will be truncated
    {"--cell-sectors=5"},
  };

  const std::vector<int> oracles = {
    -1,
    0,
    1,
    10,
    10,
    5,
  };

  for (size_t i = 0; i < tests.size(); ++i) {
    input.setArguments(tests[i]);
    EXPECT_EQ(oracles[i], input.getCellSectors())
      << "Case " << i << std::endl;
  }

#ifndef ENABLE_MPI_
  // Robust
  const std::vector<StringVec> fails = {
    {"--cell-sectors", "Undefined"},
    {"--cell-sectors", "-1"},
  };

  for (size_t i = 0; i < fails.size(); ++i) {
    input.setArguments(fails[i]);
    EXPECT_ANY_THROW(input.getCellSectors())
      << "Case " << i << ": " << fails[i][0]
      << " " << fails[i][1] << std::endl;
  }
#endif
}


/// \brief Testing --z-mesh
TEST_F(test_ConfigInputCLI_parser, parseZMesh) {
  ROOT_ONLY();

  // Object to be tested
  ConfigInputCLI input;

  // Normal functionality
  const std::vector<StringVec> tests = {
    {"--z-mesh", ""},
    {"--z-mesh", "Local"},
    {"--z-mesh", "global"},
  };

  const std::vector<bool> oracles = {
    false,
    false,
    true,
  };

  for (size_t i = 0; i < tests.size(); ++i) {
    input.setArguments(tests[i]);
    EXPECT_EQ(oracles[i], input.isGlobalZMesh())
      << "Case " << i << std::endl;
  }
}


/// \brief Testing --zones
TEST_F(test_ConfigInputCLI_parser, parseSegmentationZones) {
  ROOT_ONLY();

  ConfigInputCLI input;

  const std::vector<StringVec> tests = {
    {"--zones", ""},
    {"--zones", "auto"},
    {"--zones", "Auto"},
    {"-z", "-0.5, 3.5"},
    {"-z", "-10, -5, 0, 5.5, 10"},
    {"--zones", "-1,1"},
    {"--zones", "-7.14, 0, 7.14"},
    {"--zones=-7.14, 0, 7.14"},
  };

  const std::vector<std::vector<FP_PRECISION>> oracles = {
    std::vector<FP_PRECISION>(),  // empty
    std::vector<FP_PRECISION>(),  // empty
    std::vector<FP_PRECISION>(),  // empty
    {-0.5, 3.5},
    {-10, -5, 0, 5.5, 10},
    {-1, 1},
    {-7.14, 0, 7.14},
    {-7.14, 0, 7.14},
  };

  for (size_t i = 0; i < tests.size(); ++i) {
    input.setArguments(tests[i]);
    EXPECT_EQ(oracles[i], input.getZones())
      << "Case " << i << std::endl;
  }

#ifndef ENABLE_MPI_
  // Robust
  const std::vector<StringVec> fails = {
    {"--zones", "Undefined"},
    {"--zones", "7.14"},
    {"--zones", "3, 2, a"},
  };

  for (size_t i = 0; i < fails.size(); ++i) {
    input.setArguments(fails[i]);
    EXPECT_ANY_THROW(input.getZones())
      << "Case " << i << ": " << fails[i][0]
      << " " << fails[i][1] << std::endl;
  }
#endif
}


/// \brief Testing tally mesh type
TEST_F(test_ConfigInputCLI_parser, parseTallyMeshType) {
  ROOT_ONLY();

  // Object to be tested
  ConfigInputCLI input;

  // Normal functionality
  const std::vector<StringVec> tests = {
    {"--mesh-type", ""},
    {"--mesh-type", "Rectangle"},
    {"--mesh-type", "hexagon"},
    {"--mesh-type=hexagon"},
  };

  const std::vector<tallyMeshType> oracles = {
    +tallyMeshType::RECTANGLE,
    +tallyMeshType::RECTANGLE,
    +tallyMeshType::HEXAGON,
    +tallyMeshType::HEXAGON,
  };

  for (size_t i = 0; i < tests.size(); ++i) {
    input.setArguments(tests[i]);
    EXPECT_EQ(oracles[i], input.getTallyMeshType())
      << "Case " << i << std::endl;
  }
}


/// \brief Testing --orientation for hexagonal tally mesh
TEST_F(test_ConfigInputCLI_parser, parseTallyMeshOrientation) {
  ROOT_ONLY();

  // Object to be tested
  ConfigInputCLI input;

  // Normal functionality
  const std::vector<StringVec> tests = {
    {"--orientation", ""},
    {"--orientation", "X"},
    {"--orientation", "Y"},
  };

  const std::vector<std::string> oracles = {
    "Y",
    "X",
    "Y",
  };

  for (size_t i = 0; i < tests.size(); ++i) {
    input.setArguments(tests[i]);
    EXPECT_EQ(oracles[i], input.getTallyMeshOrientation())
      << "Case " << i << std::endl;
  }
}


/// \brief Testing --mesh
TEST_F(test_ConfigInputCLI_parser, parseTallyMesh) {
  ROOT_ONLY();

  // Object to be tested
  ConfigInputCLI input;

  const std::vector<StringVec> tests = {
    {"--mesh", ""},
    {"--mesh", "[1.26]*3, 1.26, [7.14]*3"},
    {"--mesh", "5.0, 5.0, 5"},
    {"--mesh", "5 * 5, 5.0, 5 * 5"},  // Bad style
  };

  const std::vector<std::vector<WidthVec>> oracles = {
    std::vector<WidthVec>(),  // empty
    {
      WidthVec{1.26} * 3,
      WidthVec{1.26},
      WidthVec{7.14} * 3
    },
    {
      {5.0}, {5.0}, {5.0}
    },
    {
      {25.}, {5.}, {25.}
    }
  };

  for (size_t i = 0; i < tests.size(); ++i) {
    input.setArguments(tests[i]);
    auto &oracle = oracles[i];
    auto widths = input.getTallyMesh();
    ASSERT_EQ(oracle.size(), widths.size())
      << " Case " << i << "; " << tests[i][0]
      << ' ' << tests[i][1] << '\n';

    for (size_t j = 0; j < oracle.size(); ++j)
      EXPECT_EQ(oracle[j], widths[j]);
  }

#ifndef ENABLE_MPI_
  // Robust
  const std::vector<StringVec> fails = {
    {"--mesh", "Undefined"},
    {"--mesh", "1.26"},
    {"--mesh", "1.26,"},
    {"--mesh", "1.26, 1.26"},
    {"--mesh", "1.26, 1.26,"},
    {"--mesh", "1.26 1.26"},
  };

  for (size_t i = 0; i < fails.size(); ++i) {
    input.setArguments(fails[i]);
    EXPECT_ANY_THROW(input.getTallyMesh())
      << "Case " << i << ": " << fails[i][0]
      << " " << fails[i][1] << std::endl;
  }
#endif
}


/// \brief Testing --dump-rx
/// \details Note that this test only take effects on the parser,
///          which means TallyType::All will not be expanded.
TEST_F(test_ConfigInputCLI_parser, parseDumpRXTypes) {
  ROOT_ONLY();

  // Object to be tested
  ConfigInputCLI input;

  // Normal functionality
  const std::vector<StringVec> tests = {
    {"--dump-rx", "f"},
    {"--dump-rx", "NuF, A"},
    {"--dump-rx", "t,phi"},
    {"--dump-rx", "all, none"},
  };

  const std::vector<std::set<TallyType>> oracles = {
    {TallyType::Fission_RX},
    {TallyType::NuFission_RX, TallyType::Absorption_RX},
    {TallyType::Total_RX, TallyType::Scalar_Flux},
    {TallyType::All, TallyType::None},
  };

  for (size_t i = 0; i < tests.size(); ++i) {
    input.setArguments(tests[i]);
    auto types = input.getDumpRXTypes();  // Get a set of TallyType objects
    EXPECT_EQ(oracles[i].size(), types.size());
    for (auto t : types)  // Check if two sets are equal
      EXPECT_TRUE(oracles[i].find(t) != oracles[i].end())
        << "Input RX types: " << tallyutils::join(input.getDumpRXTypes(), ", ") << '\n'
        << "Converted into: " << tallyutils::getTallyTypeName(t)
        << std::endl;
  }
}


/// \brief Testing --dump-xs
/// \details Note that this test only take effects on the parser,
///          which means TallyType::All will not be expanded.
TEST_F(test_ConfigInputCLI_parser, parseDumpXSTypes) {
  ROOT_ONLY();

  // Object to be tested
  ConfigInputCLI input;

  // Normal functionality
  const std::vector<StringVec> tests = {
    {"--dump-xs", "f"},
    {"--dump-xs", "NuF, A, t"},
    {"--dump-xs", "all, none"},
  };

  const std::vector<std::set<TallyType>> oracles = {
    {TallyType::Fission_XS},
    {TallyType::NuFission_XS, TallyType::Absorption_XS, TallyType::Total_XS},
    {TallyType::All, TallyType::None},
  };

  for (size_t i = 0; i < tests.size(); ++i) {
    input.setArguments(tests[i]);
    auto types = input.getDumpXSTypes();
    EXPECT_EQ(oracles[i].size(), types.size());
    for (auto t : types)
      EXPECT_TRUE(oracles[i].find(t) != oracles[i].end())
        << "XS types: " << tallyutils::join(input.getDumpXSTypes(), ", ") << '\n'
        << "Converted into: " << tallyutils::getTallyTypeName(t)
        << std::endl;
  }
}


/// \brief Testing --dump-tracks
TEST_F(test_ConfigInputCLI_parser, parseDumpTracksTypes) {
  ROOT_ONLY();

  // Object to be tested
  ConfigInputCLI input;

  // Normal functionality
  const std::vector<StringVec> tests = {
    {"--dump-tracks", "all, none"},
    {"--dump-tracks", "2d, 3d"},
  };

  const std::vector<std::set<TallyType>> oracles = {
    {TallyType::All, TallyType::None},
    {TallyType::Tracks_2D, TallyType::Tracks_3D},
  };

  for (size_t i = 0; i < tests.size(); ++i) {
    input.setArguments(tests[i]);
    EXPECT_EQ(input.getDumpTracksTypes(), oracles[i]);
  }
}


/// \brief Testing --tally-groups
TEST_F(test_ConfigInputCLI_parser, parseTallyGroups) {
  ROOT_ONLY();

  // Object to be tested
  ConfigInputCLI input;

  // Normal functionality
  const std::vector<StringVec> tests = {
    {"--tally-groups", "1,2,3,4,5,6,7,8"},
    {"--tally-groups", "1,2:8"},
    {"--tally-groups", "2:6:2"},
    {"--tally-groups=2:6:2"},
  };

  const std::vector<std::set<int>> oracles = {
    {1,2,3,4,5,6,7,8},
    {1,2,3,4,5,6,7,8},
    {2,4,6},
    {2,4,6},
  };

  for (size_t i = 0; i < tests.size(); ++i) {
    input.setArguments(tests[i]);
    EXPECT_EQ(input.getTallyEnergyGroups(), oracles[i])
      << "Case " << i << std::endl;
  }
}


/// \brief Testing boolean option --keff-neutron-balance
TEST_F(test_ConfigInputCLI_parser, parseKeffNeutronBalance) {
  ROOT_ONLY();

  // Object to be tested
  ConfigInputCLI input;

  // Normal functionality
  const std::vector<StringVec> tests = {
    {""},
    {"--keff-neutron-balance"},
    {"--keff-neutron-balance=true"},
    {"--keff-neutron-balance=false"},
    // expects a true because it is a boolean, which can be assigned to by '='
    {"--keff-neutron-balance", "false"},
  };

  const std::vector<bool> oracles = {
    false,
    true,
    true,
    false,
    true,
  };

  for (size_t i = 0; i < tests.size(); ++i) {
    input.setArguments(tests[i]);
    EXPECT_EQ(input.getKeffFromNeutronBalance(), oracles[i])
      << "Case " << i << std::endl;
  }
}


/// \brief Testing --required-groups
TEST_F(test_ConfigInputCLI_parser, parseRequiredGroups) {
  ROOT_ONLY();

  // Object to be tested
  ConfigInputCLI input;

  // Normal functionality
  const std::vector<StringVec> tests = {
    {"--required-groups", "ALL"},
    {"--required-groups", "all"},
    {"--required-groups", " "},
    {"--required-groups", "0"},
    {"--required-groups", "1"},
    {"--required-groups", "1000"},
  };

  const std::vector<int> oracles = {
    0,
    0,
    0,
    0,
    1,
    1000,
  };

  for (size_t i = 0; i < tests.size(); ++i) {
    input.setArguments(tests[i]);
    EXPECT_EQ(input.getRequiredEnergyGroups(), oracles[i])
      << "Case " << i << std::endl;
  }

#ifndef ENABLE_MPI_
  // Robust
  const std::vector<StringVec> fails = {
    {"--required-groups", "string"},
    {"--required-groups", "-1"},
  };

  for (size_t i = 0; i < fails.size(); ++i) {
    input.setArguments(fails[i]);
    EXPECT_ANY_THROW(input.getRequiredEnergyGroups())
      << "Case " << i << ": " << fails[i][0]
      << " " << fails[i][1] << std::endl;
  }
#endif
}


/// \brief Testing stabilization type
/// \details Type: None, Diagonal, Yamamoto, Global
TEST_F(test_ConfigInputCLI_parser, parseStabilizationType) {
  ROOT_ONLY();

  // Object to be tested
  ConfigInputCLI input;

  // Normal functionality
  const std::vector<StringVec> tests = {
    {"--stabilization", ""},
    {"--stabilization", "none"},
    {"--stabilization", "diagonal"},
    {"--stabilization", "YAMAMOTO"},
    {"--stabilization", "Global"},
    {"--stabilization=diagonal"},
  };

  const std::vector<stabilizationType> oracles = {
    +stabilizationType::NONE,
    +stabilizationType::NONE,
    +stabilizationType::DIAGONAL,
    +stabilizationType::YAMAMOTO,
    +stabilizationType::GLOBAL,
    +stabilizationType::DIAGONAL,
  };

  for (size_t i = 0; i < tests.size(); ++i) {
    input.setArguments(tests[i]);
    EXPECT_EQ(oracles[i], input.getStabilizationType());
  }

#ifndef ENABLE_MPI_
  // Robust
  const std::vector<StringVec> fails = {
    {"--stabilization", "Undefined"},
    {"--stabilization", "1.0"},
  };

  for (size_t i = 0; i < fails.size(); ++i) {
    input.setArguments(fails[i]);
    EXPECT_ANY_THROW(input.getStabilizationType())
      << "Case " << i << ": " << fails[i][0]
      << " " << fails[i][1] << std::endl;
  }
#endif
}


/// \brief Testing stabilization factor
/// \details Factor: positive floating-point number
TEST_F(test_ConfigInputCLI_parser, parseStabilizationFactor) {
  ROOT_ONLY();

  // Object to be tested
  ConfigInputCLI input;

  // Normal functionality
  const std::vector<StringVec> tests = {
    {"--stabilization-factor", "1"},
    {"--stabilization-factor", "1."},
    {"--stabilization-factor", "0.4"},
    {"--stabilization-factor", ".4"},
    {"--stabilization-factor", "1E-0"},
    {"--stabilization-factor", "1E-1"},
    {"--stabilization-factor=1e-5"},
  };

  const std::vector<double> oracles = {
    1,
    1.,
    0.4,
    .4,
    1E-0,
    1E-1,
    1E-5,
  };

  for (size_t i = 0; i < tests.size(); ++i) {
    input.setArguments(tests[i]);
    EXPECT_EQ(oracles[i], input.getStabilizationFactor());
  }

#ifndef ENABLE_MPI_
  // Robust
  const std::vector<StringVec> fails = {
    {"--stabilization-factor", "-1"},
    {"--stabilization-factor", " "},
    {"--stabilization-factor", ""},
  };

  for (size_t i = 0; i < fails.size(); ++i) {
    input.setArguments(fails[i]);
    EXPECT_ANY_THROW(input.getStabilizationFactor())
      << "Case " << i << ": " << fails[i][0]
      << " " << fails[i][1] << std::endl;
  }
#endif
}


/// \brief Testing boolean option --check-xs
TEST_F(test_ConfigInputCLI_parser, parseCheckXS) {
  ROOT_ONLY();

  // Object to be tested
  ConfigInputCLI input;

  // Normal functionality
  const std::vector<StringVec> tests = {
    {""},
    {"--check-xs"},
    {"--check-xs=true"},
    {"--check-xs=false"},
    // expects a true because it is a boolean, which can be assigned to by '='
    {"--check-xs", "false"},
  };

  const std::vector<bool> oracles = {
    false,
    true,
    true,
    false,
    true,
  };

  for (size_t i = 0; i < tests.size(); ++i) {
    input.setArguments(tests[i]);
    EXPECT_EQ(input.doesCheckXS(), oracles[i])
      << "Case " << i << std::endl;
  }
}


/// \brief Testing boolean option --correct-xs
TEST_F(test_ConfigInputCLI_parser, parseCorrectXS) {
  ROOT_ONLY();

  // Object to be tested
  ConfigInputCLI input;

  // Normal functionality
  const std::vector<StringVec> tests = {
    {""},
    {"--correct-xs"},
    {"--correct-xs=true"},
    {"--correct-xs=false"},
    // expects a true because it is a boolean, which can be assigned to by '='
    {"--correct-xs", "false"},
  };

  const std::vector<bool> oracles = {
    false,
    true,
    true,
    false,
    true,
  };

  for (size_t i = 0; i < tests.size(); ++i) {
    input.setArguments(tests[i]);
    EXPECT_EQ(input.doesCorrectXS(), oracles[i])
      << "Case " << i << std::endl;
  }
}


/// \brief Testing log level for mgxs checks
TEST_F(test_ConfigInputCLI_parser, parseCheckXSLogLevel) {
  ROOT_ONLY();

  // Object to be tested
  ConfigInputCLI input;

  // Normal functionality
  const std::vector<StringVec> tests = {
    {"--check-xs-log-level", ""},
    {"--check-xs-log-level", "debug"},
    {"--check-xs-log-level", "verbose"},
    {"--check-xs-log-level", "verbose_once"},
    {"--check-xs-log-level", "info"},
    {"--check-xs-log-level", "warn"},
    {"--check-xs-log-level", "warn_once"},
    {"--check-xs-log-level", "critical"},
    {"--check-xs-log-level", "error"},
  };

  const std::vector<log::level> oracles = {
    +log::level::error,
    +log::level::debug,
    +log::level::verbose,
    +log::level::verbose_once,
    +log::level::info,
    +log::level::warn,
    +log::level::warn_once,
    +log::level::critical,
    +log::level::error,
  };

  for (size_t i = 0; i < tests.size(); ++i) {
    input.setArguments(tests[i]);
    EXPECT_EQ(oracles[i], input.getCheckXSLogLevel());
  }

#ifndef ENABLE_MPI_
  // Robust
  const std::vector<StringVec> fails = {
    {"--check-xs-log-level", "Undefined"},
    {"--check-xs-log-level", "1.0"},
  };

  for (size_t i = 0; i < fails.size(); ++i) {
    input.setArguments(fails[i]);
    EXPECT_ANY_THROW(input.getCheckXSLogLevel())
      << "Case " << i << ": " << fails[i][0]
      << " " << fails[i][1] << std::endl;
  }
#endif
}


} /* namespace */
