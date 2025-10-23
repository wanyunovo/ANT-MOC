/// \file src/ConfigInputCLI.cpp

#include "antmoc/ConfigInputCLI.h"
#include "antmoc/log.h"
#include "antmoc/mpi_utils.h"

#include <iostream>
#include <omp.h>

namespace antmoc {

/// \brief Default constructor
ConfigInputCLI::ConfigInputCLI() {
  initializeOptions();
}


ConfigInputCLI::ConfigInputCLI(int &argc, char **argv):
  ConfigInput(argc, argv) {
  initializeOptions();
}


ConfigInputCLI::ConfigInputCLI(const StringVec &argv):
  ConfigInput(argv) {
  initializeOptions();
}

ConfigInputCLI::ConfigInputCLI(std::initializer_list<std::string> ilist):
  ConfigInput(ilist) {
  initializeOptions();
}


void ConfigInputCLI::initializeOptions() {

  // Get the underlying option parser
  auto &options = getParser();

  options.add_options()
    ("h,help",          "Display this help and exit.")
    (_l_loglevel,       "Set the minimum log level for output.\n"
                        "(default: info)\n"
                        "    error, warn, info, verbose, profile, debug"
                        , cxxopts::value<std::string>())
    (_l_logline,        "Set the length of each log line."
                        , cxxopts::value<size_t>()->default_value(std::to_string(log::get_default_line_length())))
    (_l_logranks,       "Set a subset of ranks for logging.\n(default: '')\n"
                        "This is a comma-seperated list of ranks\n"
                        "Integer range 'start:stop:step' is supported.\n"
                        "    e.g. ranks 0,2,4,6 can be written as\n"
                        "      --logging-ranks '0,2,4,6'\n"
                        "      --logging-ranks '0:6:2'"
                        , cxxopts::value<std::string>(), "LIST")
    (_a_logfile,        "Set the log path. (default: '')\n"
                        "For now, only the root can create a log file.\n"
                        , cxxopts::value<std::string>(), "FILE")
  ;

  options.add_options("Parallel computing")
    (_l_threads,        "Specify the number of omp threads.\n"
                        , cxxopts::value<int>()->default_value(std::to_string(omp_get_max_threads())))
    (_l_domains,        "Specify the number of mpi domains.\n"
                        "(default: '1,1,1,0')\n"
                        "Domains are described in 4-D\n"
                        "    nx, ny, nz, ns\n"
                        "The last domain can be detected if it is 0 or omitted. "
                        "Otherwise, the total number of domains must equal the "
                        "number of issued ranks, that is,\n"
                        "    # ranks == nx * ny * nz * ns"
                        , cxxopts::value<std::string>())
    (_l_track_mp,       "Specify the mapping algorithm for track load balance\n"
                        "(default: Auto)\n"
                        "    Block - block distribution\n"
                        "    Cyclic Track - load balancing for cyclic track decomposition\n"
                        "    Angle - load balancing for angular decomposition\n"
                        "    Auto - choose the algorithm with the highest expected efficiency"
                        , cxxopts::value<std::string>())
  ;

  options.add_options("Data input")
    (_a_conf,           "Set the path to configuration input file."
                        , cxxopts::value<std::string>(), "FILE")
    (_a_geo,            "Set the path to geometry input file."
                        , cxxopts::value<std::string>(), "FILE")
    (_l_prim,           "Set the path to global primitives file."
                        , cxxopts::value<std::string>(), "FILE")
    (_a_mat,            "Set the path to materials input file."
                        , cxxopts::value<std::string>(), "FILE")
    (_l_xsl,            "Specify the file layout of cross-section\n"
                        "data. (default: Named)\n"
                        "    Compressed, Named"
                        , cxxopts::value<std::string>())
    (_l_lat_refines,    "Specify a value to override local refines of lattices.\n"
                        "(default: 0,0,0)\n"
                        "This option overrides all of the occurrences of element 'refines' "
                        "of lattices in the primitives and geometry files.\n"
                        "It accepts 3 ints: 'x,y,z'\n"
                        "    e.g. '0,0,2' overrides only local z refines"
                        , cxxopts::value<std::string>(), "LIST")
    (_l_cell_sectors,   "Specify a value to set the number of sectors of cells.\n"
                        "(default: '')\n"
                        "This option overrides the element 'sectors' of material cells.\n"
                        "It accepts 1 int: 'n'"
                        , cxxopts::value<std::string>())
    (_l_cell_rings,     "Specify a value to set the number of rings of cells.\n"
                        "(default: '')\n"
                        "This option overrides the element 'rings' of material cells.\n"
                        "It accepts 1 int: 'n'"
                        , cxxopts::value<std::string>())
    (_l_prim_cl,        "Clean up unused primitives after building the geometry."
                        , cxxopts::value<bool>()->default_value("true"))
    (_l_mat_cl,         "Clean up unused materials after building the geometry."
                        , cxxopts::value<bool>()->default_value("true"))
  ;

  options.add_options("Data output")
    (_a_outputdir,      "Set the output directory.\n"
                        , cxxopts::value<std::string>()->default_value(getTimeStamp()), "DIR")
    (_l_dump_settings,  "Dump settings to files or not."
                        , cxxopts::value<bool>()->default_value("false"))
    (_l_dump_fsrs,      "Dump FSR data to binary files or not."
                        , cxxopts::value<bool>()->default_value("false"))
    (_l_dumptracks,     "Like --dump-rx but for tracks.\n(default: None)\n"
                        "    2D   - 2D tracks\n"
                        "    3D   - 3D tracks\n"
                        "    All  - All tracks\n"
                        "    None - No dumping"
                        , cxxopts::value<std::string>(), "LIST")
    (_l_meshtype,       "Specify the type of tally mesh.\n(default: rectangle)\n"
                        "This option affects the way that '--mesh' is interpreted\n"
                        "    rectangle - a rectangular tally mesh\n"
                        "    hexagon - a hexagonal tally mesh"
                        , cxxopts::value<std::string>())
    (_l_mesh,           "Pass a tally mesh for dumping reaction rates.\n"
                        "(default: '', dump nothing)\n"
                        "If it is a rectangular tally mesh, the argument "
                        "has three components 'widths_x, widths_y, widths_z'\n"
                        "    e.g. --mesh '[1.26]*3, [1.26]*3, 7.14'\n"
                        "If it is a hexagonal tally mesh, the argument has "
                        "three components 'num_r, width_r, widths_z'\n"
                        "    e.g. --mesh '2, 1.5, [7.14]*3'"
                        , cxxopts::value<std::string>(), "LIST")
    (_l_meshoffset,     "Specify the offset for tally mesh.\n(default: 0,0,0)\n"
                        "The offset is relative to geometry center.\n"
                        "    e.g. --mesh-offset '0.5, -0.5, 1'"
                        , cxxopts::value<std::string>(), "LIST")
    (_l_orientation,    "Specify the orientation for hexagonal tally mesh.\n"
                        "(default: y)\n"
                        "    y - a y-orientated mesh\n"
                        "    x - an x-orientated mesh"
                        , cxxopts::value<std::string>())
    (_l_dumprx,         "Specify the types of reaction rates to be dumped.\n"
                        "(default: PHI)\n"
                        "This is a comma-seperated list.\n"
                        "    e.g. --dump-rx 'f, a, phi'\n"
                        "Available types:\n"
                        "    F    - Fission\n"
                        "    NuF  - NuFission\n"
                        "    A    - Absorption\n"
                        "    T    - Total\n"
                        "    Phi  - Scalar Flux\n"
                        "    All  - All rates\n"
                        "    None - No dumping"
                        , cxxopts::value<std::string>(), "LIST")
    (_l_dumpxs,         "Like --dump-rx but for cross-sections.\n(default: None)\n"
                        "Available types:\n"
                        "    F, NuF, A, T, All, None"
                        , cxxopts::value<std::string>(), "LIST")
    (_l_tallygroups,    "Select energy groups to dump.\n(default: ALL)\n"
                        "This is a comma-seperated list of group ids or string "
                        "'ALL' which indicates all of the groups.\n"
                        "Integer range 'start:stop:step' is supported.\n"
                        "    e.g. energy groups 1,3,5,7 can be written as\n"
                        "      --tally-groups '1,3,5,7'\n"
                        "      --tally-groups '1:7:2'\n"
                        "A special group is 0, which represents the sum of each "
                        "group. This group is always dumped, but it can be the "
                        "only one to be dumped by passing it explicily.\n"
                        "    e.g. only the sum of mesh data\n"
                        "      --tally-groups '0'"
                        , cxxopts::value<std::string>(), "LIST")
  ;

  options.add_options("Quadrature")
    (_a_quad,           "Specify the quadrature type.\n(default: Equal Angle)\n"
                        "    Equal Angle\n"
                        "    Equal Weight\n"
                        "    Gauss Legendre\n"
                        "    Leonard\n"
                        "    Tabuchi Yamamoto"
                        , cxxopts::value<std::string>())
    (_a_nazims,         "Specify the number of azimuthal angles."
                        , cxxopts::value<int>()->default_value("32"))
    (_a_npolars,        "Specify the number of polar angles."
                        , cxxopts::value<int>()->default_value("8"))
    (_a_sazim,          "Specify the spacing between tracks on x-y plane."
                        , cxxopts::value<double>()->default_value("0.05"))
    (_a_spolar,         "Specify the spacing between tracks on l-z plane."
                        , cxxopts::value<double>()->default_value("0.75"))
  ;

  options.add_options("Ray tracing")
    (_a_segform,        "Specify the segmentation formation.\n(default: OTF Stacks)\n"
                        "    OTF Stacks\n"
                        "    OTF Tracks\n"
                        "    Explicit 2D"
                        , cxxopts::value<std::string>())
    (_a_segzones,       "Specify the axial segmentation zones.\n"
                        "(default: 'Auto',  means auto-generated)\n"
                        "This is a comma-separated list of axial planes. Any two "
                        "of them are supposed to define an axial zone.\n"
                        "If this argument is set to empty string, antmoc will try "
                        "to find possible zones by iterating over geometry objects.\n"
                        "    e.g. --zones '-32.13, 0, 32.13'\n"
                        "Planes below the geometry bottom or above the geometry top "
                        "will be ignored."
                        , cxxopts::value<std::string>(), "LIST")
    (_l_zmesh,          "Specify the z mesh type. (default: Local)\n"
                        "    Local, Global\n"
                        "Local z mesh is generated by vertically ray tracing."
                        , cxxopts::value<std::string>())
    (_l_modules,        "Specify the number of tracing modules.\n"
                        "(default: 1,1,1)\n"
                        "This option changes the layout of tracks."
                        , cxxopts::value<std::string>(), "LIST")
  ;

  options.add_options("Solver")

    (_l_solver,         "Specify the solver type. (default: CPU Solver)\n"
                        "    CPU Solver\n"
                        "    CPU LS Solver"
                        , cxxopts::value<std::string>())
    (_l_keff,           "Set computation method of k-eff from fission, "
                        "absorption, and leakage rates rather than from "
                        "successive fission rates.\n"
                        "    keff = fission/(absorption + leakage)"
                        , cxxopts::value<bool>()->default_value("false"))
    (_a_niters,         "Set the maximum iterations for solver."
                        , cxxopts::value<int>()->default_value("1000"))
    (_a_tol,            "Set the tolerance of convergence for solver."
                        , cxxopts::value<double>()->default_value("1e-5"))
    (_l_stab,           "Specify the stabilization option. (default: None)\n"
                        "    Diagonal, Yamamoto, Global, None"
                        , cxxopts::value<std::string>())
    (_l_stab_factor,    "Set the factor for stabilization."
                        , cxxopts::value<double>()->default_value("0.5"))
    (_l_check_xs,       "Let the solver check MGXS."
                        , cxxopts::value<bool>()->default_value("false"))
    (_l_correct_xs,     "Let the solver check and correct MGXS."
                        , cxxopts::value<bool>()->default_value("false"))
    (_l_check_xs_log_level, "Set the log level for MGXS checks."
                            , cxxopts::value<std::string>()->default_value("error"))
  ;

  options.add_options("Experimental")

    (_l_requiredgroups, "Set the number of required energy groups. "
                        "This option is used to debug the code.\n"
                        "(default: 'ALL')\n"
                        "Value 0 and string 'ALL' indicates all of the groups."
                        , cxxopts::value<std::string>())
  ;
}


void ConfigInputCLI::showHelp() {
  if (mpi::isMPIRoot()) {
    std::cout << getParser().help({"",
                                   "Parallel computing",
                                   "Data input",
                                   "Data output",
                                   "Quadrature",
                                   "Ray tracing",
                                   "Solver",
                                   "Experimental"})
              << std::endl;
  }
# ifdef ENABLE_MPI_
  MPI_Finalize();
# endif
  std::exit(0);
}


void ConfigInputCLI::showHelpAsNeeded() {

  // Print extra help which is provided by derived classes
  showExtraHelpAsNeeded();

  // Print help information if needed
  if (hasOption("help")) {
    showHelp();
  }
}

} // namespace antmoc
