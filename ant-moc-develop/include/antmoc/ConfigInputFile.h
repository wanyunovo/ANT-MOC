/// \file include/ConfigInputFile.cpp
/// \brief Parsing settings file in format XML or TOML.
/// \date June 14, 2019
/// \author An Wang, USTB (wangan@xs.ustb.edu.cn)
///         Gen Wang, USTB (17812132070@163.com)

#ifndef CONFIGINPUTFILE_H_
#define CONFIGINPUTFILE_H_

#include <string>
#include <utility>  /* std::pair */
#include <vector>

#include "antmoc/ConfigInputCLI.h"

namespace antmoc
{

///---------------------------------------------------------------------
/// \class ConfigInputFile
/// \brief A reader for settings in XML or TOML format
///---------------------------------------------------------------------
class ConfigInputFile : public ConfigInputCLI
{
public:

  /// \enum Format of the settings file.
  enum FileFormat {
    XML,
    TOML
  };

  /// \brief Default constructor.
  ConfigInputFile() = default;

  /// \brief Construct an object with standard CLI arguments.
  /// \param argc The number of arguments (at least 1).
  /// \param argv The arguments, including the executable name.
  ConfigInputFile(int &argc, char **argv);

  /// \brief Construct an object with a vector of arguments.
  /// \param argv A vector of arguments, not including the executable name
  ConfigInputFile(const StringVec &argv);

  /// \brief Construct an object with a list of arguments.
  /// \param argv A list of arguments, not including the executable name
  ConfigInputFile(std::initializer_list<std::string> ilist);

  virtual ~ConfigInputFile() = default;

  /// \brief Expand argument -c.
  void expandArguments();

  /// \brief Print help messages for settings and exit.
  void showExtraHelpAsNeeded();

private:

  /// \brief Initialize the cxxopts::Options object
  /// \details This method is supposed to be invoked after the base class
  ///          initializes its options.
  void initializeOptions();

  /// \brief Deduce file format.
  /// \details The format defaults to TOML.
  ///          TOML file: .toml
  ///          XML file: .xml
  FileFormat deduceFileFormat();

  /// \brief Read arguments from a settings file
  void readArgumentsFromFile();

  /// \brief Returns a tag of parent node
  /// \param level a integer indicates the level
  static std::string parentTag(const int level = 0);

  /// \brief Returns the level of the parent tag
  static int parentLevel(const std::string &tag);

  ///< An ordered list of options
  std::vector<std::pair<std::string, std::string>> _option_list {
    // There are two types of nodes:
    // 1. Parent nodes, which can only be defined as sub-elements.
    // 2. Option nodes, which can be defined as sub-elements or attributes.

    //-------------------------------------
    // Log options
    //-------------------------------------
    {"", parentTag()},
    {"log_level", _l_loglevel},
    {"log_line",  _l_logline},
    {"logging_ranks", _l_logranks},
    {"log_file",  _l_logfile},

    //-------------------------------------
    // Parallel computing options
    //-------------------------------------
    {"parallel", parentTag()},
    {"omp_threads",   _l_threads},
    {"domains",       _l_domains},
    {"track_mapping", _l_track_mp},

    //-------------------------------------
    // Geometry options
    // Parent: geometry
    //-------------------------------------
    {"geometry", parentTag()},
    {"path",            _l_geo},
    {"primitives",      _l_prim},
    {"lattice_refines", _l_lat_refines},
    {"cell_sectors",    _l_cell_sectors},
    {"cell_rings",      _l_cell_rings},
    {"clean_unused",    _l_prim_cl},

    //-------------------------------------
    // Materials options
    // Parent: geometry
    //-------------------------------------
    {"materials", parentTag()},
    {"path",      _l_mat},
    {"xs_layout", _l_xsl},
    {"clean_unused",  _l_mat_cl},

    //-------------------------------------
    // Output options
    // Parent: output
    //           mesh
    //-------------------------------------
    {"output", parentTag()},
    {"directory",       _l_outputdir},
    {"dump_settings",   _l_dump_settings},
    {"dump_fsrs",       _l_dump_fsrs},
    {"tracks",          _l_dumptracks},
    {"reaction_rates",  _l_dumprx},
    {"cross_sections",  _l_dumpxs},
    {"groups",          _l_tallygroups},
      // Sub options of tally
      {"mesh", parentTag(1)},
      {"type",            _l_meshtype},
      {"shape",           _l_mesh},
      {"offset",          _l_meshoffset},
      {"orientation",     _l_orientation},

    //-------------------------------------
    // Quadrature options
    // Parent: quadrature
    //-------------------------------------
    {"quadrature", parentTag()},
    {"type",          _l_quad},
    {"azims",         _l_nazims},
    {"polars",        _l_npolars},
    {"azim_spacing",  _l_sazim},
    {"z_spacing",     _l_spolar},

    //-------------------------------------
    // Ray tracing options
    // Parent: ray_tracing
    //-------------------------------------
    {"ray_tracing", parentTag()},
    {"segmentation",  _l_segform},
    {"zones",         _l_segzones},
    {"modules",       _l_modules},
    {"z_mesh",        _l_zmesh},

    //-------------------------------------
    // Solver options
    // Parent: solver
    //-------------------------------------
    {"solver", parentTag()},
    {"type",                  _l_solver},
    {"max_iterations",        _l_niters},
    {"tolerance",             _l_tol},
    {"keff_neutron_balance",  _l_keff},
    {"stabilization",         _l_stab},
    {"stabilization_factor",  _l_stab_factor},
      {"xs", parentTag(1)},
      {"check",     _l_check_xs},
      {"correct",   _l_correct_xs},
      {"log_level", _l_check_xs_log_level},

    //-------------------------------------
    // Experimental options
    // Parent: experimental
    //-------------------------------------
    {"experimental", parentTag()},
    {"required_groups", _l_requiredgroups},
  };

  ///< Keep track of boolean options because they are parsed in a different way
  std::set<std::string> _boolean_options = {
    "clean_unused",
    "keff_neutron_balance",
    "check",
    "correct",
    "dump_settings",
    "dump_fsrs",
  };
};


} // namespace antmoc

#endif  // CONFIGINPUTFILE_H_
