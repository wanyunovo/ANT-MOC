/// \file src/ConfigInputFile.cpp

#include "antmoc/ConfigInputFile.h"
#include "antmoc/file_utils.h"
#include "antmoc/log.h"
#include "antmoc/ParentNodeStack.h"  /* ParentNodeStack */
#include "antmoc/string_utils.h"

#include <iomanip>
#include <ios>
#include <memory>

namespace antmoc
{

/// \details A ConfigInputFile object will insert its arguments into
///          the front of the argument list. This is done by first
///          read a config file from the argument list. If the file
///          is not specified, This object behaves as same as a CLI
///          object.
ConfigInputFile::ConfigInputFile(int &argc, char **argv)
  :ConfigInputCLI(argc, argv) {

  initializeOptions();
  expandArguments();
}


/// \details See ConfigInputFile(int&, char**)
ConfigInputFile::ConfigInputFile(const StringVec &argv)
  :ConfigInputCLI(argv) {

  initializeOptions();
  expandArguments();
}


/// \details See ConfigInputFile(int&, char**)
ConfigInputFile::ConfigInputFile(std::initializer_list<std::string> ilist)
  :ConfigInputCLI(ilist) {

  initializeOptions();
  expandArguments();
}


void ConfigInputFile::initializeOptions() {

  // Get the underlying option parser
  auto &options = getParser();

  options.add_options()
    ("help-settings", "Display options for XML and TOML and exit.")
  ;
}


/// \details Define the way we expand arguments
void ConfigInputFile::expandArguments() {

  // Read the path to configuration file
  const auto &path = getConfInputPath();

  // Process the configuration file only if it exists
  if ( stringutils::isSpaces(path) ) {
    log::verbose("Path to configuration file is empty, skipping argument expansion");
    return;
  } else {
    readArgumentsFromFile();
  }
}


ConfigInputFile::FileFormat ConfigInputFile::deduceFileFormat() {
  const auto file = getConfInputPath();
  const auto ext = fileutils::getExtension(file);

  if (stringutils::toUpper(ext) == stringutils::toUpper(".xml"))
    return FileFormat::XML;
  else
    return FileFormat::TOML;
}


/// \details Read settings from a file. If the filename is not specified,
///          do nothing. This method iterates all of the options in the
///          pre-defined option list. Each of them will be checked. If an option
///          exists in the file, it will be put into the underlying argument
///          vector. Unknown nodes will be skipped, including the wrongly
///          positioned ones.
void ConfigInputFile::readArgumentsFromFile() {

  std::shared_ptr<ParentNodeStack> parents;

  switch (deduceFileFormat()) {
    case ConfigInputFile::FileFormat::XML:
      parents = std::make_shared<ParentNodeStackXML>();
      break;
    case ConfigInputFile::FileFormat::TOML:
      parents = std::make_shared<ParentNodeStackTOML>();
      break;
  }

  const std::string log_head = "Parsing settings";

  // Load the input file
  const auto file = getConfInputPath();

  if (!parents->initialize(file)) {
    log::warn_once("{}: failed to parse {}, using default settings", log_head, file);
    return;
  }

  // A buffer for arguments
  StringVec argv;

  // Iterate through the option list. The current level is parents.size()-2.
  // An empty node name will cause stack poping but has no impact on pushing.
  // The number of parent level starts from 0.
  for (const auto &option_pair : _option_list) {

    const auto &name_in_file = option_pair.first;
    const auto &name_in_cli = option_pair.second;

    int cur_level = (int)parents->size() - 2;
    int level = parentLevel(name_in_cli);

    // If the current node has a parent tag, it is a parent node.
    if (level >= 0) {

      log::debug("{}: saw parent node at level {}", log_head, level);

      if (level == cur_level) {
        log::debug("{}: swapping parent '{}' with current top at level {}", log_head, name_in_file, level);

        // Swap parent node at the same level
        parents->swapTop(name_in_file);
      }
      else if (level < cur_level) {
        log::debug("{}: stepping back into parent '{}' from level {} to {}", log_head, name_in_file, cur_level, level);

        // Remove elements to make the two level equal
        for (int i = 0; i < cur_level - level; ++i)
          parents->pop();
        parents->swapTop(name_in_file);
      }
      else if (level == cur_level + 1) {
        parents->pushNextParent(name_in_file);
      }
      else {
        log::error("Internal failure: wrongly formatted option list for file at "
                   "option '{}' '{}', level {}", name_in_file, name_in_cli, level);
      }
    }

    // If the current node is an actual option
    else {

      if (name_in_file.empty())
        log::error("Internal failure: the option list for file has a non-parent "
                   "empty option", name_in_file, name_in_cli);
      // Debugging
      log::debug("{}: searching for option '{}', {} parents in stack", log_head, name_in_file, parents->size());

      // Save the option if it exists
      parents->saveOption(argv, _boolean_options, name_in_file, name_in_cli);
    }

  } // end for

  // Insert the retrieved options into the underlying vector
  ConfigInputCLI::insertFrontArguments(argv);
}


std::string ConfigInputFile::parentTag(const int level) {
  if (level < 0)
    log::error("Invalid option tag level: {}", level);

  return std::string(level, ' ');
}


int ConfigInputFile::parentLevel(const std::string &tag) {
  if (stringutils::isSpaces(tag))
    return tag.size();
  else
    return -1;
}


/// \details If 'help-settings' is found in the argument list, print
///          the help and abort the program, otherwise do nothing.
void ConfigInputFile::showExtraHelpAsNeeded() {
  // Print help information if needed
  if (hasOption("help-settings")) {

    if (mpi::isMPIRoot()) {
      // Compute the maximum length of options
      size_t len_left = 0;
      size_t len_right = 0;
      int max_level = 0;
      for (const auto &entry : _option_list) {
        len_left  = std::max(len_left, entry.first.size());
        len_right = std::max(len_right, entry.second.size());
        max_level = std::max(max_level, parentLevel(entry.second));
      }

      len_left += max_level * 2;

      // Extra spaces
      std::string spaces(4, ' ');

      // Title
      fmt::print("\nNote: parent nodes must be sub-elements\n\n");
      fmt::print("{0: ^{1}}{2}{3: ^{4}}\n", "File options", len_left, spaces,
                                            "CLI options", len_right);
      fmt::print("{0:-^{1}}{2}{3:-^{4}}\n", "", len_left, spaces,
                                            "", len_right);

      // Leading spaces
      int cur_level = 0;

      // Entries
      for (const auto &entry : _option_list) {
        const auto &name_in_file = entry.first;
        const auto &name_in_cli = parentLevel(entry.second) < 0
                                  ? "--" + entry.second
                                  : "";

        // Entry indentation
        auto indent = (cur_level+1) * 2;

        // Save the level if a parent tag is encountered
        auto level = parentLevel(entry.second);
        if (level >= 0) {
          cur_level = level;
          indent = level * 2;

          if (level == 0)
            fmt::print("\n");
        }

        fmt::print("{0: <{1}}{2: <{3}}{4}{5: <{6}}\n",
                   "", indent,
                   name_in_file, len_left - indent, spaces,
                   name_in_cli, len_right);
      }
    }

#   ifdef ENABLE_MPI_
    MPI_Finalize();
#   endif
    std::exit(0);
  }
}

} // namespace antmoc
