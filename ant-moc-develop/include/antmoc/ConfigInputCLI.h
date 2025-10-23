/// \file include/ConfigInputCLI.h
/// \details Handling command line input.
/// \date June 7, 2019
/// \author An Wang, USTB (wangan@xs.ustb.edu.cn)
///         Inspired by iain from StackOverflow.

#ifndef CONFIGINPUTCLI_H_
#define CONFIGINPUTCLI_H_

#include <string>
#include <vector>

#include "antmoc/ConfigInput.h"

namespace antmoc {

///---------------------------------------------------------------------
/// \class ConfigInputCLI
/// \brief A reader for command line input
///---------------------------------------------------------------------
class ConfigInputCLI : public ConfigInput {

public:

  ConfigInputCLI();

  /// \brief Construct an object with standard CLI arguments.
  /// \param argc The number of arguments (at least 1).
  /// \param argv The arguments, including the executable name.
  ConfigInputCLI(int &argc, char **argv);

  /// \brief Construct an object with a vector of arguments.
  /// \param argv A vector of arguments, not including the executable name
  ConfigInputCLI(const StringVec &argv);

  /// \brief Construct an object with a list of arguments.
  /// \param argv A list of arguments, not including the executable name
  ConfigInputCLI(std::initializer_list<std::string> ilist);

  virtual ~ConfigInputCLI() = default;

  /// \brief Print the help message and exit.
  void showHelp();

  /// \brief Print the help message when needed and exit.
  void showHelpAsNeeded();

  // Virtual methods for derived classes
  void showExtraHelpAsNeeded() { return; }

private:

  /// \brief Initialize the cxxopts::Options object
  void initializeOptions();

};


} // namespace antmoc

#endif  // CONFIGINPUTCLI_H_
