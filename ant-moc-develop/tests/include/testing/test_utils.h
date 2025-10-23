#ifndef TEST_UTILS_H_
#define TEST_UTILS_H_

//#include "app_headers.h"
#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include "cxxopts.hpp"
#include "antmoc/mpi_utils.h"

#define ROOT_ONLY()                 \
{                                   \
  if(!antmoc::mpi::isMPIRoot())  \
    return;                         \
}


///---------------------------------------------------------------------
/// \class TestOptions
/// \brief Parses options in addition to googletest options.
///---------------------------------------------------------------------
class TestOptions {

  cxxopts::Options _options {"antmoctest", ""};
  bool _update_results;
  bool _visualize;
  std::string _log_level;

  TestOptions() = default;
  TestOptions(const TestOptions&) = delete;
  TestOptions(TestOptions&&) = delete;
  TestOptions &operator=(const TestOptions&) = delete;

public:
  
  virtual ~TestOptions() = default;

  // Returns a static instance
  static TestOptions &get() {
    static TestOptions instance;
    return instance;
  }

  /// \brief Initializes the options object and parses options.
  void parse(int argc, char **argv) {
    _options
      .allow_unrecognised_options()
      .add_options()
      ("update-results", "Update testing results"
                         , cxxopts::value<bool>(_update_results)->default_value("false"))
      ("visualize",      "Produce visualization date"
                         , cxxopts::value<bool>(_visualize)->default_value("false"))
      ("log-level",      "Set the minimum log level for output."
                         , cxxopts::value<std::string>(_log_level)->default_value("info"))
    ;

    // Creates a buffer
    char **argv_cpy = new char*[argc];
    for (int i = 0; i < argc; ++i)
      argv_cpy[i] = argv[i];

    // Parses arguments
    _options.parse(argc, argv_cpy);

    // Releases the buffer
    delete[] argv_cpy;
  }


  //--------------------------------------------------------------------
  // Get arguments
  //--------------------------------------------------------------------
  bool updateResults() { return _update_results; }
  bool visualize()     { return _visualize; }
  std::string getLogLevel() { return _log_level; }

};

#endif  // TEST_UTILS_H_
