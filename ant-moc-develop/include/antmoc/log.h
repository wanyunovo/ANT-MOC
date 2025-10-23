/// \file log.h
/// \brief A custom logger with fmt as the backend. Based on OpenMOC.
/// \date April 18, 2020
/// \author An Wang, USTB (wangan@xs.ustb.edu.cn)

#ifndef LOG_H_
#define LOG_H_

#include <set>
#include <string>

#include <fmt/core.h>
#include <fmt/format.h>
#include <fmt/chrono.h>
#include <fmt/ostream.h>
#include <fmt/printf.h>

#include "antmoc/enum_types.h"


namespace antmoc {

/// \namespace antmoc::log
namespace log {


/// \enum level
/// \brief Logging levels characterize an ordered set of message types
/// \      which may be printed to the screen.
BETTER_ENUM(level, char,
  /**< A debugging message */
  debug,

  /**< A profiling message */
  profile,

  /**< A profiling message printed by rank 0 process only */
  profile_once,

  /**< An informational but verbose message */
  verbose,

  /**< An informational verbose message printed by rank 0 process only */
  verbose_once,

  /**< A brief progress update on run progress */
  info,

  /**< A brief progress update by node on run progress */
  node,

  /**< A message of a single line of characters */
  separator,

  /**< A message centered within a line of characters */
  header,

  /**< A message sandwiched between two lines of characters */
  title,

  /**< A message to warn the user */
  warn,

  /**< A message to warn the user - to be printed by rank 0 process only */
  warn_once,

  /**< A message to warn of critical program conditions */
  critical,

  /**< A message containing program results */
  result,

  /**< A messsage for testing */
  test,

  /**< A message reporting error conditions */
  error
)


//----------------------------------------------------------------------
// Logging interfaces (fprint)
//----------------------------------------------------------------------
/// \brief Debug
template <typename... Args>
void fdebug(const char* format, const Args & ... args) {
  print_formatted_string(level::debug, fmt::sprintf(format, args...));
}

/// \brief Profile
template <typename... Args>
void fprofile(const char* format, const Args & ... args) {
  print_formatted_string(level::profile, fmt::sprintf(format, args...));
}

/// \brief Profile_once
template <typename... Args>
void fprofile_once(const char* format, const Args & ... args) {
  print_formatted_string(level::profile_once, fmt::sprintf(format, args...));
}

/// \brief Verbose
template <typename... Args>
void fverbose(const char* format, const Args & ... args) {
  print_formatted_string(level::verbose, fmt::sprintf(format, args...));
}

/// \brief Verbose_once
template <typename... Args>
void fverbose_once(const char* format, const Args & ... args) {
  print_formatted_string(level::verbose_once, fmt::sprintf(format, args...));
}

/// \brief Info
template <typename... Args>
void finfo(const char* format, const Args & ... args) {
  print_formatted_string(level::info, fmt::sprintf(format, args...));
}

/// \brief Node
template <typename... Args>
void fnode(const char* format, const Args & ... args) {
  print_formatted_string(level::node, fmt::sprintf(format, args...));
}

/// \brief Warn
template <typename... Args>
void fwarn(const char* format, const Args & ... args) {
  print_formatted_string(level::warn, fmt::sprintf(format, args...));
}

/// \brief Warn_once
template <typename... Args>
void fwarn_once(const char* format, const Args & ... args) {
  print_formatted_string(level::warn_once, fmt::sprintf(format, args...));
}

/// \brief Critical
template <typename... Args>
void fcritical(const char* format, const Args & ... args) {
  print_formatted_string(level::critical, fmt::sprintf(format, args...));
}

/// \brief Result
template <typename... Args>
void fresult(const char* format, const Args & ... args) {
  print_formatted_string(level::result, fmt::sprintf(format, args...));
}

/// \brief Test
template <typename... Args>
void ftest(const char* format, const Args & ... args) {
  print_formatted_string(level::test, fmt::sprintf(format, args...));
}

/// \brief Error
template <typename... Args>
void ferror(const char* format, const Args & ... args) {
  print_formatted_string(level::error, fmt::sprintf(format, args...));
}

/// \brief Separator
template <typename... Args>
void fseparator(const char* format, const Args & ... args) {
  print_formatted_string(level::separator, fmt::sprintf(format, args...));
}

/// \brief Header
template <typename... Args>
void fheader(const char* format, const Args & ... args) {
  print_formatted_string(level::header, fmt::sprintf(format, args...));
}

/// \brief Title
template <typename... Args>
void ftitle(const char* format, const Args & ... args) {
  print_formatted_string(level::title, fmt::sprintf(format, args...));
}

//----------------------------------------------------------------------
// Logging interfaces (print)
//----------------------------------------------------------------------
/// \brief Debug
template <typename... Args>
void debug(const char* format, const Args & ... args) {
  print_formatted_string(level::debug, fmt::format(format, args...));
}

/// \brief Profile
template <typename... Args>
void profile(const char* format, const Args & ... args) {
  print_formatted_string(level::profile, fmt::format(format, args...));
}

/// \brief Profile_once
template <typename... Args>
void profile_once(const char* format, const Args & ... args) {
  print_formatted_string(level::profile_once, fmt::format(format, args...));
}

/// \brief Verbose
template <typename... Args>
void verbose(const char* format, const Args & ... args) {
  print_formatted_string(level::verbose, fmt::format(format, args...));
}

/// \brief Verbose_once
template <typename... Args>
void verbose_once(const char* format, const Args & ... args) {
  print_formatted_string(level::verbose_once, fmt::format(format, args...));
}

/// \brief Info
template <typename... Args>
void info(const char* format, const Args & ... args) {
  print_formatted_string(level::info, fmt::format(format, args...));
}

/// \brief Node
template <typename... Args>
void node(const char* format, const Args & ... args) {
  print_formatted_string(level::node, fmt::format(format, args...));
}

/// \brief Warn
template <typename... Args>
void warn(const char* format, const Args & ... args) {
  print_formatted_string(level::warn, fmt::format(format, args...));
}

/// \brief Warn_once
template <typename... Args>
void warn_once(const char* format, const Args & ... args) {
  print_formatted_string(level::warn_once, fmt::format(format, args...));
}

/// \brief Critical
template <typename... Args>
void critical(const char* format, const Args & ... args) {
  print_formatted_string(level::critical, fmt::format(format, args...));
}

/// \brief Result
template <typename... Args>
void result(const char* format, const Args & ... args) {
  print_formatted_string(level::result, fmt::format(format, args...));
}

/// \brief Test
template <typename... Args>
void test(const char* format, const Args & ... args) {
  print_formatted_string(level::test, fmt::format(format, args...));
}

/// \brief Error
template <typename... Args>
void error(const char* format, const Args & ... args) {
  print_formatted_string(level::error, fmt::format(format, args...));
}

/// \brief Separator
template <typename... Args>
void separator(const char* format, const Args & ... args) {
  print_formatted_string(level::separator, fmt::format(format, args...));
}

/// \brief Header
template <typename... Args>
void header(const char* format, const Args & ... args) {
  print_formatted_string(level::header, fmt::format(format, args...));
}

/// \brief Title
template <typename... Args>
void title(const char* format, const Args & ... args) {
  print_formatted_string(level::title, fmt::format(format, args...));
}


//----------------------------------------------------------------------
// Legacy interfaces
//----------------------------------------------------------------------
/// \brief Interface for printing a printf-style message (legacy)
template <typename... Args>
void log_printf(level log_level, const char* format, const Args & ... args) {
  print_formatted_string(log_level, fmt::sprintf(format, args...));
}

/// \brief Interface for printing a fmt-style message (legacy)
template <typename... Args>
void log_print(level log_level, const char* format, const Args & ... args) {
  print_formatted_string(log_level, fmt::format(format, args...));
}

/// \brief The core function for creating formatted log messages.
/// \return A log message with prefix.
std::string create_log_string(level log_level, const std::string &msg);

/// \brief The core function for printing log messages.
void print_formatted_string(level log_level, const std::string &msg);

/// \brief The core function for printing error messages.
void print_formatted_error(const std::string &msg);

/// \brief The core function for printing messages to file.
void print_string_to_file(const std::string &msg);

/// \brief Breaks up a long message into multiline messages
std::string create_multiline_msg(std::string level, std::string message);


//----------------------------------------------------------------------
// Logger configuration
//----------------------------------------------------------------------
/// \brief Initializes the logger for use.
void initialize();

/// \brief Sets the rank allowed to log messages.
void set_logging_ranks(std::set<int> ranks);
std::set<int> get_logging_ranks();

/// \brief Sets and creates the log file.
/// \param path Path to the log file.
void set_path(const std::string &path);

/// \brief Returns the log path.
std::string get_path();

void set_level(const std::string &new_level);
void set_level(level new_level);
level get_level();

/// \brief Sets the maximum line length for log messages.
void set_line_length(size_t length);

/// \brief Gets the maximum line length for log messages.
size_t get_line_length();

/// \brief Gets the default line length for log messages.
size_t get_default_line_length();

void set_separator_character(char c);
char get_separator_character();
void set_header_character(char c);
char get_header_character();
void set_title_character(char c);
char get_title_character();


} // namespace log

} // namespace antmoc

#endif /* LOG_H_ */
