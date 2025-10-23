#include "antmoc/log.h"
#include "antmoc/file_utils.h"
#include "antmoc/mpi_utils.h"
#include "antmoc/string_utils.h"

#include <omp.h>

#include <cmath>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace antmoc {


namespace log {


//----------------------------------------------------------------------
// Static variables
//----------------------------------------------------------------------
/// \var _log_level
/// \brief Minimum level of logging messages printed to the screen and log file.
/// \details The default logging level is info.
static level _log_level = level::info;


/// \var _logging_ranks
/// \brief Ranks allowed to show messages (will overwirte the root process).
/// \details An empty set indicates the default behaviour.
static std::set<int> _logging_ranks;

/// \var _log_path
/// \brief The name of the output log file.
static std::string _log_path = "";


/// \var _logging
/// \brief A switch which is set to true once the first message is logged.
/// \details The logging switch is needed to indicate whether the output
///          log file has been created or not.
static bool _logging = false;


/// \var _separator_char
/// \brief The character to use for separator log messages. The default is "-".
static char _separator_char = '*';


/// \var _header_char
/// \brief The character to use for header log messages. The default is "*".
static char _header_char = '*';


/// \var _title_char
/// \brief The character to use for title log messages. The default is "*".
static char _title_char = '*';


/// \var _default_line_length
/// \brief The default line length for a log message.
constexpr size_t _default_line_length = 74;


/// \var _line_length
/// \brief The maximum line length for a log message.
/// \details The default is 74.
static size_t _line_length = 74;


/// \var _error_lock
/// \brief OpenMP mutex lock for error messages which throw exceptions
static omp_lock_t _error_lock;


//----------------------------------------------------------------------
// Functions
//----------------------------------------------------------------------
/// \details This should be immediately called when the logger is imported
///          into Python and before any of its other routines are called. The
///          routine initializes an OpenMP mutual exclusion lock which is used
///          to preclude race conditions from occurring when an error message
///          is reported and program execution is terminated.
void initialize() {

  // Initialize OpenMP mutex lock for error messages with exceptions
  omp_init_lock(&_error_lock);
}


void set_logging_ranks(std::set<int> ranks) {
  _logging_ranks = std::move(ranks);
}


std::set<int> get_logging_ranks() {
  return _logging_ranks;
}


/// \details If the directory does not exist, it creates it for the user.
void set_path(const std::string &path) {

  auto pos = path.find_last_of("/");

  if (pos != std::string::npos) {
    // Create the directory as needed
    std::string dir = path.substr(0, pos + 1);
    fileutils::createDirectory(dir);
  }

  _log_path = path;

}


std::string get_path() {
  return _log_path;
}


/// \brief Sets the character to be used when printing separator log messages.
/// \param c the character for separator log messages
void set_separator_character(char c) {
  _separator_char = c;
}


/// \brief Returns the character used to format separator log messages.
/// \return the character used for separator log messages
char get_separator_character() {
  return _separator_char;
}


/// \brief Sets the character to be used when printing header log messages.
/// \param c the character for header log messages
void set_header_character(char c) {
  _header_char = c;
}


/// \brief Returns the character used to format header type log messages.
/// \return the character used for header log messages
char get_header_character() {
  return _header_char;
}


/// \brief Sets the character to be used when printing title log messages.
/// \param c the character for title log messages
void set_title_character(char c) {
  _title_char = c;
}


/// \brief Returns the character used to format title log messages.
/// \return the character used for title log messages
char get_title_character() {
  return _title_char;
}


/// \details Messages longer than this amount will be broken up into
/// \        multiline messages.
/// \param length the maximum log message line length in characters
void set_line_length(size_t length) {
  _line_length = length;
}


size_t get_line_length() {
  return _line_length;
}


size_t get_default_line_length() {
  return _default_line_length;
}


/// \brief Sets the minimum log message level which will be printed to the
/// \      console and to the log file.
/// \param new_level the minimum logging level as a character array
void set_level(const std::string &new_level) {

  _log_level = level::_from_string_nocase(new_level.c_str());

  log::verbose_once("Logging level set to {}", new_level);
}


/// \brief Sets the minimum log message level which will be printed to the
///        console and to the log file. This is an overloaded version to handle
///        a log::level type input.
/// \param new_level the minimum logging level as an int(or enum type log::level)
void set_level(level new_level) {
  _log_level = new_level;
}


/// \brief Return the minimum level for log messages
/// \return the minimum level for log messages
level get_level() {
  return _log_level;
}


/// \param level The logging level for this message.
/// \param message The formatted message.
std::string create_log_string(level log_level, const std::string &message) {

  // If the message is too long for a line or is a multiline message,
  // split into many log messages.
  auto get_message_lines = [](const std::string &level_prefix,
                              const std::string &msg)
    {
      if (msg.length() > _line_length ||
          msg.find_last_of("\n") != std::string::npos)
        return create_multiline_msg(level_prefix, msg);

      // Puts message on single line
      else
        return level_prefix + msg + "\n";
    };

  // Append the process id to the string
  auto append_proc_id = [](std::string &level_prefix)
    {
      #ifdef ENABLE_MPI_
      level_prefix = fmt::format("{}Rank {}: ", level_prefix, mpi::getMPIRank());
      #endif
    };


  int prefix_width = 8;
  std::string msg_string;

  // Create a prefix for each line
  auto level_prefix = stringutils::toUpper(log_level._to_string());

  // Trim the prefix
  auto underscore = level_prefix.find_first_of("_");
  if (underscore != std::string::npos)
    level_prefix = level_prefix.substr(0, underscore);

  level_prefix = fmt::format("[{: ^{}}]  ", level_prefix, prefix_width);

  // Append the log level to the message
  switch (log_level) {

    // Logging messages per node
    case level::debug :
    case level::profile :
    case level::verbose :
    case level::node :
    case level::warn :
    case level::critical :
    case level::test :
    case level::error :
      {
        append_proc_id(level_prefix);
        msg_string = get_message_lines(level_prefix, message);
        break;
      }

    // Logging messages by the root
    case level::profile_once :
    case level::verbose_once :
    case level::info :
    case level::warn_once :
    case level::result :
      {
        msg_string = get_message_lines(level_prefix, message);
        break;
      }

    // Other messages
    case level::separator :
      {
        // Reset the prefix
        level_prefix = fmt::format("[{: ^{}}]  ", "SP", prefix_width);
        auto c = _separator_char;

        // Take the first character of the message as the separator if possible
        if (!message.empty())
          c = message[0];

        msg_string = level_prefix + std::string(_line_length, c) + "\n";
        break;
      }

    // Temporarily make them the same
    case (level::header):
    case (level::title):
      {
        std::string border = level_prefix + std::string(_line_length, _title_char);
        msg_string = fmt::format("{0}\n{1}{3: ^{2}}\n{0}\n",
                                 border, level_prefix, _line_length, message);
        break;
      }
  } // end of switch

  return msg_string;
}


/// \details Print a log message. Only the allowed processes could print
///          the message.
void print_formatted_string(level log_level, const std::string &message) {

  // Always show error messages
  if (log_level == +level::error) {

    auto msg = create_log_string(log_level, message);
    print_formatted_error(msg);
  }

  // Other levels
  else if (log_level >= _log_level) {

    bool is_allowed = true;

    switch (log_level) {

      // Only printed by the root
      case level::profile_once :
      case level::verbose_once :
      case level::info :
      case level::warn_once :
      case level::result :
      case level::separator :
      case level::header:
      case level::title:
        if (!mpi::isMPIRoot()) {
          is_allowed = false;
        }

      // Printed by all processes
      default :
        break;
    };

    // Overwrite the flag for allowed ranks
    if (!_logging_ranks.empty())
      is_allowed = _logging_ranks.count(mpi::getMPIRank());

    if (is_allowed) {
      auto msg = create_log_string(log_level, message);

      // Print the message to stdout
      fmt::print("{}", msg);
      fflush(stdout);

      // Print the message to file
      // FIXME: need to be adapted with multiprocessing
      if (mpi::isMPIRoot())
        print_string_to_file(msg);
    }
  }
}


/// \details The error message will be logged both to the file and stdout.
///          Then, an exception will be throwed.
void print_formatted_error(const std::string &message) {

  // Print the message to file
  print_string_to_file(message);

  // Write the log message to the shell
  omp_set_lock(&_error_lock);
#ifdef ENABLE_MPI_
  if (mpi::isDomainDecomposed()) {
    fmt::print("{}", message);
    fflush(stdout);
    MPI_Abort(mpi::getMPIComm(), 1);
  }
#endif
  omp_unset_lock(&_error_lock);
  throw std::logic_error(message);
}


/// \details Write the string to the log file if its filename is set.
///          The file will be opened before writing and be closed after
///          writing. A timestamp will be at the beginning of the log file.
/// \param msg A string to be written.
void print_string_to_file(const std::string &msg) {

  // If the log file is specified, write log messages to it
  if (!_log_path.empty()) {

    // Open the file
    std::ofstream f;
    f.open(_log_path, std::ios::app);

    // If this is our first time logging, add a header with date, time
    if (!_logging) {

      // Append date, time to the top of the log file
      std::time_t rawtime = std::time(nullptr);
      f << fmt::format("Current local date and time: {:%F %T}\n", *std::localtime(&rawtime));
      _logging = true;
    }

    f << msg;
    f.close();
  }
}


/// \brief Breaks up a message which is too long for a single line into a
///        multiline message.
/// \details This is an internal function which is called by log_printf and
///          should not be called directly by the user.
/// \param level a string containing log level prefix
/// \param message a string containing the log message
/// \return a string with a formatted multiline message
std::string create_multiline_msg(std::string level, std::string message) {

  auto len_oneline = _line_length;
  auto size = message.length();

  std::string::size_type start = 0;
  std::string::size_type end = len_oneline;

  std::string msg_string;

  // Loop over msg creating substrings for each line
  while (end < size + len_oneline) {

    // Append log level to the beginning of each line
    msg_string += level;

    // Begin multiline messages with ellipsis
    if (start != 0)
      msg_string += "... ";

    // Find the current full length substring for line
    auto substring = message.substr(start, end - start);

    // If there is a newline in the substring
    auto cur_newline = substring.find_first_of("\n");
    if (cur_newline != std::string::npos) {
      // Erase anything after the newline and leave them to the next line
      auto n_rest = substring.size() - cur_newline;
      substring.erase(cur_newline, n_rest);
      // The first character after the newline
      end = std::min(end, size) - n_rest + 1;
    }
    else {
      // Truncate substring to last complete word
      if (end < size-1) {
        int endspace = substring.find_last_of(" ");
        if (message.at(endspace+1) != ' ' &&
            endspace != int(std::string::npos)) {
          end -= len_oneline - endspace;
          substring = message.substr(start, end-start);
        }
      }
    }

    // concatenate substring to output message
    msg_string += substring + "\n";

    // Reduce line length to account for ellipsis prefix
    if (start == 0)
      len_oneline -= 4;

    // Update substring indices
    start = end;
    end += len_oneline + 1;
  }

  return msg_string;
}


} // namespace log

} // namespace antmoc
