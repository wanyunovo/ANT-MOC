#include "antmoc/Option.h"

#include <sstream>


namespace antmoc {


/// \details Use setArguments to safely set the underlying buffer.
Option::Option(int &argc, char **argv) {
  setArguments(argc, argv);
}


/// \details Use setArguments to safely set the underlying buffer.
Option::Option(const StringVec &argv) {
  setArguments(argv);
}


/// \details Use setArguments to safely set the underlying buffer.
Option::Option(std::initializer_list<std::string> ilist) {
  setArguments(ilist);
}


/// \details To parse the options, we first copy them to a new buffer.
cxxopts::ParseResult Option::parseUnderlying() {

  // Allocate memory for buffer
  int argc = _argv.size();
  char **argv = new char*[argc];

  for (int i = 0; i < argc; ++i) {
    // Since c++11, data() returns a null-terminated array just like c_str()
    argv[i] = const_cast<char*>( _argv[i].data() );
  }

  // Parse arguments
  // We are safe as long as cxxopts doesn't change the underlying data
  // of argv. In the current implementation, it has to change pointers
  // like argv[i] to skip some arguments but not the underlying data.
  auto result = _options.parse(argc, argv);

  // Free the buffer
  delete[] argv;

  return result;
}


bool Option::hasOption(const std::string &option) {

  auto parse_result = parseUnderlying();
  return parse_result.count(option);
}


std::string Option::getOptionValue(const std::string &option) {

  static const std::string empty_string("");
  auto parse_result = parseUnderlying();

  try {
    return parse_result[option].as<std::string>();
  }
  catch (...) {
    return empty_string;
  }
}


double Option::getOptionValueDouble(const std::string &option) {

  auto parse_result = parseUnderlying();
  return parse_result[option].as<double>();
}


int Option::getOptionValueInt(const std::string &option) {

  auto parse_result = parseUnderlying();
  return parse_result[option].as<int>();
}


size_t Option::getOptionValueSizet(const std::string &option) {

  auto parse_result = parseUnderlying();
  return parse_result[option].as<size_t>();
}


bool Option::getOptionValueBool(const std::string &option) {

  auto parse_result = parseUnderlying();
  return parse_result[option].as<bool>();
}


void Option::setArguments(int &argc, char **argv) {
  _argv.clear();
  for (int i = 0; i < argc; ++i)
    _argv.push_back(argv[i]);

  // Expand arguments for derived classes
  expandArguments();
}


void Option::setArguments(const StringVec &argv) {
  resetArguments();
  _argv.insert(_argv.end(), argv.begin(), argv.end());

  // Expand arguments for derived classes
  expandArguments();
}


void Option::setArguments(std::initializer_list<std::string> ilist) {
  StringVec vec {ilist};
  setArguments(vec);

  // Expand arguments for derived classes
  expandArguments();
}


void Option::appendArguments(const StringVec &argv) {
  _argv.insert(_argv.end(), argv.begin(), argv.end());
}


void Option::insertFrontArguments(const StringVec &argv) {
  auto first = _argv.front();
  _argv.pop_front();

  for (auto it = argv.rbegin(); it != argv.rend(); ++it)
    _argv.push_front(*it);

  _argv.push_front(first);
}


void Option::insertFrontArguments(std::initializer_list<std::string> ilist) {
  StringVec vec {ilist};
  insertFrontArguments(vec);
}


void Option::resetArguments() {
  _argv.clear();
  _argv.push_back("antmoc");
}


std::string Option::toString() const {
  std::stringstream ss;
  ss << "Options:";
  for (const auto &arg : _argv)
    ss << " '" << arg << "'";

  return ss.str();
}


std::ostream& operator<<(std::ostream &os, const Option &option) {
  return os << option.toString();
}


} // namespace antmoc
