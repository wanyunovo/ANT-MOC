/// \file Option.h
/// \brief An abstract class for option parsers
/// \date April 17, 2020
/// \author An Wang, USTB (wangan@xs.ustb.edu.cn)

#ifndef OPTION_H_
#define OPTION_H_

#include <deque>
#include <initializer_list>
#include <ostream>

#include <cxxopts.hpp>

#include "antmoc/string_utils.h"

namespace antmoc {


///---------------------------------------------------------------------
/// \class Option
/// \brief An abstract class for option parsers
/// \details An option parser must provide several basic interfaces
///          such as:
///            parse arguments
///            check existence of an option
///            retrieve the value of an option
///          and extra interfaces to help to manipulate the underlying
///          arguments:
///            set arguments
///            append arguments
///            reset arguments to defaults
///
///          It is not the parser's responsibility to validate the
///          arguments and print help messages.
///---------------------------------------------------------------------
class Option {

private:

  /// \brief Parses arguments with cxxopts.
  cxxopts::ParseResult parseUnderlying();

protected:

  StringDeque _argv {"antmoc"};             ///< Underlying buffer
  cxxopts::Options _options {"antmoc", ""}; ///< Underlying option parser

  StringDeque& getArgv() { return _argv; }
  cxxopts::Options& getParser() { return _options; }

  //--------------------------------------
  // Initialize options
  //--------------------------------------
  /// \brief Define options, to be implemented and called by derived classes.
  virtual void initializeOptions() = 0;

public:

  //--------------------------------------
  // Constructors
  //--------------------------------------
  /// \brief Default constructor.
  Option() = default;

  /// \brief Creates an object with standard CLI arguments.
  /// \param argc The number of arguments (at least 1).
  /// \param argv The arguments, including the executable name.
  Option(int &argc, char **argv);

  /// \brief Creates an object with a vector of arguments.
  /// \param argv A vector of arguments, not including the executable name
  Option(const StringVec &argv);

  /// \brief Creates an object with a list of arguments.
  /// \param ilist A list of arguments, not including the executable name
  Option(std::initializer_list<std::string> ilist);

  /// \brief Default destructor.
  virtual ~Option() = default;

  //--------------------------------------
  // Check existence
  //--------------------------------------
  /// \brief Check the existence of an option.
  bool hasOption(const std::string &);

  //--------------------------------------
  // Retrieve option values
  //--------------------------------------
  /// \brief Returns the value of a string option
  std::string getOptionValue(const std::string &);

  /// \brief Returns the value of a double option
  double getOptionValueDouble(const std::string &);

  /// \brief Returns the value of an int option
  int getOptionValueInt(const std::string &);

  /// \brief Returns the value of an size_t option
  size_t getOptionValueSizet(const std::string &);

  /// \brief Returns the value of a boolean option
  bool getOptionValueBool(const std::string &);

  //--------------------------------------
  // Manipulating arguments
  //--------------------------------------
  /// \brief Expand options. this is used by derived classes to expand
  ///        their own options, such as XML options.
  virtual void expandArguments() { return; }

  /// \brief Sets the arguments to new values. The interface expandArguments()
  ///        will be called after that.
  /// \param argc The number of arguments to be inserted.
  /// \param argv An array of string arguments.
  void setArguments(int &argc, char **argv);

  /// \brief Sets the arguments to new values. The interface expandArguments()
  ///        will be called after that.
  /// \param argv A vector of arguments.
  void setArguments(const StringVec &argv);

  /// \brief Sets the arguments to new values. The interface expandArguments()
  ///        will be called after that.
  /// \param ilist An initializer list of arguments.
  void setArguments(std::initializer_list<std::string> ilist);

  /// \brief Appends arguments to the end of the argument list.
  /// \param argv A vector of arguments.
  void appendArguments(const StringVec &argv);

  /// \brief Inserts arguments from the start of the argument list.
  /// \param argv A vector of arguments.
  void insertFrontArguments(const StringVec &);

  /// \brief Inserts arguments from the start of the argument list.
  /// \param ilist An initializer list of arguments.
  void insertFrontArguments(std::initializer_list<std::string> ilist);

  /// \brief Resets the argument list to default
  void resetArguments();

  //--------------------------------------
  // Output the object as a string
  //--------------------------------------
  /// \brief Converts an Option object into a string
  virtual std::string toString() const;

  /// \brief Returns the string representation of the object to operator <<
  friend std::ostream& operator<<(std::ostream &os, const Option &option);

};


} // namespace antmoc

#endif  // OPTION_H_
