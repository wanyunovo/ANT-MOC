/// \file include/string_utils.h
/// \brief String utility
/// \date Aug 26, 2019
/// \author An Wang, USTB (wangan@xs.ustb.edu.cn)

#ifndef STRING_UTILS_H_
#define STRING_UTILS_H_

#include <deque>
#include <set>
#include <string>
#include <vector>

#include <fmt/format.h>

namespace antmoc
{


// A vector of strings
using StringVec = std::vector<std::string>;

// A double-ended queue of strings
using StringDeque = std::deque<std::string>;


///---------------------------------------------------------------------
/// \namespace stringutils
/// \details Useful functions for antmoc. Most of them are intended to
///          be in-place algorithms. That is, the source string will be
///          modified. To prevent that, just prefix variable declaration
///          with const to use the overloaded version of these algorithms.
///---------------------------------------------------------------------
namespace stringutils {


  //--------------------------------------------------------------------
  // Case conversion
  //--------------------------------------------------------------------
  /// \brief Convert characters of a string to uppercase (in-place)
  /// \param s a reference to the string
  std::string &toUpper(std::string &s);

  /// \brief Convert characters of a string to uppercase
  /// \param s a constant reference to the string
  std::string toUpper(const std::string &s);

  /// \brief Convert characters of a string to lowercase (in-place)
  /// \param s a reference to the string
  std::string &toLower(std::string &s);

  /// \brief Convert characters of a string to lowercase
  /// \param s a constant reference to the string
  std::string toLower(const std::string &s);


  //--------------------------------------------------------------------
  // String spliting and concatenating
  //--------------------------------------------------------------------
  /// \brief Split a string into words in terms of the delimiter
  /// \param input a string to be split
  /// \param delimiter the delimiter used to split the string
  /// \return a vector of words
  StringVec splitString(std::string input, const std::string &delimiter);

  /// \brief Split a string into words by spaces
  /// \param input a string to be split
  /// \return a vector of words
  StringVec splitString(std::string input);

  /// \brief Join strings.
  /// \param c A container.
  /// \param delimiter The delimiter inserted between strings.
  template <typename Container>
  std::string join(const Container &c, const std::string &delimiter) {
    return fmt::format("{}", fmt::join(c, delimiter));
  }


  //--------------------------------------------------------------------
  // Trimming methods
  //--------------------------------------------------------------------
  /// \brief Trim from start (in-place)
  /// \brief s a reference to the string
  /// \return reference to the string
  std::string& ltrim(std::string &s);

  /// \brief Trim from start
  /// \param s the const lvalue or rvalue reference to a string
  std::string ltrim(const std::string &s);

  /// \brief Trim from end (in-place)
  /// \brief s a reference to the string
  /// \return reference to the string
  std::string& rtrim(std::string &s);

  /// \brief Trim from end
  /// \param s the const lvalue or rvalue reference to a string
  std::string rtrim(const std::string &s);

  /// \brief Trim from both ends (in place)
  /// \return reference to the string
  std::string& trim(std::string &s);

  /// \brief Trim from both ends
  /// \param s the const lvalue or rvalue reference to a string
  std::string trim(const std::string &s);


  //--------------------------------------------------------------------
  // String checking
  //--------------------------------------------------------------------
  /// \brief Determine whether the string is blank
  bool isSpaces(const std::string &s);


  //--------------------------------------------------------------------
  // Substitution
  //--------------------------------------------------------------------
  /// \brief Trim the string and replace spaces with underscores
  /// \param s a reference to the string
  std::string& spaceToUnderscore(std::string &s);

  /// \brief Trim the string and replace spaces with underscores
  /// \param s a const reference to the string
  std::string spaceToUnderscore(const std::string &s);

  /// \brief Trim the string and replace underscores with spaces
  /// \param s a reference to the string
  std::string& underscoreToSpace(std::string &s);

  /// \brief Trim the string and replace underscores with spaces
  /// \param s a const reference to the string
  std::string underscoreToSpace(const std::string &s);

  /// \brief Expand a interger range in string to a set of integers.
  /// \param s A range with start, stop and an optional step, e.g.
  ///          1:5 leads to 1,2,3,4,5,
  ///          5:1:-1 leads to 5,4,3,2,1
  ///          If the step is not given, it defaults to 1.
  /// \param delimiter The delimiter between start, stop and step
  ///                  (defaults to ':').
  /// \return A set of integers.
  std::set<int> toIntegerSet(std::string s, std::string delimiter = ":");

} // namespace stringutils

} // namespace antmoc

#endif  // STRING_UTILS_H_
