/// \file include/toml_utils.h
/// \date Oct 18, 2020
/// \author An Wang, USTB (wangan@xs.ustb.edu.cn)

#ifndef TOML_UTILS_H_
#define TOML_UTILS_H_

#include <string>

#include <toml.hpp>

namespace antmoc {

///---------------------------------------------------------------------
/// \namespace tomlutils
/// \details Useful TOML interfaces for antmoc.
///---------------------------------------------------------------------
namespace tomlutils {

  /// \brief Converts a toml::value to std::string.
  std::string toString(const toml::value &v);


} // namespace tomlutils

} // namespace antmoc

#endif  // TOML_UTILS_H_
