#include "antmoc/toml_utils.h"
#include "antmoc/log.h"
#include "antmoc/string_utils.h"

#include <iomanip>
#include <sstream>

namespace antmoc {

namespace tomlutils {

  std::string toString(const toml::value &v) {
    std::ostringstream ss;

    if (v.is_boolean())
      ss << std::boolalpha << toml::get<bool>(v);

    else if (v.is_integer())
      ss << toml::get<int>(v);

    else if (v.is_floating())
      ss << toml::get<double>(v);

    else if (v.is_string())
      ss << toml::get<std::string>(v);

    else
      log::error("Failed to convert TOML type to C++ type");

    return ss.str();
  }

} // namespace tomlutils

} // namespace antmoc
