/// \file include/container_utils.h
/// \brief Container utility
/// \author OpenMC

#ifndef CONTAINER_UTILS_H_
#define CONTAINER_UTILS_H_

#include <algorithm>
#include <iterator>
#include <sstream>
#include <unordered_map>
#include <vector>

namespace antmoc
{

using DoubleVec = std::vector<double>;
using IntVec    = std::vector<int>;
using IntVec2D  = std::vector<IntVec>;


/// \brief Determine if an element is in the container
/// \details This template uses standard library functions. Thus,
///          It can be instantiated with array types.
/// \param v a container of type C
/// \param x an element of type T
template<class C, class T>
inline bool contains(const C& v, const T& x)
{
  return std::end(v) != std::find(std::begin(v), std::end(v), x);
}


/// \brief Retrieve keys from a map
template<class K, class V>
inline std::vector<K>
mapExtractKeys(const std::unordered_map<K, V> &map) {
  std::vector<K> keys;
  for (const auto &e : map)
    keys.push_back(e.first);

  return keys;
}


/// \brief Convert an 1-D vector into string
/// \param vec       Vector with elements of type T.
/// \param n_oneline Number of characters per line.
template <typename T>
std::string vecToString(std::vector<T> vec, size_t n_oneline = 1000) {
  std::ostringstream ss;

  ss << "[";
  for (size_t i = 0; i < vec.size(); ++i) {
    if (i % n_oneline == 0 && i > 0) {
      ss << '\n';
    }
    ss << vec[i] << (i < vec.size() - 1 ? " " : "");
  }
  ss << "]";

  return ss.str();
}


/// \brief Convert a 2-D vector into string
template <typename T>
std::string vecToString(std::vector<std::vector<T>> vec) {
  std::ostringstream ss;

  ss << "[";
  for (size_t i = 0; i < vec.size(); ++i) {
    ss << vecToString(vec[i]) << (i < vec.size() - 1 ? "\n\n": "");
  }
  ss << "]";

  return ss.str();
}


} // namespace antmoc

#endif  // CONTAINER_UTILS_H_
