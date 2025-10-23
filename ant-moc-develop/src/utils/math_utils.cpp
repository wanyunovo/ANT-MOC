#include "antmoc/math_utils.h"
#include "antmoc/log.h"
#include "antmoc/openmp_utils.h"

#include <type_traits>

namespace antmoc {

inline namespace mathutils {


/// \details This method should work even if there are infinities.
template <typename T>
bool definitelyEqual(T x, T y) {
  return (x <= y + FLT_EPSILON) && (y <= x + FLT_EPSILON);
}

// Explicit instantiation
template bool definitelyEqual(float, float);
template bool definitelyEqual(double, double);


long gcd(long a, long b) {

  if (a <= 0 || b <= 0) {
    log::error("Can't compute GCD of non-positive numbers");
  }

  long r;
  while ((r = a % b) != 0) {
    a = b;
    b = r;
  }
  return b;
}


/// \details Different types of elements will be promoted.
template <typename U, typename V>
double dot(const std::vector<U> &vec1, const std::vector<V> &vec2) {

  if (vec1.empty() || vec2.empty()) {
    log::error("Dot production must take 2 non-empty vectors");
  }
  else if (vec1.size() != vec2.size()) {
    log::error("Dot production must take 2 vectors of the same shape: {} {}",
                      vec1.size(), vec2.size());
  }

  double result = 0;
  for (size_t i = 0; i < vec1.size(); ++i)
    result += vec1[i] * vec2[i];

  return result;
}

// Explicit instantiation
template double dot(const std::vector<double>&, const std::vector<double>&);


/// \details Empty vector will be treated as an error
template <typename T>
T sum(const std::vector<T> &vec) {
  if (vec.empty())
    log::error("Failed to call sum(...) on an empty vector");

  T sum = 0;
#pragma omp parallel for reduction(+:sum) schedule(static)
  for (size_t i = 0; i < vec.size(); ++i) {
    sum += vec[i];
  }

  return sum;
}

// Explicit instantiation
template double sum(const std::vector<double> &);


/// \brief Get the maximum element of a vector
template <typename T>
T max(const std::vector<T> &vec) {
  if (vec.empty())
    log::error("Failed to call max(...) on an empty vector");

  T maximum = vec[0];
#pragma omp parallel for reduction(max:maximum) schedule(static)
  for (size_t i = 0; i < vec.size(); ++i) {
    if (maximum < vec[i])
      maximum = vec[i];
  }

  return maximum;
}

// Explicit instantiation
template float max(const std::vector<float> &);
template double max(const std::vector<double> &);
template long max(const std::vector<long> &);


/// \details Empty vector will be treated as an error
template <typename T>
T min(const std::vector<T> &vec) {
  if (vec.empty())
    log::error("Failed to call min(...) on an empty vector");

  T minimum = vec[0];
#pragma omp parallel for reduction(min:minimum) schedule(static)
  for (size_t i = 0; i < vec.size(); ++i) {
    if (minimum > vec[i])
      minimum = vec[i];
  }

  return minimum;
}

// Explicit instantiation
template float min(const std::vector<float> &);
template double min(const std::vector<double> &);
template long min(const std::vector<long> &);


/// \details Supports floating-point numbers.
///          Empty vector will be treated as an error
template <typename T>
T mean(const std::vector<T> &vec) {
  static_assert(std::is_floating_point<T>::value,
                "mathutils::mean requires floating-point values");

  if (vec.empty())
    log::error("Failed to call mean(...) on an empty vector");

  return mathutils::sum(vec) / vec.size();
}

// Explicit instantiation
template float mean(const std::vector<float> &vec);
template double mean(const std::vector<double> &vec);


/// \details Supports floating-point numbers.
///          Empty vector will be treated as an error
template <typename T>
T stddev(const std::vector<T> &vec) {
  static_assert(std::is_floating_point<T>::value,
                "mathutils::stddev requires floating-point values");

  if (vec.empty())
    log::error("Failed to call stddev(...) on an empty vector");

  T dev = 0.;
  T maximum = max(vec);
  T mean_value = mean(vec);
  auto n = vec.size();
#pragma omp parallel for reduction(+:dev) schedule(static)
  for (size_t i = 0; i < n; ++i) {
    T tmp = (vec[i] - mean_value) / maximum;
    dev += tmp * tmp;
  }
  dev = std::sqrt(dev / n) * maximum;

  return dev;
}

// Explicit instantiation
template float stddev(const std::vector<float> &vec);
template double stddev(const std::vector<double> &vec);


template <typename T>
T rms(const std::vector<T> &v1, const std::vector<T> &v2) {
  static_assert(std::is_floating_point<T>::value,
                "mathutils::stddev requires floating-point values");

  if (v1.size() != v2.size())
    log::error("Failed to call rms(...): the vectors are of different size");
  else if (v1.empty())
    log::error("Failed to call rms(...) on empty vectors");

  T error = 0.;
  auto n = v1.size();
  for (size_t i = 0; i < n; ++i) {
    T tmp = v1[i] - v2[i];
    error += tmp * tmp;
  }
  error = std::sqrt(error / n);

  return error;
}

// Explicit instantiation
template float rms(const std::vector<float> &v1, const std::vector<float> &v2);
template double rms(const std::vector<double> &v1, const std::vector<double> &v2);


} // inline namespace mathutils
} // namespace antmoc
