/// \file include/math_utils.h
/// \brief Math utility
/// \date Sep 11, 2019
/// \author An Wang, USTB (wangan@xs.ustb.edu.cn)

#ifndef MATH_UTILS_H_
#define MATH_UTILS_H_

#define _USE_MATH_DEFINES
#include <cmath>
#include <vector>

#include "antmoc/constants.h"  /* Redefine some macros */

namespace antmoc {

inline namespace mathutils {


/// \brief Copy the sign of a number
template <typename T>
inline T copysign(T x, T y) {
  return y < 0 ? -x : x;
}


/// \brief Indicates whether two numeric numbers are equivalent.
/// \details This method should work even if there are infinities.
template <typename T> bool definitelyEqual(T x, T y);
extern template bool definitelyEqual(float, float);
extern template bool definitelyEqual(double, double);


/// \brief Computes the greatest common divisor
long gcd(long a, long b);


/// \brief Computes the dot product
template <typename U, typename V>
double dot(const std::vector<U> &vec1, const std::vector<V> &vec2);
extern template double dot(const std::vector<double>&, const std::vector<double>&);


/// \brief Computes the sum of the vector elements
template <typename T> T sum(const std::vector<T> &vec);
extern template double sum(const std::vector<double> &);


/// \brief Returns the maximum element of a vector
template <typename T> T max(const std::vector<T> &vec);
extern template float max(const std::vector<float> &);
extern template double max(const std::vector<double> &);
extern template long max(const std::vector<long> &);


/// \brief Returns the minimum element of a vector
template <typename T> T min(const std::vector<T> &vec);
extern template float min(const std::vector<float> &);
extern template double min(const std::vector<double> &);
extern template long min(const std::vector<long> &);


/// \brief Computes the mean value of the vector elements
template <typename T> T mean(const std::vector<T> &vec);
extern template float mean(const std::vector<float> &vec);
extern template double mean(const std::vector<double> &vec);


/// \brief Computes the standard deviation of the vector elements
template <typename T> T stddev(const std::vector<T> &vec);
extern template float stddev(const std::vector<float> &vec);
extern template double stddev(const std::vector<double> &vec);

/// \brief Computes the RMS error between two vectors
template <typename T> T rms(const std::vector<T> &v1, const std::vector<T> &v2);
extern template float rms(const std::vector<float> &v1, const std::vector<float> &v2);
extern template double rms(const std::vector<double> &v1, const std::vector<double> &v2);

} // inline namespace mathutils
} // namespace antmoc

#endif  // MATH_UTILS_H_
