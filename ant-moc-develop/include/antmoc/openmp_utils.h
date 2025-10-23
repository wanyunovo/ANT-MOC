/// \file openmp_utils.h
/// \brief OpenMP utilities.
/// \date May 7, 2020
/// \author An Wang, USTB (wangan@xs.ustb.edu.cn)

#ifndef OPENMP_UTILS_H_
#define OPENMP_UTILS_H_

namespace antmoc {

namespace openmp {

/// \brief Print the affinity policy, OpenMP places, and thread-place mapping.
void print_affinity_policy();

/// \brief Print OpenMP places.
void print_places();

/// \brief Print the thread-cpu mapping.
void print_cpu_bind();

} // namespace openmp

} // namespace antmoc


#endif  // OPENMP_UTILS_H_
