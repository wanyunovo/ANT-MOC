/// \file include/debug_utils.h
/// \brief Utilities for debugging
/// \author An Wang, USTB (wangan@xs.ustb.edu.cn)

#ifndef DEBUG_UTILES_H_
#define DEBUG_UTILES_H_

#include <unistd.h>
#include "antmoc/log.h"

namespace antmoc
{

/// \brief Print the pid and hostname as profiling messages.
inline void printHostnames() {
  char hostname[256];
  gethostname(hostname, sizeof(hostname));
  log::profile("PID {} is running on {}", getpid(), hostname);
  fflush(stdout);
}

} // namespace antmoc

#endif  // DEBUG_UTILES_H_
