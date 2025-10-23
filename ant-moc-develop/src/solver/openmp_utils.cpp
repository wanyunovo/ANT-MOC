#include "antmoc/openmp_utils.h"
#include "antmoc/log.h"
#include "antmoc/PyVector.h"

#include <algorithm>  /* sort */
#include <sstream>
#include <string>

#include <omp.h>
#include <sched.h>    /* sched_getcpu */

namespace antmoc {

namespace openmp {

/// \details This method only works for OpenMP >= 4.0.
void print_affinity_policy() {
#if _OPENMP >= 201307

  // Get the affinity policy
  auto policy = omp_get_proc_bind();
  std::string policy_str;

  switch (policy) {
    case 0 : policy_str = "false"; break;
    case 1 : policy_str = "true"; break;
    case 2 : policy_str = "master"; break;
    case 3 : policy_str = "close"; break;
    case 4 : policy_str = "spread"; break;
    default : policy_str = "unknown";
  }

  log::profile("Global affinity policy on current process: {}", policy_str);

#endif
}


/// \details This method only works for OpenMP >= 4.5.
void print_places() {
#if _OPENMP >= 201511

  if (omp_in_parallel()) {
    log::error("Method openmp::print_places() cannot be invoked"
               "within a parallel region");
  }

  // Number of places
  int num_places = omp_get_num_places();

  log::profile("OpenMP: number of places = {}", num_places);

  for (int i = 0; i < num_places; ++i) {
    // Number of processors
    int n = omp_get_place_num_procs(i);

    // Processor ids
    int *ids_buf = new int[n];
    omp_get_place_proc_ids(i, ids_buf);

    PyVector<int> proc_ids(ids_buf, ids_buf + n);

    delete[] ids_buf;

    log::profile("OpenMP: {} processors on place {} = [{}]", n, i, proc_ids);
  }

#endif
}


/// \details Print CPU assignments. The CPU numbers is related to the NUMA
///          nodes managed by the OS. This method will show the CPUs running
///          OpenMP threads rather than the mapping between them.
void print_cpu_bind() {

  if (omp_in_parallel()) {
    log::error("Method openmp::print_cpu_bind() cannot be invoked"
               "within a parallel region");
  }

  auto num_threads = omp_get_max_threads();

  PyVector<int> cpus;
  cpus.reserve(num_threads);

#pragma omp parallel for schedule(static)
  for (int i = 0; i < num_threads; ++i) {
#   pragma omp critical
    {
      cpus.push_back(sched_getcpu());
    }
  }
  std::sort(cpus.begin(), cpus.end());

  log::profile("CPUs occupied by {} threads: [{}]", num_threads, cpus);
}

} // namespace openmp

} // namespace antmoc
