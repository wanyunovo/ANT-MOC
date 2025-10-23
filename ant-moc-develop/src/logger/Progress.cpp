#include "antmoc/Progress.h"
#include "antmoc/Geometry.h"
#include "antmoc/log.h"
#include "antmoc/mpi_utils.h"

#include <omp.h>

namespace antmoc {


/// \brief Constructor for Progress.
Progress::Progress(long num_iterations, std::string name, double interval,
                   Geometry* geometry, bool mpi_comm) {

  if (!mpi::isDomainDecomposed())
    mpi_comm = false;

  if (geometry != nullptr && !mpi::isSpatialDecomposed())
    geometry = nullptr;

  _mpi_comm = mpi_comm;
  _geometry = geometry;
  _num_iterations = num_iterations;
  _name = name;
  _counter = 0;
  _curr_interval = 0;

  /* Based on interval fraction, set an integer number of intervals, and save 
   * each interval value in _intervals */
  int num_intervals = 1. / interval + 1;
  _intervals.resize(num_intervals);
  long interval_stride = interval * num_iterations;
  if (interval_stride == 0)
    interval_stride = 1;
  for (int i=0; i < num_intervals; i++)
    _intervals.at(i) = std::min(i * interval_stride, num_iterations-1);
  _intervals.at(num_intervals-1) = num_iterations-1;
}


/// \brief Increment the counter, print log if it has reached an interval bound.
void Progress::incrementCounter() {

  long curr_count;

  #pragma omp critical
  {
    /* Increment counter */
    curr_count = _counter++;

    /* While loop handles the case of an empty interval */
    while (_curr_interval < _intervals.size()) {

      /* Check if next interval is reached */
      if (curr_count != _intervals.at(_curr_interval))
        break;

      /* Add 1 for nicer printout */
      long count = curr_count + (curr_count != 0);
      double num_iters = _num_iterations;
      double percent = count / num_iters * 100.0;
#if defined (ENABLE_MPI_) && (ENABLE_DEBUG_)
      // Put a barrier here only if we are in DEBUG mode
      if (_mpi_comm) {
        mpi::mpiBarrier();
      }
#endif
      log::info("Progress {0}: {1:6.2f}%", _name, percent);
      _curr_interval++;
      if (_curr_interval >= _intervals.size())
        break;
    }
  }
}


/// \brief Reset the counter.
void Progress::reset() {
  _counter = 0;
  _curr_interval = 0;
}

} /* namespace antmoc */
