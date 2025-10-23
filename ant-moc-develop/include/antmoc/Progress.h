/// \file Progress.h
/// \brief A progress object
/// \date January 11, 2016
/// \author Samuel Shaner, MIT, Course 22 (shaner@mit.edu)

#ifndef PROGRESS_H_
#define PROGRESS_H_

#include <string>
#include <vector>

namespace antmoc {


// Forward declaration
class Geometry;

class Progress {

private:

  /** A list of lists representing the vector */
  std::string _name;
  long _counter;
  long _num_iterations;
  size_t _curr_interval;
  std::vector<long> _intervals;
  Geometry* _geometry;
  bool _mpi_comm;

public:
  Progress(long num_iterations, std::string name, double interval=0.1,
           Geometry* geometry=NULL, bool mpi_comm=false);
  virtual ~Progress() = default;

  /* Worker functions */
  void incrementCounter();
  void reset();
};

} /* namespace antmoc */

#endif /* PROGRESS_H_ */
