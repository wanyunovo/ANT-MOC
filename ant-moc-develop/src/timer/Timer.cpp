#include "antmoc/Timer.h"
#include "antmoc/constants.h"
#include "antmoc/log.h"

#include <unistd.h>
#include <omp.h>

#include <cmath>
#include <fstream>
#include <iomanip>
#include <ios>
#include <iostream>
#include <sstream>
#include <time.h>
#include <utility>

namespace antmoc
{

std::map<std::string, double> Timer::_timer_splits;
std::vector<double> Timer::_start_times;


/**
 * @brief Starts the Timer.
 * @details This method is similar to starting a stopwatch.
 */
void Timer::startTimer() {

  double start_time = omp_get_wtime();
  _start_times.push_back(start_time);
  _running = true;

  return;
}


/**
 * @brief Stops the Timer.
 * @details This method is similar to stopping a stopwatch.
 */
void Timer::stopTimer() {

  if (_running) {

    double end_time = omp_get_wtime();
    double start_time = _start_times.back();

    _elapsed_time = end_time - start_time;

    if (_start_times.empty())
      _running = false;

    _start_times.pop_back();
  }

  return;
}


/**
 * @brief Stops the Timer and records a message.
 * @details This method is similar to stopping a stopwatch.
 * @param msg a msg corresponding to this time split
 */
void Timer::stopTimer(const std::string &msg) {

  if (_running) {

    double end_time = omp_get_wtime();
    double start_time = _start_times.back();

    _elapsed_time = end_time - start_time;

    if (_start_times.empty())
      _running = false;

    _start_times.pop_back();
  }

  float &time = _elapsed_time;
  std::string msg_string = std::string(msg);

  if (_timer_splits.find(msg_string) != _timer_splits.end())
    _timer_splits.at(msg_string) += time;
  else
    _timer_splits.insert({msg_string, time});

  return;
}


/**
 * @brief Records a message corresponding to a time for the current split.
 * @details When this method is called it assumes that the Timer has been
 *          stopped and has the current time for the process corresponding
 *          to the message.
 * @param msg a msg corresponding to this time split
 */
void Timer::recordSplit(const std::string &msg) {

  double time = getTime();
  std::string msg_string = std::string(msg);

  if (_timer_splits.find(msg_string) != _timer_splits.end())
    _timer_splits.at(msg_string) += time;
  else
    _timer_splits.insert({msg_string, time});
}


/// \brief Append a message and a corresponding time
void Timer::appendSplit(const std::string &msg, double time) {
  if (_timer_splits.count(msg))
    log::warn("Failed to append a time record '{}': already exist", msg);
  else if (time < 0)
    log::error("Failed to append a time record '{} {}': negative time", msg, time);
  else
    _timer_splits.insert({msg, time});
}


/**
 * @brief Returns the time elapsed from startTimer() to stopTimer().
 * @return the elapsed time in seconds
 */
double Timer::getTime() {
  return _elapsed_time;
}


/**
 * @brief Returns the time associated with a particular split.
 * @details If the split does not exist, returns 0.
 * @param msg the message tag for the split
 * @return the time recorded for the split (seconds)
 */
double Timer::getSplit(const std::string &msg) {

  std::string msg_string = std::string(msg);

  if (_timer_splits.find(msg_string) == _timer_splits.end())
    return 0.0;
  else
    return _timer_splits.at(msg_string);
}


/**
 * @brief Prints the time and message for the specified split
 * @details This method print a specified message for the split. If a level
 *          is given, the message will be indented. If a floating-point
 *          value is given, a percentage will be computed and printed.
 * @param msg the message tag for the split
 * @param caption the string to be printed
 * @param level indentation level
 * @param parent_split the parent split for computing percentages
 */
void Timer::printSplit(const std::string &msg, const std::string &caption,
                       int level, const std::string &parent_split) {

  std::string cap_str = std::string(caption);
  std::string conj = "";

  // Format the string with indentation
  int width = 51 - 2 * level;

  // Conjunction symbol
  if (level > 1)
    conj = std::string(level - 1, '-');

  double split = getSplit(msg);
  double denominator = getSplit(parent_split);
  if (denominator > 0) {
    auto percent = split / denominator * 100;
    log::result("{0:{1}}{2:.<{3}} {4:1.5E} s {5}{6:5.2f}%",
                 "", 2 * level, cap_str, width, split, conj, percent);
  }
  else
    log::result("{0:{1}}{2:.<{3}} {4:1.5E} s",
                 "", 2 * level, cap_str, width, split);

}


/// \brief Prints the time and message for the specified split
/// \details The split name will be used as the caption
void Timer::printSplit(const std::string &msg, int level,
                       const std::string &parent_split) {
  printSplit(msg, msg, level, parent_split);
}


/**
 * @brief Prints the times and messages for each split to the console.
 * @details This method will loop through all of the Timer's splits and print a
 *          formatted message string (80 characters in length) to the console
 *          with the message and the time corresponding to that message.
 */
void Timer::printSplits() {

  std::string curr_msg;
  double curr_split;
  std::map<std::string, double>::iterator iter;

  for (iter = _timer_splits.begin(); iter != _timer_splits.end(); ++iter) {

    std::stringstream formatted_msg;

    curr_msg = (*iter).first;
    curr_split = (*iter).second;

    curr_msg.resize(53, '.');
    formatted_msg << curr_msg;

    log::result("{}{:1.5E} s", formatted_msg.str(), curr_split);
  }
}


/**
 * @brief Clears the time split for this message and deletes the message's
 *        entry in the Timer's splits log.
 * @param msg the message tag for the split
 */
void Timer::clearSplit(const std::string &msg) {

  std::string msg_string = std::string(msg);

  if (_timer_splits.find(msg_string) == _timer_splits.end())
    return;
  else
    _timer_splits.erase(msg_string);
}


/**
 * @brief Clears all times split messages from the Timer.
 */
void Timer::clearSplits() {
  _timer_splits.clear();
}


//TODO This function is not really about time, it could be moved elsewhere.
/**
 * @brief Read memory usage file (on a HPC installation), and process it to make
 *        it more readable. Used for profiling.
 * @param vm_usage total use of virtual memory
 * @param resident_set total use of resident memory
 */
void Timer::processMemUsage(double& vm_usage, double& resident_set) {

   vm_usage     = 0.0;
   resident_set = 0.0;

   /* Open file containing memory info */
   std::ifstream stat_stream("/proc/self/stat", std::ios_base::in);

   /* Read in dummy data */
   std::string tmp;
   for (int i=0; i < 22; i++)
     stat_stream >> tmp;

   /* Read in virtual and resident memory */
   unsigned long vsize;
   stat_stream >> vsize;
   long rss;
   stat_stream >> rss;
   stat_stream.close();

   /* Calculate memory usage */
   long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024;
   vm_usage = (double) vsize / 1024.0 / 1024.0;
   resident_set = rss * page_size_kb / 1024.0;
}


/**
 * @brief Transfer timer data across all domains.
 * @param comm a MPI communicator to transfer data
 */
#ifdef ENABLE_MPI_
void Timer::reduceTimer(MPI_Comm comm) {

  auto size = _timer_splits.size();
  double *splits = new double[size]();
  double *total_splits = new double[size]();

  int n = 0;
  for (auto &v : _timer_splits) {
    splits[n++] = v.second;
  }

  /* Collapse timing results down to one value for each category */
  MPI_Reduce(splits, total_splits, size, MPI_DOUBLE, MPI_SUM, 0, comm);


  /* On the main node, average over the number of ranks, update result */
  n = 0;
  for (auto &v : _timer_splits) {
    v.second = total_splits[n++];
    if (fabs(v.second) > FLT_EPSILON) {
      v.second /= mpi::getNumProcs();
    }
  }

  delete [] splits;
  delete [] total_splits;

}
#endif

} /* namespace antmoc */
