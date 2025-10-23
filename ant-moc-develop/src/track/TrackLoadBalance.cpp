#include "antmoc/TrackLoadBalance.h"
#include "antmoc/log.h"
#include "antmoc/math_utils.h"

#include <algorithm>
#include <numeric>

namespace antmoc {


///---------------------------------------------------------------------
/// Base class
///---------------------------------------------------------------------
/// \brief Constructor
TrackLoadBalance::TrackLoadBalance(int chain_azims,
                                   DoubleVec chain_lengths,
                                   IntVec chain_nums):
  _chain_azims(chain_azims),
  _lengths(chain_lengths),
  _chain_nums(chain_nums)
{}


/// \brief Constructor
TrackLoadBalance::TrackLoadBalance(int chain_azims,
                                   const double *chain_lengths,
                                   const long *chain_nums):
  _chain_azims(chain_azims),
  _lengths(chain_lengths, chain_lengths + chain_azims),
  _chain_nums(chain_nums, chain_nums + chain_azims)
{}


// Create a TrackLoadBalance object
TrackLoadBalancePtr TrackLoadBalance::
getTrackLoadBalancer(trackMappingType type,
                     int chain_azims,
                     DoubleVec chain_lengths,
                     IntVec chain_nums) {
  TrackLoadBalancePtr balancer;

  switch (type) {
    case trackMappingType::BLOCK :
      balancer = std::make_shared<BlockDistribution>(chain_azims,
                                                     chain_lengths,
                                                     chain_nums);
      break;
    case trackMappingType::CYCLIC_TRACK :
      balancer = std::make_shared<CyclicTrackDistribution>(chain_azims,
                                                           chain_lengths,
                                                           chain_nums);
      break;
    case trackMappingType::ANGLE :
      balancer = std::make_shared<AngleDistribution>(chain_azims,
                                                     chain_lengths,
                                                     chain_nums);
      break;
    default:
      log::error("Unrecoganized track mapping type");
  };

  return balancer;
}


// Create a TrackLoadBalance object
TrackLoadBalancePtr TrackLoadBalance::
getTrackLoadBalancer(trackMappingType type,
                     int chain_azims,
                     const double *chain_lengths,
                     const long *chain_nums) {
  return getTrackLoadBalancer(type, chain_azims,
                              DoubleVec(chain_lengths, chain_lengths + chain_azims),
                              IntVec(chain_nums, chain_nums + chain_azims));
}



///// \brief Construct an object with an object of TrackGenerator
//TrackLoadBalance::TrackLoadBalance(const TrackGenerator &track_generator):
//  _n_azims(track_generator.getNumAzims()),
//  _chain_azims(track_generator.getChainAzims()),
//  _lengths(track_generator.getLengths()),
//  _chain_nums(track_generator.getChainNums())
//{}


//----------------------------------------------------------------------
// Utility
//----------------------------------------------------------------------

/// \brief Print statistics for the loads
void TrackLoadBalance::printStatistics(const IntVec2D &load_map,
                                       int n_cpus, double p_cpu,
                                       int n_gpus, double p_gpu) {

  log::fdebug("load_map = %s", vecToString(load_map).c_str());

  // Compute the loads
  DoubleVec loads = computeLoads(load_map);

  double T_s = mathutils::sum(loads) / p_cpu;
  DoubleVec T_p = computeExpectedRunTime(loads, n_cpus, p_cpu, n_gpus, p_gpu);

  double p  = n_cpus + n_gpus * p_gpu / p_cpu;
  double speedup = T_s / mathutils::max(T_p);
  double efficiency = speedup / p * 100;

  auto max_load = mathutils::max(T_p);
  auto min_load = mathutils::min(T_p);
  log::fverbose_once("Total load of track mapping = %.2f", T_s);
  log::finfo(   "Relative error of max and min load per domain = %.2f%%",
                        (max_load - min_load)/min_load * 100);
  log::fverbose_once("Mean load of track mapping = %.2f", mathutils::mean(T_p));
  log::fverbose_once("Standard deviation of track mapping = %.2f", mathutils::stddev(T_p));
  log::finfo(   "Ideal efficiency of track mapping per domain = %.2f%%", efficiency);
}


/// \brief Compute the workloads for heterogeneous systems
DoubleVec TrackLoadBalance::computeLoads(const IntVec2D &load_map) {
  auto n_ranks = load_map.size();
  int n_azims = load_map[0].size();

  if (n_azims != _chain_azims)
    log::ferror("The number of chain azims (%d) is not the value used to "
                      "initialize the TrackLoadBalance (%d)", n_azims, _chain_azims);

  DoubleVec loads(n_ranks);
  for (size_t r = 0; r < n_ranks; ++r)
    for (int a = 0; a < _chain_azims; ++a)
      loads[r] += _lengths[a] * load_map[r][a];

  return loads;
}


/// \details For homogeneous systems, the expected time is equivalent to the
///          workload. Thus this method only affect heterogeneous systems.
DoubleVec TrackLoadBalance::
computeExpectedRunTime(const DoubleVec &loads,
                       int n_cpus, double p_cpu,
                       int n_gpus, double p_gpu) {
  DoubleVec T_p(loads);
  for (int r = 0; r < n_cpus + n_gpus; ++r) {
    // correct run times
    if (r < n_cpus)
      T_p[r] /= p_cpu;
    else
      T_p[r] /= p_gpu;
  }

  return T_p;
}


double TrackLoadBalance::
computeExpectedSpeedup(const IntVec2D &load_map,
                       int n_cpus, double p_cpu,
                       int n_gpus, double p_gpu) {
  // Compute the loads
  DoubleVec loads = computeLoads(load_map);

  double T_s = mathutils::sum(loads) / p_cpu;
  DoubleVec T_p = computeExpectedRunTime(loads, n_cpus, p_cpu, n_gpus, p_gpu);

  double speedup = T_s / mathutils::max(T_p);

  return speedup;
}


//----------------------------------------------------------------------
// Block distribution
//----------------------------------------------------------------------

/// \brief Take block distributions to parallel all of the chains
IntVec2D BlockDistribution::computeLoadMap(int n_ranks) {
  // Load map indexed by (r,a)
  IntVec2D load_map;

  for (int r = 0; r < n_ranks; ++r) {
    IntVec row;
    for (int a = 0; a < _chain_azims; ++a) {
      row.push_back(numLocal2DChainsByAzim(n_ranks, r, a));
    }
    load_map.push_back(row);
  }

  return load_map;
}


/// \brief Compute the number of 2D chains belonging to current rank
int
BlockDistribution::numLocal2DChainsByAzim(int n_ranks, int rank, int azim) {
  auto l_x = getFirst2DChainGuid(azim);
  auto u_x = getFirst2DChainGuid(azim + 1) - 1;
  auto sub1 = l_x - rank;
  auto sub2 = u_x - rank;
  int n_chains = 0;
  int a = 0, b = 0;

  if (sub2 >= 0) {
    a = (sub1 < 0)? 0 : std::ceil(sub1 / (double)n_ranks);
    b = sub2 / n_ranks;
    n_chains = b - a + 1;
  }

  return n_chains;
}


/// \brief Return the global index of the first chain of angle azim
int BlockDistribution::getFirst2DChainGuid(int azim) {
  int sum_chains = 0;
  for (int a = 0; a < azim; ++a)
    sum_chains += _chain_nums[a];

  return sum_chains;
}


//----------------------------------------------------------------------
// Cyclic track distribution
//----------------------------------------------------------------------

/// \brief Greedy strategy for cyclic track distribution
IntVec2D CyclicTrackDistribution::computeLoadMap(int n_ranks) {
  // Load map indexed by (r,a)
  IntVec2D load_map(n_ranks);
  DoubleVec deltas(n_ranks);

  // normalization
  DoubleVec norm_lengths;
  double l_max = mathutils::max(_lengths);
  for (auto e : _lengths)
    norm_lengths.push_back(e / l_max);

  for (int a = 0; a < _chain_azims; ++a) {
    // Subproblem solved by dynamic programming
    auto column = minimizeCosts(n_ranks,
                                _chain_nums[a],
                                norm_lengths[a],
                                deltas);
    for (int r = 0; r < n_ranks; ++r) {
      // Update the load map
      load_map[r].push_back(column[r]);
      // Update deltas
      auto avg_w = norm_lengths[a] * _chain_nums[a] / n_ranks;
      deltas[r] += norm_lengths[a] * column[r] - avg_w;
    }
  }

  return load_map;
}


/// \brief Exploit dynamic programming to solve subproblems
IntVec CyclicTrackDistribution::minimizeCosts(int n_ranks,
                                              int n_chains,
                                              double length,
                                              const DoubleVec &deltas) {
  // compute cost for single stage
  auto cost = [](int x, double a, double b) {
    return (a * x - b) * (a * x - b);
  };

  // compute the average workload
  auto avg_w = length * n_chains / n_ranks;

  // the coefficient in g(x) = (ax-b)^2
  DoubleVec b;
  for (auto e : deltas)
    b.push_back(avg_w - e);

  // optimal values of current and previous subproblems
  DoubleVec h(n_chains + 1);
  DoubleVec h_prev(n_chains + 1);
  // optimal choices
  IntVec2D opt_choices(n_ranks);
  for (auto &v : opt_choices)
    v.resize(n_chains + 1);

  // stage 0, allocate loads for rank 0
  for (int k = 0; k < n_chains + 1; ++k) {
    h_prev[k] = cost(k, length, b[0]);
    opt_choices[0][k] = k;
  }

  // Loop over stage 1 to n_ranks-1
  for (int r = 1; r < n_ranks; ++r) {
    // k is the total number of resouces defined for the subproblem
    // x is the number assigned to g(x)
    // k-x is the remaining number of resources for sub-subproblems
    // if r == n_ranks - 1, only k == n_chains is considered
    int k = r / (n_ranks - 1) * n_chains;
    for (; k <= n_chains; ++k) {
      h[k] = std::numeric_limits<double>::max();

      // search an optimal value
      for (int x = 0; x < k + 1; ++x) {
        auto c = h_prev[k - x] + cost(x, length, b[r]);
        if (c < h[k]) {
          h[k] = c;
          opt_choices[r][k] = x;
        }
      }
    }
    // save the current optimal values
    h_prev = h;
  }

  // Trace back for the optimal solution
  IntVec solution(n_ranks);
  int k = n_chains;
  for (int r = n_ranks - 1; r >= 0; --r) {
    solution[r] = opt_choices[r][k];
    k = k - solution[r];
  }

  return solution;
}


//----------------------------------------------------------------------
// Distribute angles over processes
//----------------------------------------------------------------------

/// \brief Constructor
AngleDistribution::AngleDistribution(int chain_azims,
                                     DoubleVec lengths,
                                     IntVec chain_nums):
  CyclicTrackDistribution(chain_azims, lengths, chain_nums),
  // Copy the vectors
  _lengths_copy(lengths),
  _chain_nums_copy(chain_nums)
{ }


/// \brief Constructor
AngleDistribution::AngleDistribution(int chain_azims,
                                     const double*lengths,
                                     const long *chain_nums):
  AngleDistribution(chain_azims,
                    DoubleVec(lengths, lengths + chain_azims),
                    IntVec(lengths, lengths + chain_azims))
{}


/// \brief Greedy strategy for angle distribution
IntVec2D AngleDistribution::computeLoadMap(int n_ranks) {

  // Accumulate values
  for (int a = 0; a < _chain_azims; ++a) {
    _lengths[a] *= _chain_nums[a];
    _chain_nums[a] = 1;
  }

  // Load map indexed by (r,a)
  IntVec2D load_map = CyclicTrackDistribution::computeLoadMap(n_ranks);

  // Restore vectors and the number of chains
  _lengths = _lengths_copy;
  _chain_nums = _chain_nums_copy;
  for (int r = 0; r < n_ranks; ++r)
    for (int a = 0; a < _chain_azims; ++a)
      load_map[r][a] *= _chain_nums[a];

  return load_map;
}


} // namespace antmoc
