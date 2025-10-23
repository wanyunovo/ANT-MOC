/// \file include/TrackLoadBalance.h
/// \brief Cyclic track mapping algorithms
/// \author An Wang, USTB (wangan.cs@gmail.com)

#ifndef TRACK_LOAD_BALANCE_H_
#define TRACK_LOAD_BALANCE_H_

#include <array>
#include <memory>

#include "antmoc/container_utils.h"
#include "antmoc/enum_types.h"
#include "antmoc/TrackGenerator.h"


namespace antmoc {


// Forward declarations
class TrackLoadBalance;

using TrackLoadBalancePtr = std::shared_ptr<TrackLoadBalance>;


///---------------------------------------------------------------------
/// \class TrackLoadBalance
/// \brief A base class for cyclic track mapping algorithms
///---------------------------------------------------------------------
class TrackLoadBalance {

protected:

  int _chain_azims;       ///< Number of angle associated with unique chains
  DoubleVec _lengths;     ///< Lengths of chains
  IntVec    _chain_nums;  ///< Number of chains with each angle

public:

  /// \brief Default constructor
  TrackLoadBalance() = default;

  /// \brief Default destructor
  virtual ~TrackLoadBalance() = default;

  /// \brief Constructor
  TrackLoadBalance(int chain_azims,
                   DoubleVec chain_lengths,
                   IntVec chain_nums);

  /// \brief Constructor
  TrackLoadBalance(int chain_azims,
                   const double *chain_lengths,
                   const long *chain_nums);

  // \brief Create a TrackLoadBalance object
  static TrackLoadBalancePtr
  getTrackLoadBalancer(trackMappingType type,
                       int chain_azims,
                       DoubleVec chain_lengths,
                       IntVec chain_nums);

  // \brief Create a TrackLoadBalance object
  static TrackLoadBalancePtr
  getTrackLoadBalancer(trackMappingType type,
                       int chain_azims,
                       const double *chain_lengths,
                       const long *chain_nums);


  ///// \brief Constructor
  //TrackLoadBalance(const TrackGenerator *track_generator);

  int getChainAzims() const    { return _chain_azims; }
  DoubleVec getLengths() const { return _lengths; }
  IntVec getChainNums() const  { return _chain_nums; }

public:

  //--------------------------------------
  // Interfaces
  //--------------------------------------
  /// \brief Compute the load map for decomposition.
  /// \param n_ranks Number of ranks.
  /// \return A 2D array with n_ranks rows.
  virtual IntVec2D computeLoadMap(int n_ranks) = 0;


  //--------------------------------------
  // Utility
  //--------------------------------------
  void printStatistics(const IntVec2D &,
                       int n_cpus = 1, double p_cpu = 1.,
                       int n_gpus = 0, double p_gpu = 1.);
  DoubleVec computeLoads(const IntVec2D &);

  /// \brief Compute the expected execution time by workloads
  DoubleVec computeExpectedRunTime(const DoubleVec &loads,
                                   int n_cpus, double p_cpu,
                                   int n_gpus, double p_gpu);

  /// \brief Compute the expected speedup
  double computeExpectedSpeedup(const IntVec2D &load_map,
                                int n_cpus = 1, double p_cpu = 1.,
                                int n_gpus = 0, double p_gpu = 1.);

};


///---------------------------------------------------------------------
/// \brief Load balancing class for block distribution
///---------------------------------------------------------------------
class BlockDistribution: public TrackLoadBalance {

public:
  using TrackLoadBalance::TrackLoadBalance;

  IntVec2D computeLoadMap(int n_ranks);

  int numLocal2DChainsByAzim(int, int, int);
  int getFirst2DChainGuid(int);

};


///---------------------------------------------------------------------
/// \brief Load balancing class for cyclic track decomposition
///---------------------------------------------------------------------
class CyclicTrackDistribution: public TrackLoadBalance {

public:
  using TrackLoadBalance::TrackLoadBalance;

  IntVec2D computeLoadMap(int n_ranks);

  IntVec minimizeCosts(int, int, double, const DoubleVec &);

};


///---------------------------------------------------------------------
/// \brief Load balancing class for angular decomposition
///---------------------------------------------------------------------
class AngleDistribution: public CyclicTrackDistribution {

private:

  DoubleVec _lengths_copy;
  IntVec    _chain_nums_copy;

public:
  /// \brief Constructor
  AngleDistribution(int chain_azims,
                    DoubleVec lengths,
                    IntVec chain_nums);

  /// \brief Constructor
  AngleDistribution(int chain_azims,
                    const double *chain_lengths,
                    const long *chain_nums);

  IntVec2D computeLoadMap(int n_ranks);
};


} // namespace

#endif  // TRACK_LOAD_BALANCE_H_
