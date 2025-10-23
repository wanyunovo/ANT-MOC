/// \file FSRDataHandler.h
/// \brief A class for fetching FSR data from a solver
/// \date March 6, 2020
/// \author An Wang, USTB (wangan@xs.ustb.edu.cn)

#ifndef FSRDATA_HANDLER_H_
#define FSRDATA_HANDLER_H_

#include <memory>
#include <set>
#include <vector>

#include "antmoc/tally_utils.h"


namespace antmoc
{

class Solver;
class Point;

///---------------------------------------------------------------------
/// \class FSRDataHandler
/// \brief Methods for fetching FSR data from a solver
///---------------------------------------------------------------------
class FSRDataHandler {

protected:

  using SolverPtr = std::shared_ptr<Solver>;

  ///< The solver from which data is extracted
  SolverPtr _solver;

public:

  FSRDataHandler(SolverPtr solver);
  ~FSRDataHandler() = default;

  int getNumEnergyGroups();
  long getNumFSRs();
  std::vector<std::vector<FP_PRECISION>> getFSRPoints();
  Point *getFSRPoint(long fsr_id);
  std::vector<std::vector<FP_PRECISION>> getFSRCentroids();
  Point *getFSRCentroid(long fsr_id);
  FP_PRECISION* getFSRVolumes();
  FP_PRECISION getFSRVolume(long fsr_id);
  FP_PRECISION* getFSRFluxes();
  FP_PRECISION getFSRFlux(long fsr_id, int group);
  FP_PRECISION getFSRXS(long fsr_id, int group, TallyType type);

  /// \brief Returns a set of energy group ids
  std::set<int> getEnergyGroupSet();
};


} // namespace antmoc

#endif  // FSRDATA_HANDLER_H_
