#include "antmoc/FSRDataHandler.h"
#include "antmoc/Geometry.h"
#include "antmoc/Material.h"
#include "antmoc/Point.h"
#include "antmoc/Solver.h"
#include "antmoc/string_utils.h"
#include "antmoc/TrackGenerator3D.h"

#include <sstream>

namespace antmoc
{


/// \brief Constructor of FSRDataHandler
FSRDataHandler::FSRDataHandler(SolverPtr solver)
  : _solver(solver) { }


/// \brief Returns the number of energy groups
int FSRDataHandler::getNumEnergyGroups() {
  return _solver->getGeometry()->getNumEnergyGroups();
}


/// \brief Returns the number of local FSRs
long FSRDataHandler::getNumFSRs() {
  return _solver->getGeometry()->getNumFSRs();
}


/// \brief Extract the characteristic point of FSRs
std::vector<std::vector<FP_PRECISION>>
FSRDataHandler::getFSRPoints() {
  auto num_fsrs = getNumFSRs();

  // Loop over all flat source regions
  std::vector<std::vector<FP_PRECISION>> fsr_points(3);
  for (auto &v : fsr_points)
    v.reserve(num_fsrs);

  for (long r = 0; r < num_fsrs; ++r) {
    Point *pt = getFSRPoint(r); // characteristic points
    fsr_points[0].push_back(pt->getX());
    fsr_points[1].push_back(pt->getY());
    fsr_points[2].push_back(pt->getZ());
  }

  return fsr_points;
}


/// \brief Returns the characteristic point of an FSR
Point* FSRDataHandler::getFSRPoint(long fsr_id) {
  return _solver->getGeometry()->getFSRPoint(fsr_id);
}


/// \brief Extract the centroids of FSRs
std::vector<std::vector<FP_PRECISION>>
FSRDataHandler::getFSRCentroids() {
  auto num_fsrs = getNumFSRs();

  // Loop over all flat source regions
  std::vector<std::vector<FP_PRECISION>> fsr_centroids(3);
  for (auto &v : fsr_centroids)
    v.reserve(num_fsrs);

  for (long r = 0; r < num_fsrs; ++r) {
    Point *pt = getFSRCentroid(r); // centroids
    fsr_centroids[0].push_back(pt->getX());
    fsr_centroids[1].push_back(pt->getY());
    fsr_centroids[2].push_back(pt->getZ());
  }

  return fsr_centroids;
}


/// \brief Returns the centroid of an FSR
Point* FSRDataHandler::getFSRCentroid(long fsr_id) {
  return _solver->getGeometry()->getFSRCentroid(fsr_id);
}


/// \brief Get the array of volumes
FP_PRECISION* FSRDataHandler::getFSRVolumes() {
  return _solver->getTrackGenerator()->getFSRVolumesBuffer();
}


/// \brief Returns the volume of an FSR
FP_PRECISION FSRDataHandler::getFSRVolume(long fsr_id) {
  return _solver->getTrackGenerator()->getFSRVolume(fsr_id);
}


/// \brief Get the array of FSR fluxes
FP_PRECISION* FSRDataHandler::getFSRFluxes() {
  return _solver->getFluxesArray();
}


/// \brief Returns the flux of an FSR
/// \param group index of energy group starting from 1
FP_PRECISION FSRDataHandler::getFSRFlux(long fsr_id, int group) {
  return _solver->getFlux(fsr_id, group);
}


/// \brief Returns the cross-section of an FSR
FP_PRECISION FSRDataHandler::getFSRXS(long fsr_id, int group, TallyType type) {

  // Get the associated material
  Material* mat = _solver->getGeometry()->findFSRMaterial(fsr_id);

  FP_PRECISION xs = 1.0;

  switch (type) {
    case TallyType::Fission_RX:
    case TallyType::Fission_XS:
      xs = mat->getSigmaFByGroup(group);
      break;

    case TallyType::NuFission_RX:
    case TallyType::NuFission_XS:
      xs = mat->getNuSigmaFByGroup(group);
      break;

    case TallyType::Total_RX:
    case TallyType::Total_XS:
      xs = mat->getSigmaTByGroup(group);
      break;

    case TallyType::Absorption_RX:
    case TallyType::Absorption_XS:
      xs = mat->getSigmaAByGroup(group);
      break;

    case TallyType::Volume:
    case TallyType::Scalar_Flux:
      // Defaults to 1.0
      break;

    default:
      log::error("Unrecognized tally type of reaction rates"
                 "/cross-sections, please make sure that ALL",
                 "/NONE have been expanded");
  }

  return xs;
}


/// \details Constructs a string for integer range and then parse it to
///          generate a set of integers.
std::set<int> FSRDataHandler::getEnergyGroupSet() {
  std::stringstream ss;
  ss << "1:" << getNumEnergyGroups();
  return stringutils::toIntegerSet(ss.str());
}


} // namespace antmoc
