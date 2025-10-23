/// \file include/Factory.h
/// \brief Object factories
/// \date Aug 6, 2019
/// \author An Wang, USTB (wangan@xs.ustb.edu.cn)

#ifndef FACTORY_H_
#define FACTORY_H_

#include <memory>
#include <string>

#include "antmoc/ConfigInputFile.h"
#include "antmoc/CPULSSolver.h"
#include "antmoc/enum_types.h"
#include "antmoc/Geometry.h"
#include "antmoc/GeoInputXml.h"
#include "antmoc/MaterialHandlerHDF5.h"
#include "antmoc/Quadrature.h"
#include "antmoc/TrackGenerator3D.h"

using std::shared_ptr;
using std::make_shared;
using std::static_pointer_cast;
using std::dynamic_pointer_cast;

namespace antmoc
{

///---------------------------------------------------------------------
/// \class Factory
/// \brief A factory for objects
/// \details This class contains class templates for object creation.
///---------------------------------------------------------------------
class Factory {

public:

  // Alias for objects
  using ConfInputPtr        = shared_ptr<ConfigInput>;
  using MaterialHandlerPtr  = shared_ptr<MaterialHandler>;
  using GeoInputPtr         = shared_ptr<GeoInput>;

  using GeometryPtr         = shared_ptr<Geometry>;
  using SolverPtr           = shared_ptr<Solver>;
  using QuadraturePtr       = shared_ptr<Quadrature>;
  using TrackGeneratorPtr   = shared_ptr<TrackGenerator>;
  using TrackGenerator3DPtr = shared_ptr<TrackGenerator3D>;

  // Factory methods
  template <typename T>
  static ConfInputPtr getConfInput(int, char **);

  template <typename T>
  static ConfInputPtr getConfInput(const StringVec &);

  template <typename T>
  static MaterialHandlerPtr getMaterialHandler(ConfInputPtr);

  template <typename T>
  static GeoInputPtr getGeoInput(Geometry *geometry = nullptr,
                                 MaterialHandlerPtr = nullptr,
                                 ConfInputPtr conf = nullptr);

  static SolverPtr getSolver(TrackGeneratorPtr, ConfInputPtr);
  static QuadraturePtr getQuadrature(ConfInputPtr);
  static TrackGeneratorPtr getTrackGenerator(Geometry *, ConfInputPtr);

};


/// \brief Create an object of ConfigInput
/// \details The child class is specified by template arguments.
/// \param argc number of CLI arguments
/// \param argv values of CLI arguments
/// \return a smart pointer to the object
template <typename T>
inline Factory::ConfInputPtr
Factory::getConfInput(int argc, char **argv) {

  return static_pointer_cast<ConfigInput>
          ( make_shared<T>(argc, argv) );

}


/// \brief Create an object of ConfigInput
/// \details The child class is specified by template arguments.
/// \param argv values of CLI arguments
/// \return a smart pointer to the object
template <typename T>
inline Factory::ConfInputPtr
Factory::getConfInput(const StringVec &argv) {

  return static_pointer_cast<ConfigInput>
          ( make_shared<T>(argv) );

}


/// \brief Create an object of MaterialHandler
/// \details The child class is specified by template arguments
/// \param conf runtime settings
/// \return a smart pointer to the object
template <typename T>
inline Factory::MaterialHandlerPtr
Factory::getMaterialHandler(ConfInputPtr conf) {

  auto path = conf->getMatInputPath();
  auto mode = HDF5Mode::ReadOnly;
  auto layout = conf->getXSFileLayout();

  return static_pointer_cast<MaterialHandler>
          ( make_shared<T>(path, mode, layout) );

}


/// \brief Create an object of GeoInput
/// \details The child class is specified by template arguments
/// \param geometry a pointer to an existed Geometry object
/// \param path geometry file to be read
/// \param matinput a pointer to an MaterialHandler object
/// \return a smart pointer to the object
template <typename T>
inline Factory::GeoInputPtr
Factory::getGeoInput(Geometry *geometry,
                     MaterialHandlerPtr matinput,
                     ConfInputPtr conf) {

  // Initialize an object of child class T
  auto geo_input = make_shared<T>(geometry, matinput);

  // Set attributes
  if (conf) {
    geo_input->setGlobalRefines(conf->getLatticeRefines());
    geo_input->setGlobalSectors(conf->getCellSectors());
    geo_input->setGlobalRings(conf->getCellRings());
  }

  return static_pointer_cast<GeoInput>(geo_input);

}


/// \brief Create an object of Solver
/// \details The child class is specified by function arguments
/// \param track_generator tracks for solver initialization
/// \param conf runtime settings
/// \return a smart pointer to the object
inline Factory::SolverPtr
Factory::getSolver(TrackGeneratorPtr track_generator,
                   ConfInputPtr conf) {

  SolverPtr solver;

  auto type = conf->getSolverType();
  switch (type) {
    case solverType::CPU_SOLVER :
      solver = make_shared<CPUSolver>(track_generator);
      break;
    case solverType::CPU_LS_SOLVER :
      solver = make_shared<CPULSSolver>(track_generator);
      break;
    default:
      solver = nullptr;
  }

  // Set the number of threads
  solver->setNumThreads(conf->getNumThreads());

  // Set tolerance
  solver->setConvergenceThreshold(conf->getTolerance());

  // Specify keff formula
  if (conf->getKeffFromNeutronBalance()) {
    solver->setKeffFromNeutronBalance();
  }

  // Set stabilization
  if (conf->getStabilizationType() != +stabilizationType::NONE) {
    solver->stabilizeTransport(conf->getStabilizationFactor(),
                               conf->getStabilizationType());
  }

  // MGXS checks and correction
  solver->setCheckXSLogLevel(conf->getCheckXSLogLevel());
  solver->setCheckXS(conf->doesCheckXS());
  solver->setCorrectXS(conf->doesCorrectXS());

  return static_pointer_cast<Solver>(solver);

}


/// \brief Create an object of Quadrature
/// \details The child class is specified by function arguments
/// \param conf runtime settings
/// \return a smart pointer to the object
inline Factory::QuadraturePtr
Factory::getQuadrature(ConfInputPtr conf) {

  QuadraturePtr quad;

  auto type = conf->getQuadratureType();
  switch (type) {
    case quadratureType::TABUCHI_YAMAMOTO :
      quad = make_shared<TYPolarQuad>();
      break;
    case quadratureType::LEONARD :
      quad = make_shared<LeonardPolarQuad>();
      break;
    case quadratureType::GAUSS_LEGENDRE :
      quad = make_shared<GLPolarQuad>();
      break;
    case quadratureType::EQUAL_WEIGHT :
      quad = make_shared<EqualWeightPolarQuad>();
      break;
    case quadratureType::EQUAL_ANGLE :
      quad = make_shared<EqualAnglePolarQuad>();
      break;
    default:
      quad = nullptr;
  }

  // Initialize quadrature attributes
  quad->setNumAzimAngles(conf->getNumAzim());
  quad->setNumPolarAngles(conf->getNumPolar());

  return static_pointer_cast<Quadrature>(quad);

}


/// \brief Create an object of TrackGenerator
/// \details The child class is specified by function arguments
/// \param geo pointer to a geometry object
/// \param conf runtime settings
/// \return a smart pointer to the object
inline Factory::TrackGeneratorPtr
Factory::getTrackGenerator(Geometry *geo, ConfInputPtr conf) {

  TrackGeneratorPtr track_gen;
  TrackGenerator3DPtr track_gen_3d;

  auto type = conf->getSegmentationType();
  switch (type) {
    case segmentationType::EXPLICIT_2D :
      track_gen = make_shared<TrackGenerator>
                    ( geo,
                      conf->getNumAzim(),
                      conf->getAzimSpacing() );
      break;

    case segmentationType::EXPLICIT_3D :
    case segmentationType::OTF_TRACKS :
    case segmentationType::OTF_STACKS :
      track_gen_3d = make_shared<TrackGenerator3D>
                      ( geo,
                        conf->getNumAzim(),
                        conf->getNumPolar(),
                        conf->getAzimSpacing(),
                        conf->getPolarSpacing() );

      // Set tracing method for 3-D cases
      track_gen_3d->setSegmentFormation(conf->getSegmentationType());
      // Set axial zones for 3-D tracing
      // If it is empty, axial zones will be auto-generated
      if (!conf->getZones().empty())
        track_gen_3d->setSegmentationZones(conf->getZones());

      if (conf->isGlobalZMesh()) {
        track_gen_3d->useGlobalZMesh();
      }

      track_gen = track_gen_3d;
      break;

    default:
      track_gen = nullptr;
  }

#ifdef ENABLE_MPI_
  // Set the track mapping algorithm
  track_gen->setTrackMappingType(conf->getTrackMappingType());
#endif

  // Set the number of threads used for ray tracing
  track_gen->setNumThreads(conf->getNumThreads());

  return static_pointer_cast<TrackGenerator>(track_gen);

}

} // namespace antmoc

#endif  // FACTORY_H_
