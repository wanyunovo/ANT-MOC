#include <array>
#include <iostream>

#include "app_headers.h"

using namespace antmoc;

int main(int argc, char* argv[]) {

#ifdef ENABLE_MPI_
  int provided;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_SERIALIZED, &provided);
  /* Initialize MPI static variables */
  mpigod::initializeMPIGod();
  /* Initialize MPI logger */
  log_set_ranks();
#endif

  /* Define simulation parameters */
  auto conf = Factory::getConfInput<ConfigInputCLI>(argc, argv);

  /* Set log level */
  conf->readLogLevel();
  conf->initializeLogLevel();

  /* Set default zones */
  conf->setSegmentationType("OTF Tracks");
  conf->setZones("-12.5 12.5");

  /* Print help information if needed */
  conf->showHelpOnNeed();

  /* Default values */
  conf->setGlobalRefines("10 10 10");
  conf->setNumModules("3 3 1");

  /* Parse conf parameters */
  conf->parseArguments();
  conf->runChecks();

  int num_threads = conf->getNumThreads();

  // Reference : Shaner-2016
  double azim_spacing = conf->getAzimSpacing();   // 0.05
  int num_azim = conf->getNumAzim();  // 32
  double polar_spacing = conf->getPolarSpacing(); // 0.05
  int num_polar = conf->getNumPolar();  // 2 to 14
  double tolerance = conf->getTolerance(); // 1E-5
  int max_iters = conf->getMaxIterations(); // 1000
  int refines = conf->getGlobalRefines()[0]; // 20
  auto nmods = conf->getNumModules();

  /* Set logging information */
  log_printf(TITLE, "Takeda Benchmark Problem - Unrodded");

  /* Report input arguments */
  conf->printArgumentsReport();

  /* Define material properties */
  log_printf(NORMAL, "Defining material properties...");

  const size_t num_groups = 2;
  std::map<std::string, std::array<double, num_groups> > sigma_a;
  std::map<std::string, std::array<double, num_groups> > sigma_t;
  std::map<std::string, std::array<double, num_groups*num_groups> > sigma_s;
  std::map<std::string, std::array<double, num_groups> > sigma_f;
  std::map<std::string, std::array<double, num_groups> > nu_sigma_f;
  std::map<std::string, std::array<double, num_groups> > chi;

  /* Define core cross-sections */
  sigma_a["Core"] = std::array<double, num_groups> {0.00852709, 0.158196};
  sigma_t["Core"] = std::array<double, num_groups> {0.223775, 1.03864};
  sigma_s["Core"] = std::array<double, num_groups*num_groups>
      {0.192423, 0.0228253,
        0.00, 0.880439};
  sigma_f["Core"] = std::array<double, num_groups> {0.0004, 0.1};
  nu_sigma_f["Core"] = std::array<double, num_groups> {0.00909319, 0.290183};
  chi["Core"] = std::array<double, num_groups> {1.0, 0.0};

  /* Define reflector cross-sections */
  sigma_a["Reflector"] = std::array<double, num_groups> {0.000416392,
                                                         0.0202999};
  sigma_t["Reflector"] = std::array<double, num_groups> {0.250367, 1.64482};
  sigma_s["Reflector"] = std::array<double, num_groups*num_groups>
      {0.193446, 0.0565042,
       0.00, 1.62452};
  sigma_f["Reflector"] = std::array<double, num_groups> {0.0, 0.0};
  nu_sigma_f["Reflector"] = std::array<double, num_groups> {0.0, 0.0};
  chi["Reflector"] = std::array<double, num_groups> {1.0, 0.0};

  /* Define void cross-sections */
  sigma_a["Void"] = std::array<double, num_groups> {0.0000465132, 0.00132890};
  sigma_t["Void"] = std::array<double, num_groups> {0.0128407, 0.0120676};
  sigma_s["Void"] = std::array<double, num_groups*num_groups>
      {0.01277, 0.0000240997,
       0.00, 0.0107387};
  sigma_f["Void"] = std::array<double, num_groups> {0.0, 0.0};
  nu_sigma_f["Void"] = std::array<double, num_groups> {0.0, 0.0};
  chi["Void"] = std::array<double, num_groups> {1.0, 0.0};

  /* Create materials */
  log_printf(NORMAL, "Creating materials...");
  std::map<std::string, Material*> materials;
  std::map<std::string, std::array<double, num_groups> >::iterator it;
  int id_num = 0;
  for (it = sigma_t.begin(); it != sigma_t.end(); it++) {

    std::string name = it->first;
    materials[name] = new Material(id_num, name.c_str());
    materials[name]->setNumEnergyGroups(num_groups);
    id_num++;

    materials[name]->setSigmaT(sigma_t[name].data(), num_groups);
    materials[name]->setSigmaS(sigma_s[name].data(), num_groups*num_groups);
    materials[name]->setSigmaF(sigma_f[name].data(), num_groups);
    materials[name]->setNuSigmaF(nu_sigma_f[name].data(), num_groups);
    materials[name]->setChi(chi[name].data(), num_groups);
  }

  /* Create boundaries of the geometry */
  log_printf(NORMAL, "Creating surfaces...");

  XPlane xmin(-12.5);
  XPlane xmax( 12.5);
  YPlane ymin(-12.5);
  YPlane ymax( 12.5);
  ZPlane zmin(-12.5);
  ZPlane zmax( 12.5);

  xmin.setBoundaryType(REFLECTIVE);
  ymin.setBoundaryType(REFLECTIVE);
  zmin.setBoundaryType(REFLECTIVE);
  xmax.setBoundaryType(VACUUM);
  ymax.setBoundaryType(VACUUM);
  zmax.setBoundaryType(VACUUM);

  /* Create cells in the geometry */
  log_printf(NORMAL, "Creating cells...");

  Cell* core_cell = new Cell();
  core_cell->setFill(materials["Core"]);

  Cell* reflector_cell = new Cell();
  reflector_cell->setFill(materials["Reflector"]);

  Cell* void_cell = new Cell();
  void_cell->setFill(materials["Void"]);

  /* Root Cell* */
  Cell* root_cell = new Cell(1, "root");
  root_cell->addSurface(+1, &xmin);
  root_cell->addSurface(-1, &xmax);
  root_cell->addSurface(+1, &ymin);
  root_cell->addSurface(-1, &ymax);
  root_cell->addSurface(+1, &zmin);
  root_cell->addSurface(-1, &zmax);

  /* Create Universes */
  log_printf(NORMAL, "Creating universes...");

  Universe* core = new Universe();
  Universe* reflector = new Universe();
  Universe* void_u = new Universe();
  Universe* root_universe = new Universe();

  core->addCell(core_cell);
  reflector->addCell(reflector_cell);
  void_u->addCell(void_cell);
  root_universe->addCell(root_cell);

  /* Setup Takeda core */
  log_printf(NORMAL, "Creating Takeda core...");

  RecLattice* lattice = new RecLattice();
  lattice->setWidth(5.0/refines, 5.0/refines, 5.0/refines);

  Universe* mold[] = {reflector, reflector, reflector, reflector, reflector,
                      reflector, reflector, reflector, reflector, reflector,
                      reflector, reflector, reflector, reflector, reflector,
                      reflector, reflector, reflector, reflector, reflector,
                      reflector, reflector, reflector, void_u,    reflector,

                      reflector, reflector, reflector, reflector, reflector,
                      reflector, reflector, reflector, reflector, reflector,
                      reflector, reflector, reflector, reflector, reflector,
                      reflector, reflector, reflector, reflector, reflector,
                      reflector, reflector, reflector, void_u,    reflector,

                      reflector, reflector, reflector, reflector, reflector,
                      reflector, reflector, reflector, reflector, reflector,
                      core,      core,      core,      reflector, reflector,
                      core,      core,      core,      reflector, reflector,
                      core,      core,      core,      void_u,    reflector,

                      reflector, reflector, reflector, reflector, reflector,
                      reflector, reflector, reflector, reflector, reflector,
                      core,      core,      core,      reflector, reflector,
                      core,      core,      core,      reflector, reflector,
                      core,      core,      core,      void_u,    reflector,

                      reflector, reflector, reflector, reflector, reflector,
                      reflector, reflector, reflector, reflector, reflector,
                      core,      core,      core,      reflector, reflector,
                      core,      core,      core,      reflector, reflector,
                      core,      core,      core,      void_u,    reflector};

  /* Refine lattice */
  Universe* refined_mold[5*5*5*refines*refines*refines];
  for (int x=0; x < refines; x++) {
    for (int i=0; i < 5; i++) {
      for (int y=0; y < refines; y++) {
        for (int j=0; j < 5; j++) {
          for (int z=0; z < refines; z++) {
            for (int k=0; k < 5; k++) {
              int index = 5*5*refines*refines*refines*k + 5*5*refines*refines*z
                + 5*refines*refines*j + 5*refines*y + i*refines + x;
              refined_mold[index] = mold[5*5*k + 5*j + i];
            }
          }
        }
      }
    }
  }
  lattice->setUniverses(5*refines, 5*refines, 5*refines, refined_mold);
  root_cell->setFill(lattice);

#ifdef ENABLECMFD
  /* Create CMFD mesh */
  log_printf(NORMAL, "Creating Cmfd mesh...");
  Cmfd* cmfd = new Cmfd();
  cmfd->setSORRelaxationFactor(1.5);
  cmfd->setLatticeStructure(5, 5, 5);
  cmfd->setKNearest(1);
#endif

  /* Create the geometry */
  log_printf(NORMAL, "Creating geometry...");

  Geometry geometry;
  geometry.setRootUniverse(root_universe);
#ifdef ENABLECMFD
  geometry.setCmfd(cmfd);
#endif
  geometry.setNumDomainModules(nmods[0], nmods[1], nmods[2]);
  geometry.initializeFlatSourceRegions();

  log_printf(NORMAL, "Number of energy groups = %d",
                     geometry.getNumEnergyGroups());
/* Create the track generator */
  log_printf(NORMAL, "Initializing the track generator...");

  auto quad = Factory::getQuadrature(conf);
  auto tg = Factory::getTrackGenerator(&geometry, conf);

  tg->setQuadrature(quad);
  tg->generateTracks();

  /* Run simulation */
  SolverPtr solver = Factory::getSolver(tg, conf);

  solver->computeEigenvalue(max_iters);
  solver->printTimerReport();

  /* Print reaction rates */
#ifdef ENABLE_MPI_
  if (mpigod::isMPIRoot()) {
#endif
    log_printf(TITLE, "Reaction Rates");
    Mesh mesh(solver);
    mesh.createLattice(5, 5, 5);
    for (int i = -1; i < num_groups; ++i) {
      mesh.printMeshDataToXML("rx", "Fission reaction rates", TallyType::Fission_RX, false, i);
      mesh.printMeshDataToXML("rx", "Scalar fluxes", TallyType::Flux_RX, false, i);
    }
#ifdef ENABLE_MPI_
  }
#endif

  log_printf(TITLE, "Finished");

#ifdef ENABLE_MPI_
  MPI_Finalize();
#endif
  return 0;
}
