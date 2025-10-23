/**
 * @file test_initializePeriodicTracks.cpp
 * @brief Test periodic track initialization
 * @date May 17, 2019
 * @author An Wang, USTB (wangan@xs.ustb.edu.cn)
 */

#ifdef USTB_

#include "testing/test_utils.h"
#include "testing/MockGeometry.h"
#include "antmoc/Geometry.h"
#include "antmoc/Point.h"
#include "antmoc/Quadrature.h"
#include "antmoc/Track3D.h"
#include "antmoc/TrackGenerator3D.h"

using ::testing::AtLeast;
using ::testing::Return;
using namespace antmoc;

namespace {
class TrackGenerator3DTest : public TrackGenerator3D {
 public:
  TrackGenerator3DTest(Geometry* geometry, int num_azim, int num_polar,
                   double azim_spacing, double z_spacing)
                   : TrackGenerator3D(geometry, num_azim, num_polar,
                     azim_spacing, z_spacing) {};

  void UTinitialize2DTracks() {
    TrackGenerator::initializeTracks();
  };

  void UTinitialize2DChains() {
    TrackGenerator::initializeTracks();
  };

  void UTinitialize3DTracks() {
    TrackGenerator::initializeTracks();
    initialize3DTracks();
  };

  void UTinitialize3DChains() {
    TrackGenerator::initializeTracks();
    initialize3DTracks();
    set3DIndexArraysAndTracks();
  };

  friend class test_initializePeriodicTracks;
  FRIEND_TEST(test_initializePeriodicTracks, project3DTrackStartingPoint);
  FRIEND_TEST(test_initializePeriodicTracks, set3DIndexArraysAndTracks);
};

// Test fixture
class test_initializePeriodicTracks : public testing::Test {
 protected:
  virtual void SetUp() {

#   ifdef ENABLE_MPI_
    mpi::setDomainDecomposition(mpi::getMPIComm(),
                                   1,1,1,2);
#   endif
    num_threads = omp_get_max_threads();

    /* Test input */
    #include "initializeTracks_input.inc"

    std::vector<FP_PRECISION> seg_zones {zmin, -5, 5, zmax};

    quad = std::make_shared<GLPolarQuad>();
    quad->setNumAzimAngles(num_azim);
    quad->setNumPolarAngles(num_polar);
    /* Initialize the quadrature set */
    quad->initialize();

    /* Initialize the mock geometry */
    /* Number of modules */
    EXPECT_CALL(mockgeo, getNumXModules())
        .WillRepeatedly(Return(1));
    EXPECT_CALL(mockgeo, getNumYModules())
        .WillRepeatedly(Return(1));
    EXPECT_CALL(mockgeo, getNumZModules())
        .WillRepeatedly(Return(1));
    /* Width */
    EXPECT_CALL(mockgeo, getWidthX())
        .WillRepeatedly(Return(xmax - xmin));
    EXPECT_CALL(mockgeo, getWidthY())
        .WillRepeatedly(Return(ymax - ymin));
    EXPECT_CALL(mockgeo, getWidthZ())
        .WillRepeatedly(Return(zmax - zmin));
    /* Boundaries */
    EXPECT_CALL(mockgeo, getMinX())
        .WillRepeatedly(Return(xmin));
    EXPECT_CALL(mockgeo, getMaxX())
        .WillRepeatedly(Return(xmax));
    EXPECT_CALL(mockgeo, getMinY())
        .WillRepeatedly(Return(ymin));
    EXPECT_CALL(mockgeo, getMaxY())
        .WillRepeatedly(Return(ymax));
    EXPECT_CALL(mockgeo, getMinZ())
        .WillRepeatedly(Return(zmin));
    EXPECT_CALL(mockgeo, getMaxZ())
        .WillRepeatedly(Return(zmax));
    EXPECT_CALL(mockgeo, getGlobalMinZ())
        .WillRepeatedly(Return(zmin));
    EXPECT_CALL(mockgeo, getGlobalMaxZ())
        .WillRepeatedly(Return(zmax));
    /* Boundary conditions */
    setMockGeometryBCs(PERIODIC);

    /* Initialize tested track generator */
    tg = std::make_shared<TrackGenerator3DTest>
              (&mockgeo, num_azim, num_polar,
               azim_spacing, polar_spacing);
    tg->setNumThreads(num_threads);
    tg->setQuadrature(quad);
    tg->setSegmentFormation(+segmentationType::OTF_STACKS);
    tg->setSegmentationZones(seg_zones);
  }

  void setMockGeometryBCs(boundaryType bc) {
    EXPECT_CALL(mockgeo, getMinXBoundaryType())
        .WillRepeatedly(Return(bc));
    EXPECT_CALL(mockgeo, getMaxXBoundaryType())
        .WillRepeatedly(Return(bc));
    EXPECT_CALL(mockgeo, getMinYBoundaryType())
        .WillRepeatedly(Return(bc));
    EXPECT_CALL(mockgeo, getMaxYBoundaryType())
        .WillRepeatedly(Return(bc));
    EXPECT_CALL(mockgeo, getMinZBoundaryType())
        .WillRepeatedly(Return(bc));
    EXPECT_CALL(mockgeo, getMaxZBoundaryType())
        .WillRepeatedly(Return(bc));
  }

  int num_threads;

  /* user input */
  double azim_spacing;
  int num_azim;
  double polar_spacing;
  int num_polar;
  double tolerance;
  int max_iters;
  int axial_refines;

  double xmin, xmax, ymin, ymax, zmin, zmax;

  MockGeometry mockgeo;
  std::shared_ptr<TrackGenerator3DTest> tg;
  QuadraturePtr quad;

};


/* Common tests for cyclic track initialization */

#undef MY_TEST_FIXTURE
#undef MY_BC

#define MY_TEST_FIXTURE test_initializePeriodicTracks
#define MY_BC PERIODIC

#include "initializeTracks_common.inc"

/* Tests for reflective track initialization */

TEST_F(test_initializePeriodicTracks, booleanPeriodic) {
  ROOT_ONLY();

  tg->UTinitialize2DTracks();

  ASSERT_FALSE(tg->isReflectiveCyclic());
}


TEST_F(test_initializePeriodicTracks, periodicCircleAfterCyclicInitialization) {
  tg->UTinitialize2DChains();

  Track ****chains = tg->get2DTrackChains();
  Point *first_start, *last_start;
  Point *first_end, *last_end;

  for (int a = 0; a < num_azim/2; a++) {
    long links = tg->getNum2DLinks(a);
    for (long x = 0; x < tg->getMyNum2DChains(a); x++) {
      // Starting point of the first link and end point of the last link
      first_start = chains[a][x][0]->getStart();  // Starting track of the cyclic track
      last_end = chains[a][x][links-1]->getEnd(); // End track of the cyclic track
      ASSERT_NEAR(first_start->getX(), last_end->getX(), azim_spacing*1E-10)
          << "  Info: Azim = " << a << ", X index = " << x
          << ", Links = " << links << std::endl
          << "  Start X of 1st link: " << first_start->getX() << std::endl
          << "    End X of nth link: " << last_end->getX() << std::endl;
      // On different boundaries
      EXPECT_EQ(first_start->getY(), ymin);
      EXPECT_EQ(last_end->getY(), ymax);

      // End point of the first link and starting point of the last link differs
      first_end = chains[a][x][0]->getEnd();
      last_start = chains[a][x][links-1]->getStart();
      ASSERT_GT(first_end->distanceToPoint(last_start), azim_spacing/2.0)
          << "  Info: Azim = " << a << ", X index = " << x
          << ", Links = " << links << std::endl
          << "  Start: (" << first_end->getX() << ","
          << first_end->getY() << ")" << std::endl
          << "    End: (" << last_start->getX() << ","
          << last_start->getY() << ")" << std::endl;
    }
  }
}

// Check if the angles of the first and last links are complementary
TEST_F(test_initializePeriodicTracks, periodicLinkAzimuthalAngles) {
  tg->UTinitialize2DChains();

  Track ****chains = tg->get2DTrackChains();
  double first_phi, last_phi;

  for (int a = 0; a < num_azim/2; a++) {
    long links = tg->getNum2DLinks(a);
    for (long x = 0; x < tg->getMyNum2DChains(a); x++) {
      // Retrieve azim indices from track objects
      first_phi = chains[a][x][0]->getPhi();
      last_phi = chains[a][x][links-1]->getPhi();
      ASSERT_DOUBLE_EQ(first_phi, last_phi)
        << "  Info: Azim = " << a << ", X index = " << x
        << ", Links = " << links << std::endl
        << "  Start angle: phi = " << first_phi << std::endl
        << "    End angle: phi = " << last_phi << std::endl;
    }
  }
}


TEST_F(test_initializePeriodicTracks, trackStartingPointsOnYAxis) {
  tg->UTinitialize2DTracks();

  Track **tracks = tg->get2DTracks();

  for (int a = 0; a < num_azim/2; a++) {
    long num_chains = tg->getMyNum2DChains(a);
    long num_links = tg->getNum2DLinks(a);
    long xy = 0;
    for (int x = 0; x < num_chains; x++) {
      Track *track = &tracks[a][xy];
      ASSERT_DOUBLE_EQ(track->getStart()->getY(), ymin)
        << "Info:" << std::endl
        << "  # azim = " << a << ", # chains = " << num_chains
        << std::endl
        << "  # links = " << num_links << std::endl
        << " x = " << x << ", xy = " << xy << std::endl;
      xy += num_links;
    }
  }
}

} /* namespace */

#endif  // USTB_
