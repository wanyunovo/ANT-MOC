/**
 * @file test_initializeReflectiveTracks.cpp
 * @brief Test track initialization
 * @date April 17, 2019
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
  }

  void UTinitialize2DChains() {
    TrackGenerator::initializeTracks();
  }

  void UTinitialize3DTracks() {
    TrackGenerator::initializeTracks();
    initialize3DTracks();
  }

  void UTinitialize3DChains() {
    TrackGenerator::initializeTracks();
    initialize3DTracks();
    set3DIndexArraysAndTracks();
  }

  friend class test_initializeReflectiveTracks;
  FRIEND_TEST(test_initializeReflectiveTracks, allPositiveNumX);
  FRIEND_TEST(test_initializeReflectiveTracks, allPositiveNumY);
  FRIEND_TEST(test_initializeReflectiveTracks, compareLocalNumXsToGlobal);
  FRIEND_TEST(test_initializeReflectiveTracks, compareLocalNumYsToGlobal);
  FRIEND_TEST(test_initializeReflectiveTracks, checkMPIRankAfter2DChainUIDtoGUID);
  FRIEND_TEST(test_initializeReflectiveTracks, checkMPIRankAfter2DTrackUIDtoGUID);
  FRIEND_TEST(test_initializeReflectiveTracks, checkMPIRankAfter2DChainIDtoGID);
  FRIEND_TEST(test_initializeReflectiveTracks, checkMPIRankAfter2DTrackIDtoGID);
  FRIEND_TEST(test_initializeReflectiveTracks, convert2DChainUIDtoGUID);
  FRIEND_TEST(test_initializeReflectiveTracks, convert2DTrackUIDtoGUID);
  FRIEND_TEST(test_initializeReflectiveTracks, convert2DChainIDtoGID);
  FRIEND_TEST(test_initializeReflectiveTracks, convert2DTrackIDtoGID);
  FRIEND_TEST(test_initializeReflectiveTracks, identicalTranslationBetween2DTracks);
  FRIEND_TEST(test_initializeReflectiveTracks, startingPoints);
  FRIEND_TEST(test_initializeReflectiveTracks, project3DTrackStartingPoint);
};

// Test fixture
class test_initializeReflectiveTracks : public testing::Test {
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
    setMockGeometryBCs(REFLECTIVE);

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

#define MY_TEST_FIXTURE test_initializeReflectiveTracks
#define MY_BC REFLECTIVE

#include "initializeTracks_common.inc"

/* Tests for reflective track initialization */


TEST_F(test_initializeReflectiveTracks, booleanReflective) {
  ROOT_ONLY();

  tg->UTinitialize2DTracks();

  ASSERT_TRUE(tg->isReflectiveCyclic());
}


TEST_F(test_initializeReflectiveTracks, reflectiveCircleAfterCyclicInitialization) {
  tg->UTinitialize2DChains();

  Track ****chains = tg->get2DTrackChains();
  Point *first_start, *last_start;
  Point *first_end, *last_end;

  for (int a = 0; a < tg->getNumChainAzims(); a++) {
    auto links = tg->getNum2DLinks(a);
    for (int x = 0; x < tg->getMyNum2DChains(a); x++) {
      // Starting point and end point must be the same
      Track *&t_start = chains[a][x][0];
      Track *&t_end = chains[a][x][links-1];

      first_start = t_start->getStart();  // Starting track of the cyclic track
      last_end = t_end->getEnd(); // End track of the cyclic track
      ASSERT_NEAR(first_start->distanceToPoint(last_end), 0, azim_spacing*1E-10)
          << "  Info: Azim = " << a << ", X index = " << x << ", Links = " << links << std::endl
          << "  Start: (" << first_start->getX() << "," << first_start->getY() << ")" << std::endl
          << "    End: (" << last_end->getX() << "," << last_end->getY() << ")" << std::endl;
      // End points must differ
      first_end = t_start->getEnd();
      last_start = t_end->getStart();
      ASSERT_GT(first_end->distanceToPoint(last_start), azim_spacing/2.0)
          << "  Info: Azim = " << a << ", X index = " << x << ", Links = " << links << std::endl
          << "  Start: (" << first_end->getX() << "," << first_end->getY() << ")" << std::endl
          << "    End: (" << last_start->getX() << "," << last_start->getY() << ")" << std::endl;

    }
  }
}

// Check if the angles of the first and last links are complementary
TEST_F(test_initializeReflectiveTracks, reflectiveLinkAzimuthalAngles) {
  tg->UTinitialize2DChains();

  Track ****chains = tg->get2DTrackChains();
  int a1, a2, a3, a4;
  double first_phi, last_phi;

  for (int a = 0; a < tg->getNumChainAzims(); a++) {
    auto links = tg->getNum2DLinks(a);
    for (int x = 0; x < tg->getMyNum2DChains(a); x++) {
      // Retrieve azim indices from track objects
      Track *&t1 = chains[a][x][0];
      Track *&t2 = chains[a][x][links-1];
      a1 = t1->getAzimIndex();  // Starting track of the cyclic track
      a2 = t2->getAzimIndex(); // End track of the cyclic track
      first_phi = t1->getPhi();
      last_phi = t2->getPhi();

      EXPECT_EQ(a1 + a2 + 1, num_azim/2);

      if (t2->getReversed()) {
        // If the track was reversed
        ASSERT_DOUBLE_EQ(first_phi + last_phi, 2 * M_PI)
          << "  Info: Azim = " << a << ", X index = " << x << ", Links = " << links << std::endl
          << "  Start angle: a = " << a1 << ", " << "phi = " << first_phi << std::endl
          << "    End angle: a = " << a2 << ", " << "phi = " << last_phi << std::endl;
      }
      else {
        ASSERT_DOUBLE_EQ(first_phi + last_phi, M_PI)
          << "  Info: Azim = " << a << ", X index = " << x << ", Links = " << links << std::endl
          << "  Start angle: a = " << a1 << ", " << "phi = " << first_phi << std::endl
          << "    End angle: a = " << a2 << ", " << "phi = " << last_phi << std::endl;
      }
      // Compute azim indices according to link indices
      a3 = tg->get2DLinkAzimIndex(a, 0);
      a4 = tg->get2DLinkAzimIndex(a, links - 1);
      ASSERT_EQ(a1, a3)
        << "Track attribute differs from computed azim index";
      ASSERT_EQ(a2, a4)
        << "Track attribute differs from computed azim index";
    }
  }
}


TEST_F(test_initializeReflectiveTracks, trackStartingPointsOnYAxis) {
  tg->UTinitialize2DTracks();

  Track **tracks = tg->get2DTracks();

  for (int a = 0; a < num_azim/2; a++) {
    long num_chains = tg->getMyNum2DChains(a);
    long num_links = tg->getNum2DLinks(a) / 2;

    long xy = 0;
    if (a >= num_azim/4) {
      xy = num_links - 1;
    }

    for (int x = 0; x < num_chains; x++) {
      Track *track = &tracks[a][xy];
      if (track->getReversed())
        ASSERT_DOUBLE_EQ(track->getEnd()->getY(), ymin)
          << "Info:" << std::endl
          << "  # azim = " << a << ", # chains = " << num_chains
          << std::endl
          << "  # links = " << num_links << std::endl
          << " x = " << x << ", xy = " << xy << std::endl;
      else
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
