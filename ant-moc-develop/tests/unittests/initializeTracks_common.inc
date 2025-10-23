/**
 * @file initializeTracks_common.cpp
 * @brief Common tests for cyclic track initialization
 * @date May 19, 2019
 * @author An Wang, USTB (wangan@xs.ustb.edu.cn)
 */

#if defined(MY_TEST_FIXTURE) && defined(MY_BC)

#define SETPARAMS()           \
{                             \
  setMockGeometryBCs(MY_BC);  \
  tg->setGeometry(&mockgeo);  \
}

#if defined (USTB_) && !defined (ENABLE_MPI_)

TEST_F(MY_TEST_FIXTURE, allPositiveNumX) {
  tg->UTinitialize2DTracks();

  for (int a=0; a < num_azim/2; a++)
    EXPECT_GT(tg->getNumX(a), 0);
}

TEST_F(MY_TEST_FIXTURE, allPositiveNumY) {
  tg->UTinitialize2DTracks();

  for (int a=0; a < num_azim/2; a++)
    EXPECT_GT(tg->getNumY(a), 0);
}


TEST_F(MY_TEST_FIXTURE, numLocal2DTracksEqualsGlobal) {
  tg->UTinitialize2DTracks();

  for (int a = 0; a < num_azim/2; a++)
    ASSERT_EQ(tg->getMyNum2DTracks(a), tg->getNum2DTracks(a));
}

TEST_F(MY_TEST_FIXTURE, numLocal2DChainsEqualsGlobal) {
  tg->UTinitialize2DTracks();

  for (int a = 0; a < tg->getNumChainAzims(); a++)
    ASSERT_EQ(tg->getMyNum2DChains(a), tg->getNum2DChains(a));
}

#endif

#if defined (USTB_) && defined (ENABLE_MPI_)

TEST_F(MY_TEST_FIXTURE, numLocal2DTracksEqualsGlobal) {
  tg->UTinitialize2DTracks();

  int n_azims = tg->getNumAzim()/2;
  long *my_nxy = new long[n_azims];

  for (int a = 0; a < n_azims; a++) {
    my_nxy[a] = tg->getMyNum2DTracks(a);
  }
  MPI_Allreduce(MPI_IN_PLACE, my_nxy, n_azims, MPI_LONG, MPI_SUM, mpi::getMPIComm());

  for (int a = 0; a < n_azims; a++) {
    EXPECT_EQ(my_nxy[a], tg->getNum2DTracks(a));
  }
  delete [] my_nxy;
}


// Check if the sum of local chains is the same as the global count
TEST_F(MY_TEST_FIXTURE, numLocal2DChainsEqualsGlobal) {
  tg->UTinitialize2DTracks();

  int n_azims = tg->getNumChainAzims();
  long *my_nc = new long[n_azims];

  for (int a = 0; a < n_azims; a++) {
    my_nc[a] = tg->getMyNum2DChains(a);
  }
  MPI_Allreduce(MPI_IN_PLACE, my_nc, n_azims, MPI_LONG, MPI_SUM, mpi::getMPIComm());

  for (int a = 0; a < n_azims; a++) {
    EXPECT_EQ(my_nc[a], tg->getNum2DChains(a));
  }
  delete [] my_nc;
}


// Check if the number of local 3D chain stacks equals the global count
TEST_F(MY_TEST_FIXTURE, numLocal3DChainStacksEqualsGlobal) {
  tg->UTinitialize3DChains();

  long my_cs = tg->getMyNum3DChainStacks();

  MPI_Allreduce(MPI_IN_PLACE, &my_cs, 1, MPI_LONG, MPI_SUM, mpi::getMPIComm());

  EXPECT_EQ(my_cs, tg->getNum3DChainStacks());
}


TEST_F(MY_TEST_FIXTURE, numLocal3DChainsEqualsGlobal) {
  tg->UTinitialize3DChains();

  long my_c = tg->getMyNum3DChains();

  MPI_Allreduce(MPI_IN_PLACE, &my_c, 1, MPI_LONG, MPI_SUM, mpi::getMPIComm());

  EXPECT_EQ(my_c, tg->getNum3DChains());
}


//TEST_F(MY_TEST_FIXTURE, getMPIRankBy2DChainGuid) {
//  tg->UTinitialize2DTracks();
//
//  for (long i = 0; i < tg->getMyNum2DChains(); i++)
//    EXPECT_EQ(tg->getRankBy2DChainGuid(tg->map2DChainUid2Guid(i)),
//              mpi::getRankSharedDomain());
//}
//
//
//TEST_F(MY_TEST_FIXTURE, getMPIRankBy2DChainAzimGid) {
//  tg->UTinitialize2DTracks();
//
//  for (int a = 0; a < tg->getNumChainAzims(); a++)
//    for (long i = 0; i < tg->getMyNum2DChains(a); i++)
//      EXPECT_EQ(tg->getRankBy2DChainGid(a, tg->map2DChainId2Gid(a, i)),
//                mpi::getRankSharedDomain());
//}


#endif /* USTBMPI */

///* Test conversion between ids */
//TEST_F(MY_TEST_FIXTURE, convert2DChainUIDtoGUID) {
//  tg->UTinitialize2DChains();
//
//  for (long i = 0; i < tg->getMyNum2DChains(); i++)
//    EXPECT_EQ(tg->map2DChainGuid2Uid(tg->map2DChainUid2Guid(i)), i);
//}


TEST_F(MY_TEST_FIXTURE, convert2DChainIDtoGID) {
  tg->UTinitialize2DChains();

  for (int a = 0; a < tg->getNumChainAzims(); a++)
    for (long i = 0; i < tg->getMyNum2DChains(a); i++)
      EXPECT_EQ(tg->map2DChainGid2Id(a, tg->map2DChainId2Gid(a, i)), i)
        << "  Info: Azim = " << a << std::endl;
}


//// This test is disabled because it only useful for rectangular geometries
// TEST_F(MY_TEST_FIXTURE, expectedNumLinksAlongCyclicTracks) {
//   tg->UTinitialize2DChains();
// 
//   for (int a = 0; a < num_azim/2; a++) {
//     if (tg->isReflectiveCyclic()) {
//       // Reflective
//       ASSERT_EQ(tg->getNum2DLinks(a), 2 * (tg->getCyclicLinksX(a) + tg->getCyclicLinksY(a)));
//     }
//     else {
//       // Periodic
//       ASSERT_EQ(tg->getNum2DLinks(a), tg->getCyclicLinksX(a) + tg->getCyclicLinksY(a));
//     }
//   }
// }


TEST_F(MY_TEST_FIXTURE, check2DLinkIndices) {
  tg->UTinitialize2DChains();

  Track ****chains = tg->get2DTrackChains();
  // Loop over links
  for (int a = 0; a < tg->getNumChainAzims(); a++) {
    for (long x = 0; x < tg->getMyNum2DChains(a); x++){
      for (long link = 0; link < tg->getNum2DLinks(a); link++) {
        Track *&track_2D = chains[a][x][link];
        EXPECT_EQ(track_2D->getChainIndex(), x)
          << "  Chain info: Azim = " << a << ", Link index = " << link << std::endl
          << "  Track info: Azim = " << track_2D->getAzimIndex()
          << ", XY index = " << track_2D->getXYIndex() << std::endl;
        EXPECT_EQ(track_2D->getLinkIndex(), link)
          << "  Chain info: Azim = " << a << ", X index = " << x << std::endl
          << "  Track info: Azim = " << track_2D->getAzimIndex()
          << ", XY index = " << track_2D->getXYIndex() << std::endl;
      }
    }
  }
}

TEST_F(MY_TEST_FIXTURE, check2DLinkAzimuthalAngles) {
  tg->UTinitialize2DChains();

  Track ****chains = tg->get2DTrackChains();
  // Loop over links
  for (int a = 0; a < tg->getNumChainAzims(); a++) {
    for (long x = 0; x < tg->getMyNum2DChains(a); x++) {
      for (long link = 0; link < tg->getNum2DLinks(a); link++) {
        int link_azim = tg->get2DLinkAzimIndex(a, link);
        Track *&track_2D = chains[a][x][link];

        ASSERT_EQ(track_2D->getAzimIndex(), link_azim);
        if (track_2D->getReversed()) {
          // If the track was reversed
          EXPECT_DOUBLE_EQ(track_2D->getPhi(), quad->getPhi(link_azim) + M_PI)
            << "  Chain info: Azim = " << a << ", Link index = " << link << std::endl
            << "  Track info: Azim = " << track_2D->getAzimIndex()
            << ", XY index = " << track_2D->getXYIndex() << std::endl;
        }
        else {
          EXPECT_DOUBLE_EQ(track_2D->getPhi(), quad->getPhi(link_azim))
            << "  Chain info: Azim = " << a << ", Link index = " << link << std::endl
            << "  Track info: Azim = " << track_2D->getAzimIndex()
            << ", XY index = " << track_2D->getXYIndex() << std::endl;
        }
      }
    }
  }
}


#if defined (USTB_) && !defined (ENABLE_MPI_)

TEST_F(MY_TEST_FIXTURE, DISABLED_check2DTrackChainLength) {
  tg->UTinitialize2DChains();

  Track ****chains = tg->get2DTrackChains();
  // Loop over links
  for (int a = 0; a < tg->getNumChainAzims(); a++) {
    for (long x = 0; x < tg->getNum2DChains(a); x++){
      double computed_length = 0;
      for (long link = 0; link < tg->getNum2DLinks(a); link++) {
        Track *&track_2D = chains[a][x][link];
        computed_length += track_2D->getLength();
      }
      ASSERT_NEAR(computed_length, tg->getCyclicLength(a), 1E-12)
        << "  Chain info: Azim = " << a << ", X index = " << x << std::endl
        << "              Links = " << tg->getNum2DLinks(a) << std::endl;
    }
  }
}


// Test two different algorithms
TEST_F(MY_TEST_FIXTURE, project3DTrackStartingPoint) {
  tg->UTinitialize3DChains();

  int num_chain_azims = tg->getNumChainAzims();
  TrackChainIndexes tci;
  int nl, nz;
  Track3D track1, track2;
  long link1, link2;

  // Loop over TCIs
  for (int a = 0; a < num_chain_azims; a++) {
    tci._azim  = a;
    for (long x = 0; x < tg->getNum2DChains(a); x++) {
      tci._x = x;
      for (int p = 0; p < num_polar; p++) {
        tci._polar = p;
        nl = tg->getNumL(a, p);
        nz = tg->getNumZ(a, p);
        // Loop over 3D chains
        for (int lz=0; lz < nl + nz; lz++) {
          tci._lz = lz;
          tci._link = -1;
          link1 = tg->getFirst2DTrackLinkIndex(&tci, &track1);
          link2 = tg->getFirst2DCyclicLinkIndex(&tci, &track2);

          /* The two algorithm produce the same result only if
           * the tracks don't connect on y-boudaries
           */
          //if (x + link1 * tg->getNumX(a) >= tg->getNumY(a))
          //  continue;

          // Check link indices
          ASSERT_EQ(link1, link2)
            << "  link1 from getFirst2DTrackLinkIndex()" << std::endl
            << "  link2 from getFirst2DCyclicLinkIndex()" << std::endl
            << "  num_chain_azims = " << num_chain_azims
            << ", num_polar = " << num_polar << std::endl
            << "  num_chains = " << tg->getNum2DChains(a)
            << ", num_links = " << tg->getNum2DLinks(a) << std::endl
            << "  dx = " << (xmax-xmin)/tg->getNumX(a)
            << ", dy = " << (ymax-ymin)/tg->getNumY(a)
            << ", dl = " << tg->getCyclicLength(a)/tg->getNumL(a, p)
            << ", dz = " << (zmax-zmin)/tg->getNumZ(a, p) << std::endl
            << "  TCI(" << a << "," << x << "," << p << "," << lz << ")";
          // Check output starting points of two algorithms
          double distance = track1.getStart()->distanceToPoint(track2.getStart());
          ASSERT_NEAR(distance, 0, 1E-06)
            << "  Distance between two starting points" << std::endl
            << "  track1 from getFirst2DTrackLinkIndex()" << std::endl
            << "  track2 from getFirst2DCyclicLinkIndex()" << std::endl
            << "  TCI(" << a << "," << x << "," << p << "," << lz << ")";
          // Check output end points of two algorithms
          ASSERT_DOUBLE_EQ(track1.getEnd()->getZ(), track2.getEnd()->getZ())
            << "  track1 from getFirst2DTrackLinkIndex()" << std::endl
            << "  track2 from getFirst2DCyclicLinkIndex()" << std::endl
            << "  TCI(" << a << "," << x << "," << p << "," << lz << ")";
        }
      }
    }
  }
}

#endif /* no USTBMPI */


TEST_F(MY_TEST_FIXTURE, conversionFromTCItoTSItoTCI) {
  tg->UTinitialize3DChains();

  TrackChainIndexes tci;
  TrackStackIndexes tsi;
  int nl, nz;

  // Loop over TCIs
  for (int a = 0; a < tg->getNumChainAzims(); a++) {
    tci._azim = a;
    for (long x = 0; x < tg->getMyNum2DChains(a); x++) {
      tci._x = x;
      for (int p = 0; p < num_polar; p++) {
        tci._polar = p;
        nl = tg->getNumL(a, p);
        nz = tg->getNumZ(a, p);
        for (long lz=0; lz < nl + nz; lz++) {
          tci._lz = lz;
          // Loop over 3D links
          for (long link = 0; link < tg->getNum3DTrackChainLinks(&tci); link++) {
            tci._link = link;
            tg->convertTCItoTSI(&tci, &tsi);
            tg->convertTSItoTCI(&tsi, &tci);
            // Check if they are the same
            EXPECT_EQ(tci._azim, a);
            EXPECT_EQ(tci._x, x);
            EXPECT_EQ(tci._polar, p);
            EXPECT_EQ(tci._lz, lz);
            EXPECT_EQ(tci._link, link);
          }
        }
      }
    }
  }
}


TEST_F(MY_TEST_FIXTURE, conversionFromTSItoTCItoTSI) {
  tg->UTinitialize3DChains();

  TrackChainIndexes tci;
  TrackStackIndexes tsi;
  int ***tracks_per_stack = tg->getTracksPerStack();

  // Loop over TSIs
  for (int a = 0; a < num_azim/2; a++) {
    tsi._azim = a;
    for (long xy = 0; xy < tg->getMyNum2DTracks(a); xy++) {
      tsi._xy = xy;
      for (int p = 0; p < num_polar; p++) {
        tsi._polar = p;
        int stack_size = tracks_per_stack[a][xy][p];
        // Loop over stack tracks
        for (long z = 0; z < stack_size; z++) {
          tsi._z = z;

          tg->convertTSItoTCI(&tsi, &tci);
          tg->convertTCItoTSI(&tci, &tsi);
          // Check if they are the same
          EXPECT_EQ(tsi._azim, a);
          EXPECT_EQ(tsi._xy, xy);
          EXPECT_EQ(tsi._polar, p);
          EXPECT_EQ(tsi._z, z);
        }
      }
    }
  }
}


TEST_F(MY_TEST_FIXTURE, set3DIndexArraysAndTracks) {
  tg->UTinitialize3DChains();
}

#endif
