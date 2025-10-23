/**
 * @file TrackGenerator3D.h
 * @brief The TrackGenerator3D class.
 * @date March 13, 2016
 * @author Geoffrey Gunow, MIT, Course 22 (geogunow@mit.edu)
 */

#ifndef TRACKGENERATOR3D_H_
#define TRACKGENERATOR3D_H_

#include "antmoc/TrackGenerator.h"
#include "antmoc/TrackIndexTypes.h"

#include <vector>

namespace antmoc
{

/** Forward declarations not in base class */
class Track;
class Track3D;
struct segment;


/**
 * @class TrackGenerator3D TrackGenerator3D.h "src/TrackGenerator3D.h"
 * @brief The TrackGenerator3D is dedicated to generating and storing Tracks
 *        which cyclically wrap across the Geometry in three dimensions.
 * @details The TrackGenerator creates Track and initializes boundary
 *          conditions (vacuum, reflective, or periodic) for each Track in
 *          three dimensions.
 */
class TrackGenerator3D : public TrackGenerator {

protected:

  /** Number of polar angles in \f$ [0, \pi] \f$ */
  int _num_polar;

  /** The requested track z-spacing (cm) */
  double _z_spacing;

  /** An array of the # of 3D tracks in each z-stack (azim, 2D track, polar) 每个轨迹堆中的轨迹数量*/
  int*** _tracks_per_stack;

  /** An array of the Track UID for the first Track and polar angle in each
    * z-stack (azim, xy, polar) */
  long*** _cum_tracks_per_stack;

  /** An array of the Track UID for the first Track in each z-stack
    * (azim, xy) */
  long** _cum_tracks_per_xy;

  /** Ann array of the first Track's l-z index for each z-stack */
  int*** _first_lz_of_stack;

  /** The total number of Tracks for all azimuthal and polar angles */
  long _num_3D_tracks;
  long _my_num_3D_tracks;

  /** Arrays of the number of Tracks starting on each axis */
  int** _num_z;
  int** _num_l;

  /* Arrays for the spacing along each axis */
  double** _dz_eff;
  double** _dl_eff;

  /** An array of 3D tracks (azim, 2D track, polar, z-stack) */
  Track3D**** _tracks_3D;

  /** The mesh defining axial heights of radial planes segmented in on-the-fly
    * calculations 存储着轴向高度的向量，这些高度用于在实时计算中划分径向平面*/
  std::vector<double> _segmentation_heights;

  /** The global axial mesh to use in on-the-fly calculations */
  std::vector<double> _global_z_mesh;

  /** Dimensions of temporary segments storage matrix */
  int _num_seg_matrix_rows;
  int _num_seg_matrix_columns;

  /** A matrix of temporary segments are created for on-the-fly methods to
    * improve efficiency */
  std::vector<segment*> _temporary_segments;
  bool _contains_temporary_segments;

  /** An array temporary Tracks are created for on-the-fly methods to
    * improve efficiency for every thread */
  std::vector<Track3D*> _temporary_3D_tracks;
  std::vector<Track**> _temporary_tracks_array;
  bool _contains_temporary_tracks;

  /** Maximum number of tracks a single 3D track stack for on-the-fly
   *  computation */
  int _max_num_tracks_per_stack;

  /** Booleans to indicate whether the Tracks and segments have been generated
   *  (true) or not (false) */
  bool _contains_3D_tracks;
  bool _contains_3D_segments;
  bool _contains_global_z_mesh;
  bool _contains_segmentation_heights;

  /** Private class methods */
  void initializeTracks();
  void initialize3DTracks();
  void set3DIndexArraysAndTracks();
  void recalibrateTracksToOrigin();
  void segmentize();
  void setContainsSegments(bool contains_segments);
  void allocateTemporarySegments();
  void allocateTemporaryTracks();
  void resetStatus();
  void initializeDefaultQuadrature();
  std::string getTestFilename(std::string directory);
  void getCycleTrackData(TrackChainIndexes* tcis, int num_cycles,
                         bool save_tracks);

  void segmentizeExtruded();
  void set3DTrackData(TrackChainIndexes* tci, Track3D* track,
                      bool create_arrays, bool save_tracks);
  double getLStart(TrackChainIndexes* tci);
  int getFirst2DTrackLinkIndex(TrackChainIndexes* tci, Track3D* track_3D = nullptr);
#ifdef USTB_
  // This method only should be used for testing.
  int getFirst2DCyclicLinkIndex(TrackChainIndexes *tci,
                                Track3D *track_3D = nullptr);
#endif

  void writeExtrudedFSRInfo(FILE* out);
  void readExtrudedFSRInfo(FILE* in);

public:

  TrackGenerator3D(Geometry* geometry, int num_azim, int num_polar,
                   double azim_spacing, double z_spacing);
  virtual ~TrackGenerator3D();

  /* Get parameters */
  int getNumPolar();
  double getDesiredZSpacing();
  long getNumTracks();
  long getNumSegments();
  long getNum3DTracks();
  long getNum3DSegments();
  Track3D**** get3DTracks();
  double getZSpacing(int azim, int polar);
  int getNumRows();
  int getNumColumns();
  segment* getTemporarySegments(int thread_id);
  Track3D* getTemporary3DTracks(int thread_id);
  Track** getTemporaryTracksArray(int thread_id);
  int getNumZ(int azim, int polar);
  int getNumL(int azim, int polar);
  int*** getTracksPerStack();
  int getMaxNumTracksPerStack();
  bool containsTracks();
  bool containsSegments();
  bool containsTemporaryTracks();
  long get3DTrackID(TrackStackIndexes* tsi);
  void convertTCItoTSI(TrackChainIndexes* tci, TrackStackIndexes* tsi);
  void convertTSItoTCI(TrackStackIndexes* tsi, TrackChainIndexes* tci);
  int get2DLinkIndex(TrackChainIndexes* tci);
  int getNum3DTrackChainLinks(TrackChainIndexes* tci);
  void getTSIByIndex(long id, TrackStackIndexes* tsi);

  void getTrackOTF(Track3D* track, TrackStackIndexes* tsi);

  /* Set parameters */
  void setNumPolar(int num_polar);
  void setDesiredZSpacing(double spacing);
  void setSegmentFormation(segmentationType segmentation_type);
  void setSegmentationZones(std::vector<double> z_mesh);
  void setLinkingTracks(TrackStackIndexes* tsi, TrackChainIndexes* tci,
                        bool outgoing, Track3D* track);
  void setTCILinkIndex(TrackChainIndexes* tci, TrackStackIndexes* tsi);
  void useGlobalZMesh();

  /* Worker functions */
  void retrieveTrackCoords(double* coords, long num_tracks);
  void retrieve3DTrackCoords(double* coords, long num_tracks);
  void retrieveGlobalZMesh(double*& z_mesh, int& num_fsrs);
  void retrieveSingle3DTrackCoords(double coords[6], long track_id);
  void retrieveSegmentCoords(double* coords, long num_segments);
  void retrieve3DSegmentCoords(double* coords, long num_segments);
  void create3DTracksArrays();
  void checkBoundaryConditions();
  long countNum3DExplicitSegments();


  //--------------------------------------------------------------------
  // Common ray tracing methods
  //--------------------------------------------------------------------

  /// \brief Get the number of my 3D tracks
  long getMyNum3DTracks() { return _my_num_3D_tracks; }

  long getNum3DChainStacks();
  long getMyNum3DChainStacks();
  long getNum3DChains();
  long getMyNum3DChains();

#ifdef ENABLE_MPI_
  void preSegmentizeExtruded(const std::vector<double> &z_coords);
#endif  // USTB_ && ENABLE_MPI_


  /* Profiling helplers */
  void printTimerReport();
#ifdef USTB_
  void printObjectStatistics();
#endif
  double countTracksMemory(char format = 'M');

  /// \brief Report the memory usage of temporary segments.
  /// \details This method requires an allreduce when MPI is enabled.
  void reportTemporarySegmentsStorage();

#ifdef ENABLE_MPI_
  /// \brief Report the range of workload by segments.
  void reportSegmentsWorkload();
#endif

  /* Dumping functions */
  std::vector<TrackMeshData> getFormattedTracks();
};

} /* namespace antmoc */

#endif /* TRACKGENERATOR3D_H_ */
