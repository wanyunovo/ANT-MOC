/**
 * @file TrackGenerator.h
 * @brief The TrackGenerator class.
 * @date January 23, 2012
 * @author William Boyd, MIT, Course 22 (wboyd@mit.edu)
 */

#ifndef TRACKGENERATOR_H_
#define TRACKGENERATOR_H_

#include "antmoc/container_utils.h"
#include "antmoc/enum_types.h"
#include "antmoc/mpi_utils.h"
#include "antmoc/Point.h"

#define _USE_MATH_DEFINES
#include <cmath>
#include <memory>
#include <omp.h>
#include <vector>

namespace antmoc
{

  /// \brief An aggregate class for track mesh dumping
  struct TrackMeshData
  {

    ///< Unique ID of the track
    long _uid;

    ///< Starting point of the track
    Point _start;

    ///< End point of the track
    Point _end;
  };

  /** Forward declarations */
  class Geometry;
  class Quadrature;
  class Timer;
  class Track;
  class TrackLoadBalance;

  using QuadraturePtr = std::shared_ptr<Quadrature>;
  using TrackLoadBalancePtr = std::shared_ptr<TrackLoadBalance>;

  /**
   * @class TrackGenerator TrackGenerator.h "src/TrackGenerator.h"
   * @brief The TrackGenerator is dedicated to generating and storing Tracks
   *        which cyclically wrap across the Geometry.
   * @details The TrackGenerator creates Track and initializes boundary
   *          conditions (vacuum, reflective, or periodic) for each Track.
   */
  class TrackGenerator
  {

  protected:
    /** The number of shared memory OpenMP threads */
    int _num_threads;

    /** Number of azimuthal angles in \f$ [0, 2 \pi] \f$ */
    int _num_azim; // 0~2

    /** The requested track azimuthal spacing (cm) */
    double _azim_spacing;

    /** The total number of Tracks for all azimuthal angles */
    long _num_2D_tracks;

    /** An integer array of the number of Tracks starting on the x-axis for each
     *  azimuthal angle */
    long *_num_x; //_num_x[a] 表示第 a 个方位角下从 x 轴（左右边）发射的轨迹条数

    /** An integer array of the number of Tracks starting on the y-axis for each
     *  azimuthal angle */
    long *_num_y; //_num_y[a] 表示从 y 轴（上下边）发射的条数。

    /** The number of cyclic tracks */
    long *_num_chains;

    /** The number of my cyclic tracks */
    long *_my_num_chains;

    /** The index of my first cyclic track */
    int *_chain_num_offsets;

    /** Length of cyclic tracks */
    double *_cyclic_length;

    /** A long integer array of the number of Tracks for each azimuthal angle */

    /** A 2D ragged array of 2D tracks (azim, track index) */
    Track **_tracks_2D; //[azim][local id] ，保存所有 2D 轨迹对象。

    /** A 1D array of Track pointers arranged by UID */
    Track **_tracks_2D_array;

    /** A 2D ragged array of 2D track chains (azim, x index, link index) */
    Track ****_tracks_2D_chains;

    /** Pointer to the Geometry */
    Geometry *_geometry;

    /** Boolean for whether to use Track input file (true) or not (false) */
    bool _use_input_file;

    /** Filename for the *.tracks input / output file */
    std::string _tracks_filename;

    /** OpenMP mutual exclusion locks for atomic FSR operations */
    omp_lock_t *_FSR_locks;

    /** Boolean indicating whether the Tracks have been generated (true) or not
     * (false) */
    bool _contains_2D_tracks;

    /** Boolean indicating whether 2D segments have been generated (true) or not
     * (false) */
    bool _contains_2D_segments;

    ///< Boolean indicating whether segments have been counted (true) or not (false)
    bool _segments_counted;

    /** Number of segments 所有 track 的 segment 数总和*/
    long _num_segments;

    /** The quadrature set */
    QuadraturePtr _quadrature;

    /** The z-coord where the 2D Tracks should be created */
    double _z_coord;

    /** Boolen to indicate whether a periodic BC exists */
    bool _periodic;

    /** Boolen to indicate whether a reflective BC exists */
    bool _reflective;

    /** Boolen to indicate whether a axial reflective BC exists */
    bool _z_reflective;

    /** Determines the type of track segmentation to use */
    segmentationType _segment_formation;

    /** Max optical segment length for Tracks before splitting */
    FP_PRECISION _max_optical_length;

    /** Maximum number of track segmenets in a single Track */
    int _max_num_segments;

    /** Boolean to indicate whether the segments should be dumped to file */
    bool _dump_segments;

    /** Boolean to indicate whether the segments have been centered around their
     * centroid or not */
    bool _segments_centered;

    /** A buffer holding the computed FSR volumes */
    FP_PRECISION *_FSR_volumes;

    /** A timer to record timing data for track generation */
    Timer *_timer;

    /** Geometry boundaries for this domain */
    double _x_min;
    double _y_min;
    double _z_min;
    double _x_max;
    double _y_max;
    double _z_max;

#ifdef ENABLE_MPI_
    /** MPI variables **/
    trackMappingType _track_mapping_type = trackMappingType::BLOCK;
#endif

    /** Private class methods */
    virtual void setContainsSegments(bool contains_segments);
    virtual void allocateTemporarySegments();
    virtual void resetStatus();
    virtual void initializeDefaultQuadrature();
    virtual void writeExtrudedFSRInfo(FILE *out);
    virtual void readExtrudedFSRInfo(FILE *in);
    virtual std::string getTestFilename(std::string directory);

  public:
    TrackGenerator(Geometry *geometry, int num_azim, double azim_spacing);
    virtual ~TrackGenerator();

    /* Get parameters */
    int getNumAzim();
    double getDesiredAzimSpacing();
    Geometry *getGeometry();
    virtual long getNumTracks();
    virtual long getNumSegments();
    long getNum2DTracks(int azim);
    long getNum2DTracks();
    long getNum2DSegments();
    bool getPeriodic();
    Track **get2DTracksArray();
    Track ****get2DTrackChains();
    virtual Track **getTracksArray();
    Track **get2DTracks();
    FP_PRECISION getMaxOpticalLength();
    int getMaxNumSegments();
    int getNumThreads();
    long getNumX(int azim);
    long getNumY(int azim);
    void exportFSRVolumes(double *out_volumes, int num_fsrs);
    void initializeVolumes();
    void initializeFSRVolumesBuffer();
    FP_PRECISION *getFSRVolumesBuffer();
    FP_PRECISION *getFSRVolumes();
    FP_PRECISION getFSRVolume(long fsr_id);
    double getZCoord();
    QuadraturePtr getQuadrature();
    FP_PRECISION retrieveMaxOpticalLength();
    omp_lock_t *getFSRLocks();
    segmentationType getSegmentFormation();
    virtual bool containsTracks();
    virtual bool containsSegments();
    long get2DTrackID(int a, long x);
    long *getTracksPerAzim();

    /* Set parameters */
    void setNumThreads(int num_threads);
    void setNumAzim(int num_azim);
    void setDesiredAzimSpacing(double spacing);
    virtual void setGeometry(Geometry *geometry);
    void setZCoord(double z_coord);
    void setQuadrature(QuadraturePtr quadrature);
    void setMaxOpticalLength(FP_PRECISION tau);
    void setMaxNumSegments(int max_num_segments);
    void setDumpSegments(bool dump_segments);

    /* Worker functions */
    virtual void segmentize();
    virtual void retrieveTrackCoords(double *coords, long num_tracks);
    void retrieve2DTrackCoords(double *coords, long num_tracks);
    virtual void retrieveSegmentCoords(double *coords, long num_segments);
    void retrieve2DSegmentCoords(double *coords, long num_segments);
    void generateFSRCentroids(FP_PRECISION *FSR_volumes);
    void generateTracks();
    void splitSegments(FP_PRECISION max_optical_length);
    double leastCommonMultiple(double a, double b);
    void dumpSegmentsToFile();
    bool readSegmentsFromFile();
    void initializeTrackFileDirectory();
    virtual void checkBoundaryConditions();
    virtual void countSegments();
    long countNumSegments();

    //--------------------------------------------------------------------
    // Common methods for track generation
    //--------------------------------------------------------------------
    virtual void initializeTracks();
    virtual void correct2DCyclicParameters(double *dx_eff, double *dy_eff);
    void correctNum2DCyclics();
    virtual void initializeTrackBCs(Track *track, int a, long i);
    virtual void initializeTracksArray();
    void set2DTrackEndpoints(Track *track, int a, long i,
                             double *dx_eff, double *dy_eff);

    /// \brief Check whether cyclic tracks are reflective
    /// \details Only ANT-MOC algorithms generate reflective cyclics. It
    ///          always returns False when called by OpenMOC algorithms.
    // 用于判断循环轨迹是否反射,由 OpenMOC 算法调用时始终返回 False。
    bool isReflectiveCyclic() { return _reflective; }

    /// \brief Check whether axial cyclic tracks are reflective
    bool isAxialReflectiveCyclic() { return _z_reflective; }

    int get2DLinkAzimIndex(int azim, long link);
    long get2DChainIndex(int azim, long xy);

    int getNumChainAzims();
    long getMyNum2DTracks(int azim);
    long getMyNum2DTracks();
    long getNum2DChains();
    long getMyNum2DChains();

    /// \brief Get the number of 2D chains per angle
    long getNum2DChains(int azim) { return _num_chains[azim]; }

    /// \brief Get the number of local 2D chains per angle
    long getMyNum2DChains(int azim) { return _my_num_chains[azim]; }

    long map2DTrackId2Gid(int azim, long xy);

    //--------------------------------------------------------------------
    // Methods for OpenMOC track generation
    //--------------------------------------------------------------------
#ifdef MIT_
    void initialize2DTracks(double *dx_eff, double *dy_eff);
    void initialize2DTrackChains();
    virtual void initializeTrackReflections(Track *track, int a, long i);
    long getNum2DLinks(int azim, long x);
#endif // MIT_

    //--------------------------------------------------------------------
    // Methods for cyclic track generation of ANT-MOC
    //--------------------------------------------------------------------
#ifdef USTB_
    void initialize2DTrackChains(double *dx_eff, double *dy_eff);
    virtual void initializeTrackReflections(Track *track, int a, long i);
    void reverseReflectiveTracks(Track *track);
    long getNum2DLinks(int azim);
#endif

    /// \brief Get the array of cyclic track lengths
    double *getCyclicLength() { return _cyclic_length; }

    /// \brief Get the length of cyclic tracks of a specified angle
    double getCyclicLength(int azim) { return _cyclic_length[azim]; }

    // Conversion between local and global ids
    long map2DChainId2Gid(int azim, long x);
    long map2DChainGid2Id(int azim, long gx);

#ifdef ENABLE_MPI_
    /// \brief Set track mapping type
    void setTrackMappingType(trackMappingType track_mapping_type);
    /// \brief Return track mapping type
    trackMappingType getTrackMappingType();

    /// \brief Create a mapping for 2D cyclic tracks
    void generate2DChainMapping();
    /// \brief Automatically choose a mapping algorithm
    /// \param load_map The load map for the selected algorithm
    TrackLoadBalancePtr generateOptimalLoadMap(IntVec2D &load_map);

    void preSegmentize2D(double z_coord);
#endif // ENABLE_MPI_

    //--------------------------------------------------------------------
    // Profiling functions
    //--------------------------------------------------------------------
    virtual void printTimerReport();
    virtual void printObjectStatistics();
    void printNumFSRs();
    double computeAzimUniformity();
    double computeStandardDeviation(const double *items, int n);

    /// \brief Report the memory usage of temporary segments.
    /// \details This method is for TrackGenerator3D.
    virtual void reportTemporarySegmentsStorage();

    //--------------------------------------------------------------------
    // Dumping functions
    //--------------------------------------------------------------------
    virtual std::vector<TrackMeshData> getFormattedTracks();
  };

} /* namespace antmoc */

#endif /* TRACKGENERATOR_H_ */
