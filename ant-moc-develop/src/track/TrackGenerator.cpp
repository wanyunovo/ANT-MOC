#include "antmoc/TrackGenerator.h"
#include "antmoc/Cmfd.h"
#include "antmoc/Cell.h"
#include "antmoc/file_utils.h"
#include "antmoc/Geometry.h"
#include "antmoc/math_utils.h"
#include "antmoc/openmp_utils.h"
#include "antmoc/Progress.h"
#include "antmoc/Quadrature.h"
#include "antmoc/string_utils.h"
#include "antmoc/Timer.h"
#include "antmoc/Track.h"
#include "antmoc/TrackLoadBalance.h"
#include "antmoc/TrackTraversingAlgorithms.h"
#include "antmoc/Lattice.h"

#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

namespace antmoc
{

  /**
   * @brief Constructor for the TrackGenerator assigns default values.
   * @param geometry a pointer to a Geometry object
   * @param num_azim number of azimuthal angles in \f$ [0, 2\pi] \f$
   * @param spacing track spacing (cm)
   */
  TrackGenerator::TrackGenerator(Geometry *geometry, int num_azim,
                                 double azim_spacing) : _segment_formation(+segmentationType::EXPLICIT_2D)
  {

    setGeometry(geometry);
    setNumThreads(1);
    setNumAzim(num_azim);
    setDesiredAzimSpacing(azim_spacing);

    _contains_2D_tracks = false;
    _contains_2D_segments = false;
    _segments_counted = false;
    _num_segments = 0;
    _quadrature = NULL;
    _z_coord = 0.0;
    _max_optical_length = std::numeric_limits<FP_PRECISION>::max();
    _max_num_segments = 0;
    _FSR_volumes = NULL;
    _dump_segments = true;
    _segments_centered = false;
    _FSR_locks = NULL;
    _tracks_2D_array = NULL;
    _tracks_2D_chains = NULL;
    _tracks_per_azim = NULL;
    _num_chains = nullptr;
    _my_num_chains = nullptr;
    _cyclic_length = nullptr;
    _timer = new Timer();
  }

  /**
   * @brief Destructor frees memory for all Tracks.
   */
  TrackGenerator::~TrackGenerator()
  {

    if (_contains_2D_tracks)
    {

      /* Delete 2D tracks and mappings */
      for (int a = 0; a < _num_azim / 2; a++)
      {
        delete[] _tracks_2D[a];
      }
      delete[] _tracks_2D;

      delete[] _tracks_2D_array;

      /* Delete 2D chains if created */
      if (_tracks_2D_chains != NULL)
      {
        for (int a = 0; a < getNumChainAzims(); a++)
        {
          for (int x = 0; x < getMyNum2DChains(a); x++)
            delete[] _tracks_2D_chains[a][x];
          delete[] _tracks_2D_chains[a];
        }
        delete[] _tracks_2D_chains;
      }

      /* Delete track laydown information */
      delete[] _num_x;
      delete[] _num_y;
      delete[] _num_chains;
      delete[] _cyclic_length;

#ifdef ENABLE_MPI_
      if (mpi::isPrdTrackDecomposed())
        delete[] _my_num_chains;
      delete[] _chain_num_offsets;
#endif
    }

    delete[] _tracks_per_azim;

    if (_FSR_locks != NULL)
      delete[] _FSR_locks;

    if (_FSR_volumes != NULL)
      delete[] _FSR_volumes;

    delete _timer;
    _timer = nullptr;
  }

  /**
   * @brief Return the number of tracks for each azimuthal angle
   * @return the number of tracks for each azimuthal angle
   */
  long *TrackGenerator::getTracksPerAzim()
  {
    return _tracks_per_azim;
  }

  /**
   * @brief Return the number of azimuthal angles in \f$ [0, 2\pi] \f$
   * @return the number of azimuthal angles in \f$ 2\pi \f$
   */
  int TrackGenerator::getNumAzim()
  {
    return _num_azim;
  }

  /**
   * @brief Return the track azimuthal spacing (cm).
   * @ditails This will return the user-specified track spacing and NOT the
   *          effective track spacing which is computed and used to generate
   *          cyclic tracks.
   * @return the track azimuthal spacing (cm)
   */
  double TrackGenerator::getDesiredAzimSpacing()
  {
    return _azim_spacing;
  }

  /**
   * @brief Return the Geometry for this TrackGenerator if one has been set.
   * @return a pointer to the Geometry
   */
  Geometry *TrackGenerator::getGeometry()
  {
    if (_geometry == NULL)
      log::ferror("Unable to return the TrackGenerator's Geometry "
                  "since it has not yet been set");

    return _geometry;
  }

  /**
   * @brief Return the array of FSR locks for atomic FSR operations.
   * @return an array of FSR locks
   */
  omp_lock_t *TrackGenerator::getFSRLocks()
  {
    if (_FSR_locks == NULL)
      log::ferror("Unable to return the TrackGenerator's FSR locks "
                  "since they have not yet been created");

    return _FSR_locks;
  }

  /**
   * @brief Initialize an array to contain the FSR volumes.
   * 为存储 FSR 体积的数组分配内存
   */
  void TrackGenerator::initializeFSRVolumesBuffer()
  {

    if (_FSR_volumes != NULL)
      delete[] _FSR_volumes;

#pragma omp critical // OpenMP 临界区，保证多线程下只有一条线程执行此代码块，防止并发分配冲突
    {
      long num_FSRs = _geometry->getNumFSRs();
      _FSR_volumes = new FP_PRECISION[num_FSRs](); // 分配新数组并初始化为0
    }
  }

  /**
   * @brief Return the array used to store the FSR volumes
   * @return _FSR_volumes the FSR volumes array indexed by FSR ID
   */
  FP_PRECISION *TrackGenerator::getFSRVolumesBuffer()
  {

    return _FSR_volumes;
  }

  /**
   * @brief Return the total number of Tracks across the Geometry.
   * @return the total number of Tracks
   */
  long TrackGenerator::getNumTracks()
  {
    return getNum2DTracks();
  }

  /**
   * @brief Return the number of 2D Tracks across the Geometry.
   * @return the total number of 2D Tracks
   */
  long TrackGenerator::getNum2DTracks(int azim)
  {
    return _num_x[azim] + _num_y[azim];
  }

  /**
   * @brief Return the total number of 2D Tracks across the Geometry.
   * @return the total number of 2D Tracks
   */
  long TrackGenerator::getNum2DTracks()
  {

    long num_2D_tracks = 0;

    for (int a = 0; a < _num_azim / 2; a++)
      num_2D_tracks += _num_x[a] + _num_y[a];

    return num_2D_tracks;
  }

  /**
   * @brief Return the total number of Track segments across the Geometry.
   * @return the total number of Track segments
   */
  long TrackGenerator::getNumSegments()
  {
    return getNum2DSegments();
  }

  /**
   * @brief Return the total number of 2D Track segments across the Geometry.
   * @return the total number of 2D Track segments
   */
  long TrackGenerator::getNum2DSegments()
  {

    if (!TrackGenerator::containsSegments())
      log::error("Cannot get the number of 2D segments since they have not been generated.");

    long num_2D_segments = 0;

    long num_2D_tracks = getMyNum2DTracks();
#pragma omp parallel for reduction(+ : num_2D_segments)
    for (long i = 0; i < num_2D_tracks; i++)
    {
      num_2D_segments += _tracks_2D_array[i]->getNumSegments();
    }

    return num_2D_segments;
  }

  /**
   * @brief Returns an array of the Track pointers by increasing UID
   * @details An array of pointers to all 2D Track objects in the Geometry is
   *          returned, arranged by increasing unique identifier (UID).
   * @return the array of Track pointers
   */
  Track **TrackGenerator::get2DTracksArray()
  {

    if (!TrackGenerator::containsTracks())
      log::ferror("Unable to return the 1D array of Tracks "
                  "since Tracks have not yet been generated.");

    return _tracks_2D_array;
  }

  /**
   * @brief Returns an array of the Track pointers by increasing UID.
   * @details Calls TrackGenerator::get2DTracksArray to return all 2D Tracks
   * @return the array of Track pointers
   */
  Track **TrackGenerator::getTracksArray()
  {
    return get2DTracksArray();
  }

  /**
   * @brief Get the array of 2D Chains.
   * @return a pointer to the array of 2D Chains
   */
  Track ****TrackGenerator::get2DTrackChains()
  {
    return _tracks_2D_chains;
  }

  /**
   * @brief Returns a 2D jagged array of the 2D Tracks.
   * @details The first index into the array is the azimuthal angle and the
   *          second index is the Track number.
   * @return the 2D jagged array of 2D Tracks
   */
  Track **TrackGenerator::get2DTracks()
  {

    if (!TrackGenerator::containsTracks())
      log::ferror("Unable to return the 3D ragged array of the 2D Tracks "
                  "since Tracks have not yet been generated.");

    return _tracks_2D;
  }

  /**
   * @brief Calculates and returns the maximum optcial length for any segment
   *        in the Geomtry.
   * @details The _max_optical_length value is recomputed, updated, and returned.
   *          This value determines the when segments must be split during ray
   *          tracing.
   * @return _max_optical_length the maximum optical length of any segment in the
   *         Geometry
   */
  FP_PRECISION TrackGenerator::getMaxOpticalLength()
  {
    MaxOpticalLength update_max_optical_length(this);
    update_max_optical_length.execute();

    // Indicate that segments should be counted again.
    _segments_counted = false;

    return _max_optical_length;
  }

  /**
   * @brief Returns the maximum number of segments along a single track
   * @details The TrackGenerator::countSegments routine must be called before
   *          this function will return a correct value
   * @return the maximum number of segments
   */
  int TrackGenerator::getMaxNumSegments()
  {
    return _max_num_segments;
  }

  /**
   * @brief Returns the number of shared memory OpenMP threads in use.
   * @return the number of threads
   */
  int TrackGenerator::getNumThreads()
  {
    return _num_threads;
  }

  /**
   * @brief Returns the number of 2D Tracks in the x-direction for a given
   *        azimuthal angle index
   * @param azim the azimuthal angle index
   * @return the number of 2D Tracks in the x-direction of the Geometry
   */
  long TrackGenerator::getNumX(int azim)
  {
    return _num_x[azim];
  }

  /**
   * @brief Returns the number of 2D Tracks in the y-direction for a given
   *        azimuthal angle index
   * @param azim the azimuthal angle index
   * @return the number of 2D Tracks in the y-direction of the Geometry
   */
  long TrackGenerator::getNumY(int azim)
  {
    return _num_y[azim];
  }

  /**
   * @brief FSR volumes are coppied to an array input by the user
   * @param out_volumes The array to which FSR volumes are coppied
   * @param num_fsrs The number of FSR volumes to copy. The first num_fsrs
   *        volumes stored in the FSR volumes array are coppied.
   */
  void TrackGenerator::exportFSRVolumes(double *out_volumes, int num_fsrs)
  {

    for (int i = 0; i < num_fsrs; i++)
      out_volumes[i] = _FSR_volumes[i];
  }

  /**
   * @brief Computes and returns an array of volumes indexed by FSR.
   * @details Note: The memory is stored in the FSR volumes buffer of the
   *          TrackGenerator and is freed during deconstruction.
   * @return a pointer to the array of FSR volumes
   * 计算并返回 FSR 体积数组
   * _FSR_volumes[fsr_id] =该 FSR 的体积
   * 在 MOC 中，区域的体积通常是通过统计穿过该区域的所有特征线段的长度及权重累加得到的（数值积分）。它内部会调用一个 VolumeCalculator 类来执行这个计算
   */
  FP_PRECISION *TrackGenerator::getFSRVolumes()
  {

    log::finfo("Computing volumes of FSRs...");

    /* Reset FSR volumes to zero */
    long num_FSRs = _geometry->getNumFSRs();

    if (_FSR_volumes != NULL)
      memset(_FSR_volumes, 0., num_FSRs * sizeof(FP_PRECISION));

    /* Create volume calculator and calculate new FSR volumes */
    VolumeCalculator volume_calculator(this);
    volume_calculator.execute();

#ifdef ENABLE_MPI_
    if (mpi::isPrdTrackDecomposed())
    {
      // Synchronize volumes across processes
      _timer->startTimer();

      float size_mb = static_cast<float>(num_FSRs) * sizeof(FP_PRECISION) / (1 << 20);
      log::verbose_once("Reducing FSR volumes...");
      log::verbose_once("Buffer size for reducing FSR volumes = {:.3f} MiB", size_mb);

      MPI_Allreduce(MPI_IN_PLACE,
                    _FSR_volumes,
                    num_FSRs,
                    mpi::getDatatype<FP_PRECISION>(),
                    MPI_SUM,
                    mpi::getCommSharedDomain());
      _timer->stopTimer("Transfer FSRs");
    }
#endif

    double tot_volume = 0;
    /* Check to ensure all FSRs are crossed by at least one track */
    for (long i = 0; i < num_FSRs; i++)
    {
      tot_volume += _FSR_volumes[i];

      // Found zero-volume FSR
      if (fabs(_FSR_volumes[i]) < FLT_EPSILON)
      {
        auto point_xyz = _geometry->getFSRPoint(i)->getXYZ();
        log::fwarn("Zero volume calculated for FSR %d at Point (%.2f %.2f %.2f). Perhaps there is no track "
                   "traversed this FSR. The FSR is contained in Cell: %d",
                   i,
                   point_xyz[0], point_xyz[1], point_xyz[2], _geometry->findCellContainingFSR(i)->getId());
      }
    }

#ifdef ENABLE_MPI_
    if (mpi::isSpatialDecomposed())
    {
      MPI_Allreduce(MPI_IN_PLACE, &tot_volume, 1, mpi::getDatatype<double>(),
                    MPI_SUM, mpi::getMPIComm());
    }
#endif

    log::info("Total FSR volume       = {:.10e}", tot_volume);
    log::info("Average volume per FSR = {:.10e}", tot_volume / _geometry->getNumTotalFSRs());

    return _FSR_volumes;
  }

  /**
   * @brief Returns the volume of an FSR.
   * @param fsr_id the ID for the FSR of interest
   * @return the FSR volume
   */
  FP_PRECISION TrackGenerator::getFSRVolume(long fsr_id)
  {

    if (_FSR_volumes == NULL)
      log::ferror("Unable to get the FSR volume since FSR volumes "
                  "have not yet been generated");

    else if (fsr_id < 0 || fsr_id > _geometry->getNumFSRs())
      log::ferror("Unable to get the volume for FSR %d since the FSR IDs "
                  "lie in the range (0, %d)",
                  fsr_id, _geometry->getNumFSRs());
    return _FSR_volumes[fsr_id];
  }

  /**
   * @brief Returns the z-coord of the radial plane used in 2D calcualtions
   * @return the z-coord of the 2D calculation
   */
  double TrackGenerator::getZCoord()
  {
    return _z_coord;
  }

  /**
   * @brief Returns the Quadrature object
   * @return the Quadrature object
   */
  QuadraturePtr TrackGenerator::getQuadrature()
  {
    return _quadrature;
  }

  /**
   * @brief Sets the number of shared memory OpenMP threads to use (>0).
   * @param num_threads the number of threads
   */
  void TrackGenerator::setNumThreads(int num_threads)
  {

    if (num_threads <= 0)
      log::error("Unable to set the number of threads for the "
                 "TrackGenerator to {} since it is less than or equal to 0",
                 num_threads);

#ifdef ENABLE_MPI_
    /* Check the MPI library has enough thread support */
    int provided;
    MPI_Query_thread(&provided);
    if (num_threads > 1 && provided < MPI_THREAD_SERIALIZED)
      log::warn("Not enough thread support level in the MPI library, re-compile "
                "with another library. Thread support level should be at least "
                "MPI_THREAD_SERIALIZED.");
#endif

    _num_threads = num_threads;

    /* Set the number of threads for OpenMP */
    omp_set_num_threads(_num_threads);
    if (_geometry != NULL)
      _geometry->reserveKeyStrings(num_threads);

    // Verbose messages about OpenMP
    if (num_threads > 1)
    {
      // Print thread affinity
      openmp::print_affinity_policy();

      // Print OpenMP places
      openmp::print_places();

      // Print CPU assignments
      openmp::print_cpu_bind();
    }

#ifdef ENABLE_MPI_
    mpi::mpiBarrier();
#endif
  }

  /**
   * @brief Set the number of azimuthal angles in \f$ [0, 2\pi] \f$.
   * @param num_azim the number of azimuthal angles in \f$ 2\pi \f$
   */
  void TrackGenerator::setNumAzim(int num_azim)
  {

    if (num_azim < 0)
      log::ferror("Unable to set a negative number of azimuthal angles "
                  "%d for the TrackGenerator.",
                  num_azim);

    if (num_azim % 4 != 0)
      log::ferror("Unable to set the number of azimuthal angles to %d for "
                  "the TrackGenerator since it is not a multiple of 4",
                  num_azim);

    _num_azim = num_azim;
    resetStatus();
  }

  /**
   * @brief Set the suggested azimuthal track spacing (cm).
   * @param spacing the suggested track azimuthal spacing
   */
  void TrackGenerator::setDesiredAzimSpacing(double spacing)
  {
    if (spacing < 0)
      log::ferror("Unable to set a negative track azimuthal spacing "
                  "%f for the TrackGenerator.",
                  spacing);

    _azim_spacing = spacing;
    resetStatus();
  }

  /**
   * @brief Set a pointer to the Geometry to use for track generation.
   * @param geometry a pointer to the Geometry
   */
  void TrackGenerator::setGeometry(Geometry *geometry)
  {
    _geometry = geometry;
    _x_min = geometry->getMinX();
    _y_min = geometry->getMinY();
    _z_min = geometry->getMinZ();
    _x_max = geometry->getMaxX();
    _y_max = geometry->getMaxY();
    _z_max = geometry->getMaxZ();
    resetStatus();

#ifdef USTB_
    // Set BC for cyclic tracing
    if (_geometry->getMaxXBoundaryType() == REFLECTIVE ||
        _geometry->getMinXBoundaryType() == REFLECTIVE ||
        _geometry->getMaxYBoundaryType() == REFLECTIVE ||
        _geometry->getMinYBoundaryType() == REFLECTIVE)
    {
      log::info("Boundary condition (XY) for cyclic tracks = REFLECTIVE");
      _reflective = true;
    }
    else
      _reflective = false;

#elif MIT_
    _reflective = false;
#endif

    if (_geometry->getMaxZBoundaryType() == REFLECTIVE ||
        _geometry->getMinZBoundaryType() == REFLECTIVE)
    {
      log::info("Boundary condition (Z) for cyclic tracks = REFLECTIVE");
      _z_reflective = true;
    }
    else
      _z_reflective = false;
  }

  /**
   * @brief Sets the z-coord of the radial plane used in 2D calculations.
   * @param z_coord the z-coord of the radial plane
   */
  void TrackGenerator::setZCoord(double z_coord)
  {
    _z_coord = z_coord;

    /* Move the CMFD lattice near the plane of interest */
    if (_geometry->getCmfd() != NULL)
      _geometry->getCmfd()->getLattice()->getOffset()->setZ(z_coord - 0.5);
  }

  /**
   * @brief Sets the Quadrature used for integrating the MOC equations.
   * @param quadrature a pointer to the Quadrature object used in calculation
   */
  void TrackGenerator::setQuadrature(QuadraturePtr quadrature)
  {
    _quadrature = quadrature;
  }

  /**
   * @brief Returns whether or not the TrackGenerator contains Tracks
   *        for its current number of azimuthal angles, track spacing and
   *        geometry.
   * @return true if the TrackGenerator conatains Tracks; false otherwise
   */
  bool TrackGenerator::containsTracks()
  {
    return _contains_2D_tracks;
  }

  /**
   * @brief Returns whether or not the TrackGenerator contains segments
   *        for its current number of azimuthal angles, track spacing and
   *        geometry for it's current segmentation type.
   * @return true if the TrackGenerator conatains segments; false otherwise
   */
  bool TrackGenerator::containsSegments()
  {
    return _contains_2D_segments;
  }

  /**
   * @brief Fills an array with the x,y coordinates for each Track.
   * @details This class method is intended to be called by the OpenMOC
   *          Python "plotter" module as a utility to assist in plotting
   *          tracks. Although this method appears to require two arguments,
   *          in reality it only requires one due to SWIG and would be called
   *          from within Python as follows:
   *
   * @code
   *          num_tracks = track_generator.getNumTracks()
   *          coords = track_generator.retrieveTrackCoords(num_tracks*4)
   * @endcode
   *
   * @param coords an array of coords of length 4 times the number of Tracks
   * @param num_tracks the total number of Tracks
   */
  void TrackGenerator::retrieveTrackCoords(double *coords, long num_tracks)
  {
    retrieve2DTrackCoords(coords, num_tracks);
  }

  /**
   * @brief Fills an array with the x,y coordinates for each Track.
   * @details This class method is intended to be called by the OpenMOC
   *          Python "plotter" module as a utility to assist in plotting
   *          tracks. Although this method appears to require two arguments,
   *          in reality it only requires one due to SWIG and would be called
   *          from within Python as follows:
   *
   * @code
   *          num_tracks = track_generator.getNum2DTracks()
   *          coords = track_generator.retrieve2DTrackCoords(num_tracks*4)
   * @endcode
   *
   * @param coords an array of coords of length 4 times the number of Tracks
   * @param num_tracks the total number of Tracks
   */
  void TrackGenerator::retrieve2DTrackCoords(double *coords, long num_tracks)
  {

    if (num_tracks != NUM_VALUES_PER_RETRIEVED_TRACK * getNum2DTracks())
      log::ferror("Unable to retrieve the Track coordinates since the "
                  "TrackGenerator contains %d Tracks with %d coordinates but an "
                  "array of length %d was input",
                  getNum2DTracks(), NUM_VALUES_PER_RETRIEVED_TRACK * getNum2DTracks(), num_tracks);

    /* Fill the array of coordinates with the Track start and end points */
    int counter = 0;
    for (int a = 0; a < _num_azim / 2; a++)
    {
      for (long i = 0; i < _num_x[a] + _num_y[a]; i++)
      {
        coords[counter] = _tracks_2D[a][i].getStart()->getX();
        coords[counter + 1] = _tracks_2D[a][i].getStart()->getY();
        coords[counter + 2] = _tracks_2D[a][i].getStart()->getZ();
        coords[counter + 3] = _tracks_2D[a][i].getEnd()->getX();
        coords[counter + 4] = _tracks_2D[a][i].getEnd()->getY();
        coords[counter + 5] = _tracks_2D[a][i].getEnd()->getZ();

        counter += NUM_VALUES_PER_RETRIEVED_TRACK;
      }
    }
  }

  /**
   * @brief Fills an array with the x,y coordinates for each Track segment.
   * @details This class method is intended to be called by the OpenMOC
   *          Python "plotter" module as a utility to assist in plotting
   *          segments. Although this method appears to require two arguments,
   *          in reality it only requires one due to SWIG and would be called
   *          from within Python as follows:
   *
   * @code
   *          num_segments = track_generator.getNumSegments()
   *          coords = track_generator.retrieveSegmentCoords(num_segments*5)
   * @endcode
   *
   * @param coords an array of coords of length 5 times the number of segments
   * @param num_segments the total number of Track segments
   */
  void TrackGenerator::retrieveSegmentCoords(double *coords, long num_segments)
  {
    retrieve2DSegmentCoords(coords, num_segments);
  }

  /**
   * @brief Fills an array with the x,y coordinates for each Track segment.
   * @details This class method is intended to be called by the OpenMOC
   *          Python "plotter" module as a utility to assist in plotting
   *          segments. Although this method appears to require two arguments,
   *          in reality it only requires one due to SWIG and would be called
   *          from within Python as follows:
   *
   * @code
   *          num_segments = track_generator.getNum2DSegments()
   *          coords = track_generator.retrieve2DSegmentCoords(num_segments*5)
   * @endcode
   *
   * @param coords an array of coords of length 5 times the number of segments
   * @param num_segments the total number of Track segments
   */
  void TrackGenerator::retrieve2DSegmentCoords(double *coords, long num_segments)
  {

    if (num_segments != NUM_VALUES_PER_RETRIEVED_SEGMENT * getNum2DSegments())
      log::ferror("Unable to retrieve the Track segment coordinates since "
                  "the TrackGenerator contains %d segments with %d coordinates "
                  "but an array of length %d was input",
                  getNum2DSegments(), NUM_VALUES_PER_RETRIEVED_SEGMENT * getNum2DSegments(), num_segments);

    segment *curr_segment = NULL;
    double x0, x1, y0, y1, z0, z1;
    double phi;
    segment *segments;

    int counter = 0;

    /* Loop over Track segments and populate array with their FSR ID and *
     * start/end points */
    for (int a = 0; a < _num_azim / 2; a++)
    {
      for (long i = 0; i < _num_x[a] + _num_y[a]; i++)
      {

        x0 = _tracks_2D[a][i].getStart()->getX();
        y0 = _tracks_2D[a][i].getStart()->getY();
        z0 = _tracks_2D[a][i].getStart()->getZ();
        phi = _tracks_2D[a][i].getPhi();

        segments = _tracks_2D[a][i].getSegments();

        for (int s = 0; s < _tracks_2D[a][i].getNumSegments(); s++)
        {
          curr_segment = &segments[s];

          coords[counter] = curr_segment->_region_id;

          coords[counter + 1] = x0;
          coords[counter + 2] = y0;
          coords[counter + 3] = z0;

          x1 = x0 + cos(phi) * curr_segment->_length;
          y1 = y0 + sin(phi) * curr_segment->_length;
          z1 = z0;

          coords[counter + 4] = x1;
          coords[counter + 5] = y1;
          coords[counter + 6] = z1;

          x0 = x1;
          y0 = y1;
          z0 = z1;

          counter += NUM_VALUES_PER_RETRIEVED_SEGMENT;
        }
      }
    }
  }

  /**
   * @brief Checks the boundary conditions for all 2D surfaces for inconsistent
   *        periodic boundary conditions
   */
  void TrackGenerator::checkBoundaryConditions()
  {

    /* Extract the X and Y boundaries for whole Geometry */
    Universe *root_universe = _geometry->getRootUniverse();
    boundaryType min_x_bound = root_universe->getMinXBoundaryType();
    boundaryType max_x_bound = root_universe->getMaxXBoundaryType();
    boundaryType min_y_bound = root_universe->getMinYBoundaryType();
    boundaryType max_y_bound = root_universe->getMaxYBoundaryType();

    /* Check X and Y boundaries for consistency */
    /* 边界条件一致指的是，周期性边界条件要成对出现。换句话说，成对的边界要么都是
     * 周期性边界，要么都不是。这一点无论对四边形还是六边形都是成立的。
     */
    if ((min_x_bound == PERIODIC && max_x_bound != PERIODIC) ||
        (min_x_bound != PERIODIC && max_x_bound == PERIODIC))
      log::ferror("Cannot create tracks with only one x boundary"
                  " set to PERIODIC");

    else if ((min_y_bound == PERIODIC && max_y_bound != PERIODIC) ||
             (min_y_bound != PERIODIC && max_y_bound == PERIODIC))
      log::ferror("Cannot create tracks with only one y boundary"
                  " set to PERIODIC");

    /* Check that there are no periodic boundaries if domain decomposed */
    if (mpi::isSpatialDecomposed())
      if (min_x_bound == PERIODIC || min_y_bound == PERIODIC)
        log::ferror("Periodic boundaries are not supported for domain "
                    "decomposition");

    /* Check for correct track method if a PERIODIC bc is present */
    if (mpi::isSpatialDecomposed() || min_x_bound == PERIODIC ||
        min_y_bound == PERIODIC)
      _periodic = true;
    else
      _periodic = false;
  }

  /**
   * @brief Generates tracks for some number of azimuthal angles and track spacing
   * @details Computes the effective angles and track spacing. Computes the
   *          number of Tracks for each azimuthal angle, allocates memory for
   *          all Tracks at each angle and sets each Track's starting and ending
   *          Points, azimuthal angle, and azimuthal angle quadrature weight.
   *
   * @brief 为一定数量的方位角和轨迹间距生成轨迹
   * @details 计算有效角度和轨迹间距。计算每个方位角的轨迹数量，为每个角度的所有轨迹分配内存，
   *          并设置每个轨迹的起点和终点、方位角以及方位角求积权重。
   */
  void TrackGenerator::generateTracks()
  {

    /* Start recording track generation time */
    /* 开始记录轨迹生成时间 */
    _timer->startTimer();

    /* Check for valid quadrature */
    /* 检查求积组是否有效 */
    if (_quadrature != nullptr)
    {
      if (_quadrature->getNumAzimAngles() != (size_t)_num_azim)
      {
        _quadrature = nullptr;
      }
    }

    /* Generate Tracks, perform ray tracing across the geometry, and store
     * the data to a Track file */
    /* 生成轨迹，执行几何体的射线追踪，并将数据存储到轨迹文件中 */
    try
    {

      _timer->startTimer();
      /* Create default quadrature set if user one has not been set */
      /* 如果用户未设置求积组，则创建默认求积组 */
      if (_quadrature == NULL)
        initializeDefaultQuadrature();

      /* Initialize the quadrature set */
      /* 初始化求积组 */
      _quadrature->initialize();
      _timer->stopTimer("Initialize Quadrature");

      /* Check periodic BCs for symmetry */
      /* 检查边界条件 */
      checkBoundaryConditions();

      /* Lay down Tracks accross the Geometry */
      /* 在几何体上布置轨迹 */
      if (_geometry == NULL)
        log::ferror("Unable to lay down Tracks since no Geometry "
                    "has been set for the TrackGenerator");
      /* 无法布置轨迹，因为没有为轨迹生成器设置几何体 */

      /* Initialize the Tracks */
      /* 初始化轨迹 */
      _timer->startTimer();
      // 调用TrackGenerator3D中的initializeTracks()来生成2D和3D轨迹
      // 这一部分的3D轨迹链和轨迹堆没搞懂，后续对着文档再过一遍
      initializeTracks();
      _timer->stopTimer("Initialize Tracks");

      /* Initialize the track file directory and read in tracks if they exist */
      /* 初始化轨迹文件目录，如果轨迹文件存在则读取 */
      // FIXME initializeTrackFileDirectory();

      /* If track file not present, generate segments */
      /* 如果轨迹文件不存在，则生成轨迹段 */
      if (_use_input_file == false)
      {

        /* Segmentize the tracks */
        _timer->startTimer();
        // 调用TrackGenerator3D中的segmentize()来进行3D轨迹分段
        segmentize(); // CMFD只与轨迹分段相关，与轨迹生成无关，CMFD主要在3D轨迹分段中使用
        _timer->stopTimer("Segmenting");
        // FIXME HERE dumpSegmentsToFile();
      }

      /* Allocate array of mutex locks for each FSR
      在 MOC 求解过程中，计算是高度并行的：
并行维度：程序会同时在多个 CPU 核心上处理不同的特征线（Tracks）。
共享资源：虽然每条线是独立的，但它们都会穿过同一个几何体，进而穿过相同的 FSR（平源区）。
冲突场景：
当一条特征线穿过某个 FSR 时，它会计算该 FSR 的中子通量贡献，并试图累加到该 FSR 的总通量变量中。
如果线程 A 正在处理穿过 FSR #10 的特征线，同时线程 B 也在处理另一条穿过 FSR #10 的特征线。
如果两个线程同时尝试修改 FSR #10 的通量值（flux += delta），就会发生写冲突，导致数据错误。

_FSR_locks 数组为每一个 FSR 提供了一把独立的“锁”。
加锁：当某个线程想要更新 FSR #i 的数据时，它必须先拿到 _FSR_locks[i]。
互斥：如果锁已经被别的线程拿走了，当前线程就必须等待（阻塞），直到锁被释放。
解锁：更新完成后，线程释放锁，其他线程才能继续操作该 FSR。

3. 为什么是“细粒度”锁？
代码中是为每个 FSR 分配一个锁（new omp_lock_t[num_FSRs]），而不是用一个全局大锁。
全局锁：如果所有 FSR 共用一把锁，那么线程 A 更新 FSR #1 时，线程 B 连 FSR #2 都不能更新，这会严重降低并行效率，变成串行程序。
细粒度锁（Fine-grained Locking）：每个 FSR 有自己的锁。线程 A 更新 FSR #1 时，完全不影响线程 B 更新 FSR #2。只有当它们真的撞在同一个 FSR 上时才需要排队。这最大化了并行效率。
      */
      /* 为每个FSR（平源区）分配互斥锁数组 */
      long num_FSRs = _geometry->getNumFSRs();
      _FSR_locks = new omp_lock_t[num_FSRs];

/* Loop over all FSRs to initialize OpenMP locks */
/* 遍历所有FSR以初始化OpenMP锁 */
#pragma omp parallel for schedule(guided)
      for (long r = 0; r < num_FSRs; r++)
        omp_init_lock(&_FSR_locks[r]);

      /* Precompute the quadrature weights */
      /* 预计算求积权重 */
      _quadrature->precomputeWeights(_segment_formation != +segmentationType::EXPLICIT_2D); // 计算权重，传入true
    }
    catch (std::exception &e)
    {
      log::error("Unable to allocate memory needed to generate "
                 "Tracks. Backtrace:\n{}",
                 e.what());
    }

#ifdef ENABLE_MPI_
    mpi::mpiBarrier();
    // MPI同步点：确保所有进程在继续之前完成轨迹生成
#endif

    /* Stop recording track generation time and print */
    /* 停止记录轨迹生成时间并打印 */
    _timer->stopTimer("Track Generation Time");
    // 记 generateTracks() 花费的时间，方便做性能分析。

    printTimerReport();
    // 打印计时器报告，显示各个步骤的耗时情况
  }

  void TrackGenerator::printTimerReport()
  {

#ifdef ENABLE_MPI_
    _timer->reduceTimer(mpi::getMPIComm());
#endif

    log::header("Timing Report - 2D Ray Tracing (average)");

    const std::string tot_string = "Track Generation Time";
    _timer->printSplit(tot_string, "Total Track Generation & Segmentation Time");

    _timer->printSplit("Initialize Quadrature", "Quadrature Initialization", 1, tot_string);

    const std::string init_track_string = "Initialize Tracks";
    _timer->printSplit(init_track_string, "Tracks Initialization", 1, tot_string);

#ifdef ENABLE_MPI_
    if (mpi::isPrdTrackDecomposed())
    {
      _timer->printSplit("2D Track Mapping", "Time to mapping 2D tracks",
                         2, init_track_string);
    }
#endif

    const std::string segment_string = "Segmenting";
    _timer->printSplit(segment_string, "Segmenting Tracks", 1, tot_string);

#ifdef ENABLE_MPI_
    if (mpi::isPrdTrackDecomposed())
    {
      _timer->printSplit("Pre-segmenting", "Time to create FSRs", 2, segment_string);
    }
#endif

    _timer->printSplit("2D Segmenting", "2D Tracks Segmenting", 2, segment_string);

    _timer->printSplit("Initialize 2D FSRs", "2D FSRs Initialization", 2, segment_string);

    /* Show the statistics of objects */
    printObjectStatistics();
  }

  void TrackGenerator::printObjectStatistics()
  {

    auto part_length = log::get_line_length() / 3;

    long n1 = getNum2DChains();
    long n2 = getNum2DTracks();
    long n3 = _geometry->getNumTotalFSRs();

    log::separator("-");

    log::result("{1: ^{0}}{2: ^{0}}{3: ^{0}}",
                part_length, "# 2D chains", "# 2D tracks", "# FSRs");

    log::separator("-");

    log::result("{1: ^{0}}{2: ^{0}}{3: ^{0}}",
                part_length, n1, n2, n3);

    log::separator("-");
  }

  /**
   * @brief Allocates a new Quadrature with the default Quadrature
   * @details The defualt quadrature for 2D calculations is the TY quadrature
   */
  void TrackGenerator::initializeDefaultQuadrature()
  {
    _quadrature = std::make_shared<TYPolarQuad>();
    _quadrature->setNumAzimAngles(_num_azim);
    _quadrature->setNumPolarAngles(6);
  }

  /**
   * @brief calcualtes the least common multiple of two numbers a and b
   * @param first number a
   * @param second nuber b (order does not matter)
   * @return the least common multiple of a and b
   */
  double TrackGenerator::leastCommonMultiple(double a, double b)
  {

    bool _found = false;
    int lcm_a = 1;
    int lcm_b;
    double residual;

    /* For efficiency, make a the longer length */
    if (a < b)
    {
      double a_temp = a;
      a = b;
      b = a_temp;
    }

    while (!_found)
    {

      lcm_b = (int)round((lcm_a * a) / b);
      residual = fabs(lcm_a * a - lcm_b * b);

      if (residual < LCM_TOLERANCE)
        _found = true;
      else
        lcm_a++;
    }

    return lcm_a * a;
  }

  /**
   * @brief Returns the type of ray tracing used for segment formation
   * @return the segmentation type
   */
  segmentationType TrackGenerator::getSegmentFormation()
  {
    return _segment_formation;
  }

  /**
   * @brief Initializes Track azimuthal angles, start and end Points.
   * @details This method computes the azimuthal angles and effective track
   *          spacing to use to guarantee cyclic Track wrapping. Based on the
   *          angles and spacing, the number of Tracks per angle and the start
   *          and end Points for each Track are computed.
   */
  void TrackGenerator::initializeTracks()
  {

    /* Make sure that the width and height of the Geometry are nonzero */
    if (_geometry->getWidthX() <= 0 || _geometry->getWidthY() <= 0)
      log::ferror("The total height and width of the Geometry must"
                  " be nonzero for Track generation. Create a CellFill which "
                  "is filled by the entire geometry and bounded by XPlanes "
                  "and YPlanes to enable the Geometry to determine the "
                  "total width and height of the model.");

    log::finfo("Initializing 2D tracks...");

    /* Allocate memory for arrays */
    _tracks_2D = new Track *[_num_azim / 2];
    _num_x = new long[_num_azim / 2];
    _num_y = new long[_num_azim / 2];
    _num_chains = new long[_num_azim / 2];
    _my_num_chains = _num_chains;
    _tracks_per_azim = new long[_num_azim / 2];
    _cyclic_length = new double[_num_azim / 2];
    _chain_num_offsets = new int[_num_azim / 2]();

    double *dx_eff = new double[_num_azim / 2];
    double *dy_eff = new double[_num_azim / 2];

    /* Step 1: Compute nx, ny, dx, dy, phi, length, etc. */
    correct2DCyclicParameters(dx_eff, dy_eff);

#ifdef ENABLE_MPI_
    /* Compute my number of 2D chains */
    if (mpi::isPrdTrackDecomposed())
      generate2DChainMapping();
#endif

    /* Step 2: Initialize the 2D Tracks array along track chains */
#ifdef USTB_
    initialize2DTrackChains(dx_eff, dy_eff);
#elif MIT_
    initialize2DTracks(dx_eff, dy_eff);
#endif

    /* Set the flag indicating 2D tracks have been generated */
    _contains_2D_tracks = true;

    /* Step 3: Initialize the 1D array of Tracks for all Tracks */
    initializeTracksArray();

#ifdef MIT_
    // In this diagram, chains are supposed to be initialized
    // after _tracks_2D_array
    initialize2DTrackChains();
#endif

    log::fresult("Total number of 2D Tracks = %ld", getNum2DTracks());

    delete[] dx_eff;
    delete[] dy_eff;
  }

  /**
   * @brief Correct parameters for 2D cyclic tracks
   * @details Compute nx and ny according to user-defined args.
   *          Correct phi, dx and dy, and compute the length of
   *          2D cyclics.
   * @param dx_eff dx_eff to be computed
   * @param dy_eff dy_eff to be computed
   */
  void TrackGenerator::correct2DCyclicParameters(double *dx_eff, double *dy_eff)
  {

    double phi;
    double width = _geometry->getWidthX();
    double height = _geometry->getWidthY();
    long num_tracks = 0;

    /* Determine number of cyclics */
    for (int a = 0; a < _num_azim / 4; a++)
    {

      /* Get the desired azimuthal angle */
      phi = _quadrature->getPhi(a);

      /* The number of intersections with x,y-axes */
      log::fdebug("azim is %f, num X modules is %d, num Y modules is %d", phi, _geometry->getNumXModules(), _geometry->getNumYModules());
      double module_width = width / _geometry->getNumXModules();
      double module_height = height / _geometry->getNumYModules();

      /* The number of intersections with x,y-axes */
      _num_x[a] = (long)(fabs(module_width / _azim_spacing * sin(phi))) + 1;
      _num_y[a] = (long)(fabs(module_height / _azim_spacing * cos(phi))) + 1;
      log::fdebug("The number of intersections with x = %ld,y = %ld", _num_x[a], _num_y[a]);

      /* Extend the number of intersections with x,y axes for modules */
      _num_x[a] *= _geometry->getNumXModules();
      _num_y[a] *= _geometry->getNumYModules();
      num_tracks += _num_x[a] + _num_y[a];
#ifdef USTB_
      _num_chains[a] = antmoc::gcd(_num_x[a], _num_y[a]);
#elif MIT_
      _num_chains[a] = _num_x[a];
#endif
    }

    /* Correct the number of cyclics */
    // if (mpi::getNumProcs() > 1) {
    //   correctNum2DCyclics();
    // }

    /* Determine angular quadrature and track spacing */
    long extra_tracks = -num_tracks;
    double max_length = 0;
    double min_length = std::numeric_limits<double>::infinity();
    for (int a = 0; a < _num_azim / 4; a++)
    {
      /* Save number of intersections for supplementary angles */
      _num_x[_num_azim / 2 - a - 1] = _num_x[a];
      _num_y[_num_azim / 2 - a - 1] = _num_y[a];
      _num_chains[_num_azim / 2 - a - 1] = _num_chains[a];

      /* Effective/actual angle (not the angle we desire, but close) */
      phi = atan((height * _num_x[a]) / (width * _num_y[a]));
      _quadrature->setPhi(phi, a);

      /* Effective Track spacing (not spacing we desire, but close) */
      dx_eff[a] = (width / _num_x[a]);
      dy_eff[a] = (height / _num_y[a]);
      double azim_spacing = dx_eff[a] * sin(phi); // 修正后的d
      _quadrature->setAzimSpacing(azim_spacing, a);

      /* Save spacings for supplementary angles */
      dx_eff[_num_azim / 2 - a - 1] = dx_eff[a];
      dy_eff[_num_azim / 2 - a - 1] = dy_eff[a];

      /* Length of periodic tracks */
      _cyclic_length[a] = (_num_x[a] * _num_y[a]) / _num_chains[a] *
                          dx_eff[a] / cos(phi);

      if (isReflectiveCyclic())
      {
        _cyclic_length[a] *= 2;
      }
      _cyclic_length[_num_azim / 2 - a - 1] = _cyclic_length[a];

      extra_tracks += _num_x[a] + _num_y[a];

      /* Compute the max length of cyclics */
      max_length = std::max(max_length, _cyclic_length[a]);
      min_length = std::min(min_length, _cyclic_length[a]);
    }

    // Check after correction
    log::finfo("SD of phis after correction  = %.4f",
               computeAzimUniformity());

    log::finfo("Extra tracks due to correction = %ld, %.4f%%",
               extra_tracks, 100.0 * extra_tracks / num_tracks);

    log::finfo("Max/Min length of cyclics = %.2f, %.2f",
               max_length, min_length);

#ifdef ENABLE_MPI_
    if (mpi::isPrdTrackDecomposed())
    {
      long nc = getNum2DChains();
      int np = mpi::getNumProcsSharedDomain();
      if (nc < np && mpi::isRootSharedDomain())
      {
        log::fwarn("Too few cyclics: %d cyclics for %d processes",
                   nc, np);
      }
    }
#endif
  }

  /*
   * @brief Correct the number of 2D cyclics for load balancing
   */
  void TrackGenerator::correctNum2DCyclics()
  {

    double expect_div = 5.0;              // expected divisor
    double ulimit = 1 + 0.5 / expect_div; // upper limit of skipped ratio
    double llimit = 1 / ulimit;           // lower limit of skipped ratio
    long factor;

    for (int a = 0; a < _num_azim / 4; a++)
    {
      /* Skip special cases */
      double ratio = 1.0 * _num_y[a] / _num_x[a];
      if (_num_chains[a] == _num_x[a] ||
          _num_chains[a] == _num_y[a] ||
          (ratio < ulimit && ratio > llimit))
      {
        continue;
      }

      /* Correct nx and ny to balance the load */
      long *d1, *d2;
      d1 = &_num_x[a];
      d2 = &_num_y[a];

      if (ratio < 1.0)
      {
        d1 = &_num_y[a];
        d2 = &_num_x[a];
        ratio = 1.0 / ratio;
      }

      // Compute nxy
      double div = expect_div;
      if (*d1 < 2 * div)
      { // *d1/div <= 1
        div = 1;
      }
      factor = std::round(*d1 / div);         // determine the factor
      *d1 = div * factor;                     // recompute nx or ny
      *d2 = std::round(div * ratio) * factor; // recompute nx or ny

      /* Save the number of cyclic tracks */
      _num_chains[a] = antmoc::gcd(_num_x[a], _num_y[a]);
    }
  }

#ifdef USTB_

  //--------------------------------------------------------------------
  // Methods for cyclic track generation of ANT-MOC
  //--------------------------------------------------------------------

  /// \brief Initializes 2D Tracks array along cyclic tracks
  /// \details This method creates an array of 2D Tracks ordered by azimuthal
  ///          angle, x index, and link index. It then loops over chain links
  ///          to initialize them. If tracks are decomposed, their global ids
  ///          will be mapped to local ones.
  /// \param dx_eff the array of corrected dx
  /// \param dy_eff the array of corrected dy
  void TrackGenerator::initialize2DTrackChains(double *dx_eff, double *dy_eff)
  {

    log::finfo("Initializing 2D tracks along cyclic tracks...");

    // Allocate memory for the 2D tracks array
    for (int a = 0; a < _num_azim / 2; a++)
    {
      _tracks_2D[a] = new Track[getMyNum2DTracks(a)];
    }

    // Pointer to the current track
    Track *track;

    // Reflective and periodic tracks have different numbers
    int num_azim = getNumChainAzims();

    _tracks_2D_chains = new Track ***[num_azim];
    // Loop over azims of chains
    for (int a = 0; a < num_azim; a++)
    {
      long num_chains = getMyNum2DChains(a);
      long num_links = getNum2DLinks(a);
      _tracks_2D_chains[a] = new Track **[num_chains];

      // The local id of the current track
      long local_id = 0;

      // Loop over track chains
      for (long x = 0; x < num_chains; x++)
      {
        // The actual azim index assigned to tracks
        int ai = a;
        // Cyclic linking track ids
        long next_id = 0;
        long cur_gid = map2DChainId2Gid(a, x);
        long prev_gid = -1;

        /* Allocate memory for the track chains */
        _tracks_2D_chains[a][x] = new Track *[num_links];

        /* Cycle through the tracks in the 2D chain */
        for (long link_index = 0; link_index < num_links; link_index++)
        {
          /* Initialize the mapping */
          if (isReflectiveCyclic())
            next_id = local_id++ / 2;
          else
            next_id = local_id++;

          // Get the next link
          track = &_tracks_2D[ai][next_id];

          /* Initialize the current track and its reflections */
          double phi = _quadrature->getPhi(ai);
          track->setPhi(phi);
          track->setAzimIndex(ai);
          track->setXYIndex(next_id); // Set local index

          // Set track endpoints and reflections according to track gid
          set2DTrackEndpoints(track, ai, cur_gid, dx_eff, dy_eff);
          TrackGenerator::initializeTrackBCs(track, ai, cur_gid);
          TrackGenerator::initializeTrackReflections(track, ai, cur_gid);

          // Set linking info
          track->setLinkIndex(link_index);
          track->setChainIndex(x);

          // Set the current chain link
          _tracks_2D_chains[a][x][link_index] = track;

          /* Get the uid of the next link along the cyclic track */
          if (isReflectiveCyclic())
          {
            // Find the correct uid of the next track
            long next_gid = track->getTrackReflFwd();
            // If the direction is wrong, the track should be reversed
            if (next_gid == prev_gid)
            {
              next_gid = track->getTrackReflBwd();
              reverseReflectiveTracks(track);
            }
            ai = _num_azim / 2 - ai - 1; // complement
            prev_gid = cur_gid;
            cur_gid = next_gid;
          }
          else
          {
            cur_gid = track->getTrackPrdcFwd();
          }
        } // end for links
      } // end for chains
    } // end for azims

    log::fdebug("Finish initialization of 2D track chains array...");

    // Statistics
    log::fresult("Total number of 2D Chains = %ld", getNum2DChains());
#ifdef ENABLE_MPI_
    log::fverbose("Number of local 2D Chains = %ld", getMyNum2DChains());
#endif
  }

  /// \brief Initializes 2D Track reflections
  /// \details This method computes the connecting Tracks for a 2D Track in
  ///          the TrackGenerator analytically, handling both reflective and
  ///          periodic boundaries.
  ///          Track IDs used in this method are local IDs.
  void
  TrackGenerator::initializeTrackReflections(Track *track, int a, long i)
  {

    /* Set connecting tracks in forward direction */
    if (i < _num_y[a])
    { // 从右边/左边出射的轨迹
      track->setNextFwdFwd(true);
      track->setTrackPrdcFwd(i + _num_x[a]); // 周期边界，相接的轨迹的id
      track->setTrackReflFwd(i + _num_x[a]); // 反射边界，相接的轨迹的id

      if (!isReflectiveCyclic())
        track->setTrackNextFwd(i + _num_x[a]);
      else
        track->setTrackNextFwd(i + _num_x[a]);
    }
    else
    { // 从顶边出射的轨迹
      track->setTrackPrdcFwd(i - _num_y[a]);
      track->setTrackReflFwd((_num_x[a] + _num_y[a]) - (i - _num_y[a]) - 1);

      if (!isReflectiveCyclic())
      {
        track->setNextFwdFwd(true);
        track->setTrackNextFwd(i - _num_y[a]);
      }
      else
      {
        track->setNextFwdFwd(false);
        track->setTrackNextFwd((_num_x[a] + _num_y[a]) - (i - _num_y[a]) - 1);
      }
    }

    /* Set connecting tracks in backward direction */
    if (i < _num_x[a])
    { // 从底边出射的轨迹
      track->setTrackPrdcBwd(i + _num_y[a]);
      track->setTrackReflBwd(_num_x[a] - i - 1);

      if (!isReflectiveCyclic())
      {
        track->setNextBwdFwd(false);
        track->setTrackNextBwd(i + _num_y[a]);
      }
      else
      {
        track->setTrackNextBwd(_num_x[a] - i - 1);
        track->setNextBwdFwd(true);
      }
    }
    else
    { // 从左边/右边出射的轨迹
      track->setNextBwdFwd(false);
      track->setTrackPrdcBwd(i - _num_x[a]);
      track->setTrackReflBwd(i - _num_x[a]);

      if (!isReflectiveCyclic())
        track->setTrackNextBwd(i - _num_x[a]);
      else
        track->setTrackNextBwd(i - _num_x[a]);
    }
  }

  /// \brief Reverse the 'forward' and 'backward' direction
  /// \details A track in 'forward' direction will be traced along the
  ///          cyclic track it is on. For periodic tracks, this method
  ///          should change nothing. For reflective tracks, some of the
  ///          links will be reverse to subject to the direction of cyclic
  ///          tracks.
  void TrackGenerator::reverseReflectiveTracks(Track *track)
  {
    track->setReversed(true);

    // Reverse endpoints and direction
    Point *start = track->getStart();
    Point *end = track->getEnd();
    track->setValues(end->getX(), end->getY(),
                     start->getX(), start->getY(),
                     track->getPhi() + M_PI);

    // Reverse linking info
    auto bc = track->getBCFwd();
    track->setBCFwd(track->getBCBwd());
    track->setBCBwd(bc);
    auto surf = track->getSurfaceOut();
    track->setSurfaceOut(track->getSurfaceIn());
    track->setSurfaceIn(surf);
    auto nextfwd = track->getNextFwdFwd();
    track->setNextFwdFwd(track->getNextBwdFwd());
    track->setNextBwdFwd(nextfwd);
    auto prdc = track->getTrackPrdcFwd();
    track->setTrackPrdcFwd(track->getTrackPrdcBwd());
    track->setTrackPrdcBwd(prdc);
    auto refl = track->getTrackReflFwd();
    track->setTrackReflFwd(track->getTrackReflBwd());
    track->setTrackReflBwd(refl);
    auto next = track->getTrackNextFwd();
    track->setTrackNextFwd(track->getTrackNextBwd());
    track->setTrackNextBwd(next);
  }

#elif MIT_

  //--------------------------------------------------------------------
  // Methods for OpenMOC track generation
  //--------------------------------------------------------------------

  /// \brief Initializes 2D Tracks array
  /// \details See OpenMOC
  /// \param dx_eff the array of corrected dx
  /// \param dy_eff the array of corrected dy
  void TrackGenerator::initialize2DTracks(double *dx_eff, double *dy_eff)
  {

    for (int a = 0; a < _num_azim / 2; a++)
    {
      // Allocate memory for the 2D tracks array
      _tracks_2D[a] = new Track[getNum2DTracks(a)];
      for (int i = 0; i < getNum2DTracks(a); i++)
      {

        /* Get track and set angle and track indices */
        Track *track = &_tracks_2D[a][i];
        track->setPhi(_quadrature->getPhi(a));
        track->setAzimIndex(a);
        track->setXYIndex(i);

        // Set the starting point and end point of the track
        set2DTrackEndpoints(track, a, i, dx_eff, dy_eff);

        // Initialize 2D track reflections
        TrackGenerator::initializeTrackBCs(track, a, i);
        TrackGenerator::initializeTrackReflections(track, a, i);
      }
    }
  }

  /// \brief Initializes 2D Track reflections
  /// \details See OpenMOC.
  ///          Track IDs used in this method are global IDs.
  void
  TrackGenerator::initializeTrackReflections(Track *track, int a, long i)
  {

    // Supplementary angle
    int ac = _num_azim / 2 - a - 1;

    /* Set connecting tracks in forward direction */
    if (i < _num_y[a])
    { // 从右边/左边出射的轨迹
      track->setNextFwdFwd(true);

      auto prdc_fwd = get2DTrackID(a, i + _num_x[a]);
      auto refl_fwd = get2DTrackID(ac, i + _num_x[a]);
      track->setTrackPrdcFwd(prdc_fwd);
      track->setTrackReflFwd(refl_fwd);

      if (track->getBCFwd() == PERIODIC ||
          track->getBCFwd() == INTERFACE)
        track->setTrackNextFwd(prdc_fwd);
      else
        track->setTrackNextFwd(refl_fwd);
    }
    else
    { // 从顶边出射的轨迹
      auto prdc_fwd = get2DTrackID(a, i - _num_y[a]);
      auto refl_fwd = get2DTrackID(ac, (_num_x[a] + _num_y[a]) - (i - _num_y[a]) - 1);
      track->setTrackPrdcFwd(prdc_fwd);
      track->setTrackReflFwd(refl_fwd);

      if (_geometry->getMaxYBoundaryType() == PERIODIC ||
          _geometry->getMaxYBoundaryType() == INTERFACE)
      {
        track->setNextFwdFwd(true);
        track->setTrackNextFwd(prdc_fwd);
      }
      else
      {
        track->setNextFwdFwd(false);
        track->setTrackNextFwd(refl_fwd);
      }
    }

    /* Set connecting tracks in backward direction */
    if (i < _num_x[a])
    { // 从底边出射的轨迹
      auto prdc_fwd = get2DTrackID(a, i + _num_y[a]);
      auto refl_fwd = get2DTrackID(ac, _num_x[a] - i - 1);
      track->setTrackPrdcBwd(prdc_fwd);
      track->setTrackReflBwd(refl_fwd);

      if (_geometry->getMinYBoundaryType() == PERIODIC ||
          _geometry->getMinYBoundaryType() == INTERFACE)
      {
        track->setNextBwdFwd(false);
        track->setTrackNextBwd(prdc_fwd);
      }
      else
      {
        track->setNextBwdFwd(true);
        track->setTrackNextBwd(refl_fwd);
      }
    }
    else
    { // 从左边/右边出射的轨迹
      track->setNextBwdFwd(false);

      auto prdc_fwd = get2DTrackID(a, i - _num_x[a]);
      auto refl_fwd = get2DTrackID(ac, i - _num_x[a]);
      track->setTrackPrdcBwd(prdc_fwd);
      track->setTrackReflBwd(refl_fwd);

      if (track->getBCBwd() == PERIODIC ||
          track->getBCBwd() == INTERFACE)
        track->setTrackNextBwd(prdc_fwd);
      else
        track->setTrackNextBwd(refl_fwd);
    }
  }

  /**
   * @brief Initializes 2D Track chains array.
   * @details This method creates an array of 2D Tracks ordered by azimuthal
   *          angle, x index, and link index. A track chain in OpenMOC is not
   *          a cyclic track but a part of it which is cut off by the YMAX
   *          surface.
   */
  void TrackGenerator::initialize2DTrackChains()
  {

    Track *track;
    int link_index;

    _tracks_2D_chains = new Track ***[_num_azim / 2];
    for (int a = 0; a < _num_azim / 2; a++)
    {
      long num_chains = getNum2DChains(a);
      _tracks_2D_chains[a] = new Track **[num_chains];
      for (long x = 0; x < num_chains; x++)
      {

        /* Get the first track in the 2D chain */
        link_index = 0;
        track = &_tracks_2D[a][x];
        track->setLinkIndex(link_index);

        /* Cycle through 2D chain's tracks, set their index, get chain length */
        while (track->getXYIndex() < _num_y[a])
        {
          link_index++;
          track = _tracks_2D_array[track->getTrackPrdcFwd()];
          track->setLinkIndex(link_index);
        }

        /* Allocate memory for the track chains */
        _tracks_2D_chains[a][x] = new Track *[link_index + 1];

        /* Assign tracks to the chains array */
        link_index = 0;
        track = &_tracks_2D[a][x];
        _tracks_2D_chains[a][x][link_index] = track;

        while (track->getXYIndex() < _num_y[a])
        {
          link_index++;
          track = _tracks_2D_array[track->getTrackPrdcFwd()];
          _tracks_2D_chains[a][x][link_index] = track;
        }
      }
    }
  }

#endif // USTB_ or MIT_

  /// \brief Set the starting point and end point of a 2D track
  /// \param track a pointer to the 2D track
  /// \param a azim index of the track
  /// \param i xy index of the track
  /// \param dx_eff an array of spacings
  /// \param dy_eff an array of spacings
  void
  TrackGenerator::set2DTrackEndpoints(Track *track, int a, long i,
                                      double *dx_eff, double *dy_eff)
  {
    if (track == nullptr)
      log::ferror("I can't forge a NULL pointer!");

    double width = _geometry->getWidthX();
    double height = _geometry->getWidthY();
    double x_min = _geometry->getMinX();
    double y_min = _geometry->getMinY();
    double dx = dx_eff[a];
    double dy = dy_eff[a];

    /* Set start point */
    if (a < _num_azim / 4)
    {                    // 锐角
      if (i < _num_x[a]) // 从底边入射的轨迹起点
        track->getStart()->setCoords(x_min + width - dx * (i + 0.5),
                                     y_min);
      else // 从左边入射的轨迹起点
        track->getStart()->setCoords(x_min, y_min + dy *
                                                        (i - _num_x[a] + 0.5));
    }
    else
    {                    // 钝角
      if (i < _num_x[a]) // 从底边入射的轨迹起点
        track->getStart()->setCoords(x_min + dx * (i + 0.5), y_min);
      else // 从右边入射的轨迹起点
        track->getStart()->setCoords(x_min + width, y_min + dy *
                                                                (i - _num_x[a] + 0.5));
    }

    /* Set end point */
    if (a < _num_azim / 4)
    {                    // 锐角
      if (i < _num_y[a]) // 从右边出射的轨迹终点
        track->getEnd()->setCoords(x_min + width, y_min + dy *
                                                              (i + 0.5));
      else // 从顶边出射的轨迹终点
        track->getEnd()->setCoords(x_min + width - dx * ((i - _num_y[a]) + 0.5), y_min + height);
    }
    else
    {                    // 钝角
      if (i < _num_y[a]) // 从左边出射的轨迹终点
        track->getEnd()->setCoords(x_min, y_min + dy * (i + 0.5));
      else // 从顶边出射的轨迹终点
        track->getEnd()->setCoords(x_min + dx * (i - _num_y[a] + 0.5),
                                   y_min + height);
    }
  }

  /// \brief Initializes boundary conditions for 2D tracks
  /// \details This method also sets domain connections for 2D tracks.
  void
  TrackGenerator::initializeTrackBCs(Track *track, int a, long i)
  {

    /* Set the foward boundary conditions */
    if (a < _num_azim / 4)
    { // 锐角
      if (i < _num_y[a])
      { // 从右边出射的轨迹
        track->setBCFwd(_geometry->getMaxXBoundaryType());
        track->setSurfaceOut(SURFACE_X_MAX);
#ifdef ENABLE_MPI_
        track->setDomainFwd(_geometry->getNeighborDomain(1, 0, 0));
#endif
      }
      else
      { // 从顶边出射的轨迹
        track->setBCFwd(_geometry->getMaxYBoundaryType());
        track->setSurfaceOut(SURFACE_Y_MAX);
#ifdef ENABLE_MPI_
        track->setDomainFwd(_geometry->getNeighborDomain(0, 1, 0));
#endif
      }

      if (i < _num_x[a])
      { // 从底边入射的轨迹
        track->setBCBwd(_geometry->getMinYBoundaryType());
        track->setSurfaceIn(SURFACE_Y_MIN);
#ifdef ENABLE_MPI_
        track->setDomainBwd(_geometry->getNeighborDomain(0, -1, 0));
#endif
      }
      else
      { // 从左边入射的轨迹
        track->setBCBwd(_geometry->getMinXBoundaryType());
        track->setSurfaceIn(SURFACE_X_MIN);
#ifdef ENABLE_MPI_
        track->setDomainBwd(_geometry->getNeighborDomain(-1, 0, 0));
#endif
      }
    }

    /* Set the backward boundary conditions */
    else
    { // 钝角
      if (i < _num_y[a])
      { // 从左边出射的轨迹
        track->setBCFwd(_geometry->getMinXBoundaryType());
        track->setSurfaceOut(SURFACE_X_MIN);
#ifdef ENABLE_MPI_
        track->setDomainFwd(_geometry->getNeighborDomain(-1, 0, 0));
#endif
      }
      else
      { // 从顶边出射的轨迹
        track->setBCFwd(_geometry->getMaxYBoundaryType());
        track->setSurfaceOut(SURFACE_Y_MAX);
#ifdef ENABLE_MPI_
        track->setDomainFwd(_geometry->getNeighborDomain(0, 1, 0));
#endif
      }

      if (i < _num_x[a])
      { // 从底边入射的轨迹
        track->setBCBwd(_geometry->getMinYBoundaryType());
        track->setSurfaceIn(SURFACE_Y_MIN);
#ifdef ENABLE_MPI_
        track->setDomainBwd(_geometry->getNeighborDomain(0, -1, 0));
#endif
      }
      else
      { // 从右边入射的轨迹
        track->setBCBwd(_geometry->getMaxXBoundaryType());
        track->setSurfaceIn(SURFACE_X_MAX);
#ifdef ENABLE_MPI_
        track->setDomainBwd(_geometry->getNeighborDomain(1, 0, 0));
#endif
      }
    }
  }

  /**
   * @brief Generate segments for each Track across the Geometry.
   */
  {
    void TrackGenerator::segmentize()

        log::finfo("Ray tracing for 2D track segmentation...");

    /* Check to ensure the Geometry is infinite in axial direction */
    double max_z = _geometry->getGlobalMaxZ();
    double min_z = _geometry->getGlobalMinZ();
    if (max_z - min_z < FLT_INFINITY)
    {
      log::fwarn_once("The Geometry was set with non-infinite "
                      "z-boundaries and supplied to a 2D TrackGenerator. The min-z "
                      "boundary was set to %5.2f and the max-z boundary was set to "
                      "%5.2f. Z-boundaries are assumed to be infinite in 2D "
                      "TrackGenerators.",
                      min_z, max_z);

#ifdef MIT_
      if (mpi::isSpatialDecomposed())
      {
        /* Check that the geometry is not domain decomposed in Z */
        int domains_xyz[3];
        _geometry->getDomainStructure(domains_xyz);
        if (domains_xyz[2] > 1)
          log::ferror("A geometry with an axial domain domain decomposition "
                      "has been supplied to a 2D ray tracer.");
      }
#endif

      Cmfd *cmfd = _geometry->getCmfd();
      if (cmfd != NULL)
      {
        cmfd->setWidthZ(std::numeric_limits<double>::infinity());
        Point offset;
        offset.setX(cmfd->getLattice()->getOffset()->getX());
        offset.setY(cmfd->getLattice()->getOffset()->getY());
        offset.setZ(0.0);
        cmfd->initializeLattice(&offset);
      }
    }

#ifdef ENABLE_MPI_
    long num_fsrs = 0;
    if (mpi::isPrdTrackDecomposed())
    {
      // pre-segmenting for FSR generation
      log::finfo("Pre-segmenting tracks to generate 2D FSRs...");

      preSegmentize2D(_z_coord);

      // number of fsrs after pre-segmentation
      num_fsrs = _geometry->getNumFSRs();
      log::fverbose_once("Number of FSRs created before segmenting "
                         "= %ld",
                         num_fsrs);
    }
#endif

    /* Make sure CMFD lattice is initialized and has right offset for 2D */
    Cmfd *cmfd = _geometry->getCmfd();
    if (cmfd != NULL)
    {

      /* Check that CMFD has been initialized */
      if (cmfd->getLattice() == NULL)
        log::error("CMFD has not been initialized before generating tracks."
                   "A call to geometry.initializeFlatSourceRegions() may be missing.");

      /* Re-initialize CMFD lattice with 2D dimensions */
      Point offset;
      offset.setX(cmfd->getLattice()->getOffset()->getX());
      offset.setY(cmfd->getLattice()->getOffset()->getY());
      offset.setZ(_z_coord);
      cmfd->initializeLattice(&offset);
    }

    _timer->startTimer();

    /* Loop over all Tracks */
    long num_2D_tracks = getMyNum2DTracks();
    Progress progress(num_2D_tracks, "Segmenting 2D Tracks", 0.1, _geometry,
                      true);

    /* Loop over all Tracks */
#pragma omp parallel for schedule(dynamic)
    for (long uid = 0; uid < num_2D_tracks; uid++)
    {
      progress.incrementCounter();
      _geometry->segmentize2D(_tracks_2D_array[uid], _z_coord);
    }

    _timer->stopTimer("2D Segmenting");

#ifdef ENABLE_MPI_
    if (mpi::isPrdTrackDecomposed())
    {
      long tot_num_fsrs = _geometry->getNumFSRs();
      if (num_fsrs < tot_num_fsrs)
      {
        log::fwarn("Pre-segmentation failed, %ld out of %ld, please "
                   "increase the number or azimuthals or the spacing "
                   "of 2D tracks",
                   num_fsrs, tot_num_fsrs);
      }

      long global_fsrs;
      MPI_Allreduce(&tot_num_fsrs, &global_fsrs, 1, mpi::getDatatype<long>(), MPI_MAX,
                    mpi::getCommSharedDomain());

      if (tot_num_fsrs != global_fsrs)
      {
        log::ferror("The number of FSRs is inconsistant across processes: "
                    "Local = %ld, Global = %ld",
                    tot_num_fsrs, global_fsrs);
      }
    }
#endif

    /* Initialize FSRs and their associated vectors */
    _timer->startTimer();
    _geometry->initializeFSRVectors();
    _timer->stopTimer("Initialize 2D FSRs");

    _contains_2D_segments = true;
  }

  /**
   * @brief This method creates a directory to store Track files, and reads
   *        in ray tracing data for Tracks and segments from a Track file
   *        if one exists.
   * @details This method is called by the TrackGenerator::generateTracks()
   *          class method. If a Track file exists for this Geometry, number
   *          of azimuthal angles, and track spacing, then this method will
   *          import the ray tracing Track and segment data to fill the
   *          appropriate data structures.
   */
  void TrackGenerator::initializeTrackFileDirectory()
  {

    struct stat buffer;
    std::stringstream directory;

    /** Create directory to store Track files with pre-generated ray tracing data
     *  if the directory does not yet exist */
    directory << fileutils::getOutputDirectory() << "/tracks";
    struct stat st;
    if (!(stat(directory.str().c_str(), &st) == 0))
      mkdir(directory.str().c_str(), S_IRWXU);

    /* Check to see if a Track file exists for this geometry, number of azimuthal
     * angles, and track spacing, and if so, import the ray tracing data */
    _tracks_filename = getTestFilename(directory.str());
    if (!stat(_tracks_filename.c_str(), &buffer))
    {
      if (readSegmentsFromFile())
      {
        _use_input_file = true;
        setContainsSegments(true);
      }
    }
  }

  /**
   * @brief Returns the filename for writing tracking data
   */
  std::string TrackGenerator::getTestFilename(std::string directory)
  {

    std::stringstream test_filename;

    if (_geometry->getCmfd() != NULL)
      test_filename << directory << "/2D_"
                    << _num_azim << "_azim_"
                    << _azim_spacing << "_cm_spacing_cmfd_"
                    << _geometry->getCmfd()->getNumX()
                    << "x" << _geometry->getCmfd()->getNumY();
    else
      test_filename << directory << "/2D_"
                    << _num_azim << "_angles_"
                    << _azim_spacing << "_cm_spacing";

    test_filename << ".data";

    return test_filename.str();
  }

  /**
   * @brief Updates whether the TrackGenerator contains segments
   * @param contains_segments whether the TrackGenerator contains segments
   */
  void TrackGenerator::setContainsSegments(bool contains_segments)
  {
    _contains_2D_segments = contains_segments;
  }

  /**
   * @brief Writes all Track and segment data to a "*.tracks" binary file.
   * @details Storing Tracks in a binary file saves time by eliminating ray
   *          tracing for Track segmentation in commonly simulated geometries.
   */
  void TrackGenerator::dumpSegmentsToFile()
  {

    /* Check whether the segments should be dumped */
    if (!_dump_segments)
      return;

    log::finfo("Dumping segments to file...");

    if (!containsSegments())
      log::ferror("Unable to dump Segments to a file since no Segments "
                  "have been generated for %d azimuthal angles and %f track "
                  "spacing",
                  _num_azim, _azim_spacing);

    FILE *out;
    out = fopen(_tracks_filename.c_str(), "w");

    /* Get a string representation of the Geometry's attributes. This is used to
     * check whether or not ray tracing has been performed for this Geometry */
    std::string geometry_to_string = _geometry->toString();
    int string_length = geometry_to_string.length() + 1;

    /* Write geometry metadata to the Track file */
    fwrite(&string_length, sizeof(int), 1, out);
    fwrite(geometry_to_string.c_str(), sizeof(char) * string_length, 1, out);

    /* Write segment data to Track file */
    DumpSegments dump_segments(this);
    dump_segments.setOutputFile(out);
    if (_segment_formation == +segmentationType::EXPLICIT_2D ||
        _segment_formation == +segmentationType::EXPLICIT_3D)
      dump_segments.execute();

    /* Get FSR vector maps */
    auto &FSR_keys_map = _geometry->getFSRKeysMap();
    auto &FSRs_to_keys = _geometry->getFSRsToKeys();
    auto &FSRs_to_material_IDs = _geometry->getFSRsToMaterialIDs();

    std::string fsr_key;
    long fsr_id;
    double x, y, z;

    /* Write number of FSRs */
    long num_FSRs = _geometry->getNumFSRs();
    fwrite(&num_FSRs, sizeof(long), 1, out);

    /* Write FSR vector maps to file */
    auto fsr_key_list = FSR_keys_map.keys();
    auto fsr_data_list = FSR_keys_map.values();
    Cmfd *cmfd = _geometry->getCmfd();
    for (long i = 0; i < num_FSRs; i++)
    {

      /* Write data to file from FSR_keys_map */
      fsr_key = fsr_key_list[i];
      string_length = fsr_key.length() + 1;
      fwrite(&string_length, sizeof(int), 1, out);
      fwrite(fsr_key.c_str(), sizeof(char) * string_length, 1, out);

      fsr_id = fsr_data_list[i]->_fsr_id;
      x = fsr_data_list[i]->_point->getX();
      y = fsr_data_list[i]->_point->getY();
      z = fsr_data_list[i]->_point->getZ();
      fwrite(&fsr_id, sizeof(long), 1, out);
      fwrite(&x, sizeof(double), 1, out);
      fwrite(&y, sizeof(double), 1, out);
      fwrite(&z, sizeof(double), 1, out);

      /* Write data to file from FSRs_to_material_IDs */
      fwrite(&(FSRs_to_material_IDs.at(i)), sizeof(int), 1, out);

      /* Write data to file from FSRs_to_keys */
      fsr_key = FSRs_to_keys.at(i);
      string_length = fsr_key.length() + 1;
      fwrite(&string_length, sizeof(int), 1, out);
      fwrite(fsr_key.c_str(), sizeof(char) * string_length, 1, out);
    }

    /* Write cmfd_fsrs vector of vectors to file */
    if (cmfd != NULL)
    {
      std::vector<std::vector<long>> *cell_fsrs = cmfd->getCellFSRs();
      std::vector<long>::iterator iter;
      int num_cells = cmfd->getNumCells();
      fwrite(&num_cells, sizeof(int), 1, out);

      /* Loop over CMFD cells */
      for (int cell = 0; cell < num_cells; cell++)
      {
        int num_cell_FSRs = cell_fsrs->at(cell).size();
        fwrite(&num_cell_FSRs, sizeof(int), 1, out);

        /* Loop over FSRs within cell */
        for (iter = cell_fsrs->at(cell).begin();
             iter != cell_fsrs->at(cell).end(); ++iter)
          fwrite(&(*iter), sizeof(long), 1, out);
      }
    }

    /* Delete key and value lists */
    delete[] fsr_key_list;
    delete[] fsr_data_list;

    /* Write 2D basis information for 3D solvers */
    if (_segment_formation != +segmentationType::EXPLICIT_2D &&
        _segment_formation != +segmentationType::EXPLICIT_3D)
      writeExtrudedFSRInfo(out);

    /* Close the Track file */
    fclose(out);

    /* Inform other the TrackGenerator::generateTracks() method that it may
     * import ray tracing data from this file if it is called and the ray
     * tracing parameters have not changed */
    _use_input_file = true;
  }

  /**
   * @brief Write information of all Extruded FSRs to a file
   //TODO Use implementation in 3D track generator
   * @param out file to write to
   */
  void TrackGenerator::writeExtrudedFSRInfo(FILE *out) {}

  /**
   * @brief Reads Tracks in from a "*.tracks" binary file.
   * @details Storing Tracks in a binary file saves time by eliminating ray
   *          tracing for Track segmentation in commonly simulated geometries.
   * @return true if able to read Tracks in from a file; false otherwise
   */
  bool TrackGenerator::readSegmentsFromFile()
  {

    int ret;
    FILE *in = fopen(_tracks_filename.c_str(), "r");

    int string_length;

    /* Import Geometry metadata from the Track file */
    ret = _geometry->twiddleRead(&string_length, sizeof(int), 1, in);
    char *geometry_to_string = new char[string_length];
    ret = _geometry->twiddleRead(geometry_to_string, sizeof(char) * string_length, 1, in);

    /* Check if our Geometry is exactly the same as the Geometry in the
     * Track file for this number of azimuthal angles and track spacing */
    if (_geometry->toString().compare(std::string(geometry_to_string)) != 0)
      return false;

    delete[] geometry_to_string;

    log::finfo("Importing ray tracing data from file...");

    /* Load all segment data into Tracks */
    ReadSegments read_segments(this);
    read_segments.setInputFile(in);
    if (_segment_formation == +segmentationType::EXPLICIT_2D ||
        _segment_formation == +segmentationType::EXPLICIT_3D)
      read_segments.execute();

    /* Create FSR vector maps */
    auto &FSR_keys_map = _geometry->getFSRKeysMap();
    auto &FSRs_to_material_IDs = _geometry->getFSRsToMaterialIDs();
    auto &FSRs_to_keys = _geometry->getFSRsToKeys();
    auto &FSRs_to_centroids = _geometry->getFSRsToCentroids();
    auto &FSRs_to_CMFD_cells = _geometry->getFSRsToCMFDCells();

    long num_FSRs;
    std::string fsr_key;
    long fsr_key_id;
    double x, y, z;

    /* Get number of FSRs */
    ret = _geometry->twiddleRead(&num_FSRs, sizeof(long), 1, in);

    /* Resize vectors */
    FSRs_to_centroids.resize(num_FSRs);
    FSRs_to_CMFD_cells.resize(num_FSRs);

    /* Read FSR vector maps from file */
    for (long fsr_id = 0; fsr_id < num_FSRs; fsr_id++)
    {

      /* Read key for FSR_keys_map */
      ret = _geometry->twiddleRead(&string_length, sizeof(int), 1, in);
      char *char_buffer1 = new char[string_length];
      ret = _geometry->twiddleRead(char_buffer1, sizeof(char) * string_length, 1, in);
      fsr_key = std::string(char_buffer1);

      /* Read data from file for FSR_keys_map */
      ret = _geometry->twiddleRead(&fsr_key_id, sizeof(long), 1, in);
      ret = _geometry->twiddleRead(&x, sizeof(double), 1, in);
      ret = _geometry->twiddleRead(&y, sizeof(double), 1, in);
      ret = _geometry->twiddleRead(&z, sizeof(double), 1, in);
      auto fsr = new FSRData();
      fsr->_fsr_id = fsr_key_id;
      Point *point = new Point();
      point->setCoords(x, y, z);
      fsr->_point = point;
      FSR_keys_map.insert(fsr_key, fsr);

      /* Read data from file for FSR_to_materials_IDs */
      int material_id;
      ret = _geometry->twiddleRead(&material_id, sizeof(int), 1, in);
      FSRs_to_material_IDs.push_back(material_id);

      /* Read data from file for FSR_to_keys */
      ret = _geometry->twiddleRead(&string_length, sizeof(int), 1, in);
      char *char_buffer2 = new char[string_length];
      ret = _geometry->twiddleRead(char_buffer2, sizeof(char) * string_length, 1, in);
      fsr_key = std::string(char_buffer2);
      FSRs_to_keys.push_back(fsr_key);
    }

    /* Read cmfd cell_fsrs vector of vectors from file */
    Cmfd *cmfd = _geometry->getCmfd();
    if (cmfd != NULL)
    {
      std::vector<std::vector<long>> cell_fsrs;
      int num_cells;
      ret = _geometry->twiddleRead(&num_cells, sizeof(int), 1, in);
      long fsr_id;

      /* Loop over CMFD cells */
      for (int cell = 0; cell < num_cells; cell++)
      {
        std::vector<long> *fsrs = new std::vector<long>;
        cell_fsrs.push_back(*fsrs);
        int num_cell_FSRs;
        ret = _geometry->twiddleRead(&num_cell_FSRs, sizeof(int), 1, in);

        /* Loop over FSRs within cell */
        for (int fsr = 0; fsr < num_cell_FSRs; fsr++)
        {
          ret = _geometry->twiddleRead(&fsr_id, sizeof(long), 1, in);
          cell_fsrs.at(cell).push_back(fsr_id);
          FSRs_to_CMFD_cells.at(fsr_id) = cell;
        }
      }

      /* Set CMFD cell_fsrs vector of vectors */
      cmfd->setCellFSRs(&cell_fsrs);
    }

    /* Read 2D basis information for OTF 3D solvers */
    if (_segment_formation != +segmentationType::EXPLICIT_2D &&
        _segment_formation != +segmentationType::EXPLICIT_3D)
      readExtrudedFSRInfo(in);

    log::fverbose("Status of reading segments from file: %d", ret);

    /* Close the Track file */
    fclose(in);

    return true;
  }

  /**
   * @brief Read information of all Extruded FSRs from a file.
   * @param in file to read from
   */
  void TrackGenerator::readExtrudedFSRInfo(FILE *in) {}

  /**
   * @brief Splits Track segments into sub-segments for a user-defined
   *        maximum optical length for the problem.
   * @details This routine is needed so that all segment lengths fit
   *          within the exponential interpolation table used in the MOC
   *          transport sweep.
   * @param max_optical_length the maximum optical length
   */
  void TrackGenerator::splitSegments(FP_PRECISION max_optical_length)
  {

    if (!containsSegments())
      log::error("Unable to split segments since segments have not yet been generated");

    if (_segment_formation != +segmentationType::EXPLICIT_3D &&
        _segment_formation != +segmentationType::EXPLICIT_2D)
      log::error("Segments cannot be split for on-the-fly ray tracing");

    /* Split all segments along all Tracks */
    _max_optical_length = max_optical_length;
    SegmentSplitter segment_splitter(this);
    segment_splitter.execute();

    // Indicate that segments should be counted again.
    _segments_counted = false;
  }

  /**
   * @brief Generates the numerical centroids of the FSRs.
   * @details This routine generates the numerical centroids of the FSRs
   *          by weighting the average x and y values of each segment in the
   *          FSR by the segment's length and azimuthal weight. The numerical
   *          centroid fomula can be found in R. Ferrer et. al. "Linear Source
   *          Approximation in CASMO 5", PHYSOR 2012.
   * @param FSR_volumes An array of FSR volumes.
   * 计算每个 FSR 的几何中心（形心）。
   */
  void TrackGenerator::generateFSRCentroids(FP_PRECISION *FSR_volumes)
  {

    log::finfo("Computing centroids of FSRs...");

    long num_FSRs = _geometry->getNumFSRs();

    /* Create temporary array of centroids and initialize to origin */
    Point **centroids = new Point *[num_FSRs];
    for (long r = 0; r < num_FSRs; r++)
    {
      centroids[r] = new Point();
      centroids[r]->setCoords(0.0, 0.0, 0.0);
    }

    /* Generate FSR centroids by looping over all Tracks */
    CentroidGenerator centroid_generator(this);
    centroid_generator.setCentroids(centroids);
    centroid_generator.execute();

#ifdef ENABLE_MPI_
    // If there are processes sharing FSRs, reduce FSR centroids
    /**
     * 如果开启了 MPI，不同进程可能只负责一部分径迹。因此，每个进程算出来的形心只是“局部”的。
     * MPI_Allreduce: 这是一个 MPI 通信函数，将所有进程计算出的坐标加在一起（MPI_SUM），然后把结果分发给所有进程。这样每个进程都能得到完整的、正确的形心坐标。
     */
    if (mpi::isPrdTrackDecomposed())
    {
      _timer->startTimer();

      log::fverbose_once("Reducing FSR centroids...");

      // Number of coordinates
      auto size = 3 * num_FSRs;
      log::verbose_once("Buffer size for reducing FSR centroids = {:.3f} MiB",
                        size * sizeof(double) / (double)(1 << 20));

      // Pack the buffer
      double *buffer = new double[size];
      for (long r = 0; r < num_FSRs; r++)
      {
        buffer[3 * r] = centroids[r]->getX();
        buffer[3 * r + 1] = centroids[r]->getY();
        buffer[3 * r + 2] = centroids[r]->getZ();
      }

      MPI_Allreduce(MPI_IN_PLACE, buffer, size,
                    MPI_DOUBLE, MPI_SUM,
                    mpi::getCommSharedDomain());

      _timer->stopTimer("Transfer FSRs");

      // Unpack the buffer and set the centroid for FSRs
      for (long r = 0; r < num_FSRs; r++)
      {
        centroids[r]->setCoords(buffer[3 * r],
                                buffer[3 * r + 1],
                                buffer[3 * r + 2]);
      }

      delete[] buffer;
    }
#endif

    /* Set the centroid for the FSR */
    for (long r = 0; r < num_FSRs; r++)
    {
      _geometry->setFSRCentroid(r, centroids[r]);
    }

    /* Recenter the segments around FSR centroid */
    if ((_segment_formation == +segmentationType::EXPLICIT_2D ||
         _segment_formation == +segmentationType::EXPLICIT_3D) &&
        _segments_centered == false)
    {
      log::finfo("Centering segments around FSR centroid...");
      RecenterSegments rs(this);
      rs.execute();
      _segments_centered = true;
    }

    delete[] centroids;
  }

  /**
   * @brief Sets the max optical path length of 3D segments for use in
   *        on-the-fly computation
   * @param tau maximum optical path length
   */
  void TrackGenerator::setMaxOpticalLength(FP_PRECISION tau)
  {
    if (_max_optical_length != tau)
    {
      _max_optical_length = tau;

      // Indicate that segments should be counted again.
      _segments_counted = false;
    }
  }

  /**
   * @brief Sets the maximum number of segments per Track
   * @param max_num_segments the maximum number of segments per Track
   */
  void TrackGenerator::setMaxNumSegments(int max_num_segments)
  {
    _max_num_segments = max_num_segments;
  }

  /**
   * @brief Retrieves the max optical path length of 3D segments for use in
   *        on-the-fly computation
   * @return maximum optical path length
   */
  FP_PRECISION TrackGenerator::retrieveMaxOpticalLength()
  {
    return _max_optical_length;
  }

  /**
   * @brief Counts the number of segments for each Track in the Geomtery
   * @details All segments are subject to the max optical path length to
   *          determine the number of segments for each track as well as the
   *          maximum number of segments per Track in the Geometry. For
   *          on-the-fly calculations, the temporary segment buffer is expanded
   *          to fit the calculated maximum number of segments per Track.
   */
  void TrackGenerator::countSegments()
  {

    log::info("Counting segments and allocate temporary segments...");

    // Count the number of segments on each track and update the maximum
    countNumSegments();

    // 内存分配（如果是 OTF 模式）
    // 如果不是显式存储所有线段（EXPLICIT_2D/3D），说明是 On-The-Fly (OTF) 模式。
    // OTF 模式下，线段是算的时候临时生成的，所以需要预先分配一块足够大的“临时缓冲区”来放这些临时线段。
    if (_segment_formation != +segmentationType::EXPLICIT_3D &&
        _segment_formation != +segmentationType::EXPLICIT_2D)
      allocateTemporarySegments();
  }

  /// \brief Count the number of segments on each track and update the maximum.
  /// \details This behaviour is affected by the maximum optical length.
  /// \return Total number of segments.
  // 统计每条 track 的分段数并更新最大值。最大合法段数会受到“最大光学长度”配置的影响，所以这一步通常在分段完成后立即执行
  long TrackGenerator::countNumSegments()
  {

    // 1. 安全检查
    // After segmentize() or segmentizeExtruded(), this flag will be set to true.
    // 如果还没生成线段（segmentize 还没跑），就没法数，报错。
    if (!containsSegments())
      log::error("Cannot count segments since they have not been generated.");

    log::verbose_once("Counting segments and update the maximum number...");

    // 2. 使用辅助类 SegmentCounter 进行统计
    // SegmentCounter 是一个专门用来数数的工具类。
    // 为什么要专门搞个类？因为数数的过程可能很复杂（涉及并行、最大光学长度拆分等逻辑）。
    SegmentCounter counter(this);
    counter.countTotalNumSegments();
    counter.execute();

    // Update the number of segments
#pragma omp critical
    _num_segments = counter.getTotalNumSegments();

    // Set the flag
    _segments_counted = true;

    return _num_segments; // 返回的是所有 track 的 segment 数总和
  }

  /**
   * @brief Creates a Track array by increasing uid
   * @details An array is created which indexes Tracks by increasing uid.
   *          Parallel groups are also initialized -- groups of Tracks that can
   *          be computed in parallel without the potential of overwriting
   *          angular fluxes of connecting tracks prematurely.
   *  重新构建 _tracks_2D_array，让它按 uid 顺序指向本进程负责的每一条 2D 轨迹。这样后续求解、调度时只需要顺着一个连续数组扫描，不必再走原来的“按方位角→链→轨迹”的复杂嵌套结构，每条轨迹也已经知道自己的 uid，方便快速定位。
   */
  void TrackGenerator::initializeTracksArray()
  {

    log::finfo("Initializing 2D tracks array...");

    /* Allocate memory for tracks array */
    if (_tracks_2D_array != NULL)
      delete[] _tracks_2D_array;

    long num_2D_tracks = getMyNum2DTracks();

    _tracks_2D_array = new Track *[num_2D_tracks];

    /* Loop over all 2D tracks */
    long uid = 0;
    for (int a = 0; a < _num_azim / 2; a++)
    {
      for (long i = 0; i < getMyNum2DTracks(a); i++)
      {

        /* Get current track and azim group ids */
        Track *track = &_tracks_2D[a][i];

        track->setUid(uid); // local id
        _tracks_2D_array[uid] = track;
        uid++;
      }
    }
  }

  /**
   * @brief returns whether periodic boundaries are present in Track generation
   * @return a boolean value - true if periodic; false otherwise
   */
  bool TrackGenerator::getPeriodic()
  {
    return _periodic;
  }

  /**
   * @brief Sets a flag to record all segment information in the tracking file
   * @param A boolean value to determine whether or not to record segment
   *        information in the tracking file: true to record, false not to record
   */
  void TrackGenerator::setDumpSegments(bool dump_segments)
  {
    _dump_segments = dump_segments;
  }

  /**
   * @brief Resets the TrackGenerator to not contain tracks or segments
   */
  void TrackGenerator::resetStatus()
  {
    _contains_2D_tracks = false;
    _contains_2D_segments = false;
    _use_input_file = false;
    _tracks_filename = "";
  }

  /**
   * @brief Allocates memory for temporary segment storage if necessary
   * @details Temporary segments are not allocated for 2D calculations
   */
  void TrackGenerator::allocateTemporarySegments() {}

  /**
   * @brief Get the id of a 2D Track based on its azimuthal angle and index in the
   *        azimuthal stack
   * @param a azimuthal angle of the Track
   * @param x index in azimuthal stack
   * @return Track unique id
   */
  long TrackGenerator::get2DTrackID(int a, long x)
  {

    long uid = 0;

    for (int ai = 0; ai < a; ai++)
      uid += getNumX(ai) + getNumY(ai);

    uid += x;
    return uid;
  }

  //--------------------------------------------------------------------
  // Common methods for track generation
  //--------------------------------------------------------------------

  /// \brief Returns azim index of links
  /// \details Along a reflective track, links have differenct azims.
  ///          Every time the cyclic track reach a boundary, its direction
  ///          changes.
  /// \param azim original azim index
  /// \param link link index
  /// \return azim index of current link
  int TrackGenerator::get2DLinkAzimIndex(int azim, long link)
  {
    if (!isReflectiveCyclic())
      return azim;
    else
    {
      // The azim will change every reflection
      return link % 2 == 0 ? azim : _num_azim / 2 - azim - 1;
    }
  }

  /// \brief Get the index of the chain a track resides
  /// \details In OpenMOC, the index can be computed by simply translating
  ///          the track, while in ANT-MOC we have a different way.
  long TrackGenerator::get2DChainIndex(int azim, long xy)
  {
#ifdef USTB_
    return _tracks_2D[azim][xy].getChainIndex();
#elif MIT_
    return xy % _num_x[azim];
#endif
  }

  /// \brief Get the number of azims of unique chains
  /// \details Reflective chains only have na/4 azims because they are
  ///          circles. Only ANT-MOC algorithms generate reflective cyclic
  ///          tracks so that it always returns na/2 in OpenMOC algorithms.
  int TrackGenerator::getNumChainAzims()
  {
    if (isReflectiveCyclic())
      return _num_azim / 4;
    else
      return _num_azim / 2;
  }

  /// \brief get the total number of 2D chains
  /// \details Reflective tracks are different from periodic ones
  long TrackGenerator::getNum2DChains()
  {
    long num_chains = 0;
    for (int a = 0; a < getNumChainAzims(); a++)
    {
      num_chains += getNum2DChains(a);
    }
    return num_chains;
  }

  /// \brief Get the number of 2D chains belonging to me
  long TrackGenerator::getMyNum2DChains()
  {
    long num_chains = 0;
    for (int a = 0; a < getNumChainAzims(); a++)
    {
      num_chains += getMyNum2DChains(a);
    }
    return num_chains;
  }

#ifdef MIT_
  //--------------------------------------------------------------------
  // Methods for OpenMOC track generation
  //--------------------------------------------------------------------

  /// \brief Compute the number of links of 2D chain tracks
  /// \details A track chain is a part of a cyclic track truncated by
  ///          the YMAX surface. The number of links varies among
  ///          track chains of the same azimuthal angle.
  long TrackGenerator::getNum2DLinks(int azim, long x)
  {
    long link_index = 0;
    Track *track = &_tracks_2D[azim][x];
    while (track->getXYIndex() < _num_y[azim])
    {
      link_index++;
      track = _tracks_2D_array[track->getTrackPrdcFwd()];
    }
    return link_index + 1;
  }

#elif USTB_
  //----------------------------------------------------------------------
  // Cyclic ray tracing methods for ANT-MOC
  //----------------------------------------------------------------------

  /// \brief Compute the number of links of cyclic tracks
  /// \details Cyclic tracks of the same azim have same number of links.
  ///          That is, n = (nx+ny)/gcd(nx,ny).
  ///          If cyclic tracks are reflective, double the number.
  long TrackGenerator::getNum2DLinks(int azim)
  {
    long num_prdc = (_num_x[azim] + _num_y[azim]) / getNum2DChains(azim);
    return (1 + 1 * isReflectiveCyclic()) * num_prdc;
  }

#endif // MIT_ or USTB_

  //----------------------------------------------------------------------
  // ID conversion methods for ANT-MOC
  //----------------------------------------------------------------------

  /// \brief Convert the local id of a chain to global id
  /// \details This method support conversion between local ids
  ///          and global ids. The local id is converted to global
  ///          uid, which is then converted to global id.
  long TrackGenerator::map2DChainId2Gid(int azim, long x)
  {

#ifdef ENABLE_MPI_
    return x + _chain_num_offsets[azim];
#else
    return x;
#endif
  }

  /// \brief Convert the global id of a chain to local id
  /// \details This method support conversion between local ids
  ///          and global ids. The global id is converted to local
  ///         uid, which is then converted to local id.
  long TrackGenerator::map2DChainGid2Id(int azim, long gx)
  {

#ifdef ENABLE_MPI_
    return gx - _chain_num_offsets[azim];
#else
    return gx;
#endif
  }

#ifdef ENABLE_MPI_
  //----------------------------------------------------------------------
  // Cyclic ray tracing methods for parallel ANT-MOC
  //
  // The following methods are supposed to be activated
  // when periodic tracks are decomposed
  //----------------------------------------------------------------------

  void TrackGenerator::setTrackMappingType(trackMappingType track_mapping_type)
  {
    _track_mapping_type = track_mapping_type;
  }

  trackMappingType TrackGenerator::getTrackMappingType()
  {
    return _track_mapping_type;
  }

  void TrackGenerator::generate2DChainMapping()
  {

    log::finfo("Initializing 2D track mapping...");

    _timer->startTimer();

    auto rank = mpi::getRankAuto();
    auto num_procs = mpi::getNumProcsAuto();
    int chain_azims = getNumChainAzims();
    TrackLoadBalancePtr balancer;
    IntVec2D load_map;

    // Automatically choose the proper algorithm or specify one manually
    if (_track_mapping_type == +trackMappingType::AUTO)
    {
      balancer = generateOptimalLoadMap(load_map);
    }
    else
    {
      balancer = TrackLoadBalance::getTrackLoadBalancer(_track_mapping_type,
                                                        chain_azims,
                                                        _cyclic_length,
                                                        _num_chains);
      load_map = balancer->computeLoadMap(num_procs);
    }

    balancer->printStatistics(load_map, num_procs);

    _my_num_chains = new long[_num_azim / 2]();

    for (int a = 0; a < chain_azims; a++)
    {
      // Compute the number of local chains
      _my_num_chains[a] = load_map[rank][a];

      // Compute the offset
      _chain_num_offsets[a] = 0;
      for (int r = 0; r < rank; ++r)
        _chain_num_offsets[a] += load_map[r][a];

      if (isReflectiveCyclic())
      {
        auto ca = _num_azim / 2 - a - 1;
        _my_num_chains[ca] = _my_num_chains[a];
        _chain_num_offsets[ca] = _chain_num_offsets[a];
      }
    }

    _timer->stopTimer("2D Track Mapping");
  }

  /// \details Automatically choose the proper algorithm.
  TrackLoadBalancePtr TrackGenerator::generateOptimalLoadMap(IntVec2D &load_map)
  {

    auto num_procs = mpi::getNumProcsAuto();
    int chain_azims = getNumChainAzims();
    std::vector<trackMappingType> types = {
        trackMappingType::BLOCK,
        trackMappingType::CYCLIC_TRACK,
        trackMappingType::ANGLE};

    trackMappingType opt_type = trackMappingType::BLOCK;
    TrackLoadBalancePtr balancer;
    double speedup = -1;

    // Try each of the algorithms to determine the best one of them
    for (auto type : types)
    {
      balancer = TrackLoadBalance::getTrackLoadBalancer(type,
                                                        chain_azims,
                                                        _cyclic_length,
                                                        _num_chains);
      auto cur_load_map = balancer->computeLoadMap(num_procs);
      auto cur_speedup = balancer->computeExpectedSpeedup(cur_load_map, num_procs);

      if (cur_speedup > speedup)
      {
        opt_type = type;
        speedup = cur_speedup;
        load_map = std::move(cur_load_map);
      }
    }

    log::info("Auto-selected track mapping algorithm = {}",
              stringutils::underscoreToSpace(opt_type._to_string()));

    return balancer;
  }

  /// \brief Generate all of the FSRs before segmentation
  /// \details This method is supposed to ensure the order of FSR IDs.
  void TrackGenerator::preSegmentize2D(double z_coord)
  {

    std::vector<double> z_coords = {z_coord};
    long num_2D_tracks = getMyNum2DTracks();
    Progress progress(num_2D_tracks, "Pre-segmenting 2D Tracks", 0.1, _geometry, true);

    _timer->startTimer();

#pragma omp parallel for schedule(dynamic)
    for (long uid = 0; uid < num_2D_tracks; uid++)
    {
      progress.incrementCounter();
      _geometry->segmentize2D(_tracks_2D_array[uid], z_coord, true);
    }

    _timer->stopTimer("Pre-segmenting");

    // Reduce all of the FSRs
    _timer->startTimer();
    _geometry->allReduceFSRs(z_coords, false);
    _timer->stopTimer("Reduce FSRs");
  }

#endif // ENABLE_MPI_

  /**
   * @brief Return the number of local 2D tracks of a specified azim
   * @details The number of 2D tracks under an specified azim is
   *          num_chains * num_links
   * 这个方位角方向下的 2D track 有多少条
   */
  long TrackGenerator::getMyNum2DTracks(int azim)
  {
    // #ifdef USTB_ / #elif MIT_ / #endif
    // 通过编译宏选择不同学校版本的实现：
#ifdef USTB_
    long num_tracks = getMyNum2DChains(azim) * getNum2DLinks(azim);
    return num_tracks / (1 + 1 * isReflectiveCyclic());
#elif MIT_
    return getNum2DTracks(azim);
#endif
  }

  /**
   * @brief Return the number of 2D tracks belonging to current rank
   * @details The number of 2D tracks under an specified azim is
   *          num_chains * num_links
   *所有方位角方向上的总 track 数
   */
  long TrackGenerator::getMyNum2DTracks()
  {

    long num_2D_tracks = 0;

    for (int a = 0; a < _num_azim / 2; a++)
      num_2D_tracks += getMyNum2DTracks(a);

    return num_2D_tracks;
  }

  /**
   * @brief Convert the local id of a track to global id
   * @details This method support conversion between local ids
   *          and global ids.
   */
  long TrackGenerator::map2DTrackId2Gid(int azim, long xy)
  {
    long gxy = xy;
#ifdef ENABLE_MPI_
    if (mpi::isPrdTrackDecomposed())
      gxy = _tracks_2D[azim][xy].getXYIndex();
#endif
    return gxy;
  }

  //----------------------------------------------------------------------
  // Profiling utils
  //----------------------------------------------------------------------

  double TrackGenerator::computeAzimUniformity()
  {
    int n = _num_azim / 2 - 1;
    double *intervals = new double[n]();

    for (int i = 0; i < n; i++)
      intervals[i] = _quadrature->getPhi(i + 1) - _quadrature->getPhi(i);

    double sd = computeStandardDeviation(intervals, n);
    delete[] intervals;
    return sd;
  }

  double TrackGenerator::computeStandardDeviation(const double *items, int n)
  {
    double tot_value = 0;
    for (int i = 0; i < n; i++)
      tot_value += items[i];

    double avg = tot_value / n; // average

    tot_value = 0;
    for (int i = 0; i < n; i++)
      tot_value += std::pow(items[i] - avg, 2);

    return std::sqrt(tot_value / (n - 1));
  }

  void TrackGenerator::reportTemporarySegmentsStorage() {}

  //----------------------------------------------------------------------
  // Dumping Tracks
  //----------------------------------------------------------------------

  /// \brief Returns a vector containing endpoints of each 2D track
  /// \details Tracks are traversed by indices of track chains
  std::vector<TrackMeshData> TrackGenerator::getFormattedTracks()
  {

    std::vector<TrackMeshData> cyclic_coords; ///< A vector of track mesh data

    // Loop over tracks along cyclic tracks
    for (int a = 0; a < getNumChainAzims(); ++a)
    {
      for (long x = 0; x < getMyNum2DChains(a); ++x)
      {
        long link = 0;
        while (true)
        {
          Track *track_2D = _tracks_2D_chains[a][x][link++];

          // Push back an object of TrackMeshData
          cyclic_coords.push_back(
              {track_2D->getUid(), // FIXME, supposed to be guid
               *track_2D->getStart(),
               *track_2D->getEnd()});

#ifdef USTB_
          if (link >= getNum2DLinks(a))
            break;
#elif MIT_
          if (track_2D->getXYIndex() >= _num_y[a])
            break;
#endif
        }
      }
    }

    return cyclic_coords;
  }

} /* namespace antmoc */
