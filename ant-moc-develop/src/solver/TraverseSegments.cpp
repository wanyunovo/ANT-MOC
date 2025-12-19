#include "antmoc/TraverseSegments.h"
#include "antmoc/mpi_utils.h"
#include "antmoc/Cmfd.h"
#include "antmoc/Geometry.h"
#include "antmoc/Track3D.h"
#include "antmoc/TrackGenerator3D.h"

namespace antmoc
{

  /**
   * @brief Constructor for the TraverseSegments class assigns the TrackGenerator
   *        and pulls relevant information from it.
   */
  TraverseSegments::TraverseSegments(TrackGenerator *track_generator) : _segment_formation(track_generator->getSegmentFormation())
  {

    /* Save the track generator */
    _track_generator = track_generator;

    /* Determine if a global z-mesh is used for 3D calculations */
    _track_generator_3D = dynamic_cast<TrackGenerator3D *>(track_generator);
    if (_track_generator_3D != nullptr)
    {
      _track_generator_3D->retrieveGlobalZMesh(_global_z_mesh, _mesh_size);
    }
  }

  /**
   * @brief Destructor for TraverseSegments.
   */
  TraverseSegments::~TraverseSegments()
  {
  }

  /**
   * @brief Loops over Tracks, applying the provided kernel to all segments and
   *        the functionality described in onTrack(...) to all Tracks.
   * @details The segment formation method imported from the TrackGenerator
   *          during construction is used to redirect to the appropriate looping
   *          scheme. If a kernel is provided (not nullptr) then it is deleted at
   *          the end of the looping scheme.
   * @param kernel MOCKernel to apply to all segments
   */
  void TraverseSegments::loopOverTracks(MOCKernel *kernel)
  {

    switch (_segment_formation)
    {
    case segmentationType::EXPLICIT_2D:
      loopOverTracks2D(kernel);
      break;
    case segmentationType::EXPLICIT_3D:
      loopOverTracksExplicit(kernel);
      break;
    case segmentationType::OTF_TRACKS:
      loopOverTracksByTrackOTF(kernel);
      break;
    case segmentationType::OTF_STACKS:
      loopOverTracksByStackOTF(kernel);
      break;
    }

    if (kernel != nullptr)
      delete kernel;
  }

  /**
   * @brief Loops over all explicit 2D Tracks.
   * @details The onTrack(...) function is applied to all 2D Tracks and the
   *          specified kernel is applied to all segments. If nullptr is provided
   *          for the kernel, only the onTrack(...) functionality is applied.
   * @param kernel The MOCKernel dictating the functionality to apply to
   *        segments
   */
  void TraverseSegments::loopOverTracks2D(MOCKernel *kernel)
  {

    /* Loop over all parallel tracks for each azimuthal angle */
    Track **tracks_2D = _track_generator->get2DTracksArray();
    long num_tracks = _track_generator->getMyNum2DTracks();

#pragma omp for schedule(dynamic)
    for (long t = 0; t < num_tracks; t++)
    {

      Track *track_2D = tracks_2D[t];
      segment *segments = track_2D->getSegments();

      /* Operate on segments if necessary */
      if (kernel != nullptr)
      {
        kernel->newTrack(track_2D);
        traceSegmentsExplicit(track_2D, kernel);
      }

      /* Operate on the Track */
      onTrack(track_2D, segments);
    }
  }

  /**
   * @brief Loops over all explicit 3D Tracks.
   * @details The onTrack(...) function is applied to all 3D Tracks and the
   *          specified kernel is applied to all segments. If NULL is provided
   *          for the kernel, only the onTrack(...) functionality is applied.
   * @param kernel The MOCKernel dictating the functionality to apply to
   *        segments
   */
  void TraverseSegments::loopOverTracksExplicit(MOCKernel *kernel)
  {

    Track3D ****tracks_3D = _track_generator_3D->get3DTracks();
    int num_azim = _track_generator_3D->getNumAzim();
    int num_polar = _track_generator_3D->getNumPolar();
    int ***tracks_per_stack = _track_generator_3D->getTracksPerStack();

    /* Loop over all tracks, parallelizing over parallel 2D tracks */
    for (int a = 0; a < num_azim / 2; a++)
    {
      int num_xy = _track_generator->getMyNum2DTracks(a);
#pragma omp for schedule(dynamic) collapse(2)
      for (int i = 0; i < num_xy; i++)
      {

        /* Loop over polar angles */
        for (int p = 0; p < num_polar; p++)
        {

          /* Loop over tracks in the z-stack */
          for (int z = 0; z < tracks_per_stack[a][i][p]; z++)
          {

            /* Extract 3D track */
            Track *track_3D = &tracks_3D[a][i][p][z];

            /* Operate on segments if necessary */
            if (kernel != nullptr)
            {

              /* Reset kernel for a new Track */
              kernel->newTrack(track_3D);

              /* Trace the segments on the track */
              traceSegmentsExplicit(track_3D, kernel);
            }

            /* Operate on the Track */
            segment *segments = track_3D->getSegments();
            onTrack(track_3D, segments);
          }
        }
      }
    }
  }

  /**
   * @brief Loops over all 3D Tracks using axial on-the-fly ray tracing by Track.
   * @details The onTrack(...) function is applied to all 3D Tracks and the
   *          specified kernel is applied to all segments. If nullptr is provided
   *          for the kernel, only the onTrack(...) functionality is applied.
   * @param kernel The MOCKernel dictating the functionality to apply to
   *        segments
   */
  void TraverseSegments::loopOverTracksByTrackOTF(MOCKernel *kernel)
  {

    int num_2D_tracks = _track_generator_3D->getMyNum2DTracks();
    Track **tracks_2D = _track_generator_3D->get2DTracksArray();
    int ***tracks_per_stack = _track_generator_3D->getTracksPerStack();
    // unused
    // int num_azim = _track_generator->getNumAzim();
    int num_polar = _track_generator_3D->getNumPolar();
    int tid = omp_get_thread_num();

    /* Loop over flattened 2D tracks */
#pragma omp for schedule(dynamic)
    for (int ext_id = 0; ext_id < num_2D_tracks; ext_id++)
    {

      /* Extract indices of 3D tracks associated with the flattened track */
      Track *flattened_track = tracks_2D[ext_id];
      TrackStackIndexes tsi;
      tsi._azim = flattened_track->getAzimIndex();
      tsi._xy = flattened_track->getXYIndex();

      /* Loop over polar angles */
      for (int p = 0; p < num_polar; p++)
      {

        /* Loop over tracks in the z-stack */
        for (int z = 0; z < tracks_per_stack[tsi._azim][tsi._xy][p]; z++)
        {

          /* Extract 3D track and retrieve its information */
          Track3D track_3D;
          tsi._polar = p;
          tsi._z = z;
          _track_generator_3D->getTrackOTF(&track_3D, &tsi);

          /* Operate on segments if necessary */
          if (kernel != nullptr)
          {

            /* Reset kernel for a new Track */
            kernel->newTrack(&track_3D);
            double theta = track_3D.getTheta();
            Point *start = track_3D.getStart();

            /* Trace the segments on the track */
            traceSegmentsOTF(flattened_track, start, theta, kernel);
            track_3D.setNumSegments(kernel->getCount());
          }

          /* Operate on the Track */
          segment *segments = _track_generator_3D->getTemporarySegments(tid);
          onTrack(&track_3D, segments);
        }
      }
    }
  }

  /**
   * @brief Loops over all 3D Tracks using axial on-the-fly ray tracing by
   *        z-stack.
   * @details The onTrack(...) function is applied to all 3D Tracks and the
   *          specified kernel is applied to all segments. If NULL is provided
   *          for the kernel, only the onTrack(...) functionality is applied.
   * @param kernel The MOCKernel dictating the functionality to apply to
   *        segments
   */
  void TraverseSegments::loopOverTracksByStackOTF(MOCKernel *kernel)
  {

    int num_2D_tracks = _track_generator_3D->getMyNum2DTracks();
    Track **flattened_tracks = _track_generator_3D->get2DTracksArray();
    int ***tracks_per_stack = _track_generator_3D->getTracksPerStack();
    int num_polar = _track_generator_3D->getNumPolar();
    int tid = omp_get_thread_num();

    /* Allocate array of current Tracks */
    Track3D *current_stack = _track_generator_3D->getTemporary3DTracks(tid);

    /* Loop over flattened 2D tracks */
#pragma omp for schedule(dynamic)
    for (int ext_id = 0; ext_id < num_2D_tracks; ext_id++)
    {

      /* Extract indices of 3D tracks associated with the flattened track */
      TrackStackIndexes tsi;
      Track *flattened_track = flattened_tracks[ext_id];
      tsi._azim = flattened_track->getAzimIndex();
      tsi._xy = flattened_track->getXYIndex();

      /* Loop over polar angles */
      for (int p = 0; p < num_polar; p++)
      {

        /* Retrieve information for the first 3D Track in the z-stack */
        tsi._polar = p;
        int stack_size = tracks_per_stack[tsi._azim][tsi._xy][tsi._polar];
        for (int z = 0; z < stack_size; z++)
        {
          tsi._z = z;
          _track_generator_3D->getTrackOTF(&current_stack[z], &tsi);
        }

        if (kernel != nullptr)
        {

          /* Reset kernel to for the new base Track */
          kernel->newTrack(&current_stack[0]);

          /* Trace all segments in the z-stack */
          traceStackOTF(flattened_track, p, kernel);
          current_stack[0].setNumSegments(kernel->getCount());
        }

        /* Operate on the Track */
        segment *segments = _track_generator_3D->getTemporarySegments(tid);
        onTrack(&current_stack[0], segments);
      }
    }
  }

  /**
   * @brief Loops over segments in a Track when segments are explicitly generated.
   * @details All segments in the provided Track are looped over and the provided
   *          MOCKernel is applied to them.
   * @param track The Track whose segments will be traversed
   * @param kernel The kernel to apply to all segments
   */
  void TraverseSegments::traceSegmentsExplicit(Track *track, MOCKernel *kernel)
  {

    /* Get direction of the track */
    double phi = track->getPhi();
    double theta = M_PI_2;
    Track3D *track_3D = dynamic_cast<Track3D *>(track);
    if (track_3D != NULL)
      theta = track_3D->getTheta();

    for (int s = 0; s < track->getNumSegments(); s++)
    {
      segment *seg = track->getSegment(s);
      kernel->execute(seg->_length, seg->_material, seg->_region_id, 0,
                      seg->_cmfd_surface_fwd, seg->_cmfd_surface_bwd,
                      seg->_starting_position[0], seg->_starting_position[1],
                      seg->_starting_position[2], phi, theta);
    }
  }

  /**
   * @brief Computes 3D segment lengths on-the-fly for a single 3D track given an
   *        associated 2D Track with a starting point and a polar angle. The
   *        computed segments are passed to the provided kernel.
   * @details Segment lengths are computed on-the-fly using 2D segment lengths
   *          stored in a 2D Track object and 1D meshes from the extruded
   *          FSRs. Note: before calling this function with a SegmentationKernel,
   *          the memory for the segments should be allocated and referenced by
   *          the kernel using the setSegments routine in the kernel.
   * @param flattened_track the 2D track associated with the 3D track for which
   *        3D segments are computed
   * @param start the starting coordinates of the 3D track
   * @param theta the polar angle of the 3D track
   * @param kernel An MOCKernel object to apply to the calculated 3D segments
   */
  void TraverseSegments::traceSegmentsOTF(Track *flattened_track, Point *start,
                                          double theta, MOCKernel *kernel)
  {

    /* Create unit vector */
    double phi = flattened_track->getPhi();
    double cos_phi = cos(phi);
    double sin_phi = sin(phi);
    double cos_theta = cos(theta);
    double sin_theta = sin(theta);
    int sign = (cos_theta > 0) - (cos_theta < 0);

    /* Extract starting coordinates */
    double x_start_3D = start->getX();
    double x_start_2D = flattened_track->getStart()->getX();
    double x_coord = x_start_3D;
    double y_coord = start->getY();
    double z_coord = start->getZ();

    /* Find 2D distance from 2D edge to start of track */
    double start_dist_2D = (x_start_3D - x_start_2D) / cos_phi;

    /* Find starting 2D segment */
    int seg_start = 0;
    segment *segments_2D = flattened_track->getSegments();
    for (int s = 0; s < flattened_track->getNumSegments(); s++)
    {

      /* Determine if start point of track is beyond current 2D segment */
      double seg_len_2D = segments_2D[s]._length;
      if (start_dist_2D > seg_len_2D)
      {
        // 扣除，最终的值是3D轨迹起点离它所在2D线段起点的距离
        start_dist_2D -= seg_len_2D;
        seg_start++;
      }
      else
      {
        break;
      }
    }

    // 后面要用到ExtrudedFSR
    Geometry *geometry = _track_generator_3D->getGeometry();
    Cmfd *cmfd = geometry->getCmfd();

    /* For very short tracks, it's possible no significant segments will be
     * traversed */
    if (seg_start == flattened_track->getNumSegments())
    {
      log::fwarn("Track of zero length encountered at starting point "
                 "%s traveling on 2D Track: %s at polar angle cos %3.2f "
                 "degrees on domain with z-bounds %3.2f and %3.2f",
                 start->toString().c_str(), flattened_track->toString().c_str(),
                 cos(theta), geometry->getMinZ(), geometry->getMaxZ());
      return;
    }

    /* Extract the appropriate starting mesh */
    int num_fsrs;
    double *axial_mesh;
    bool contains_global_z_mesh;
    // 有全局轴向网就用它，没有就用局部轴向网
    if (_global_z_mesh != nullptr)
    {
      contains_global_z_mesh = true;
      num_fsrs = _mesh_size;
      axial_mesh = _global_z_mesh;
    }
    else
    {
      contains_global_z_mesh = false;
      int extruded_fsr_id = segments_2D[seg_start]._region_id; // 每条2D线段对应的extruded fsr
      ExtrudedFSR *extruded_FSR = geometry->getExtrudedFSR(extruded_fsr_id);
      num_fsrs = extruded_FSR->_num_fsrs;
      axial_mesh = extruded_FSR->_mesh; // 轴向网包括ZMAX和ZMIN，因此比fsr数量多1
    }

    /* Get the starting z index */
    int z_ind = findMeshIndex(axial_mesh, num_fsrs + 1, z_coord, sign);

    /* Loop over 2D segments */
    bool first_segment = true;
    bool segments_complete = false;
    for (int s = seg_start; s < flattened_track->getNumSegments(); s++)
    {

      /* Extract extruded FSR */
      int extruded_fsr_id = segments_2D[s]._region_id;
      ExtrudedFSR *extruded_FSR = geometry->getExtrudedFSR(extruded_fsr_id);

      /* Determine new mesh and z index */
      if (first_segment || contains_global_z_mesh)
      {
        first_segment = false;
      }
      else
      {
        /* Determine the axial region */
        num_fsrs = extruded_FSR->_num_fsrs;
        axial_mesh = extruded_FSR->_mesh;
        z_ind = findMeshIndex(axial_mesh, num_fsrs + 1, z_coord, sign);
      }

      /* Extract 2D segment length */
      double remaining_length_2D = segments_2D[s]._length - start_dist_2D;
      start_dist_2D = 0;

      /* Transport along the 2D segment until it is completed */
      while (remaining_length_2D > 0)
      {
        /* 因为外层循环在遍历线段，所以这里的while结束时，保证一条2D线段所对的
         * 所有3D轨迹都被处理完了。
         */
        /* Calculate 3D distance to z intersection */
        // 即，轨迹起点到它和某一个轴向网交点的距离
        double z_dist_3D;
        if (sign > 0)
          z_dist_3D = (axial_mesh[z_ind + 1] - z_coord) / cos_theta;
        else
          z_dist_3D = (axial_mesh[z_ind] - z_coord) / cos_theta;

        /* Calculate 3D distance to end of segment */
        double seg_dist_3D = remaining_length_2D / sin_theta;

        /* 下面的分支是个分类，seg_dist_3D比较大意味着轨迹从mesh的
         * 顶边/底边出去了，或者说先到达顶边/底边；反之，说明轨迹
         * 先到达左边/右边。
         */
        /* Calcualte shortest distance to intersection */
        double dist_2D;
        double dist_3D;
        int z_move;
        if (z_dist_3D <= seg_dist_3D)
        {
          dist_2D = z_dist_3D * sin_theta;
          dist_3D = z_dist_3D;
          z_move = sign;
        }
        else
        {
          dist_2D = remaining_length_2D;
          dist_3D = seg_dist_3D;
          z_move = 0;
        }

        /* Get the 3D FSR */
        long fsr_id = extruded_FSR->_fsr_ids[z_ind];

        /* Calculate CMFD surface */
        int cmfd_surface_bwd = -1;
        int cmfd_surface_fwd = -1;
        if (cmfd != nullptr && dist_3D > TINY_MOVE)
        {

          /* Determine if this is the first 3D segment handled for the flattened
             2D segment. If so, get the 2D cmfd surface. */
          if (segments_2D[s]._length - remaining_length_2D <= TINY_MOVE)
            cmfd_surface_bwd = segments_2D[s]._cmfd_surface_bwd;

          /* Determine if this is the last 3D segment handled for the flattened
             2D segment. If so, get the 2D cmfd surface. */
          double next_dist_3D = (remaining_length_2D - dist_2D) / sin_theta;
          if (z_move == 0 || next_dist_3D <= TINY_MOVE)
            cmfd_surface_fwd = segments_2D[s]._cmfd_surface_fwd;

          /* Get CMFD cell */
          int cmfd_cell = geometry->getCmfdCell(fsr_id);

          /* Find the backwards surface */
          cmfd_surface_bwd = cmfd->findCmfdSurfaceOTF(cmfd_cell, z_coord,
                                                      cmfd_surface_bwd);

          /* Find forward surface */
          double z_coord_end = z_coord + dist_3D * cos_theta;
          cmfd_surface_fwd = cmfd->findCmfdSurfaceOTF(cmfd_cell, z_coord_end,
                                                      cmfd_surface_fwd);
        }

        /* Operate on segment */
        if (dist_3D > TINY_MOVE)
        {
          double x_centroid = 0;
          double y_centroid = 0;
          double z_centroid = 0;
          if (geometry->containsFSRCentroids())
          {
            Point *centroid = geometry->getFSRCentroid(fsr_id);
            x_centroid = centroid->getX();
            y_centroid = centroid->getY();
            z_centroid = centroid->getZ();
          }

          kernel->execute(dist_3D, extruded_FSR->_materials[z_ind], fsr_id, 0,
                          cmfd_surface_fwd, cmfd_surface_bwd,
                          x_coord - x_centroid, y_coord - y_centroid,
                          z_coord - z_centroid, phi, theta);
        }

        /* Move axial height to end of segment */
        x_coord += dist_3D * sin_theta * cos_phi;
        y_coord += dist_3D * sin_theta * sin_phi;
        z_coord += dist_3D * cos_theta;

        /* Shorten remaining 2D segment length and move axial level */
        remaining_length_2D -= dist_2D;
        // 全局轴向网的情况下会使用这个值
        z_ind += z_move;

        /* Check if the track has crossed a Z boundary */
        if (z_ind < 0 or z_ind >= num_fsrs)
        {

          /* Reset z index */
          if (z_ind < 0)
            z_ind = 0;
          else
            z_ind = num_fsrs - 1;

          /* Mark the 2D segment as complete */
          segments_complete = true;
          break;
        }
      }

      /* Check if the track is completed due to an axial boundary */
      if (segments_complete)
        break;
    }
  }

  /** 使用 OTF 方法追踪整个轨迹堆中的轨迹
   * 该函数针对给定的2D轨迹及其关联的极角索引，计算相应的3D段长度，并将这些段长度应用于指定的MOCKernel。
   * @brief Computes 3D segment lengths on-the-fly for all tracks in a z-stack
   *        for a given associated 2D Track and a polar index on-the-fly and
   *        passes the computed segments to the provided kernel.
   * @details Segment lengths are computed on-the-fly using 2D segment lengths
   *          stored in a 2D Track object and 1D meshes from the extruded
   *          FSRs. Note: before calling this function with SegmentationKernels,
   *          the memory for the segments should be allocated and referenced by
   *          the kernel using the setSegments routine.
   * @param flattened_track the 2D track associated with the z-stack for which
   *        3D segments are computed   与z堆栈关联的2D轨迹对象，计算3D段长度的基础
   * @param polar_index the index into the polar angles which is associated with
   *        the polar angle of the z-stack
   * @param kernel The MOCKernel to apply to the calculated 3D segments
   */
  void TraverseSegments::traceStackOTF(Track *flattened_track, int polar_index,
                                       MOCKernel *kernel)
  {

    /* Extract information about the z-stack */
    int azim_index = flattened_track->getAzimIndex();
    long track_index = flattened_track->getXYIndex();
    int ***tracks_per_stack = _track_generator_3D->getTracksPerStack();       // 轨迹堆数组
    int num_z_stack = tracks_per_stack[azim_index][track_index][polar_index]; // 获取轨迹堆中有多少轨迹
    double z_spacing = _track_generator_3D->getZSpacing(azim_index, polar_index);

    /* Get information for the first Track in the z-stack */
    TrackStackIndexes tsi;
    Track3D first;
    tsi._azim = azim_index;
    tsi._xy = track_index;
    tsi._polar = polar_index;
    tsi._z = 0; // 将 z 轴索引初始化为 0
    /*调用 getTrackOTF 方法，使用 tsi 结构体中的索引信息获取与 tsi 对应的第一条3D轨迹，
    并存储在 first 对象中。这个方法会动态生成3D轨迹的相关信息。*/
    _track_generator_3D->getTrackOTF(&first, &tsi); // 即时生成3D轨迹相关信息,仅一条(第一条)轨迹&track 3D
    double theta = first.getTheta();                // 3D轨迹的极角角度

    /* Create unit vector */
    double phi = flattened_track->getPhi();
    double cos_phi = cos(phi);
    double sin_phi = sin(phi);
    double cos_theta = cos(theta);
    double sin_theta = sin(theta);
    double tan_theta = sin_theta / cos_theta;
    int sign = (cos_theta > 0) - (cos_theta < 0);
    double track_spacing_3D = z_spacing / std::abs(cos_theta);

    /* Find 2D distance from 2D edge to start of track 查找从二维边到轨迹起点的二维距离*/
    double x_start_3D = first.getStart()->getX();
    double x_start_2D = flattened_track->getStart()->getX();
    double y_start_2D = flattened_track->getStart()->getY();
    double start_dist_2D = (x_start_3D - x_start_2D) / cos_phi; // 从2D轨迹起点到3D轨迹起点在2D平面上的距离

    /* Calculate starting intersection of lowest track with z-axis 计算最低轨迹与z轴的起始交点*/
    double z0 = first.getStart()->getZ();
    double start_z = z0 - start_dist_2D / tan_theta; // 第0条轨迹（也就是第一条轨迹）与z轴交点的纵坐标  start_dist_2D / tan_theta 计算了沿z方向的距离

    /* Adapt for traceStackTwoWay reverse direction  kernel 是一个 TransportKernel 类型的对象，并且其方向为反向（逆向）。*/
    // NOTE If more applications for this arise, make 'reverse' an argument
    if (dynamic_cast<TransportKernel *>(kernel) &&
        !dynamic_cast<TransportKernel *>(kernel)->getDirection())
    {
      // 修改方位角角度和方向
      phi += M_PI;
      cos_phi *= -1;
      sin_phi *= -1;

      // 更新极角角度相关参数
      theta = M_PI - theta;
      cos_theta = cos(theta);
      sin_theta = sin(theta);
      tan_theta = sin_theta / cos_theta;
      sign = (cos_theta > 0) - (cos_theta < 0); // sign值为1时，代表轨迹沿Z轴正方向 为-1时代表负方向

      x_start_3D = first.getEnd()->getX();
      x_start_2D = flattened_track->getEnd()->getX();
      y_start_2D = flattened_track->getEnd()->getY();
      start_dist_2D = (x_start_3D - x_start_2D) / cos_phi;

      /* 计算第0条轨迹与z轴交点的纵坐标，即z0(0) */
      z0 = first.getEnd()->getZ();
      start_z = z0 - start_dist_2D / tan_theta;
    }

    /* Get the Geometry and CMFD mesh */
    Geometry *geometry = _track_generator_3D->getGeometry();
    Cmfd *cmfd = geometry->getCmfd();

    /* Extract the appropriate starting mesh */
    size_t num_fsrs = 0;
    double *axial_mesh;
    if (_global_z_mesh != nullptr)
    { // 如果有全局轴向网
      num_fsrs = _mesh_size;
      axial_mesh = _global_z_mesh;
    }

    /* Set the current x and y coordinates */
    double x_curr = x_start_2D;
    double y_curr = y_start_2D;

    /* Loop over 2D segments */
    double first_start_z = start_z;
    segment *segments_2D = flattened_track->getSegments();
    for (int s = 0; s < flattened_track->getNumSegments(); s++)
    { // 遍历2D线段，一次处理一个轴向挤出FSR内的线段

      /* Get segment length and extruded FSR */
      double seg_length_2D = segments_2D[s]._length;
      int extruded_fsr_id = segments_2D[s]._region_id; // 轴向挤出FSR_ID
      ExtrudedFSR *extruded_FSR = geometry->getExtrudedFSR(extruded_fsr_id);

      /* Determine new mesh and z index */
      if (_global_z_mesh == nullptr)
      {
        num_fsrs = extruded_FSR->_num_fsrs;
        axial_mesh = extruded_FSR->_mesh;
      }

      /* Calculate the end z coordinate of the first track */
      // 如果是钝角，正切为负
      double first_end_z = first_start_z + seg_length_2D / tan_theta; /* 锐角和钝角都适用 */

      /* 处理锐角和钝角的情况，即（z0,lower）, （z0,upper） */
      /* Find the upper and lower z coordinates of the first track 找到第一个轨迹的上下z坐标*/
      double first_track_lower_z;
      double first_track_upper_z;
      if (sign > 0)
      { // 锐角
        first_track_lower_z = first_start_z;
        first_track_upper_z = first_end_z; // 确保upper值一定是较大的z轴值
      }
      else
      {
        first_track_lower_z = first_end_z;
        first_track_upper_z = first_start_z;
      }

      /* Loop over all 3D FSRs in the Extruded FSR to find intersections */
      double first_seg_len_3D;

      /* 遍历当前轴向挤出FSR中的所有FSR */
      for (size_t z_iter = 0; z_iter < num_fsrs; z_iter++)
      {

        /* If traveling in negative-z direction, loop through FSRs from top */
        int z_ind = z_iter;
        /* 若为钝角，反向遍历，即从顶部开始*/
        if (sign < 0)
          z_ind = num_fsrs - z_iter - 1;

        /* Extract the FSR ID and Material ID of this 3D FSR */
        long fsr_id = extruded_FSR->_fsr_ids[z_ind];
        Material *material = extruded_FSR->_materials[z_ind];

        /* Find CMFD cell if necessary */
        int cmfd_cell;
        if (cmfd != nullptr)
          cmfd_cell = geometry->getCmfdCell(fsr_id); // 找到该轴向挤出FSR所在的CMFD网格（域中的局部一维ID）

        /* Get boundaries of the current mesh cell 获取当前轴向网格单元的Z边界*/
        double z_min = axial_mesh[z_ind];
        double z_max = axial_mesh[z_ind + 1];

        /* Calculate the local x and y centroid of the Extruded FSR 计算当前轴向挤出FSR的局部x和y质心*/
        double fsr_x_start = 0;
        double fsr_y_start = 0;
        double z_cent = 0.0;
        if (geometry->containsFSRCentroids())
        {
          Point *centroid = geometry->getFSRCentroid(fsr_id);
          fsr_x_start = x_curr - centroid->getX(); // 全局坐标转换为相对于FSR质心的局部坐标,相对于FSR质心的局部 x 坐标
          fsr_y_start = y_curr - centroid->getY();
          z_cent = geometry->getFSRCentroid(fsr_id)->getZ(); // 得到该轴向挤出FSR的Z轴质心坐标点的值
        }

        /* Calculate z-stack track indexes that cross the 3D FSR */
        int start_track = std::ceil((z_min - first_track_upper_z) / z_spacing); // 第一条轨迹的索引
        int start_full = std::ceil((z_min - first_track_lower_z) / z_spacing);  // 第一条H-track，也就是横向穿过整个FSR的轨迹
        int end_full = std::ceil((z_max - first_track_upper_z) / z_spacing);    // 最后一条H-track的后一条。在有V-track的FSR中，也是第一条V-track
        int end_track = std::ceil((z_max - first_track_lower_z) / z_spacing);   // 最后一条轨迹的索引再加1

        /* Check track bounds */
        start_track = std::max(start_track, 0);
        end_track = std::min(end_track, num_z_stack);

        /* Treat lower tracks that do not cross the entire 2D length   M-tracks case A(istart ≤ i < min(start_full, end_full))*/
        int min_lower = std::min(start_full, end_full);
        first_seg_len_3D = (first_track_upper_z - z_min) / std::abs(cos_theta); // 算了一部分，没算完
        for (int i = start_track; i < min_lower; i++)
        { // i为轨迹索引

          /* Calculate distance traveled in 3D FSR */
          double seg_len_3D = first_seg_len_3D + i * track_spacing_3D; // 和文档对应起来，算出穿过FSR的长度 P178

          /* Determine if segment length is large enough to operate on 太小的不必处理*/
          if (seg_len_3D > TINY_MOVE)
          {

            /* Initialize CMFD surfaces to none (-1) */
            int cmfd_surface_fwd = -1;
            int cmfd_surface_bwd = -1;

            /* Get CMFD surface if necessary */
            double lower_z = first_track_lower_z + i * z_spacing; // 找到该轨迹穿过轴向挤出FSR的上下高度
            double upper_z = first_track_upper_z + i * z_spacing;
            /*
            如果 dist_to_corner 小于或等于某个极小值 TINY_MOVE，表示轨迹段非常接近某个二维段的拐角。
            这种情况下，可能需要特殊处理，因为在拐角处的物理现象（如反射和折射）会比较复杂，需要更精确的处理
             */
            double dist_to_corner = std::abs((z_min - lower_z) / cos_theta); // 未进入FSR的一段长度
            if (cmfd != nullptr)
            { // CMFD
              if (sign > 0)
              {                                                      // 锐角
                cmfd_surface_fwd = segments_2D[s]._cmfd_surface_fwd; // 二维CMFD表面

                cmfd_surface_fwd = cmfd->findCmfdSurfaceOTF(cmfd_cell, upper_z,
                                                            cmfd_surface_fwd); // 求轨迹穿过轴向挤出FSR的upper_z点所在的CMFD网格面的ID
                if (dist_to_corner <= TINY_MOVE)                               // 如果z_min与lower_z重合
                  cmfd_surface_bwd = segments_2D[s]._cmfd_surface_bwd;
                cmfd_surface_bwd = cmfd->findCmfdSurfaceOTF(cmfd_cell, z_min,
                                                            cmfd_surface_bwd); // 求轨迹穿过轴向挤出FSR的z_min点所在的CMFD网格面的ID
                else
                { // 钝角
                  if (dist_to_corner <= TINY_MOVE)
                    cmfd_surface_fwd = segments_2D[s]._cmfd_surface_fwd;
                  cmfd_surface_fwd = cmfd->findCmfdSurfaceOTF(cmfd_cell, z_min,
                                                              cmfd_surface_fwd);
                  cmfd_surface_bwd = segments_2D[s]._cmfd_surface_bwd;
                  cmfd_surface_bwd = cmfd->findCmfdSurfaceOTF(cmfd_cell, upper_z,
                                                              cmfd_surface_bwd);
                }
              }

              /* Calculate the entry point of the segment into the FSR 计算该段进入FSR的入口点*/
              double x_entry = fsr_x_start;
              double y_entry = fsr_y_start;
              double z_entry = 0;
              if (sign > 0)
              { // 锐角
                double partial_2D = dist_to_corner * sin_theta;
                x_entry += partial_2D * cos_phi;
                y_entry += partial_2D * sin_phi;
                z_entry = z_min - z_cent; // 获取相对Z坐标
              }
              else
              {
                z_entry = upper_z - z_cent;
              }

              /* Operate on segment */
              kernel->execute(seg_len_3D, material, fsr_id, i,
                              cmfd_surface_fwd, cmfd_surface_bwd,
                              x_entry, y_entry, z_entry, phi, theta); // 涉及CMFD网格面中子流计算
            }
          }

          /* Find if there are tracks that traverse the entire 2D length */
          if (end_full > start_full)
          { // 处理H-tracks的情况，此时没有 V-tracks

            /* Calculate distance traveled in 3D FSR */
            double seg_len_3D = seg_length_2D / sin_theta;

            /* Determine if segment length is large enough to operate on */
            if (seg_len_3D > TINY_MOVE)
            {

              /* Treat tracks that do cross the entire 2D length */
              for (int i = start_full; i < end_full; i++)
              {

                /* Initialize CMFD surfaces to 2D CMFD surfaces */
                int cmfd_surface_fwd = segments_2D[s]._cmfd_surface_fwd;
                int cmfd_surface_bwd = segments_2D[s]._cmfd_surface_bwd;

                /* Get CMFD surfaces if necessary */
                double start_z = first_start_z + i * z_spacing;
                if (cmfd != nullptr)
                {

                  /* Calculate start and end z */
                  double end_z = first_end_z + i * z_spacing;
                  cmfd_surface_fwd = cmfd->findCmfdSurfaceOTF(cmfd_cell, end_z,
                                                              cmfd_surface_fwd);
                  cmfd_surface_bwd = cmfd->findCmfdSurfaceOTF(cmfd_cell, start_z,
                                                              cmfd_surface_bwd);
                }

                /* Calculate the entry point of the segment into the FSR */
                double x_entry = fsr_x_start;
                double y_entry = fsr_y_start;
                double z_entry = start_z - z_cent;

                /* Operate on segment */
                kernel->execute(seg_len_3D, material, fsr_id, i,
                                cmfd_surface_fwd, cmfd_surface_bwd,
                                x_entry, y_entry, z_entry, phi, theta); // 涉及CMFD网格面中子流计算
              }
            }
          }

          /* Find if there are tracks that cross both upper and lower boundaries
             NOTE: this will only be true if there are no tracks that cross the
             entire 2D length in the FSR */
          else if (start_full > end_full)
          { // V-tracks的情况，此时没有 H-tracks

            /* Calculate distance traveled in 3D FSR */
            double seg_len_3D = (z_max - z_min) / std::abs(cos_theta);

            /* Determine if segment length is large enough to operate on */
            if (seg_len_3D > TINY_MOVE)
            {

              /* Treat tracks that cross through both the upper and lower axial
                 boundaries */
              for (int i = end_full; i < start_full; i++)
              {

                /* Initialize CMFD surfaces to none (-1) */
                int cmfd_surface_bwd = -1;
                int cmfd_surface_fwd = -1;

                /* Determine start and end z */
                double enter_z;
                double exit_z;
                if (sign > 0)
                {
                  enter_z = z_min;
                  exit_z = z_max;
                }
                else
                {
                  enter_z = z_max;
                  exit_z = z_min;
                }

                /* Get CMFD surfaces if necessary */
                double track_start_z = first_start_z + i * z_spacing;
                double dist_to_corner_bwd = (enter_z - track_start_z) / cos_theta;
                if (cmfd != nullptr)
                {

                  /* Determine if any corners in the s-z plane are hit */
                  if (dist_to_corner_bwd <= TINY_MOVE)
                    cmfd_surface_bwd = segments_2D[s]._cmfd_surface_bwd;

                  double track_end_z = first_end_z + i * z_spacing;
                  double dist_to_corner_fwd = (track_end_z - exit_z) / cos_theta;
                  if (dist_to_corner_fwd <= TINY_MOVE)
                    cmfd_surface_fwd = segments_2D[s]._cmfd_surface_fwd;

                  /* Find CMFD surfaces */
                  cmfd_surface_fwd = cmfd->findCmfdSurfaceOTF(cmfd_cell, exit_z,
                                                              cmfd_surface_fwd);
                  cmfd_surface_bwd = cmfd->findCmfdSurfaceOTF(cmfd_cell, enter_z,
                                                              cmfd_surface_bwd);
                }

                /* Calculate the entry point of the segment into the FSR */
                double partial_2D = dist_to_corner_bwd * sin_theta;
                double x_entry = fsr_x_start + partial_2D * cos_phi;
                double y_entry = fsr_y_start + partial_2D * sin_phi;
                double z_entry = enter_z - z_cent;

                /* Operate on segment */
                kernel->execute(seg_len_3D, material, fsr_id, i,
                                cmfd_surface_fwd, cmfd_surface_bwd,
                                x_entry, y_entry, z_entry, phi, theta);
              }
            }
          }

          /* Treat upper tracks that do not cross the entire 2D length   处理M-tracks case B*/
          int min_upper = std::max(start_full, end_full);
          first_seg_len_3D = (z_max - first_track_lower_z) / std::abs(cos_theta);
          for (int i = min_upper; i < end_track; i++)
          {

            /* Calculate distance traveled in 3D FSR */
            double seg_len_3D = first_seg_len_3D - i * track_spacing_3D;

            /* Determine if segment length is large enough to operate on */
            if (seg_len_3D > TINY_MOVE)
            {

              /* Initialize CMFD surfaces to none (-1) */
              int cmfd_surface_fwd = -1;
              int cmfd_surface_bwd = -1;

              /* Get CMFD surface if necessary */
              double lower_z = first_track_lower_z + i * z_spacing;
              double upper_z = first_track_upper_z + i * z_spacing;
              double dist_to_corner = (upper_z - z_max) / std::abs(cos_theta);
              if (cmfd != nullptr)
              {
                if (sign > 0)
                {
                  if (dist_to_corner <= TINY_MOVE)
                    cmfd_surface_fwd = segments_2D[s]._cmfd_surface_fwd;
                  cmfd_surface_fwd = cmfd->findCmfdSurfaceOTF(cmfd_cell, z_max,
                                                              cmfd_surface_fwd);
                  cmfd_surface_bwd = segments_2D[s]._cmfd_surface_bwd;
                  cmfd_surface_bwd = cmfd->findCmfdSurfaceOTF(cmfd_cell, lower_z,
                                                              cmfd_surface_bwd);
                }
                else
                {
                  cmfd_surface_fwd = segments_2D[s]._cmfd_surface_fwd;
                  cmfd_surface_fwd = cmfd->findCmfdSurfaceOTF(cmfd_cell, lower_z,
                                                              cmfd_surface_fwd);
                  if (dist_to_corner <= TINY_MOVE)
                    cmfd_surface_bwd = segments_2D[s]._cmfd_surface_bwd;
                  cmfd_surface_bwd = cmfd->findCmfdSurfaceOTF(cmfd_cell, z_max,
                                                              cmfd_surface_bwd);
                }
              }

              /* Calculate the entry point of the segment into the FSR */
              double x_entry = fsr_x_start;
              double y_entry = fsr_y_start;
              double z_entry = 0;
              if (sign < 0)
              {
                double partial_2D = dist_to_corner * sin_theta;
                x_entry += partial_2D * cos_phi;
                y_entry += partial_2D * sin_phi;
                z_entry = z_max - z_cent;
              }
              else
              {
                z_entry = lower_z - z_cent;
              }

              /* Operate on segment */
              kernel->execute(seg_len_3D, material, fsr_id, i,
                              cmfd_surface_fwd, cmfd_surface_bwd,
                              x_entry, y_entry, z_entry, phi, theta);
            }
          }
        }
        /* Traverse segment on first track */
        first_start_z = first_end_z;
        x_curr += seg_length_2D * cos_phi;
        y_curr += seg_length_2D * sin_phi;
      }
    }

    /**
     * @brief A function that searches for the index into a values mesh using a
     *        binary search.
     * @details A binary search is used to calculate the index into a mesh of where
     *          the value val resides. If a mesh boundary is hit, the upper region
     *          is selected for positive-z traversing rays and the lower region is
     *          selected for negative-z traversing rays.
     * @param values an array of monotonically increasing values
     * @param size the size of the values array
     * @param val the level to be searched for in the mesh
     * @param sign the direction of the ray in the z-direction
     */
    int TraverseSegments::findMeshIndex(double *values, int size,
                                        double val, int sign)
    {

      /* Initialize indexes into the values array */
      int imin = 0;
      int imax = size - 1;

      /* Check if val is outside the range */
      if (val < values[imin] || val > values[imax])
      {
        log::ferror("Value out of the mesh range in binary search");
        return -1;
      }

      /* Search for interval containing val */
      while (imax - imin > 1)
      {

        int imid = (imin + imax) / 2;

        if (val > values[imid])
          imin = imid;
        else if (val < values[imid])
          imax = imid;
        else
        {
          if (sign > 0)
            return imid;
          else
            return imid - 1;
        }
      }
      return imin;
    }

    /**
     * @brief Loops over all 3D Tracks using axial on-the-fly ray tracking by
     *        z-stack, going forward then backward on each 3D Track.
     * @details The onTrack(...) function is applied to all 3D Tracks and the
     *          specified kernel is applied to all segments. If nullptr is provided
     *          for the kernel, only the onTrack(...) functionality is applied.
     * @param kernel The TransportKernel dictating the functionality to apply to
     *        segments
     */
    void TraverseSegments::loopOverTracksByStackTwoWay(TransportKernel * kernel)
    { // 使用轨迹堆Z-stack 遍历所有3D轨迹

      if (_segment_formation != +segmentationType::OTF_STACKS)
        log::ferror("Two way on-the-fly transport has only been implemented "
                    "for ray tracing by z-stack");

      int num_2D_tracks = _track_generator_3D->getMyNum2DTracks(); // 2D轨迹数量
      Track **flattened_tracks = _track_generator_3D->get2DTracksArray();
      // unused
      // int*** tracks_per_stack = _track_generator_3D->getTracksPerStack();
      int num_polar = _track_generator_3D->getNumPolar();
      int tid = omp_get_thread_num();

      /* Loop over flattened 2D tracks */
#pragma omp for schedule(dynamic)
      for (int ext_id = 0; ext_id < num_2D_tracks; ext_id++)
      {

        /* Extract indices of 3D tracks associated with the flattened track 提取与二维轨迹关联的三维轨迹的索引*/
        TrackStackIndexes tsi;
        Track *flattened_track = flattened_tracks[ext_id];
        tsi._azim = flattened_track->getAzimIndex();
        tsi._xy = flattened_track->getXYIndex();

        /* Loop over polar angles */
        for (int p = 0; p < num_polar; p++)
        {

          /* Retrieve information for the first 3D Track in the z-stack */
          tsi._polar = p;
          tsi._z = 0;
          Track3D track_3D;
          _track_generator_3D->getTrackOTF(&track_3D, &tsi); // 即时生成3D轨迹相关信息,实例化一个轨迹堆&track_3D

          if (kernel != nullptr)
          {

            /* Reset kernel for a new base Track */
            kernel->newTrack(&track_3D);

            /* Trace all segments in the z-stack 遍历一个确定极角的轨迹中的所有分段*/
            // 建立一个 2D 线段的临时数组，先使用按轨迹堆的追踪沿正向追踪（极角为锐角），随后把线段反过来再追踪，再把线段反回去
            traceStackTwoWay(flattened_track, p, kernel); // 追踪一个2D轨迹中的所有分段
            track_3D.setNumSegments(kernel->getCount());
          }

          /* Operate on the Track */
          segment *segments = _track_generator_3D->getTemporarySegments(tid);
          onTrack(&track_3D, segments); // 遍历轨迹上的线段，应用 MOC 方程
        }
      }
    }

    /**
     * @brief Traces the 3D segments of 3D Tracks in a z-stack both forward and
     *        backward across the geometry, applying the kernel provided by the
     *        user when the segment information is calculated.
     * @details This function copies information of the 3D z-stack, ray traces the
     *          z-stack forward using TrackGenerator::traceStackOTF, then reverses
     *          the tracks so that they point backwards, and ray traces in the
     *          reverse direction. This allows segments to be applied to
     *          TransportKernels during the on-the-fly ray tracing process.
     * @param flattened_track the 2D track associated with the z-stack for which
     *        3D segments are computed
     * @param polar_index the polar index of the 3D Track z-stack
     * @param kernel The TransportKernel applied to the calculated 3D segments
     *
     * @brief 在几何体中正向和反向追踪z轨迹堆中3D轨迹的3D分段，
     *        在计算分段信息时应用用户提供的核函数。
     * @details 此函数复制3D z轨迹堆的信息，使用TrackGenerator::traceStackOTF
     *          正向射线追踪z轨迹堆，然后反转轨迹使其指向反方向，并在反方向
     *          进行射线追踪。这允许在即时(on-the-fly)射线追踪过程中将分段
     *          应用于传输核函数。
     * @param flattened_track 与计算3D分段的z轨迹堆相关联的2D轨迹
     * @param polar_index 3D轨迹z轨迹堆的极角索引
     * @param kernel 应用于计算出的3D分段的传输核函数
     */
    void TraverseSegments::traceStackTwoWay(Track * flattened_track, int polar_index,
                                            TransportKernel *kernel)
    {

      /* Get segments from flattened track 取片段*/
      segment *segments = flattened_track->getSegments();
      MOCKernel *moc_kernel = dynamic_cast<MOCKernel *>(kernel);

      /* Trace stack forwards 设置轨迹方向为forward方向*/
      kernel->setDirection(true);
      // 注：z-stacks 指的是对应于某条 2D 轨迹的，某个极角的所有轴向轨迹，也称为轨迹堆
      traceStackOTF(flattened_track, polar_index, moc_kernel); // 使用 OTF 方法追踪整个轨迹堆中的轨迹
      kernel->post();

      /* Reverse segments in flattened track */
      /* 反转二维轨迹中的分段 */
      int num_segments = flattened_track->getNumSegments();
      for (int s = 0; s < num_segments / 2; s++)
      {
        segment tmp_segment = segments[num_segments - s - 1];
        segments[num_segments - s - 1] = segments[s];
        segments[s] = tmp_segment;
      }

      /* Flip CMFD surfaces on segments in flattened track 在二维轨道的分段上翻转CMFD surface*/
      /* 翻转二维轨迹分段上的CMFD表面 */
      for (int s = 0; s < num_segments; s++)
      {
        int tmp_surface = segments[s]._cmfd_surface_fwd;
        segments[s]._cmfd_surface_fwd = segments[s]._cmfd_surface_bwd;
        segments[s]._cmfd_surface_bwd = tmp_surface;
      }

      /* Trace stack backwards */
      /* 反向追踪轨迹堆 */
      kernel->setDirection(false);
      traceStackOTF(flattened_track, polar_index, moc_kernel);
      // 使用相同的OTF方法，但方向设置为反向
      kernel->post();

      /* Reverse segments in flattened track */
      /* 再次反转二维轨迹中的分段 将分段顺序恢复到原始状态 */
      for (int s = 0; s < num_segments / 2; s++)
      {
        segment tmp_segment = segments[num_segments - s - 1];
        segments[num_segments - s - 1] = segments[s];
        segments[s] = tmp_segment;
      }

      /* Flip CMFD surfaces on segments in flattened track 再反转回来*/
      for (int s = 0; s < num_segments; s++)
      {
        int tmp_surface = segments[s]._cmfd_surface_fwd;
        segments[s]._cmfd_surface_fwd = segments[s]._cmfd_surface_bwd;
        segments[s]._cmfd_surface_bwd = tmp_surface;
      }
    }

  } /* namespace antmoc */
