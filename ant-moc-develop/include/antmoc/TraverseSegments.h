/**
 * @file TraverseSegments.h
 * @brief A TraverseSegments object
 * @date February 15, 2016
 * @author Geoffrey Gunow, MIT, Course 22 (geogunow@mit.edu)
 */

#ifndef TRAVERSE_SEGMENTS_H_
#define TRAVERSE_SEGMENTS_H_
#include <typeinfo>

#include "antmoc/enum_types.h"
#include "antmoc/MOCKernel.h"

namespace antmoc
{

  /** Forward declarations */
  class Point;
  class Track;
  class TrackGenerator;
  class TrackGenerator3D;
  struct segment;

  /**
   * @class TraverseSegments TraverseSegments.h "src/TraverseSegments.h"
   * @brief TraverseSegments 对象定义了在给定各种不同分段方案的情况下如何遍历 Tracks，
   *        并如何将算法应用于 Tracks 和相关的段。
   * @details TraverseSegments 对象描述了如何针对各种不同的分段方案（如 2D 显式、3D 显式和
   *          动态光线追踪）遍历 Tracks。TraverseSegments 对象是一个抽象类，旨在由
   *          TrackTraversingAlgorithms.h 中定义的类扩展。这个父类的主要目的是抽象出遍历过程，
   *          并对每个 Track 应用 onTrack(...) 函数，同时对每个段应用提供的 MOCKernel。
   *          如果为 MOCKernel 提供了 nullptr，则仅对每个 Track 应用 onTrack(...) 中定义的功能。
   */
  class TraverseSegments
  {

  private:
    /* Functions defining how to loop over Tracks */
    void loopOverTracks2D(MOCKernel *kernel);
    void loopOverTracksExplicit(MOCKernel *kernel);
    void loopOverTracksByTrackOTF(MOCKernel *kernel);
    void loopOverTracksByStackOTF(MOCKernel *kernel);

    /* Functions defining how to traverse segments */
    void traceSegmentsExplicit(Track *track, MOCKernel *kernel);
    void traceSegmentsOTF(Track *flattened_track, Point *start,
                          double theta, MOCKernel *kernel);
    void traceStackOTF(Track *flattened_track, int polar_index,
                       MOCKernel *kernel);

    void traceStackTwoWay(Track *flattened_track, int polar_index,
                          TransportKernel *kernel);

    int findMeshIndex(double *values, int size, double val, int sign);

  protected:
    /** Pointer to the associated TrackGenerator */
    TrackGenerator *_track_generator;
    TrackGenerator3D *_track_generator_3D;

    /** A pointer to the associated global z-mesh (if applicable) */
    double *_global_z_mesh;

    /** The size of the global z-mesh */
    int _mesh_size;

    /** The type of segmentation used for segment formation */
    segmentationType _segment_formation;

    TraverseSegments(TrackGenerator *track_generator);
    virtual ~TraverseSegments();

    /* Functions defining how to loop over and operate on Tracks */
    void loopOverTracks(MOCKernel *kernel);
    virtual void onTrack(Track *track, segment *segments) = 0;

    // FIXME Rework function calls to make this private
    void loopOverTracksByStackTwoWay(TransportKernel *kernel);

    /* Returns a kernel of the requested type */
    template <class KernelType>
    MOCKernel *getKernel()
    {

      /* Check for segmentation kernel in explicit methods */
      if ((typeid(KernelType) != typeid(SegmentationKernel)) ||
          ((_segment_formation != +segmentationType::EXPLICIT_2D) &&
           (_segment_formation != +segmentationType::EXPLICIT_3D)))
      {

        /* Allocate kernel */
        MOCKernel *kernel = new KernelType(_track_generator);
        return kernel;
      }
      else
        return nullptr;
    }

  public:
    virtual void execute() = 0;
  };

} /* namespace antmoc */

#endif
