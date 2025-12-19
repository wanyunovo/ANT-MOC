/// \file   FSRDataTypes.h
/// \brief  FSR data structures.
/// \date   April 18, 2019
/// \author William Boyd, MIT, Course 22 (wboyd@mit.edu)
/// \author An Wang (wangan@xs.ustb.edu.cn)

#ifndef FSR_DATA_TYPES_H_
#define FSR_DATA_TYPES_H_

#include <cstddef>

namespace antmoc
{

  /** Forward declarations */
  class Material;
  class LocalCoords;
  class Point;

  /// \struct FSRData
  /// \brief  A FSRData struct represents an FSR with a unique FSR ID
  ///         and a characteristic point that lies within the FSR that
  ///         can be used to recompute the hierarchical LocalCoords
  ///         linked list.
  struct FSRData
  {

    /** The FSR ID */
    long _fsr_id;

    /** The CMFD Cell */
    int _cmfd_cell; // 所在的cmfd cell

    /** The Material ID */
    int _mat_id; // 包含的材料的编号

    /** Characteristic point in Root Universe that lies in FSR */
    Point *_point; // FSR中的一个特征点,坐标相对于Root Universe

    /** Global numerical centroid in Root Universe */
    Point *_centroid; // 数值中心点,坐标相对于Root Universe

    /// \brief Constructor
    FSRData();

    /// \brief Destructor
    ~FSRData();
  };

  /// \struct ExtrudedFSR
  /// \brief An ExtrudedFSR struct represents a FSR region in the superposition
  ///        plane for axial on-the-fly ray tracing. It contains a characteristic
  ///        point that lies within the FSR, an axial mesh, and an array of 3D
  ///        FSR IDs contained within the extruded region along with their
  ///        corresponding materials.
  struct ExtrudedFSR
  {

    /** Array defining the axial mesh */
    double *_mesh; // 轴向网

    /** Axial extruded FSR ID */
    int _fsr_id; // FSR的编号,按轨迹分段的顺序编号

    /** Array of 3D FSR IDs */
    long *_fsr_ids; // 包含的3D FSR的编号

    /** Array of material pointers for each FSR */
    Material **_materials; // 每个3D FSR的材料指针

    /** Number of FSRs in the axially extruded FSR */
    int _num_fsrs; // 这个轴向挤出FSR中的FSR数量

    /** Coordinates inside the FSR */
    LocalCoords *_coords; // ExtrudedFSR中的一个点,创建时分配

    /*
            extruded_FSR->_num_fsrs = num_regions;
            extruded_FSR->_materials = new Material *[num_regions];
            extruded_FSR->_fsr_ids = new long[num_regions];
    */
    /// \brief Constructor
    ExtrudedFSR();

    /// \brief Destructor
    ~ExtrudedFSR();
  };

} /* namespace antmoc */

#endif /* FSR_DATA_TYPES_H_ */
