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
struct FSRData {

  /** The FSR ID */
  long _fsr_id;

  /** The CMFD Cell */
  int _cmfd_cell;

  /** The Material ID */
  int _mat_id;

  /** Characteristic point in Root Universe that lies in FSR 一个局部点链表对应了一个 FSR（它也就是 FSR 的特征点），一个 FSR 被这个局部点
   * 链表的 CSG 层次信息唯一标识（CSG 层次的交集）。*/
  Point* _point;

  /** Global numerical centroid in Root Universe */
  Point* _centroid;

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
struct ExtrudedFSR {

  /** Array defining the axial mesh */
  double* _mesh;

  /** Axial extruded FSR ID */
  int _fsr_id;

  /** Array of 3D FSR IDs */
  long* _fsr_ids;

  /** Array of material pointers for each FSR */
  Material** _materials;

  /** Number of FSRs in the axially extruded FSR */
  int _num_fsrs;

  /** Coordinates inside the FSR */
  LocalCoords* _coords;

  /// \brief Constructor
  ExtrudedFSR();

  /// \brief Destructor
  ~ExtrudedFSR();
};


} /* namespace antmoc */

#endif /* FSR_DATA_TYPES_H_ */
