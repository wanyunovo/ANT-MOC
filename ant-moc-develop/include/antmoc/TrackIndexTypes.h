/// \file   TrackIndexTypes.h
/// \brief  Classes for indexing tracks.
/// \date   March 13, 2016
/// \author Geoffrey Gunow, MIT, Course 22 (geogunow@mit.edu)

#ifndef TRACK_INDEX_TYPES_H_
#define TRACK_INDEX_TYPES_H_

namespace antmoc {

/**
 * @struct TrackChainIndexes
 * @brief A Track chain represents a list of Tracks connected via periodic BCs
 *        that extend in the l-z plane and this struct contains the azim, x, polar,
 *        lz, and link indexes of a Track.
 */
struct TrackChainIndexes {

  /** The azimuthal index of the first link (in 0 to _num_azim / 4) */
  int _azim = -1;

  /** The x index (in 0 to _num_x[_azim]) */
  long _x = -1;

  /** The polar index (in 0 to _num_polar) */
  int _polar = -1;

  /** The lz index (in 0 to _num_l[_azim][_polar] + _num_z[_azim][_polar]) */
  int _lz = -1;

  /** The link index of the 3D chain */
  int _link = -1;
};


/**
 * @struct TrackStackIndexes
 * @brief A stack track represents a track in a z-stack of tracks and this
 *        struct contains the azim, xy, polar, and z indices of a track.
 */
struct TrackStackIndexes {

  /** The azimuthal index (in 0 to _num_azim / 2) */
  int _azim = -1;

  /** The xy index (in 0 to _num_x[_azim] + _num_y[_azim]) */
  long _xy = -1;

  /** The polar index (in 0 to _num_polar) */
  int _polar = -1;

  /** The z index in the z-stack (in 0 to _tracks_per_stack[_azim][_xy][_polar]) */
  int _z = -1;
};

} // namespace antmoc

#endif  // TRACK_INDEX_TYPES_H_
