/// \file tally_utils.h
/// \brief Utilities for tallying
/// \date Aug 26, 2019
/// \author An Wang, USTB (wangan@xs.ustb.edu.cn)

#ifndef TALLY_UTILS_H_
#define TALLY_UTILS_H_

#include <set>
#include <string>
#include <unordered_map>

#include "antmoc/enum_types.h"

namespace antmoc
{

/// \enum TallyType
/// \brief The type of tallied mesh data
BETTER_ENUM(TallyType, char,

  Fission_RX,     /**< Fission reaction rates */
  NuFission_RX,   /**< NuF reaction rates */
  Absorption_RX,  /**< Absorption reaction rates */
  Total_RX,       /**< Total reaction rates */
  Scalar_Flux,    /**< Scalar flux */

  Fission_XS,     /**< Fission cross-sections */
  NuFission_XS,   /**< NuF cross-sections */
  Absorption_XS,  /**< Absorption cross-sections */
  Total_XS,       /**< Total cross-sections */

  Volume,         /**< Volumes of mesh or FSRs */

  Tracks_2D,      /**< 2D tracks */
  Tracks_3D,      /**< 3D tracks */

  All,            /**< All types */
  None            /**< Nothing */

)

//----------------------------------------------------------------------
/// \namespace tallyutils
/// \brief Common funtions for tallying
//----------------------------------------------------------------------
namespace tallyutils {

  std::string getTallyTypeName(TallyType);
  TallyType codeToTallyType(std::string, std::string);
  std::string tallyTypeToCode(TallyType type);
  std::string join(const std::set<TallyType> &types, const std::string &delimiter);

  bool isXSTallyType(TallyType);
  bool isRXTallyType(TallyType);
  bool isFieldTallyType(TallyType);
  bool isTracksTallyType(TallyType);

  void removeTallyTypeNone(std::set<TallyType> &types);
  void removeTallyTypeAll(std::set<TallyType> &types);
  void expandTallyTypeForRX(std::set<TallyType> &types);
  void expandTallyTypeForXS(std::set<TallyType> &types);
  void expandTallyTypeForVolume(std::set<TallyType> &types);
  void expandTallyTypeForTrack(std::set<TallyType> &types);

  std::set<TallyType> combineRXAndXS(std::set<TallyType>, std::set<TallyType>);

} // namespace tallyutils

} // namespace antmoc

#endif  // TALLY_UTILS_H_
