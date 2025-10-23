#include "antmoc/tally_utils.h"
#include "antmoc/string_utils.h"

#include <sstream>

namespace antmoc {

namespace tallyutils {

  /// \brief Get the name of the tally type
  /// \details This function is intended to convert a user-defined
  ///          enumeration into a string
  std::string getTallyTypeName(TallyType t) {
    std::string name = t._to_string();
    return stringutils::underscoreToSpace(name);
  }


  /// \brief Determine if a given type is a TallyType of XS
  bool isXSTallyType(TallyType t) {
    bool result;
    switch (t) {
      case TallyType::Fission_XS:
      case TallyType::NuFission_XS:
      case TallyType::Total_XS:
      case TallyType::Absorption_XS:
      case TallyType::All:
      case TallyType::None:
        result = true;
        break;

      default:
        result = false;
    }
    return result;
  }


  /// \brief Determine if a given type is a TallyType of reaction rates
  bool isRXTallyType(TallyType t) {
    bool result;
    switch (t) {
      case TallyType::Fission_RX:
      case TallyType::NuFission_RX:
      case TallyType::Total_RX:
      case TallyType::Absorption_RX:
      case TallyType::Scalar_Flux:
      case TallyType::All:
      case TallyType::None:
        result = true;
        break;

      default:
        result = false;
    }
    return result;
  }


  /// \brief Determine if a given type is a field type
  bool isFieldTallyType(TallyType type) {
    bool result = isRXTallyType(type) ||
                  isXSTallyType(type) ||
                  (type == +TallyType::Volume);
    return result;
  }


  /// \brief Determine if a given type is a TallyType of track data
  bool isTracksTallyType(TallyType t) {
    bool result;
    switch (t) {
      case TallyType::Tracks_2D:
      case TallyType::Tracks_3D:
      case TallyType::All:
      case TallyType::None:
        result = true;
        break;

      default:
        result = false;
    }
    return result;
  }


  /// \brief Simply remove None from the TallyType set
  void removeTallyTypeNone(std::set<TallyType> &types) {
    if (types.find(TallyType::None) != types.end()) {
      types.erase(TallyType::None);
    }
  }


  /// \brief Simply remove All from the TallyType set
  void removeTallyTypeAll(std::set<TallyType> &types) {
    if (types.find(TallyType::All) != types.end()) {
      types.erase(TallyType::All);
    }
  }


  /// \brief Remove TallyType::None and expand TallyType::All
  /// \details TallyType::All will be expanded to RX TallyTypes
  void expandTallyTypeForRX(std::set<TallyType> &types) {
    // If no TallyType is specified, dump nothing
    removeTallyTypeNone(types);

    if (types.find(TallyType::All) != types.end()) {
      types.clear();
      types.insert(TallyType::Fission_RX);
      types.insert(TallyType::NuFission_RX);
      types.insert(TallyType::Absorption_RX);
      types.insert(TallyType::Total_RX);
      types.insert(TallyType::Scalar_Flux);
    }
  }


  /// \brief Remove TallyType::None and expand TallyType::All
  /// \details TallyType::All will be expanded to XS TallyTypes
  void expandTallyTypeForXS(std::set<TallyType> &types) {
    // If no TallyType is specified, dump nothing
    removeTallyTypeNone(types);

    if (types.find(TallyType::All) != types.end()) {
      types.clear();
      types.insert(TallyType::Fission_XS);
      types.insert(TallyType::NuFission_XS);
      types.insert(TallyType::Absorption_XS);
      types.insert(TallyType::Total_XS);
    }
  }


  /// \brief Remove TallyType::None and expand TallyType::All
  /// \details TallyType::All will be expanded to Volume TallyTypes
  void expandTallyTypeForVolume(std::set<TallyType> &types) {
    // If no TallyType is specified, dump nothing
    removeTallyTypeNone(types);

    if (types.find(TallyType::All) != types.end()) {
      types.clear();
      types.insert(TallyType::Volume);
    }
  }


  /// \brief Remove TallyType::None and expand TallyType::All
  /// \details TallyType::All will be expanded to Track TallyTypes
  void expandTallyTypeForTrack(std::set<TallyType> &types) {
    // If no TallyType is specified, dump nothing
    removeTallyTypeNone(types);

    if (types.find(TallyType::All) != types.end()) {
      types.clear();
      types.insert(TallyType::Tracks_2D);
      types.insert(TallyType::Tracks_3D);
    }
  }


  /// \brief A helper method
  /// \details For example, "RX" and "PHI" will be combined to "RX_PHI"
  /// \param type a code of a tally type
  /// \param category RX, XS or VOL
  TallyType codeToTallyType(std::string type, std::string category) {
    std::unordered_map<std::string, TallyType> case_map {
      // reaction rates
      {"RX_F",   TallyType::Fission_RX},
      {"RX_NUF", TallyType::NuFission_RX},
      {"RX_A",   TallyType::Absorption_RX},
      {"RX_T",   TallyType::Total_RX},
      {"RX_PHI", TallyType::Scalar_Flux},

      // cross sections
      {"XS_F",   TallyType::Fission_XS},
      {"XS_NUF", TallyType::NuFission_XS},
      {"XS_A",   TallyType::Absorption_XS},
      {"XS_T",   TallyType::Total_XS},

      // tracks
      {"TRACKS_2D", TallyType::Tracks_2D},
      {"TRACKS_3D", TallyType::Tracks_3D},
    };

    type = stringutils::toUpper(type);
    category = stringutils::toUpper(category);
    std::string case_string = category + "_" + type;

    // Defaults to none
    TallyType obj = TallyType::None;

    if ( case_map.count(case_string) )
      obj = case_map[case_string];

    if (type == "ALL")
      obj = TallyType::All;

    return obj;
  }


  /// \brief Convert a tally type to its code
  std::string tallyTypeToCode(TallyType type) {

    switch (type) {
      // reaction rates
      case TallyType::Fission_RX:     return "F";
      case TallyType::NuFission_RX:   return "NuF";
      case TallyType::Absorption_RX:  return "A";
      case TallyType::Total_RX:       return "T";
      case TallyType::Scalar_Flux:    return "PHI";

      // cross sections
      case TallyType::Fission_XS:     return "F";
      case TallyType::NuFission_XS:   return "NuF";
      case TallyType::Absorption_XS:  return "A";
      case TallyType::Total_XS:       return "T";

      // tracks
      case TallyType::Tracks_2D:      return "2D";
      case TallyType::Tracks_3D:      return "3D";

      case TallyType::Volume:         return "VOL";

      case TallyType::All:            return "ALL";
      case TallyType::None:           return "NONE";
    }

    return "";
  }


  /// \brief Concatenate tally types
  std::string join(const std::set<TallyType> &types,
                   const std::string &delimiter) {
    std::stringstream ss;

    for (auto type : types)
      ss << tallyTypeToCode(type) << delimiter;

    auto s = ss.str();

    return s.substr(0, s.find_last_not_of(delimiter) + 1);
  }


  // Combine rx, xs and volume types
  std::set<TallyType> combineRXAndXS(std::set<TallyType> rx_types,
                                     std::set<TallyType> xs_types) {

    expandTallyTypeForRX(rx_types);
    expandTallyTypeForXS(xs_types);

    rx_types.insert(xs_types.begin(), xs_types.end());

    return rx_types;
  }

} // namespace tallyutils

} // namespace antmoc
