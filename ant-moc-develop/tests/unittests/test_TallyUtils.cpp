/// \file test_TallyUtils.cpp
/// \brief Test tally utility
/// \date March 26, 2019
/// \author An Wang, USTB (wangan@xs.ustb.edu.cn)

#include "testing/test_utils.h"
#include "antmoc/tally_utils.h"
#include "antmoc/string_utils.h"

using namespace antmoc;

namespace {

/// Testing fixture
class test_TallyUtils : public testing::Test {

protected:

  const std::vector<TallyType> all_types = {
    TallyType::Fission_RX,
    TallyType::NuFission_RX,
    TallyType::Absorption_RX,
    TallyType::Total_RX,
    TallyType::Scalar_Flux,

    TallyType::Fission_XS,
    TallyType::NuFission_XS,
    TallyType::Absorption_XS,
    TallyType::Total_XS,

    TallyType::Volume,
    TallyType::Tracks_2D,
    TallyType::Tracks_3D,

    TallyType::All,
    TallyType::None
  };

  const StringVec all_names = {
    "Fission RX",
    "NuFission RX",
    "Absorption RX",
    "Total RX",
    "Scalar Flux",

    "Fission XS",
    "NuFission XS",
    "Absorption XS",
    "Total XS",

    "Volume",
    "Tracks 2D",
    "Tracks 3D",

    "All",
    "None"
  };

};


TEST_F(test_TallyUtils, getTallyTypeName) {
  ROOT_ONLY();

  auto &tests = all_types;
  auto &oracles = all_names;

  for (size_t i = 0; i < tests.size(); ++i) {
    EXPECT_EQ(oracles[i], tallyutils::getTallyTypeName(tests[i]));
  }
}


TEST_F(test_TallyUtils, isRXTallyType) {
  ROOT_ONLY();

  auto &tests = all_types;

  const std::vector<bool> oracles = {
    true,
    true,
    true,
    true,
    true,

    false,
    false,
    false,
    false,

    false,
    false,
    false,

    true,
    true
  };

  for (size_t i = 0; i < tests.size(); ++i) {
    EXPECT_EQ(oracles[i], tallyutils::isRXTallyType(tests[i]));
  }
}


TEST_F(test_TallyUtils, isXSTallyType) {
  ROOT_ONLY();

  auto &tests = all_types;

  const std::vector<bool> oracles = {
    false,
    false,
    false,
    false,
    false,

    true,
    true,
    true,
    true,

    false,
    false,
    false,

    true,
    true
  };

  for (size_t i = 0; i < tests.size(); ++i) {
    EXPECT_EQ(oracles[i], tallyutils::isXSTallyType(tests[i]));
  }
}


TEST_F(test_TallyUtils, isFieldTallyType) {
  ROOT_ONLY();

  auto &tests = all_types;

  const std::vector<bool> oracles = {
    true,
    true,
    true,
    true,
    true,

    true,
    true,
    true,
    true,

    true,
    false,
    false,

    true,
    true
  };

  for (size_t i = 0; i < tests.size(); ++i) {
    EXPECT_EQ(oracles[i], tallyutils::isFieldTallyType(tests[i]));
  }
}


TEST_F(test_TallyUtils, isTracksTallyType) {
  ROOT_ONLY();

  auto &tests = all_types;

  const std::vector<bool> oracles = {
    false,
    false,
    false,
    false,
    false,

    false,
    false,
    false,
    false,

    false,
    true,
    true,

    true,
    true
  };

  for (size_t i = 0; i < tests.size(); ++i) {
    EXPECT_EQ(oracles[i], tallyutils::isTracksTallyType(tests[i]));
  }
}


TEST_F(test_TallyUtils, removeTallyTypeNone) {
  ROOT_ONLY();

  std::set<TallyType> types = {
    TallyType::None,
    TallyType::Fission_RX,
    TallyType::Fission_XS
  };

  const std::set<TallyType> oracles = {
    TallyType::Fission_RX,
    TallyType::Fission_XS
  };

  tallyutils::removeTallyTypeNone(types);
  EXPECT_EQ(types, oracles);
}


TEST_F(test_TallyUtils, removeTallyTypeAll) {
  ROOT_ONLY();

  std::set<TallyType> types = {
    TallyType::All,
    TallyType::Fission_RX,
    TallyType::Fission_XS
  };

  const std::set<TallyType> oracles = {
    TallyType::Fission_RX,
    TallyType::Fission_XS
  };

  tallyutils::removeTallyTypeAll(types);
  EXPECT_EQ(types, oracles);
}


TEST_F(test_TallyUtils, expandTallyTypeForRX) {
  ROOT_ONLY();

  std::set<TallyType> types = {
    TallyType::All,
    TallyType::Fission_RX
  };

  const std::set<TallyType> oracles = {
    TallyType::Fission_RX,
    TallyType::NuFission_RX,
    TallyType::Absorption_RX,
    TallyType::Total_RX,
    TallyType::Scalar_Flux
  };

  tallyutils::expandTallyTypeForRX(types);
  EXPECT_EQ(types, oracles);
}


TEST_F(test_TallyUtils, expandTallyTypeForXS) {
  ROOT_ONLY();

  std::set<TallyType> types = {
    TallyType::All,
    TallyType::Fission_XS
  };

  const std::set<TallyType> oracles = {
    TallyType::Fission_XS,
    TallyType::NuFission_XS,
    TallyType::Absorption_XS,
    TallyType::Total_XS
  };

  tallyutils::expandTallyTypeForXS(types);
  EXPECT_EQ(types, oracles);
}


TEST_F(test_TallyUtils, expandTallyTypeForVolume) {
  ROOT_ONLY();

  std::set<TallyType> types = {
    TallyType::All
  };

  const std::set<TallyType> oracles = {
    TallyType::Volume
  };

  tallyutils::expandTallyTypeForVolume(types);
  EXPECT_EQ(types, oracles);
}


TEST_F(test_TallyUtils, expandTallyTypeForTrack) {
  ROOT_ONLY();

  std::set<TallyType> types = {
    TallyType::All,
    TallyType::Fission_XS
  };

  const std::set<TallyType> oracles = {
    TallyType::Tracks_2D,
    TallyType::Tracks_3D
  };

  tallyutils::expandTallyTypeForTrack(types);
  EXPECT_EQ(types, oracles);
}


/// \brief Testing shorthands of tally types
/// \details This function is used to handle user input
TEST_F(test_TallyUtils, codeToTallyType) {
  ROOT_ONLY();

  std::vector<std::pair<std::string, std::string>> tests = {
    {"f", "rx"},
    {"NUF", "rx"},
    {"a", "RX"},
    {"t", "rx"},
    {"phi", "rx"},

    {"f", "xs"},
    {"nuf", "xs"},
    {"a", "xs"},
    {"t", "xs"},

    {"2d", "tracks"},
    {"3d", "tracks"},
  };

  const std::vector<TallyType> oracles = {
    TallyType::Fission_RX,
    TallyType::NuFission_RX,
    TallyType::Absorption_RX,
    TallyType::Total_RX,
    TallyType::Scalar_Flux,

    TallyType::Fission_XS,
    TallyType::NuFission_XS,
    TallyType::Absorption_XS,
    TallyType::Total_XS,

    TallyType::Tracks_2D,
    TallyType::Tracks_3D
  };

  for (size_t i = 0; i < tests.size(); ++i) {
    auto &test = tests[i];
    auto &oracle = oracles[i];
    EXPECT_EQ(oracle, tallyutils::codeToTallyType(test.first, test.second));
  }
}


/// \brief Testing names of tally types
/// \details This function is used by data output
TEST_F(test_TallyUtils, tallyTypeOutputNames) {
  ROOT_ONLY();

  const std::vector<TallyType> tests = {
    TallyType::Fission_RX,
    TallyType::NuFission_RX,
    TallyType::Absorption_RX,
    TallyType::Total_RX,
    TallyType::Scalar_Flux,

    TallyType::Fission_XS,
    TallyType::NuFission_XS,
    TallyType::Absorption_XS,
    TallyType::Total_XS,

    TallyType::Tracks_2D,
    TallyType::Tracks_3D
  };

  std::vector<std::string> oracles = {
    "Fission RX",
    "NuFission RX",
    "Absorption RX",
    "Total RX",
    "Scalar Flux",

    "Fission XS",
    "NuFission XS",
    "Absorption XS",
    "Total XS",

    "Tracks 2D",
    "Tracks 3D"
  };

  for (size_t i = 0; i < tests.size(); ++i) {
    auto &test = tests[i];
    auto &oracle = oracles[i];
    EXPECT_EQ(oracle, tallyutils::getTallyTypeName(test));
  }
}

} // namespace
