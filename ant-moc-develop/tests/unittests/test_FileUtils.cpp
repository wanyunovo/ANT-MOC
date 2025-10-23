/// \file test_FileUtils.cpp
/// \brief Test file utility
/// \date Oct 18, 2020
/// \author An Wang, USTB (wangan@xs.ustb.edu.cn)

#include "testing/test_utils.h"
#include "antmoc/file_utils.h"

using namespace antmoc;

namespace {

/// Testing fixture
class test_FileUtils : public testing::Test {
};

TEST_F(test_FileUtils, getExtension) {
  ROOT_ONLY();

  EXPECT_EQ(fileutils::getExtension("settings.xml"), ".xml");
  EXPECT_EQ(fileutils::getExtension("settings.toml"), ".toml");
  EXPECT_EQ(fileutils::getExtension("settings."), ".");
  EXPECT_EQ(fileutils::getExtension("settings"), "");
}

} // namespace
