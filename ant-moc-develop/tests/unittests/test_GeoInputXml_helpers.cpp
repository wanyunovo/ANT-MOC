/**
 * @file test_GeoInputXml_Helpers.cpp
 * @brief Test helper functions
 * @date April 19, 2019
 * @author An Wang, USTB (wangan@xs.ustb.edu.cn)
 */

#include "testing/geometry_testharness.h"

#include <vector>

using namespace antmoc;

namespace {

// Test fixture
class test_GeoInputXml_Helpers : public GeometryTestHarness {
  // do nothing
};


TEST_F(test_GeoInputXml_Helpers, replicateANumberOnce) {
  ROOT_ONLY();

  int times = 1;
  int number = 1;
  std::string oracle = "1";
  EXPECT_EQ(oracle, _geo_input->replicateANumber(number, times));
  // Two digits
  number = 11;
  oracle = "11";
  EXPECT_EQ(oracle, _geo_input->replicateANumber(number, times));
}

TEST_F(test_GeoInputXml_Helpers, replicateANumberManyTimes) {
  ROOT_ONLY();

  int times = 2;
  int number = 1;
  std::string oracle = "1 1";
  EXPECT_EQ(oracle, _geo_input->replicateANumber(number, times));
  number = 11;
  oracle = "11 11";
  EXPECT_EQ(oracle, _geo_input->replicateANumber(number, times));
}

TEST_F(test_GeoInputXml_Helpers, extendOnlyDigits) {
  ROOT_ONLY();

  EXPECT_EQ("1", _geo_input->extendR("1"));
  EXPECT_EQ("1 2", _geo_input->extendR("1 2"));
  EXPECT_EQ("1 2", _geo_input->extendR(" 1 2"));
  EXPECT_EQ("1 2", _geo_input->extendR("1 2 "));
  EXPECT_EQ("10 2 34", _geo_input->extendR("10 2 34"));
  // Handle tabs
  EXPECT_EQ("10 2 34", _geo_input->extendR("10	2 34"));
}

TEST_F(test_GeoInputXml_Helpers, extendOneSymbolWithOneDigits) {
  ROOT_ONLY();

  std::string line = "1 1R";
  std::string oracle = "1 1";
  EXPECT_EQ(oracle, _geo_input->extendR(line));

  line = "2 2More";
  oracle = "2 2 2";
  EXPECT_EQ(oracle, _geo_input->extendR(line));
}

TEST_F(test_GeoInputXml_Helpers, extendOneSymbolWithTwoDigits) {
  ROOT_ONLY();

  std::string line = "1 11R";
  std::string oracle = "1 1 1 1 1 1 1 1 1 1 1 1";
  EXPECT_EQ(oracle, _geo_input->extendR(line));

  line = "2 12*";
  oracle = "2 2 2 2 2 2 2 2 2 2 2 2 2";
  EXPECT_EQ(oracle, _geo_input->extendR(line));
}

TEST_F(test_GeoInputXml_Helpers, extendTwoSymbolWithOneDigits) {
  ROOT_ONLY();

  std::string line = "1 1R 2 1R";
  std::string oracle = "1 1 2 2";
  EXPECT_EQ(oracle, _geo_input->extendR(line));

  line = "10 1R 20 2R";
  oracle = "10 10 20 20 20";
  EXPECT_EQ(oracle, _geo_input->extendR(line));

  line = "300 2R 400 1R";
  oracle = "300 300 300 400 400";
  EXPECT_EQ(oracle, _geo_input->extendR(line));
}

TEST_F(test_GeoInputXml_Helpers, extendTwoSymbolWithTwoDigits) {
  ROOT_ONLY();

  std::string line = "1 11R 2 11R";
  std::string oracle = "1 1 1 1 1 1 1 1 1 1 1 1"
                       " 2 2 2 2 2 2 2 2 2 2 2 2";
  EXPECT_EQ(oracle, _geo_input->extendR(line));

  line = "10 11R 20 12R";
  oracle = "10 10 10 10 10 10 10 10 10 10 10 10"
           " 20 20 20 20 20 20 20 20 20 20 20 20 20";
  EXPECT_EQ(oracle, _geo_input->extendR(line));

  line = "30 12R 40 1R";
  oracle = "30 30 30 30 30 30 30 30 30 30 30 30 30"
           " 40 40";
  EXPECT_EQ(oracle, _geo_input->extendR(line));
}

TEST_F(test_GeoInputXml_Helpers, extendAdjacentNumbersFollowedSymbols) {
  ROOT_ONLY();

  std::string line = "1 2 1A 3 4 2B 5";
  std::string oracle = "1 2 2 3 4 4 4 5";
  EXPECT_EQ(oracle, _geo_input->extendR(line));
}

TEST_F(test_GeoInputXml_Helpers, extendAdjacentSymbols) {
  ROOT_ONLY();

  std::string line = "1 2A 2B 2";
  std::string oracle = "1 1 1 1 1 2";
  EXPECT_EQ(oracle, _geo_input->extendR(line));
}

TEST_F(test_GeoInputXml_Helpers, extendLeadingSymbols) {
  ROOT_ONLY();

  std::string line = "2R 2";
  std::string oracle = "";
  EXPECT_EQ(oracle, _geo_input->extendR(line));
}

} /* namespace */
