/// \file test_XMLUtils.cpp
/// \brief Test xml utility
/// \date Dec 22, 2019
/// \author An Wang, USTB (wangan@xs.ustb.edu.cn)

#include "testing/test_utils.h"
#include "antmoc/string_utils.h"
#include "antmoc/xml_utils.h"

using namespace antmoc;

namespace {

/// Testing fixture
class test_XMLUtils : public testing::Test {
 protected:
  tinyxml2::XMLDocument doc;  ///< XML document object
};

#ifndef ENABLE_MPI_
/// Query for an nonexisting element
TEST_F(test_XMLUtils, queryFirstChildRobust) {
  std::string test = "<head/>";

  doc.Parse(test.c_str());
  auto element = doc.RootElement();
  EXPECT_ANY_THROW(antmoc::queryFirstChild(element, "tail", true));
}


/// Function queryStringAttr() is used to extract a string from
/// attributes.
TEST_F(test_XMLUtils, queryStringAttr) {
  std::string test = "<head name='Alice'/>";

  std::vector<bool> enforces = {
    false,
    true,
  };

  doc.Parse(test.c_str());
  auto element = doc.RootElement();
  for (auto enforce : enforces) {
    EXPECT_STREQ("Alice", antmoc::queryStringAttr(element, "name", "head", enforce));
  }
}


/// Function queryStringAttr() is used to extract a string from
/// attributes.
/// This will test its robustness.
TEST_F(test_XMLUtils, queryStringAttrRobust) {
  std::string test = "<head/>";

  doc.Parse(test.c_str());
  auto element = doc.RootElement();
  // No exceptions because enforced == false
  EXPECT_NO_THROW(antmoc::queryStringAttr(element, "name", "head", false));
  // Exception throwed because enforced == true
  EXPECT_ANY_THROW(antmoc::queryStringAttr(element, "name", "head", true));
}


/// Function queryDoubleAttr() is used to extract a floating-point
/// number from attributes.
TEST_F(test_XMLUtils, queryDoubleAttr) {
  std::string test = "<head value='1.0'/>";

  std::vector<bool> enforces = {
    false,
    true,
  };

  doc.Parse(test.c_str());
  auto element = doc.RootElement();
  for (auto enforce : enforces) {
    EXPECT_EQ(1.0, antmoc::queryDoubleAttr(element, "value", "head", enforce));
  }
}


/// Function queryDoubleAttr() is used to extract a floating-point
/// number from attributes.
/// This will test its robustness.
TEST_F(test_XMLUtils, queryDoubleAttrRobust) {
  std::string test = "<head/>";

  doc.Parse(test.c_str());
  auto element = doc.RootElement();
  // No exceptions because enforced == false
  EXPECT_NO_THROW(antmoc::queryDoubleAttr(element, "value", "head", false));
  // Exception throwed because enforced == true
  EXPECT_ANY_THROW(antmoc::queryDoubleAttr(element, "value", "head", true));
}


/// Function queryIntAttr() is used to extract an integer number
/// from attributes.
TEST_F(test_XMLUtils, queryIntAttr) {
  std::string test = "<head id='1'/>";

  std::vector<bool> enforces = {
    false,
    true,
  };

  doc.Parse(test.c_str());
  auto element = doc.RootElement();
  for (auto enforce : enforces) {
    EXPECT_EQ(1, antmoc::queryIntAttr(element, "id", "head", enforce));
  }
}


/// Function queryIntAttr() is used to extract an integer number
/// from attributes.
/// This will test its robustness.
TEST_F(test_XMLUtils, queryIntAttrRobust) {
  std::string test = "<head/>";

  doc.Parse(test.c_str());
  auto element = doc.RootElement();
  // No exceptions because enforced == false
  EXPECT_NO_THROW(antmoc::queryIntAttr(element, "id", "head", false));
  // Exception throwed because enforced == true
  EXPECT_ANY_THROW(antmoc::queryIntAttr(element, "id", "head", true));
}


/// Function queryNodeDouble() is used to extract a floating-point
/// number either from attributes or elements but not both.
/// This will test if the query succeeds or if there are any duplicated
/// elements or attributes.
TEST_F(test_XMLUtils, queryNodeDoubleRobust) {
  StringVec tests = {
    // Nonexisting attribute/sub-element
    "<head/>",
    // Conflicts between attributes and elements
    "<head value='1.0'>"
    "  <value> 2.0 </value>"
    "</head>",
    // Conflicts between elements
    "<head>"
    "  <value> 1.0 </value>"
    "  <value> 2.0 </value>"
    "</head>"
  };

  for (auto &s : tests) {
    doc.Parse(s.c_str());
    auto element = doc.RootElement();
    EXPECT_ANY_THROW(antmoc::queryNodeDouble(element, "value", true));
  }
}


/// Function queryNodeInt() is used to extract an integer
/// number either from attributes or elements but not both.
/// This will test if the query succeeds or if there are any duplicated
/// elements or attributes.
TEST_F(test_XMLUtils, queryNodeIntRobust) {
  StringVec tests = {
    // Nonexisting attribute/sub-element
    "<head/>",
    // Conflicts between attributes and elements
    "<head id='1'>"
    "  <id> 2 </id>"
    "</head>",
    // Conflicts between elements
    "<head>"
    "  <id> 1 </id>"
    "  <id> 2 </id>"
    "</head>"
  };

  for (auto &s : tests) {
    doc.Parse(s.c_str());
    auto element = doc.RootElement();
    EXPECT_ANY_THROW(antmoc::queryNodeInt(element, "id", true));
  }
}
#endif  // ENABLE_MPI_


} // namespace antmoc
