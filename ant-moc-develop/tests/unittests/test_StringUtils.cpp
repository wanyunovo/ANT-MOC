/**
 * @date Aug 26, 2019
 * @author An Wang, USTB (wangan@xs.ustb.edu.cn)
 */

#include "testing/test_utils.h"
#include "antmoc/string_utils.h"

using namespace antmoc;

namespace {

/// Test fixture
class test_StringUtils : public testing::Test { };


TEST_F(test_StringUtils, ltrim) {
  ROOT_ONLY();

  StringVec tests = {
    "a",
    "1 2",
    " 1 2 ",
    "  1 2   ",
    "\t\v\f\n1 2\t\v\f\n",
  };

  StringVec oracles = {
    "a",
    "1 2",
    "1 2 ",
    "1 2   ",
    "1 2\t\v\f\n",
  };

  for (size_t i = 0; i < tests.size(); ++i) {
    stringutils::ltrim(tests[i]);
    EXPECT_EQ(oracles[i], tests[i]);
  }
}


TEST_F(test_StringUtils, rtrim) {
  ROOT_ONLY();

  StringVec tests = {
    "a",
    "1 2",
    " 1 2 ",
    "  1 2   ",
    "\t\v\f\n1 2\t\v\f\n",
  };

  StringVec oracles = {
    "a",
    "1 2",
    " 1 2",
    "  1 2",
    "\t\v\f\n1 2",
  };

  for (size_t i = 0; i < tests.size(); ++i) {
    stringutils::rtrim(tests[i]);
    EXPECT_EQ(oracles[i], tests[i]);
  }
}


TEST_F(test_StringUtils, trim) {
  ROOT_ONLY();

  StringVec tests = {
    "a",
    "1 2",
    " 1 2 ",
    "  1 2   ",
    "\t\v\f\n1 2\t\v\f\n",
  };

  std::vector<const char*> tests2 = {
    "a",
    "1 2",
    " 1 2 ",
    "  1 2   ",
    "\t\v\f\n1 2\t\v\f\n",
  };

  StringVec oracles = {
    "a",
    "1 2",
    "1 2",
    "1 2",
    "1 2",
  };

  for (size_t i = 0; i < tests.size(); ++i) {
    stringutils::trim(tests[i]);
    EXPECT_EQ(oracles[i], tests[i]);
    EXPECT_EQ(oracles[i], stringutils::trim(std::string(tests2[i])));
  }
}


/// \brief Testing antmoc::toUpper
TEST_F(test_StringUtils, toupper) {
  ROOT_ONLY();

  std::string str("aBcDe");
  stringutils::toUpper(str);
  EXPECT_STREQ("ABCDE", str.c_str());

  auto str2 = stringutils::toUpper(std::string("FgHiJ"));
  EXPECT_STREQ("FGHIJ", str2.c_str());
}


/// \brief Testing the word spliting function
TEST_F(test_StringUtils, splitSingleWordString) {
  ROOT_ONLY();

  constexpr size_t num_words = 1;
  std::string input = "Hello";
  std::string oracle = "Hello";

  auto parsed = stringutils::splitString(input);
  ASSERT_EQ(num_words, parsed.size());
  EXPECT_EQ(oracle, parsed.at(0));
}


/// \brief Testing the word spliting function
TEST_F(test_StringUtils, splitNormalString) {
  ROOT_ONLY();

  constexpr size_t num_words = 4;
  std::string input = "Hello USTB HPCer !";
  std::string oracle[num_words] = {"Hello", "USTB", "HPCer", "!"};

  auto parsed = stringutils::splitString(input);
  ASSERT_EQ(num_words, parsed.size());
  for (size_t i = 0; i < num_words; i++)
    EXPECT_EQ(oracle[i], parsed.at(i));
}


/// \brief Testing the word spliting function
TEST_F(test_StringUtils, splitNormalStringWithTabs) {
  ROOT_ONLY();

  constexpr size_t num_words = 4;
  std::string input = "Hello	USTB	HPCer !";
  std::string oracle[num_words] = {"Hello", "USTB", "HPCer", "!"};

  auto parsed = stringutils::splitString(input);
  ASSERT_EQ(num_words, parsed.size());
  for (size_t i = 0; i < num_words; i++)
    EXPECT_EQ(oracle[i], parsed.at(i));
}


/// \brief Testing the word spliting function
TEST_F(test_StringUtils, splitStringWithLeadingSpaces) {
  ROOT_ONLY();

  constexpr size_t num_words = 4;
  std::string input = " Hello USTB HPCer !";
  std::string oracle[num_words] = {"Hello", "USTB", "HPCer", "!"};

  auto parsed = stringutils::splitString(input);
  ASSERT_EQ(num_words, parsed.size());
  for (size_t i = 0; i < num_words; i++)
    EXPECT_EQ(oracle[i], parsed.at(i));

  // If there is a leading tab
  input = "	Hello USTB HPCer !";
  parsed = stringutils::splitString(input);
  ASSERT_EQ(num_words, parsed.size());
  for (size_t i = 0; i < num_words; i++)
    EXPECT_EQ(oracle[i], parsed.at(i));
}


/// \brief Testing the word spliting function
TEST_F(test_StringUtils, splitStringWithTrailingSpaces) {
  ROOT_ONLY();

  constexpr size_t num_words = 4;
  std::string input = "Hello USTB HPCer ! ";
  std::string oracle[num_words] = {"Hello", "USTB", "HPCer", "!"};

  auto parsed = stringutils::splitString(input);
  ASSERT_EQ(num_words, parsed.size());
  for (size_t i = 0; i < num_words; i++)
    EXPECT_EQ(oracle[i], parsed.at(i));

  // If there is a trailing tab
  input = "Hello USTB HPCer !	";
  parsed = stringutils::splitString(input);
  ASSERT_EQ(num_words, parsed.size());
  for (size_t i = 0; i < num_words; i++)
    EXPECT_EQ(oracle[i], parsed.at(i));
}


/// \brief Testing the word spliting function
TEST_F(test_StringUtils, splitStringWithAdjacentDelimeters) {
  ROOT_ONLY();

  constexpr size_t num_words = 4;
  std::string input = "Hello     USTB    HPCer ! ";
  std::string oracle[num_words] = {"Hello", "USTB", "HPCer", "!"};

  auto parsed = stringutils::splitString(input);
  ASSERT_EQ(num_words, parsed.size());
  for (size_t i = 0; i < num_words; i++)
    EXPECT_EQ(oracle[i], parsed.at(i));
}


/// \brief Testing the string concatenation function
TEST_F(test_StringUtils, joinStrings) {
  ROOT_ONLY();

  const std::vector<StringVec> tests = {
    {""},
    {"A", "", "B"}, // the mid empty string will also introduce a delimiter
    {"A", "B C", "D"},
  };

  const StringVec delimiter = {
    "", " ", ", ",
  };
  const std::vector<StringVec> oracles = {
    {"", "AB", "AB CD"},
    {"", "A  B", "A B C D"},
    {"", "A, , B", "A, B C, D"},
  };

  for (size_t i = 0; i < delimiter.size(); ++i) {
    for (size_t k = 0; k < tests.size(); ++k)
      EXPECT_EQ(oracles[i][k], stringutils::join(tests[k], delimiter[i]));
  }
}


/// \brief Testing the vector concatenation function
TEST_F(test_StringUtils, joinVectorElements) {
  ROOT_ONLY();

  const std::vector<std::vector<int>> tests = {
    {0},
    {-1, 1, 100},
  };

  const StringVec delimiter = {
    "", " ", ", ",
  };

  const std::vector<StringVec> oracles = {
    {"0", "0", "0"},
    {"-11100", "-1 1 100", "-1, 1, 100"},
  };

  for (size_t k = 0; k < tests.size(); ++k)
    for (size_t i = 0; i < delimiter.size(); ++i)
      EXPECT_EQ(oracles[k][i], stringutils::join(tests[k], delimiter[i]));
}


/// \brief Testing the string concatenation function
/// \details Taking an empty vector as input results in
///          an empty string
TEST_F(test_StringUtils, joinStringsRobust) {
  ROOT_ONLY();

  const StringVec test;
  const std::string delimiter;

  EXPECT_EQ(std::string(), stringutils::join(test, delimiter));
}


/// \brief Testing antmoc::stringutils::isSpaces
TEST_F(test_StringUtils, isSpaces) {
  ROOT_ONLY();

  const StringVec strvec = {
    "Fail", " \t\n\v\f\r  "
  };

  const bool oracle[2] = {
    false, true
  };

  for (size_t i = 0; i < strvec.size(); ++i)
    EXPECT_EQ(oracle[i], stringutils::isSpaces(strvec[i]));
}


/// \brief Testing string sustitution
TEST_F(test_StringUtils, spaceToUnderscore) {
  ROOT_ONLY();

  const StringVec tests = {
    "CPU LS Solver",
    " CPU Solver ",
    " Gauss Legendre",
    "Tabuchi Yamamoto ",
  };

  const StringVec oracles = {
    "CPU_LS_Solver",
    "CPU_Solver",
    "Gauss_Legendre",
    "Tabuchi_Yamamoto",
  };

  for (size_t i = 0; i < tests.size(); ++i) {
    EXPECT_EQ(oracles[i], stringutils::spaceToUnderscore(tests[i]));
  }
}


/// \brief Testing string sustitution
TEST_F(test_StringUtils, underscoreToSpace) {
  ROOT_ONLY();

  const StringVec tests = {
    "Fission_RX",
    " Fission_RX ",
    " CPU_LS_Solver",
    "Tabuchi Yamamoto ",
  };

  const StringVec oracles = {
    "Fission RX",
    "Fission RX",
    "CPU LS Solver",
    "Tabuchi Yamamoto",
  };

  for (size_t i = 0; i < tests.size(); ++i) {
    EXPECT_EQ(oracles[i], stringutils::underscoreToSpace(tests[i]));
  }
}


/// \brief Testing integer range expansion
TEST_F(test_StringUtils, toIntegerList) {
  ROOT_ONLY();

  const StringVec tests = {
    "1",

    "1:5",
    " 1 : 5",
    "1:5:2",
    "1:5:5",

    "5:1:-1",
    "5:1:-5",

    "5:5",
    "5:5:0",
    "5:5:1",
    "5:5:-1",
  };

  const std::vector<std::set<int>> oracles = {
    {1},

    {1,2,3,4,5},
    {1,2,3,4,5},
    {1,3,5},
    {1},

    {5,4,3,2,1},
    {5},

    {5},
    {5},
    {5},
    {5},
  };

  for (size_t i = 0; i < tests.size(); ++i)
    EXPECT_EQ(oracles[i], stringutils::toIntegerSet(tests[i], ":"));
}


#ifndef ENABLE_MPI_
TEST_F(test_StringUtils, toIntegerListRobust) {
  ROOT_ONLY();

  // Fail
  const StringVec tests_fail = {
    "1:5:-1",
    "5:1",
    "5:1:2",
    " ",
    ":",
  };

  for (auto test : tests_fail)
    EXPECT_ANY_THROW(stringutils::toIntegerSet(test))
      << "test: " << test << std::endl;

  // Pass
  const StringVec tests_pass = {
    "",
  };

  for (auto test : tests_pass)
    EXPECT_NO_THROW(stringutils::toIntegerSet(test))
      << "test: " << test << std::endl;
}
#endif


} // namespace
