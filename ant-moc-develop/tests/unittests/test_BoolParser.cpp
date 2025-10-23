/// \file test_BoolParser.cpp
/// \brief Test class BoolParser

#include "testing/test_utils.h"
#include "antmoc/BoolParser.h"
#include "antmoc/string_utils.h"

using namespace antmoc;
namespace {

/// Derived class
class BoolParserTest : public BoolParser {
  BoolParserTest() = default;

  // Friend classed for testing private members
  friend class test_BoolParser;
  FRIEND_TEST(test_BoolParser, getPriority);
  FRIEND_TEST(test_BoolParser, getPriorityRobust);
  FRIEND_TEST(test_BoolParser, hasOperator);
};


/// Test fixture
class test_BoolParser : public testing::Test {
protected:

  BoolParserTest parser;

  // Shorthands
  int OP_U = BoolParser::OP_UNION;
  int OP_I = BoolParser::OP_INTERSECTION;
  int OP_C = BoolParser::OP_COMPLEMENT;
  int OP_LP = BoolParser::OP_LEFT_PAREN;
  int OP_RP = BoolParser::OP_RIGHT_PAREN;
};


#ifndef ENABLE_MPI_
/// Test the robustness of BoolParser
TEST_F(test_BoolParser, boolParserRobust){
  StringVec tests = {
    // mismatched parenthesis
    "(", "(1|2", "1|2)",
    "~(1|2))",
    // duplicated signs
    "--1", "-+1", "+-1",
    "---1", "+-+1", "-+-1",
    "-", "+",
    // missing operands
    "1|", "1|2~"
  };

  for(auto &test : tests)
    EXPECT_ANY_THROW(parser.parse(test))
      << " Test: " << test << '\n';
}
#endif  // ENABLE_MPI_


/// Test function getPriority()
TEST_F(test_BoolParser, getPriority) {
  ROOT_ONLY();

  const std::vector<char> operators = {
    '(', ')', '~', '&', '|',
  };

  const std::vector<int> oracles = {
    OP_LP,
    OP_RP,
    OP_C,
    OP_I,
    OP_U,
  };

  for (size_t i = 0; i < operators.size(); ++i) {
    EXPECT_EQ(oracles[i], parser.getPriority(operators[i]));
  }
}


#ifndef ENABLE_MPI_
/// Test the robustness of function getPriority()
TEST_F(test_BoolParser, getPriorityRobust) {
  const std::vector<char> operators = {
    '+', '-', '*', '/'
  };

  for (auto op : operators) {
    EXPECT_ANY_THROW(parser.getPriority(op));
  }
}
#endif  // ENABLE_MPI_


/// Test function hasOperator()
TEST_F(test_BoolParser, hasOperator) {
  ROOT_ONLY();

  EXPECT_TRUE(parser.hasOperator('&'));
  EXPECT_FALSE(parser.hasOperator('+'));
}


#ifndef ENABLE_MPI_
/// Test the robustness of function tokenize()
/// It returns an empty vector if there is nothing to be tokenized.
/// A exception will be throwed if there is any undefined characters.
TEST_F(test_BoolParser, tokenizeRobust) {
  // Tokenize an empty string
  std::string test = "";
  PyVector<int> result;
  EXPECT_NO_THROW(result = parser.tokenize(test));

  // Tokenize a string with undefined characters
  test = "1*2";
  EXPECT_ANY_THROW(parser.tokenize(test));
}
#endif  // ENABLE_MPI_


/**
 * test the tokenize function:
 *  convert the bool_expression to the vector element  whithout space.
 */
TEST_F(test_BoolParser, tokenize){
  ROOT_ONLY();

  StringVec bool_exprs = {
    "-1",
    "2 3",
    "2  3    -4",

    "1|2",
    "1|2 3",
    "~1 -2",
    "~(-2)",

    "(1|2|-3)",
    "(1 2(~-3))|(2 3)"
  };

  std::vector<StrPyVec> oracles = {
   {"-1"},
   {"+2", "+3"},
   {"+2", "+3", "-4"},

   {"+1", "|", "+2"},
   {"+1", "|", "+2", "+3"},
   {"~", "+1", "-2"},
   {"~", "(", "-2", ")"},

   {"(", "+1", "|", "+2", "|", "-3", ")"},
   {"(", "+1", "+2", "(", "~", "-3", ")", ")", "|", "(", "+2", "+3", ")"}

  };

  for(size_t i = 0; i < bool_exprs.size(); i++) {

    PyVector<int> tokens = parser.tokenize(bool_exprs[i]);

    EXPECT_EQ(oracles[i].toString(), parser.toString(tokens));
  }
}


/*
 * test the add_insetction function:
 *  add the operator intersection where a missing operator is needed
 */
TEST_F(test_BoolParser, addIntersection) {
  ROOT_ONLY();

  StringVec bool_exprs = {
    "-1",
    "2 3",
    "2  3    -4",

    "1|2",
    "1|2 3",
    "~1 -2",
    "~(-2)",

    "(1|2|-3)",
    "(1 2(~-3))|(2 3)"
  };

  std::vector<StrPyVec> oracles = {
    {"-1"},
    {"+2", "&","+3"},
    {"+2", "&", "+3", "&", "-4"},

    {"+1", "|", "+2"},
    {"+1", "|", "+2", "&", "+3"},
    {"~", "+1", "&", "-2"},
    {"~", "(", "-2", ")"},

    {"(", "+1", "|", "+2", "|", "-3", ")"},
    {"(", "+1","&", "+2", "&", "(", "~", "-3", ")", ")", "|", "(", "+2", "&", "+3", ")"}

  };


  for(size_t i = 0; i<bool_exprs.size();i++){

    auto infix = parser.tokenize(bool_exprs[i]);
    parser.addIntersection(infix);

    EXPECT_EQ(oracles[i].toString(), parser.toString(infix));
  }
}

/// Test function BoolParser::parse()
/// This is the main interface for converting expressions to RPNs
TEST_F(test_BoolParser, complicatedExpressions){
  ROOT_ONLY();

  StringVec bool_exprs = {
   // one operator
   "1 -211",
   "1|-3",
   "~3",
   "~(-3)",
   //complicated operator
   "1 2|3",
   "1  2 3",
   "1 2 ~3",
   "~1|2  &3",
   // Associativity
   "~(1 2)",
   "~(1|2)",
   "~(1 2|-3|2) -4",
   "~(1 2|3|2) &4(~10 -122)|3"
  };

  std::vector<StrPyVec> oracles = {

   {"+1", "-211", "&"},
   {"+1", "-3", "|"},
   {"+3", "~"},
   {"-3", "~"},

   {"+1", "+2", "&", "+3", "|"},
   {"+1", "+2", "&", "+3", "&"},
   {"+1", "+2", "&", "+3", "~", "&"},
   {"+1", "~", "+2", "+3", "&", "|"},

   {"+1", "+2", "&", "~"},
   {"+1", "+2", "|", "~"},
   {"+1", "+2", "&", "-3", "|", "+2", "|", "~", "-4", "&"},
   {"+1", "+2", "&", "+3", "|", "+2", "|", "~", "+4", "&", "+10", "~", "-122", "&", "&", "+3", "|"}
  };

  for(size_t i = 0; i<bool_exprs.size(); i++){
    auto rpn = parser.parse(bool_exprs[i]);

    EXPECT_EQ(oracles[i].toString(), parser.toString(rpn));
  }

}


/// Test function BoolParser::parseToRPN()
/// This function first tokenizes the string and then parse it.
/// Results will be saved into a vector of integers.
TEST_F(test_BoolParser, parseToRPN) {
  ROOT_ONLY();

  const StringVec bool_exprs = {
   "1|-3",
   "1 2 ~3",
   "1 2 ~(3 4)",
  };

  const std::vector<IntPyVec> oracles = {
   {1, -3, OP_U},
   {1, 2, OP_I, 3, OP_C, OP_I},
   {1, 2, OP_I, 3, 4, OP_I, OP_C, OP_I},
  };

  for(size_t i = 0; i<bool_exprs.size(); ++i) {
    auto rpn = parser.parseToRPN(bool_exprs[i]);
    EXPECT_EQ(oracles[i].toString(), rpn.toString());
  }
}


#ifndef ENABLE_MPI_
/// Test the robustness of BoolParser::parseToRPN
TEST_F(test_BoolParser, parseToRPNRobust){
  StringVec tests = {
    // mismatched parenthesis
    "(",
    "(1|2",
    "1|2)",
    "~(1|2))",
    // duplicated signs
    "--1", "-+1", "+-1",
    "---1", "+-+1", "-+-1",
    "-", "+",
  };

  for(auto &test : tests)
    EXPECT_ANY_THROW(parser.parseToRPN(test))
      << " Test: " << test << '\n';
}
#endif  // ENABLE_MPI_


/// Test function BoolParser::parseToAST()
/// This function first tokenizes the string and then parse it.
/// In fact, this test makes use of toPostfixVector() to compare the results,
/// and the conversion function relies on a custom iterator.
TEST_F(test_BoolParser, parseToAST) {
  ROOT_ONLY();

  const StringVec bool_exprs = {
   "1|-3",
   "1 2 ~3",
  };

  const std::vector<IntPyVec> oracles = {
   {1, -3, OP_U},
   {1, 2, OP_I, 3, OP_C, OP_I},
  };

  for(size_t i = 0; i<bool_exprs.size(); ++i) {
    auto ast = parser.parseToAST(bool_exprs[i]);
    EXPECT_EQ(oracles[i].toString(), ast->toPostfixVector().toString());
  }
}


/// Test BoolParser::removeOperators(PyVector<int>&)
/// This function takes a boolean expression and remove all of the valid operators
/// from it. After that, the vector will be resized and shrinked.
/// However, shrinking is a non-binding request so that we do not check the capacity
/// of the vector.
TEST_F(test_BoolParser, removeOperators) {
  ROOT_ONLY();

  std::vector<IntPyVec> bool_exprs = {
   {1, -3, OP_U},
   {1, 2, OP_I, 3, OP_C, OP_I},
  };

  const std::vector<IntPyVec> oracles = {
   {1, -3},
   {1, 2, 3},
  };

  for(size_t i = 0; i < bool_exprs.size(); ++i) {
    // Pass the reference to the expression
    parser.removeOperators(bool_exprs[i]);
    EXPECT_EQ(oracles[i].size(), bool_exprs[i].size());
    EXPECT_EQ(oracles[i].toString(), bool_exprs[i].toString());
  }
}

} /* namespace*/
