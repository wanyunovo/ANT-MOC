/**
 * @date Aug 17, 2019
 * @author An Wang, USTB (wangan@xs.ustb.edu.cn)
 */

#include "testing/test_utils.h"
#include "antmoc/VecParser.h"
#include "antmoc/string_utils.h"

using namespace antmoc;

namespace {

// Test fixture
class test_VecParser : public testing::Test {
 protected:
  VecParser parser;
};


TEST_F(test_VecParser, singleOperand) {
  ROOT_ONLY();

  std::vector<std::string> exprs = {
    "0", "1.", "21.42", "5", "[0]", "[1]"
  };

  std::vector<double> oracles = {
    0, 1., 21.42, 5., 0., 1.
  };

 for (size_t i = 0; i < exprs.size(); ++i) {
    auto result = parser.parse(exprs[i]);
    ASSERT_EQ(1, result.size())
      << "  Case: " << i << '\n';
    EXPECT_EQ(oracles.at(i), result.at(0))
      << "  Case: " << i << '\n';
  }
}


TEST_F(test_VecParser, vectorExpressions) {
  ROOT_ONLY();

  std::vector<std::string> exprs = {
    // Simple expressions
    "[1.0] + 2.0",
    " 2.0 + [1.0]",
    "[1.0] * 1",
    "[1.0] * 2 ",
    " 3 * [1.0] ",
    "3*[2]",
    "[3]*2",

    "[2.0] + (1 + 2)",
    "[2.0] + 1 + 2",  // equals ([2.0]+1)+2
    "[2.0] + [1] + 2",
    "1 + [2.0] + 2",
    // Commutativity
    "0.2 + [0.3] * 2",
    "2 * [0.3] + 0.2",
    // Associativity
    "(0.1 + 2 * ([0.2] + [0.3])) * 2",
    // Nested square brackets
    "2.0 + [[[1.2]]]",
    "2.0 + [[1.2] * 2 + 0.3] * 2",
  };

  std::vector<DblPyVec> oracles = {
    {1., 2.},
    {2., 1.},
    {1.},
    {1., 1.},
    {1., 1., 1.},
    {2., 2., 2.},
    {3., 3.},

    {2.0, 3.0},
    {2.0, 1.0, 2.0},
    {2.0, 1.0, 2.0},
    {1.0, 2.0, 2.0},

    {0.2, 0.3, 0.3},
    {0.3, 0.3, 0.2},

    {0.1, 0.2, 0.3, 0.2, 0.3, 0.1, 0.2, 0.3, 0.2, 0.3},

    {2.0, 1.2},
    {2.0, 1.2, 1.2, 0.3, 1.2, 1.2, 0.3},
  };
  ASSERT_EQ(exprs.size(), oracles.size());

  for (size_t i = 0; i < exprs.size(); ++i) {
    auto result = parser.parse(exprs[i]);
    EXPECT_EQ(oracles[i].toString(), result.toString())
      << "  Case: " << i << '\n';
  }
}


TEST_F(test_VecParser, normalExpressions) {
  ROOT_ONLY();

  std::vector<std::string> exprs = {
    "1 + 2",
    "2 * 3",
    // Commutativity
    "2 + 3 * 2",
    "2 * 3.1 + 2",
    // Associativity
    "(1 + 2 * (2 + 3)) * 2",
  };

  std::vector<DblPyVec> oracles = {
    {3.},
    {6.},
    {8.},
    {8.2},
    {22},
  };

  for (size_t i = 0; i < exprs.size(); ++i) {
    auto result = parser.parse(exprs[i]);
    EXPECT_EQ(oracles[i].toString(), result.toString())
      << "  Case: " << i << '\n';
  }
}


#ifndef ENABLE_MPI_
TEST_F(test_VecParser, vecParserRobust) {
  ROOT_ONLY();

  std::vector<std::string> exprs = {
    // Non-exist operators
    "1. - 1",
    "1. % 1",
    // Missing operands
    "1. + 1. **",
    "1. ++ 1.",
    // Ill-formed
    "(2.0 * 2) + 1.0)",
    "(2.0 * 2 + (1.0)",
    "2.0 * [2] + [1.0",
    "2.0 * [2] + 1.0]",
    "[a]",
    "[]",
    // Invalid mulplication
    "[2.0] * [3.0]",
    "[2.0] * 3.0",
    "2.0 * [3.0]",
  };

  for (size_t i = 0; i < exprs.size(); ++i) {
    EXPECT_ANY_THROW(parser.parse(exprs[i]))
      << "  Case " << i << ' ' << exprs[i] << '\n';
  }
}
#endif  // ENABLE_MPI_


/// Test function toString()
TEST_F(test_VecParser, toString) {
  ROOT_ONLY();

  const std::vector<Operand> operands = {
    Operand(size_t(1)),
    Operand(double(2.2)),
    Operand(DblPyVec({1.1, 2.2})),
  };

  const StringVec oracles = {
    "1",
    "2.2",
    "1.1 2.2",
  };

  for (size_t i = 0; i < operands.size(); ++i) {
    EXPECT_EQ(oracles[i], operands[i].toString());
  }
}


/// Test function toIntVector()
/// This is a helpler function. Any floating-point numbers
/// will be truncated by the function.
TEST_F(test_VecParser, toIntVector) {
  ROOT_ONLY();

  const std::vector<Operand> operands = {
    Operand(size_t(1)),
    Operand(double(2.2)),
    Operand(DblPyVec({1.1, 2.2})),
  };

  const std::vector<IntPyVec> oracles = {
    {1},
    {2},
    {1, 2}
  };

  for (size_t i = 0; i < operands.size(); ++i) {
    EXPECT_EQ(oracles[i].toString(), operands[i].toIntVector().toString());
  }
}

} // namespace
