/// \file include/BoolParser.h
/// \brief A parser for boolean expressions
/// \date Aug 17, 2019
/// \author An Wang, USTB (wangan@xs.ustb.edu.cn)
///         Gen Wang, USTB (17812132070@163.com)

#ifndef BOOLPARSER_H_
#define BOOLPARSER_H_

#include <antmoc/PyVector.h>

#include <algorithm>
#include <limits>
#include <memory>
#include <stack>
#include <string>

namespace antmoc
{

// Forward decleration
struct BoolNode;
class ReverseNodeIter;

/// A pointer to BoolNode
using NodePtr = BoolNode*;

///---------------------------------------------------------------------
/// \class BoolNode
/// \brief The root node of an abstract syntax tree
///---------------------------------------------------------------------
struct BoolNode {
  friend class ReverseNodeIter;
  ReverseNodeIter rbegin();
  ReverseNodeIter rend();

  BoolNode(const int t, const NodePtr l = nullptr,
    const NodePtr r = nullptr)
    : token(t), left(l), right(r) { }

  BoolNode(const BoolNode&) = delete;
  BoolNode(BoolNode&&) = delete;
  ~BoolNode() {
    delete left;
    delete right;
  }

  NodePtr clone();
  void remove(int id);
  PyVector<int> toPostfixVector();
  std::string toPostfixString();

  int token;
  NodePtr left;
  NodePtr right;
};


///---------------------------------------------------------------------
/// \class BoolParser
/// \brief A boolean expression parser
/// \details Infix expressions could be converted into ASTs or RPNs.
///          The default behaviour of the public interface parse(...) is
///          to generate an RPN for the expression.
///---------------------------------------------------------------------
class BoolParser {

  // Operator:Insection、Union、Complement
  const std::string op_list = "&|~";

protected:
  int getPriority(const char token);
  bool hasOperator(const char token);

public:
  // Available operators for boolean expressions
  static const int OP_LEFT_PAREN;
  static const int OP_RIGHT_PAREN;
  static const int OP_COMPLEMENT;
  static const int OP_INTERSECTION;
  static const int OP_UNION;

  // Default constructor and destructor
  BoolParser() = default;
  ~BoolParser() = default;

  PyVector<int> tokenize(const std::string &expression);
  void addIntersection(PyVector<int> &tokens);

  /// \brief Parse an infix expression to an RPN
  /// \param expression an infix expression
  /// \return a vector storing a postfix expression
  PyVector<int> parse(const std::string &expression);

  /// \brief Parse an infix expression to an RPN
  /// \param expression an infix expression
  /// \return a vector storing a postfix expression
  PyVector<int> parse(const PyVector<int> &infix);

  // Parse expressions to ASTs
  void shuntAST(std::stack<NodePtr>&, std::stack<int>&);
  NodePtr parseToAST(const std::string &expression);
  NodePtr parseToAST(const PyVector<int> &infix);

  // Parse expressions to RPNs
  PyVector<int> parseToRPN(const std::string &expression);
  PyVector<int> parseToRPN(const PyVector<int> &infix);

  // Helpers
  bool isOperator(int token) { return token >= OP_UNION; }
  void removeOperators(PyVector<int> &expression);
  bool isSimple(const PyVector<int> &tokens);

  static std::string tokenToString(int token);
  static std::string toString(PyVector<int> tokens);

};


///---------------------------------------------------------------------
/// \class ReverseNodeIter
/// \brief An iterator for ASTs
/// \details The iterator iterates symbols of an AST by first
///          generating an RPN from it.
///---------------------------------------------------------------------
class ReverseNodeIter {
public:
  ///< Position of the past-the-last element
  static constexpr size_t end_pos = std::numeric_limits<size_t>::max();

  ReverseNodeIter(BoolNode &root, size_t pos)
    : _root(root), _pos(pos) {

    // Initialize data members
    initStack();
    if (_pos == end_pos)
      _pos = _rpn.size();
  }

  bool operator==(const ReverseNodeIter &rhs)
    { return (_pos == rhs._pos); }

  bool operator!=(const ReverseNodeIter &rhs)
    { return !(*this == rhs); }

  // Return a token
  int operator*() { return _rpn[_pos]->token; }

  /// \brief Increment operator
  ReverseNodeIter& operator++() {
    if (_pos < _rpn.size())
      ++_pos;
    else
      // The index past the last element
      _pos = _rpn.size();
    return *this;
  }

  void initStack();

  /// \brief Return the current position of the iterator
  size_t getPos() const { return _pos; }

protected:
  BoolNode &_root;  ///< Reference to the root node
  size_t _pos;      ///< Current position of the iterator
  std::vector<NodePtr> _rpn;  ///< Underlying RPN
};


/// \brief Converts an AST to an RPN stored in a stack
/// \details This implementation is a version of two-stack postorder
///          traversals.
inline
void ReverseNodeIter::initStack() {
  std::stack<NodePtr> aux_stack;
  aux_stack.push(&_root);

  NodePtr curr;
  while (!aux_stack.empty()) {
    // Push the current node
    curr = aux_stack.top();
    aux_stack.pop();
    _rpn.push_back(curr);

    // Push the childs of the current node
    if (curr->left)
      aux_stack.push(curr->left);
    if (curr->right)
      aux_stack.push(curr->right);
  }
  std::reverse(_rpn.begin(), _rpn.end());
}


} // namespace antmoc

#endif  // BOOLPARSER_H_
