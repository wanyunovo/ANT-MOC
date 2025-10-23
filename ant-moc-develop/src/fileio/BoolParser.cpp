/// \file src/BoolParser.cpp

#include "antmoc/BoolParser.h"
#include "antmoc/constants.h"
#include "antmoc/log.h"
#include "antmoc/string_utils.h"

#include <cmath>
#include <sstream>
#include <stdexcept>

namespace antmoc
{

//----------------------------------------------------------------------
//  BoolNode
//----------------------------------------------------------------------

/// \brief Returns an iterator pointing to the first node
ReverseNodeIter BoolNode::rbegin() {
  return ReverseNodeIter(*this, 0);
}


/// \brief Returns an iterator pointing to the past-the-end node
ReverseNodeIter BoolNode::rend() {
  return ReverseNodeIter(*this, ReverseNodeIter::end_pos);
}


/// \brief Removes all the nodes containing a specified id.
/// \details Signs are ignored. This is an simple implementation
///          of the Preoder Traversal algorithm.
///          FIXME: The current implementation of remove() doesn't clean
///          superfluous operators.
void BoolNode::remove(int id) {
  if (id >= BoolParser::OP_UNION)
    log::error("The value of id is too large: {}", id);

  std::stack<NodePtr> stack;
  stack.push(this);

  while (!stack.empty()) {
    auto curr = stack.top();
    stack.pop();
    // Determine if it has a child to be removed.
    // If there is nothing to be removed, childs will be pushed
    // back into the stack.
    if (curr->right) {
      if (id == std::abs(curr->right->token))
        curr->right = nullptr;  // memory is not freed
      else
        stack.push(curr->right);
    }
    if (curr->left) {
      if (id == std::abs(curr->left->token))
        curr->left = nullptr;
      else
        stack.push(curr->left);
    }
  }
}


/// \brief Deep copy
/// \details This is a recursive preorder traversal.
NodePtr BoolNode::clone() {
  auto parent = new BoolNode(token);
  if (left)
    parent->left = left->clone();
  if (right)
    parent->right = right->clone();

  return parent;
}


/// \brief Converts an AST to an RPN vector
PyVector<int> BoolNode::toPostfixVector() {
  PyVector<int> rpn;

  for (auto it = this->rbegin(); it != this->rend(); ++it)
    rpn.push_back(*it);

  return rpn;
}


/// \brief Converts an AST to an RPN string
std::string BoolNode::toPostfixString() {
  std::ostringstream os;

  for (auto it = this->rbegin(); it != this->rend(); ++it) {
    if (*it < BoolParser::OP_UNION && *it > 0)  // positive operands
      os << '+';
    os << BoolParser::tokenToString(*it) << ' ';
  }

  return stringutils::rtrim(os.str());
}


//----------------------------------------------------------------------
//  BoolParser
//----------------------------------------------------------------------

/// Operator left parenthesis '('
const int BoolParser::OP_LEFT_PAREN   = std::numeric_limits<int>::max();
/// Operator right parenthesis ')'
const int BoolParser::OP_RIGHT_PAREN  = BoolParser::OP_LEFT_PAREN-1;
/// Operator complement '~'
const int BoolParser::OP_COMPLEMENT   = BoolParser::OP_LEFT_PAREN-2;
/// Operator intersection '&'
const int BoolParser::OP_INTERSECTION = BoolParser::OP_LEFT_PAREN-3;
/// Operator union '|'
const int BoolParser::OP_UNION        = BoolParser::OP_LEFT_PAREN-4;


/// \details First parse the expression to an AST to validate the
///          expression, then convert the AST to an RPN.
PyVector<int> BoolParser::parse(const std::string &expression) {
  auto ast = parseToAST(expression);
  auto vec = ast->toPostfixVector();
  delete ast;
  return vec;
}


/// \details First parse the expression to an AST to validate the
///          expression, then convert the AST to an RPN.
PyVector<int> BoolParser::parse(const PyVector<int> &infix) {
  auto ast = parseToAST(infix);
  auto vec = ast->toPostfixVector();
  delete ast;
  return vec;
}


/// \brief Get the priority of operators
/// \param token the operator
/// \return priority of the operator
int BoolParser::getPriority(const char token){
  int prio = 0;
  switch(token)
  {
    case '(': prio = OP_LEFT_PAREN; break;
    case ')': prio = OP_RIGHT_PAREN; break;
    case '~': prio = OP_COMPLEMENT; break;
    case '&': prio = OP_INTERSECTION; break;
    case '|': prio = OP_UNION; break;
    default:
      throw std::logic_error("Nonexist operator");
  }
  return prio;
}


/// \brief Determin whether the operator is defined
bool BoolParser::hasOperator(const char token){
  return (op_list.find(token) != std::string::npos);
}


/// \brief Tokenize a boolean expression
/// \details The original expression will be tokenized and stored in
///          a PyVector.
/// \return a vector of tokens
PyVector<int> BoolParser::tokenize(const std::string &expression) {

  PyVector<int> tokens;
  if(expression.empty())
    return tokens;

  // Loop over characters of the expression
  for (size_t i = 0; i < expression.size(); ) {
    // Spaces will be skipped.
    // Any necessary intersection operator will be added later
    // in method addIntersection().
    if (isspace(expression[i])) {
      i++;
    }
    else if (expression[i] == '(') {
      tokens.push_back(OP_LEFT_PAREN);
      i++;
    }
    else if (expression[i] == ')') {
      tokens.push_back(OP_RIGHT_PAREN);
      i++;
    }
    else if (expression[i] == '~') {
      tokens.push_back(OP_COMPLEMENT);
      i++;
    }
    else if (expression[i] == '&') {
      tokens.push_back(OP_INTERSECTION);
      i++;
    }
    else if (expression[i] == '|') {
      tokens.push_back(OP_UNION);
      i++;
    }
    else if (expression[i] == '-' ||
             expression[i] == '+' ||
             isdigit(expression[i])) {
      auto j = i + 1;
      // Find an id in the expression
      while ( j < expression.size() &&
              isdigit(expression[j]) ) {
        j++;
      }
      tokens.push_back(std::stoi(expression.substr(i, j - i)));
      i = j;
    }
    else{
       log::error("Invalid character found in '{}'", expression);
    }

  }

  return tokens;
}


/// \brief Add intersection operators on demand
/// \details Intersection operators may be omitted in the original
///          boolean expression. This method detects all of the
///          possible positions to insert intersections.
void BoolParser::addIntersection(PyVector<int> &tokens) {

  size_t i = 0;
  while (i < tokens.size() - 1) {

    // Indicate there is an id or ')'
    bool curr_compat = (tokens[i] < OP_UNION) ||
                       (tokens[i] == OP_RIGHT_PAREN);
    // Indicate there is an id or '(' or '~'
    bool next_compat = (tokens[i + 1] < OP_UNION) ||
                       (tokens[i + 1] == OP_LEFT_PAREN) ||
                       (tokens[i + 1] == OP_COMPLEMENT);

    if (curr_compat && next_compat) {
      tokens.insert(tokens.begin() + i + 1, OP_INTERSECTION);
    }
    i++;
  }
}


/// \brief One-step shunting
/// \param operands a stack of AST nodes
/// \param operators a stack of operators
void BoolParser::shuntAST(std::stack<NodePtr> &operands,
                          std::stack<int> &operators) {

  // Make a new node to be the parent
  int op = operators.top();
  auto parent = new BoolNode(op);

  if (operands.empty()) {
    log::error("Missing operands in the boolean expression");
  }
  parent->right = operands.top();
  operands.pop();

   // binary operators
  if (op != OP_COMPLEMENT) {
    if (operands.empty()) {
      log::error("Missing operands in the boolean expression");
    }
    parent->left = operands.top();
    operands.pop();
  }

  // Push the result back to the stack
  operands.push(parent);
  operators.pop();
}


/// \brief Parse an infix expression to an AST
/// \param expression an infix expression
/// \return the root of an AST
NodePtr BoolParser::parseToAST(const std::string &expression){

  PyVector<int> infix = tokenize(expression);
  // Add intersection operators on demand
  addIntersection(infix);

  return parseToAST(infix);
}


/// \brief Parse an infix expression to an AST
/// \param expression an infix expression
/// \return the root of an AST
NodePtr BoolParser::parseToAST(const PyVector<int> &infix) {

  if (infix.empty())
    log::error("Cannot parse an empty infix expression");

  std::stack<NodePtr> ast;
  std::stack<int> operators;

  for (auto token: infix) {

    if (token < OP_UNION) {
      // Push any operands to the ast stack
      ast.push(new BoolNode(token));
    }
    else if (token == OP_LEFT_PAREN) {
      // Push left parentheses onto the stack
      operators.push(token);
    }
    else if (token < OP_RIGHT_PAREN) { // token: &,|,~
      while (!operators.empty()) {
        auto op = operators.top();
        if ( op < OP_RIGHT_PAREN &&
             (token != OP_COMPLEMENT && token <= op) ) {
          // Shunting when the token is of lower precedence and
          // left-associative
          shuntAST(ast, operators);
        }
        else
          break;
      }
      // There is no operator which has higher precedence than '~'.
      // Thus, operator '~' is always pushed
      operators.push(token);
    }
    else {  // token == OP_RIGHT_PAREN
      while ( !operators.empty() &&
              (operators.top() != OP_LEFT_PAREN) ) {
        // shunting
        shuntAST(ast, operators);
      }
      if (operators.empty()) {
        log::error("Mismatched parentheses in '{}'", infix);
      }
      else {
        operators.pop(); // Pop the left parenthesis.
      }
    }

  }

  // Add operators to the rpn in reverse order
  while (!operators.empty()) {
    auto op = operators.top();

    if (op >= OP_RIGHT_PAREN) {
      log::error("Mismatched parentheses in '{}'", infix);
    }
    // shunting
    shuntAST(ast, operators);
  }

  return ast.top();
}


/// \brief Parse an infix expression to an RPN
/// \param expression an infix expression
/// \return a vector storing an RPN
PyVector<int> BoolParser::parseToRPN(const std::string &expression){

  PyVector<int> infix = tokenize(expression);
  // Add intersection operators on demand
  addIntersection(infix);

  return parseToRPN(infix);
}


/// \brief Parse an infix expression to an RPN
/// \param infix an infix expression
/// \return a vector storing an RPN
PyVector<int> BoolParser::parseToRPN(const PyVector<int> &infix) {

  if (infix.empty())
    log::error("Cannot parse an empty infix expression");

  PyVector<int> rpn;
  std::stack<int> stack;

  for (auto token: infix) {

    if (token < OP_UNION) {
      // Add non-operator tokens immediately
      rpn.push_back(token);
    }
    else if (token == OP_LEFT_PAREN) {
      // Push left parentheses onto the stack
      stack.push(token);
    }
    else if (token < OP_RIGHT_PAREN) { // token: &,|,~

      while(!stack.empty()) {
        auto op = stack.top();
        if ( op < OP_RIGHT_PAREN &&
             (token != OP_COMPLEMENT && token <= op) ) {
          // Shunting when the token is of lower precedence
          // and left-associative
          rpn.push_back(op);
          stack.pop();
        }
        else
          break;
      }
      // There is no operator which has higher precedence than '~'.
      // Thus, operator '~' is always pushed
      stack.push(token);
    }
    else {  // token == OP_RIGHT_PAREN
      while ( !stack.empty() &&
              (stack.top() != OP_LEFT_PAREN) ) {
        // shunting
        rpn.push_back(stack.top());
        stack.pop();
      }
      if (stack.empty()) {
        log::error("Mismatched parentheses in '{}'", infix);
      }
      else {
        stack.pop(); // Pop the left parenthesis.
      }
    }

  }

  // Add operators to the rpn in reverse order
  while(!stack.empty()) {
    auto op = stack.top();

    if (op >= OP_RIGHT_PAREN) {
      log::error("Mismatched parentheses in '{}'", infix);
    }
    rpn.push_back(stack.top());
    stack.pop();
  }

  return rpn;
}


/// \brief Remove all the superfluous operators from a boolean expression
void BoolParser::removeOperators(PyVector<int> &expression) {
  size_t id = 0;
  size_t old_id = 0;
  while (old_id < expression.size()) {
    if (expression[old_id] < OP_UNION) {
      expression[id] = expression[old_id];
      ++id;
    }
    ++old_id;
  }
  expression.resize(id);
  expression.shrink_to_fit();
}


/// \brief Determine whether there is operator '|' or '~'
bool BoolParser::isSimple(const PyVector<int> &tokens) {
  bool simple = true;
  for (auto token : tokens) {
    if ( (token == BoolParser::OP_UNION) ||
         (token == BoolParser::OP_COMPLEMENT) ) {
      simple = false;
      break;
    }
  }
  return simple;
}


/// \brief Converts a token to string
std::string BoolParser::tokenToString(int token) {
  switch (token) {
    case OP_LEFT_PAREN :   return "(";
    case OP_RIGHT_PAREN :  return ")";
    case OP_COMPLEMENT :   return "~";
    case OP_INTERSECTION : return "&";
    case OP_UNION :        return "|";
    default :              return std::to_string(token);
  }
}


/// \brief Converts a set of tokens to string
std::string BoolParser::toString(PyVector<int> tokens){

  std::ostringstream os;
  for (auto token : tokens) {
    if (token < OP_UNION && token > 0)  // positive operands
      os << '+';
    os << tokenToString(token) << ' ';
  }

  return stringutils::rtrim(os.str());
}


} // namespace antmoc
