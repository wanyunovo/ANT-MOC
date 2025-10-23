#include "antmoc/VecParser.h"
#include "antmoc/log.h"
#include "antmoc/string_utils.h"

#include <algorithm>
#include <sstream>
#include <stdexcept>

namespace antmoc
{

/// \briefs Return flags to indicate the contents
/// \details 1 - vector, 2 - double, 4 - u
uint8_t Operand::getFlags() const {
  uint8_t flags = 0;
  if     (isVec())    flags += V;
  else if(isDouble()) flags += D;
  else              flags += U;
  return flags;
}


/// \brief Convert the underlying data to a vector of ints
/// \details Note that this algorithm may cause information loss
///          If the object holds a Python-like vector, it will be returned
///          with each element being truncated.
///          If the object holds a double or an unsigned int, a vector of
///          a single element will be returned.
PyVector<int> Operand::toIntVector() const {
  PyVector<int> v;
  if (isVec())
    std::transform(_v.begin(), _v.end(), std::back_inserter(v),
                   [](double d){ return (int)std::lround(d); }
                  );
  else if (isDouble())
    v.push_back((int)_d);
  else if (isUint())
    v.push_back(_u);
  return v;
}


/// \brief Convert the underlying data to a Python-like vector of doubles
/// \details If the object holds a Python-like vector, it will be returned.
///          If the object holds a double or an unsigned int, a vector of
///          a single element will be returned.
DblPyVec Operand::toDblPyVector() const {
  DblPyVec v(_v);
  if (isDouble())
    v.push_back(_d);
  else if (isUint())
    v.push_back((double)_u);
  return v;
}


/// \brief Represent the Operand object as a string
std::string Operand::toString() const {
  std::stringstream ss;
  if (isVec()) {
    for (size_t i = 0; i < _v.size(); ++i)
      ss << _v[i] << (i != _v.size()-1 ? " " : "");
  }
  else if (isDouble())
    ss << _d;
  else
    ss << _u;
  return ss.str();
}


int VecParser::getPrio(const char token) {
  int prio = -1;
  switch (token) {
    case '(': prio = -2; break;
    case '+': case '-': prio = 1; break;
    case '*': case '/': prio = 2; break;
    default:
      log::ferror("Nonexist operator: '%c'", token);
  }
  return prio;
}

bool VecParser::hasOperator(const char token) {
  return (op_list.find(token) != std::string::npos);
}

Operand VecParser::evaluate(const Operand &lhs, const Operand &rhs,
                 const char op) {
  Operand result;
  switch (op) {
    case '+': result = lhs + rhs; break;
    case '*': result = lhs * rhs; break;
    default:
      log::ferror("Nonexist operator: '%c'", op);
  }
  return result;
}

void VecParser::shunt() {
  if (operands.empty())
    log::ferror("Missing operands in vector expression");
  auto rhs = operands.top();
  operands.pop();

  if (operands.empty())
    log::ferror("Missing operands in vector expression");
  auto lhs = operands.top();
  operands.pop();

  operands.push(evaluate(lhs, rhs, operators.top()));
  operators.pop();
}

void VecParser::takeOperand(std::string &digits) {
  if (digits.find('.') != std::string::npos)
    operands.push(std::stod(digits));
  else
    operands.push((size_t)std::stoul(digits));
  digits.clear();
}


/// \brief Take a vector from the string stream
/// \details The input string stream is assumed to be a vector expression,
///          that is, something enclosed by a pair of square brackets '[]'.
///          However, the leading '[' was taken from it.
///          As a result, this method parses an expression in the form of
///
///           '...something...]'
/// \param ss a vector expression without leading '['
void VecParser::takeVector(std::istringstream &ss) {
  std::ostringstream expr;
  std::stack<char> brackets;  ///< A stack of brackets '[' and ']'

  char token;
  while (ss >> token) {
    if (token == '[')
      brackets.push(token);
    else if (token == ']') {
      if (brackets.empty())
        break;  // this is the end
      else
        brackets.pop();
    }
    expr << token;  // Push everything back except brackets
  }

  if (!brackets.empty())
    log::ferror("No matching '[' found: '%s'", ss.str());
  if (stringutils::isSpaces(expr.str()))
    log::ferror("Square brackets [] cannot be empty: '%s'", ss.str());

  // Now the expression is 'unpacked' and we could recursively
  // parse the contents between enclosed brackets.
  VecParser p;
  operands.push(p.parse(expr.str()));
}


void VecParser::clearStacks() {
  while (!operands.empty())
    operands.pop();

  while (!operators.empty())
    operators.pop();
}


Operand VecParser::getOperand() const {
  if (operands.empty())
    return Operand();
  else
    return operands.top();
}


DblPyVec VecParser::getVector() const {
  if (operands.empty())
    return DblPyVec();
  else
    return operands.top().getVector();
}


std::string VecParser::getString() const {
  if (operands.empty())
    return std::string("");
  else
    return operands.top().toString();
}


/// \brief Parse a vector expression
/// \details This method implements the standard shunting-yard
///          algorithm.
DblPyVec
VecParser::parse(const std::string &expression) {
  clearStacks();

  std::istringstream ss(expression);

  char token;
  bool hit_decimal = false;
  std::string digit_buf;
  while (ss >> token) {
    if (token == '[') { // parse the contents
      takeVector(ss);
    }
    else if (isdigit(token) || token == '.') {
      if (isdigit(token))
        digit_buf += token; // add to buffer
      else if (token == '.' && !hit_decimal) {
        hit_decimal = true;
        digit_buf += token; // add to buffer
      }
      else
        log::ferror("Ill-formed floating point: '%s'", expression);
    }
    else {
      hit_decimal = false;
      if (!digit_buf.empty())
        takeOperand(digit_buf);

      if (token == '(')
        operators.push(token);
      else if (token == ')') {
        while (!operators.empty() && operators.top() != '(') { // left associative
          shunt();
        }

        if (operators.empty())
          log::ferror("No matching '(': '%s'", expression);
        else
          operators.pop(); // pop '('
      }
      else if (hasOperator(token))
      {
        while (!operators.empty()) {
          if (getPrio(token) <= getPrio(operators.top())) {
            shunt();
          }
          else
            break;
        }
        operators.push(token);
      }
      else
        log::ferror("Syntax error: '%s'", expression);
    }
  }

  if (!digit_buf.empty())
    takeOperand(digit_buf);

  while (!operators.empty()) {
    shunt();
  }

  return operands.top().toDblPyVector();
}

} /* namespace antmoc */
