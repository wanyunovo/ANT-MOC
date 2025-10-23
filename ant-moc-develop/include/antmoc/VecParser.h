/**
 * @file VecParser.h
 * @details User-defined parser to facilitate processing
 * @date Aug 17, 2019
 * @author An Wang, USTB (wangan@xs.ustb.edu.cn)
 */

#ifndef VECPARSER_H_
#define VECPARSER_H_

#include <cstdint>
#include <stack>
#include <stdexcept>
#include <string>

#include "antmoc/PyVector.h"

namespace antmoc
{

///< \class A helpler class for the expression parser
///< \brief An Operand acts like a Union
class Operand {
public:
  Operand() = default;
  Operand(const size_t i): _u(i) { };
  Operand(const double d): _d(d) { };
  Operand(const DblPyVec &vec): _v(vec) { };

  // Copy operations
  Operand(const Operand &op): _u(op._u), _d(op._d), _v(op._v) { };
  Operand &operator=(const Operand &) = default;
  Operand &operator=(const size_t &);
  Operand &operator=(const double &);
  Operand &operator=(const DblPyVec &);
  // Move operations
  Operand(Operand &&op) noexcept: _u(op._u), _d(op._d), _v(std::move(op._v)) { }
  Operand &operator=(Operand &&) noexcept;

  ~Operand() = default;

  size_t getUint() const { return _u; }
  double getDouble() const { return _d; }
  DblPyVec getVector() const { return _v; }

  bool isVec() const { return !_v.empty(); }
  bool isDouble() const { return !isVec() && _d >= 0; }
  bool isUint() const { return !isVec() && !isDouble(); }
  uint8_t getFlags() const;

  // Conversions
  PyVector<int> toIntVector() const;
  DblPyVec toDblPyVector() const;
  std::string toString() const;

  /* Computed assignment */
  friend Operand operator*(const Operand &lhs, const Operand &rhs);
  friend Operand operator+(const Operand &lhs, const Operand &rhs);

private:
  size_t _u = 0;    //< Holds an unsigned integer
  double _d = -1;   //< Holds a double (make it positive to be valid)
  DblPyVec _v;      //< Holds a Python-like vector

  // The following table takes 6 bits to determine the types of two
  // objects of Operand.
  static constexpr uint8_t V = 1;
  static constexpr uint8_t D = 2;
  static constexpr uint8_t U = 4;

  static constexpr uint8_t VV = V+(V<<3);
  static constexpr uint8_t VD = V+(D<<3);
  static constexpr uint8_t VU = V+(U<<3);
  static constexpr uint8_t DV = D+(V<<3);
  static constexpr uint8_t DD = D+(D<<3);
  static constexpr uint8_t DU = D+(U<<3);
  static constexpr uint8_t UV = U+(V<<3);
  static constexpr uint8_t UD = U+(D<<3);
  static constexpr uint8_t UU = U+(U<<3);
};

inline Operand &
Operand::operator=(const size_t &u) {
  _u = u;
  _d = -1;
  _v.clear();
  return *this;
}

inline Operand &
Operand::operator=(const double &d) {
  _d = d;
  _v.clear();
  return *this;
}

inline Operand &
Operand::operator=(const DblPyVec &v) {
  _v = v;
  return *this;
}


// Move-assignment operator
inline Operand&
Operand::operator=(Operand &&rhs) noexcept {
  if (this != &rhs) {
    _u = rhs._u;
    _d = rhs._d;
    _v = std::move(rhs._v);
  }
  return *this;
}


///< \brief Multiply two objects of Operand
///< \details Only 6 out of 9 cases are supported.
///<          Vector-double, double-vector and vector-vector multiplication
///<          are not intended to be implemented in class Operand.
///<          The behaviour of multiplication between a PyVector and an int
///<          is defined by PyVector itself, which means it could be just
///<          duplication or something else.
///< \see PyVector
///< \param lhs left-hand side
///< \param rhs right-hand side
inline Operand
operator*(const Operand &lhs, const Operand &rhs) {
  auto flags = lhs.getFlags() + (rhs.getFlags() << 3);
  Operand op;
  switch (flags) {
    case Operand::VU: op = lhs._v * rhs._u; break;
    case Operand::DD: op = double(lhs._d * rhs._d); break;
    case Operand::DU: op = double(lhs._d * rhs._u); break;
    case Operand::UV: op = lhs._u * rhs._v; break;
    case Operand::UD: op = double(lhs._u * rhs._d); break;
    case Operand::UU: op = size_t(lhs._u * rhs._u); break;
    default:
      throw std::logic_error("Cannot multiply vector by vector or double");
  }
  return op;
}

///< \brief Add up two objects of Operand
///< \details The behaviour of addition involving objects of PyVector is
///<          defined by PyVector. Normally, it never changes existing
///<          elements of the vector but expands or duplicates the vector.
///< \see PyVector
///< \param lhs left-hand side
///< \param rhs right-hand side
inline Operand
operator+(const Operand &lhs, const Operand &rhs) {
  auto flags = lhs.getFlags() + (rhs.getFlags() << 3);
  Operand op;
  switch (flags) {
    case Operand::VV: op = lhs._v + rhs._v; break;
    case Operand::VD: op = lhs._v + rhs._d; break;
    case Operand::VU: op = lhs._v + rhs._u; break;
    case Operand::DV: op = lhs._d + rhs._v; break;
    case Operand::DD: op = double(lhs._d + rhs._d); break;
    case Operand::DU: op = double(lhs._d + rhs._u); break;
    case Operand::UV: op = lhs._u + rhs._v; break;
    case Operand::UD: op = double(lhs._u + rhs._d); break;
    case Operand::UU: op = size_t(lhs._u + rhs._u); break;
    default:
      throw std::logic_error("Undefined error");
  }
  return op;
}


////////////////////////////////////////////////////////////////////////
/// \class An expression parser
/// A vector expression is a string including brackets ('[' and ']').
/// So, basically, this class is designed to imitate Python lists.
/// Expressions between brackets will be considered as objects of
/// PyVector. For example, "[0.43]" will be treated as a vector of a
/// single element.
/// Some kinds of arithmetic are supported, including multiplication
/// between a vector and an int, addition between a vector and a single
/// number, and arithmetic of doubles and ints.
/// For example, "([0.43] + [1.26])*2" gives a vector of size 4
///   0.43, 1.26, 0.43, 1.26
////////////////////////////////////////////////////////////////////////
class VecParser {
  const std::string op_list = "+-*/";   ///< A string of supported operators
  std::stack<char> operators;           ///< A stack of operators (shunting-yard)
  std::stack<Operand> operands;         ///< A stack of operands (shunting-yard)

  int getPrio(const char token);
  bool hasOperator(const char token);
  void shunt();
  Operand evaluate(const Operand &lhs, const Operand &rhs,
                   const char op);
  void takeOperand(std::string &digits);
  void takeVector(std::istringstream &ss);
  void clearStacks();

public:
  VecParser() = default;

  // Copy operations
  VecParser(const VecParser &p)
    : operators(p.operators), operands(p.operands) { };
  VecParser &operator=(const VecParser &);

  // Move operations
  VecParser(VecParser &&) noexcept ;
  VecParser &operator=(VecParser &&) noexcept;

  ~VecParser() = default;

  Operand getOperand() const;
  DblPyVec getVector() const;
  std::string getString() const;

  DblPyVec parse(const std::string &expression);

};

/// \brief Copy-assignment operator
inline
VecParser &VecParser::operator=(const VecParser &other) {
  operators = other.operators;
  operands = other.operands;
  return *this;
}

/// \brief Move constructor
inline VecParser::VecParser(VecParser &&other) noexcept {
  operators = std::move(other.operators);
  operands = std::move(other.operands);
}

/// \brief Move-assignment operator
inline VecParser &
VecParser::operator=(VecParser &&other) noexcept {
  if (this != &other) {
    operators = std::move(other.operators);
    operands = std::move(other.operands);
  }
  return *this;
}

} /* namespace antmoc */

#endif /* VECPARSER_H_ */
