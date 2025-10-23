/// \file include/PyVector.h
/// \brief User-defined vector to facilitate processing
/// \date Aug 8, 2019
/// \author An Wang, USTB (wangan@xs.ustb.edu.cn)

#ifndef PYVECTOR_H_
#define PYVECTOR_H_

#define _USE_MATH_DEFINES
#include <cmath>
#include <sstream>
#include <string>
#include <vector>

namespace antmoc
{

using std::vector;

///---------------------------------------------------------------------
/// \class PyVector
/// \brief A derived vector imitating Python list
///---------------------------------------------------------------------
template <typename T>
class PyVector: public vector<T>
{
public:
  using vector<T>::vector;

  PyVector() = default;

  vector<T> get() const;
  std::string toString() const;

  template <typename U>
  friend std::ostream& operator<<(std::ostream &, const PyVector<U> &);

  // conversion constructor
  PyVector(const std::vector<T> &vec);

  // computed assignment
  PyVector<T>& operator+=(const PyVector<T>& other);
  PyVector<T>& operator+=(const T &rhs);
  PyVector<T>& operator+=(std::initializer_list<T> ilist);
  PyVector<T>& operator*=(size_t rep);

  // binary operators
  template <typename U>
  friend PyVector<U> operator*(const PyVector<U> &lhs, size_t rep);

  template <typename U>
  friend PyVector<U> operator*(size_t rep, const PyVector<U> &rhs);

  template <typename U>
  friend PyVector<U> operator+(const PyVector<U> &lhs,
                               const PyVector<U> &rhs);

  template <typename U, typename V>
  friend PyVector<U> operator+(const PyVector<U> &lhs, const V &rhs);

  template <typename U, typename V>
  friend PyVector<U> operator+(const V &lhs, const PyVector<U> &rhs);
};


/// \brief Get the equivalent vector
/// \return a vector equivalent to the PyVector
template <typename T>
inline vector<T> PyVector<T>::get() const {
  vector<T> v;
  v.reserve(this->size());
  v.insert(v.end(), this->begin(), this->end());
  return v;
}

/// \brief Convert the vector to a string
template <typename T>
inline std::string PyVector<T>::toString() const {
  std::stringstream ss;
  for (auto &e : *this)
    ss << e << ' ';
  auto s = ss.str();
  return s.substr(0, s.find_last_not_of(" ") + 1);
}

/// \brief Return the vector as a string by operator <<
template <typename T>
std::ostream& operator<<(std::ostream &os, const PyVector<T> &vec) {
  return os << vec.toString();
}

/// \brief Convert a normal vector to a Python-like vector
/// \param a vector of elements
template <typename T>
PyVector<T>::PyVector(const std::vector<T> &vec) {
  this->reserve(this->size() + vec.size());
  this->insert(this->end(), vec.begin(), vec.end());
}

/// \brief Addition assignment operator +=
/// \details Append another vector to this one.
template <typename T>
inline PyVector<T>&
PyVector<T>::operator+=(const PyVector<T>& other) {
  this->reserve(this->size() + other.size());
  this->insert(this->end(), other.begin(), other.end());
  return *this;
}

/// \brief Addition assignment operator +=
/// \param rhs an element of type T
template <typename T>
inline PyVector<T>&
PyVector<T>::operator+=(const T& rhs) {
  this->push_back(rhs);
  return *this;
}

/// \brief Addition assignment operator +=
/// \param ilist a list of initializers
template <typename T>
inline PyVector<T>&
PyVector<T>::operator+=(std::initializer_list<T> ilist) {
  return (*this) += PyVector<T>(ilist);
}


/// \brief Duplicate the vector many times
/// \details Duplication is performed by doubling the current
///          vector each time.
/// \param rep the number of copies
template <typename T>
inline PyVector<T>&
PyVector<T>::operator*=(size_t rep) {
  if (!rep) {
    this->clear();
  }
  else {
    // Backup
    PyVector<T> v(*this);
    this->reserve(this->size() * rep);

    // Rewrite rep to 2^k + r
    // where k = log2(rep) and r = rep^(1<<k)
    size_t k = static_cast<size_t>(std::log2(rep));

    for (size_t i=0; i < k; ++i)
      this->insert(this->end(), this->begin(), this->end());

    for (size_t i=0; i < (rep^(1<<k)); ++i)
      this->insert(this->end(), v.begin(), v.end());
  }
  return *this;
}


/// \brief Duplicate a given vector many times
/// \param lhs a given vector
/// \param rep the number of copies
template <typename T>
inline PyVector<T> operator*(const PyVector<T> &lhs, size_t rep) {
  PyVector<T> v(lhs);
  return (v *= rep);
}

/// \brief Duplicate a given vector many times
/// \param rep the number of copies
/// \param lhs a given vector
template <typename T>
inline PyVector<T> operator*(size_t rep, const PyVector<T> &rhs) {
  PyVector<T> v(rhs);
  return (v *= rep);
}

/// \brief Concatenate two vectors
template <typename T>
inline PyVector<T> operator+(const PyVector<T> &lhs,
                             const PyVector<T> &rhs) {
  PyVector<T> v(lhs);
  return (v += rhs);
}


/// \brief Append an element of type V to a vector of type U
template <typename U, typename V>
inline PyVector<U> operator+(const PyVector<U> &lhs, const V &rhs) {
  PyVector<U> v(lhs);
  return v += static_cast<U>(rhs);
}


/// \brief Insert an element at the front of a vector
/// \details Addition against an element and PyVector is not
///          commutative. Thus, this implementation trys to
///          insert the whole vector into another vector which
///          has only 1 element. The complexity is O(N).
template <typename U, typename V>
inline PyVector<U> operator+(const V &lhs, const PyVector<U> &rhs) {
  PyVector<U> v = {static_cast<U>(lhs)};
  return (v += rhs);
}


/// Type alias
using DblPyVec   = PyVector<double>;      ///< A Python-like vector of doubles
using DblPyVec2D = PyVector<DblPyVec>;
using DblPyVec3D = PyVector<DblPyVec2D>;
using IntPyVec   = PyVector<int>;         ///< A Python-like vector of ints
using IntPyVec2D = PyVector<IntPyVec>;
using IntPyVec3D = PyVector<IntPyVec2D>;
using StrPyVec   = PyVector<std::string>;   ///< A PyVector-like vector of strings


} // namespace antmoc

#endif  // PYVECTOR_H_
