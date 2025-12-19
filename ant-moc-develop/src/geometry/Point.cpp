#include "antmoc/Point.h"
#include "antmoc/math_utils.h"

#include <sstream>

namespace antmoc {


/// \brief 用二维坐标初始化一个点。
/// \param x x 坐标
/// \param y y 坐标
/// \param z z 坐标
void Point::setCoords(const double x, const double y, const double z) {
  _xyz[0] = x;
  _xyz[1] = y;
  _xyz[2] = z;
}


/// \brief Returns this Point's x-coordinate.
/// \return the x-coordinate
double Point::getX() const {
  return _xyz[0];
}


/// \brief Returns this Point's y-coordinate.
/// \return the y-coordinate
double Point::getY() const {
  return _xyz[1];
}


/// \brief Returns this Point's z-coordinate.
/// \return the z-coordinate
double Point::getZ() const {
  return _xyz[2];
}


/// \brief Returns this Point's coordinates.
/// \return the xyz coordinates
double* Point::getXYZ() {
  return _xyz;
}


/// \brief Set the Point's x-coordinate.
/// \param x the new x-coordinate
void Point::setX(const double x) {
  _xyz[0] = x;
}


/// \brief Set the Point's y-coordinate.
/// \param y the new y-coordinate
void Point::setY(const double y) {
  _xyz[1] = y;
}


/// \brief Set the Point's z-coordinate.
/// \param z the new z-coordinate
void Point::setZ(const double z) {
  _xyz[2] = z;
}


/// \brief Set the Point's x, y and z-coordinates.
/// \param xyz array with the three coordinates
void Point::setXYZ(double* xyz) {
  _xyz[0] = xyz[0];
  _xyz[1] = xyz[1];
  _xyz[2] = xyz[2];
}


/// \brief Compute the distance from this Point to another Point of interest.
/// \param x the x-coordinate of the Point of interest
/// \param y the y-coordinate of the Point of interest
/// \param z the z-coordinate of the Point of interest
/// \return distance to the Point of interest
double Point::distance(const double x, const double y, const double z) const {
  double deltax = _xyz[0] - x;
  double deltay = _xyz[1] - y;
  double deltaz = _xyz[2] - z;
  return sqrt(deltax*deltax + deltay*deltay + deltaz*deltaz);
}


/// \brief Compute the distance from this Point to another Point of interest.
/// \param point a pointer to the Point of interest
/// \return distance to the Point of interest
double Point::distanceToPoint(const Point* point) {
  return distance(point->getX(), point->getY(), point->getZ());
}


/// \brief Copy the coordinates from another point.
/// \param point a pointer to the Point that has the coordinates of interest
void Point::copyCoords(Point* point) {
  _xyz[0] = point->getX();
  _xyz[1] = point->getY();
  _xyz[2] = point->getZ();
}


//----------------------------------------------------------------------
// Arithmetic
//----------------------------------------------------------------------

/// \brief Addition assignment operator +=
/// \details Add another point to this one.
Point& Point::operator+=(const Point& other) {
  this->_xyz[0] += other._xyz[0];
  this->_xyz[1] += other._xyz[1];
  this->_xyz[2] += other._xyz[2];
  return *this;
}


/// \brief Subtraction assignment operator +=
/// \details Subtract another point from this one.
Point& Point::operator-=(const Point& other) {
  this->_xyz[0] -= other._xyz[0];
  this->_xyz[1] -= other._xyz[1];
  this->_xyz[2] -= other._xyz[2];
  return *this;
}


/// \brief Add up two points
Point operator+(const Point &lhs, const Point &rhs) {
  Point pt(lhs);
  return (pt += rhs);
}


/// \brief Subtract one point from another
Point operator-(const Point &lhs, const Point &rhs) {
  Point pt(lhs);
  return (pt -= rhs);
}


/// \brief Equality operator
bool operator==(const Point &lhs, const Point &rhs) {
  return definitelyEqual(lhs._xyz[0], rhs._xyz[0]) &&
         definitelyEqual(lhs._xyz[1], rhs._xyz[1]) &&
         definitelyEqual(lhs._xyz[2], rhs._xyz[2]);
}


/// \brief Equality operator
bool operator!=(const Point &lhs, const Point &rhs) {
  return !(lhs == rhs);
}


/// \brief Converts this Point to a character representation of its attributes.
/// \details The character array includes the x-coordinate, y-coordinate, 
///          and z-coordinate
/// \return a character array of this Point's attributes
std::string Point::toString() const {
  std::stringstream string;
  string << "Point: x = " << _xyz[0] << ", y = " << _xyz[1] << ", z = " << _xyz[2];
  return string.str();
}


std::ostream& operator<<(std::ostream &os, const Point &point) {
  return os << point.toString();
}

} // namespace antmoc
