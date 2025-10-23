/// \file Point.h
/// \brief The Point class.
/// \date January 18, 2012
/// \author William Boyd, MIT, Course 22 (wboyd@mit.edu)
/// \author An Wang, USTB (wangan@xs.ustb.edu.cn)

#ifndef POINT_H_
#define POINT_H_

#include <array>
#include <ostream>
#include <string>


namespace antmoc {


/// \class Point Point.h "include/Point.h"
/// \brief Class to represent a 2D/3D point in space.
class Point {

private:

  // The Point's xyz coordinates
  double _xyz[3];

public:

  /// \brief Construct a point by three values
  Point(const double x=0., const double y=0., const double z=0.)
    : _xyz{x, y, z} { };

  /// \brief Construct a point by an array
  Point(const std::array<double, 3> &xyz)
    : _xyz{xyz[0], xyz[1], xyz[2]} { };

  /// \brief Copy constructor
  Point(const Point& point)
    : _xyz{point._xyz[0], point._xyz[1], point._xyz[2]} { };

  Point &operator=(const Point& point) {
    _xyz[0] = point._xyz[0];
    _xyz[1] = point._xyz[1];
    _xyz[2] = point._xyz[2];
    return *this;
  }

  virtual ~Point() = default;

  void setCoords(const double x, const double y, const double z=0.0);
  double getX() const;
  double getY() const;
  double getZ() const;
  double* getXYZ();
  void setX(const double x);
  void setY(const double y);
  void setZ(const double z);
  double distance(const double x, const double y, const double z) const;
  double distanceToPoint(const Point* point);
  void copyCoords(Point* point);
  void setXYZ(double* xyz);

  // computed assignment
  Point& operator+=(const Point& other);
  Point& operator-=(const Point& other);

  // binary operators
  friend Point operator+(const Point &lhs, const Point &rhs);
  friend Point operator-(const Point &lhs, const Point &rhs);
  friend bool operator==(const Point &lhs, const Point &rhs);
  friend bool operator!=(const Point &lhs, const Point &rhs);

  // convert to string
  std::string toString() const;
  friend std::ostream& operator<<(std::ostream &, const Point &);
};


} // namespace antmoc

#endif /* POINT_H_ */
