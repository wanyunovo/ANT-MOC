/**
 * @file LocalCoords.h
 * @brief The LocalCoords class.
 * @date January 25, 2012
 * @author William Boyd, MIT, Course 22 (wboyd@mit.edu)
 */

#ifndef LOCALCOORDS_H_
#define LOCALCOORDS_H_

#include <string>

#include "antmoc/Point.h"

namespace antmoc
{

/* Forward declarationss to resolve circular dependencies */
class Cell;
class Lattice;
class Universe;

/**
 * @enum coordType
 * @brief The type of Universe level on which the LocalCoords reside
 */
enum coordType {
  /** A Universe level coordinate type */
  UNIV,

  /** A Lattice level coordinate type */
  LAT
};


/**
 * @class LocalCoords LocalCoords.h“openmoc/src/host/LocalCoords.h”
 * @brief LocalCoords 表示某些上的一组局部坐标
 *构成几何体的嵌套宇宙的级别。
 *_next 和 _prev 允许使用 LocalCoords 作为链表
 *但 _next_array 也可用于访问坐标。
 */
class LocalCoords {

private:
  /** 本地坐标类型（UNIV 或 LAT） */
  coordType _type;

  /** LocalCoords 所在的宇宙 */
  Universe* _universe;

  /** LocalCoords 所在的单元格 */
  Cell* _cell;

  /** LocalCoords 所在的格子 */
  Lattice* _lattice;

  /** LocalCoords 所在的晶格单元的第一个索引
   *  驻留*/
  int _lattice_x;

  /** LocalCoords 所在的晶格单元的第二个索引
   *  驻留*/
  int _lattice_y;

  /** LocalCoords 所在的晶格单元的第三个索引
   *  驻留*/
  int _lattice_z;

  /** 表示此 LocalCoords 的 3D 坐标的 Point */
  Point _coords;

  /** 相对于 x 轴的方向角（以弧度表示） */
  double _phi;

  /** 相对于 z 轴的方向角（以弧度表示） */
  double _polar;

  /** 指向下一个较低嵌套 Universe 级别的 LocalCoords 的指针 */
  LocalCoords* _next;

  /** 指向下一个更高嵌套 Universe 级别的 LocalCoords 的指针 */
  LocalCoords* _prev;

  /** 包含指向所有下一个 LocalCoords 的指针的数组
   * // 该局部点所在的链表表头 */
  LocalCoords* _next_array;

  /** _next_array 坐标中的位置 */
  int _position;

  /** _next_array 坐标的大小 */
  int _array_size;

  /** 一个整数，用于区分其他匹配的坐标 FSR 键 */
  int _version_num;

  void setArrayPosition(LocalCoords* array, int position, int array_size);

public:
  LocalCoords(double x=0.0, double y=0.0, double z=0.0, bool first=false);
  virtual ~LocalCoords();
  coordType getType();
  Universe* getUniverse() const;
  Cell* getCell() const;
  Lattice* getLattice() const;
  int getLatticeX() const;
  int getLatticeY() const;
  int getLatticeZ() const;
  double getX() const;
  double getY() const;
  double getZ() const;
  double getPhi() const;
  double getPolar() const;
  Point* getPoint();
  LocalCoords* getNext() const;
  LocalCoords* getNextCreate(double x, double y, double z);
  LocalCoords* getPrev() const;
  int getVersionNum();
  int getPosition();

  void setType(coordType type);
  void setUniverse(Universe* universe);
  void setCell(Cell* cell);
  void setLattice(Lattice* lattice);
  void setLatticeX(int lattice_x);
  void setLatticeY(int lattice_y);
  void setLatticeZ(int lattice_z);
  void setX(double x);
  void setY(double y);
  void setZ(double z);
  void setPhi(double phi);
  void setPolar(double polar);
  void setNext(LocalCoords *next);
  void setPrev(LocalCoords* coords);
  void setVersionNum(int version_num);

  LocalCoords* getLowestLevel();
  LocalCoords* getHighestLevel();
  void adjustCoords(double delta_x, double delta_y, double delta_z=0.0);
  void updateMostLocal(Point* point);
  void prune();
  void deleteArray();
  void copyCoords(LocalCoords* coords);
  std::string toString();
  void detectLoop();
};

} /* namespace antmoc */

#endif /* LOCALCOORDS_H_ */
