/// \file Lattice.h
/// \brief The Lattice derivative class.
/// \details Refactored on April 18, 2019, in HPC&DE, USTB
/// \date January 9, 2012
/// \author William Boyd, MIT, Course 22 (wboyd@mit.edu)

#ifndef LATTICE_H_
#define LATTICE_H_

#include <limits>
#include <map>
#include <set>
#include <vector>

#include "antmoc/constants.h"
#include "antmoc/enum_types.h"
#include "antmoc/container_utils.h"
#include "antmoc/math_utils.h"
#include "antmoc/Point.h"
#include "antmoc/Universe.h"

namespace antmoc
{

  /** Forward declarations */
  class Cell;
  class LocalCoords;
  class LatticeIter;

  ///---------------------------------------------------------------------
  /// \class Lattice Lattice.h "include/Lattice.h"
  /// \brief Represents a repeating 3D Lattice of Universes.
  ///---------------------------------------------------------------------
  class Lattice : public Universe
  {

  protected:
    /** Lattice的布局类型（枚举类型 latticeType，可能值：RECTANGULAR矩形、HEXAGONAL六边形等）
     *  用于标识当前晶格是矩形排列还是六边形排列 */
    latticeType _lattice_type;

    /** Lattice在父几何体中的偏移坐标（Point类型包含x、y、z三个坐标值）
     *  表示此Lattice的几何中心在父Universe坐标系中的位置
     *  例如：_offset = (1.5, 2.0, 0.0) 表示lattice中心位于父坐标系的(1.5, 2.0, 0.0)处
     *  作用：坐标转换时使用，全局坐标 = 局部坐标 + _offset
     * */
    Point _offset;

    /** 标识晶格是否为非均匀网格（布尔类型）
     *  true  = 非均匀晶格（每个格子可以有不同的尺寸，需要用_widths数组存储）
     *  false = 均匀晶格（所有格子尺寸相同，只需要_width单个值） */
    bool _non_uniform;

    /** 沿z轴（高度/轴向）的Lattice层数
     *  例如：_num_z = 10 表示z方向堆叠了10层 */
    int _num_z;

    /** z方向每层的统一高度（单位：cm），仅用于均匀晶格
     *  例如：_width_z = 5.0 表示每层高5.0cm
     *  特殊值：infinity（无穷大）表示2D晶格（不沿z方向分层）
     * */
    double _width_z;

    /** z方向各层的高度数组，用于非均匀晶格
     *  长度为 _num_z，_widths_z[i] 存储第i层的高度
     *  例如：_widths_z = {5.0, 8.0, 5.0} 表示3层分别高5.0、8.0、5.0cm
     * */
    DoubleVec _widths_z;

    /** z方向的累积高度数组（前缀和），用于快速定位在哪一层
     *  长度为 _num_z+1，_accumulate_z[i] 表示从底部到第i层底面的累积高度
     *  工作原理与x/y方向的_accumulate完全相同
     * */
    DoubleVec _accumulate_z;

    /** 存储所有子Universe指针的容器
     *  _universes[i] 指向第i个lattice单元中填充的Universe对象
     *  Universe*：指针类型，指向Universe对象
     *  长度等于晶格单元总数，对于RecLattice是 _num_x * _num_y * _num_z */
    std::vector<Universe *> _universes;

  public:
    //-----------------------------------------------------
    // Iterators
    //-----------------------------------------------------
    friend class LatticeIter; // LatticeIter为友元类，可访问Lattice的保护/私有成员

    virtual LatticeIter begin() = 0; // 纯虚函数：返回Lattice迭代器的起始位置（子类必须实现）
    LatticeIter end();

    //-----------------------------------------------------
    // Copy-control members
    //-----------------------------------------------------
    Lattice(const int id = -1, const char *name = "");
    virtual ~Lattice() = default;

    //-----------------------------------------------------
    // Non-virtual functions implemented in the base class
    //-----------------------------------------------------
    latticeType getLatticeType() const { return _lattice_type; }
    void setLatticeType(latticeType type) { _lattice_type = type; }

    Point *getOffset() { return &_offset; }
    void setOffset(double x, double y, double z = 0.0);
    void setOffset(Point offset);

    /// \brief Determine if this is a 3-D lattice
    bool is3D() const
    {
      return _width_z != std::numeric_limits<double>::infinity();
      // 只要 _width_z 不等于无穷大，就代表每层有确定宽度，整个 Lattice 需要沿 z 轴重复，从而视为 3D。
      // Lattice 既支持纯 2D（只在 x-y 平面铺展）也支持 3D（沿 z 方向分层）。
    }

    // Utils for CSG construction
    Cell *findCell(LocalCoords *coords);               // 找到局部点所在的 Cell
    void subdivideCells(double max_radius = INFINITY); // 将 Lattice 中所有以材料填充的 Cell,按照垂直于 z 轴方向,进行环状和扇区分割
    void buildNeighbors();                             // 为了优化特征线跟踪,对 Lattice 中 Universe 相关的全部 Cell,  均为其建立邻居 Cell 集合

    // Control child universes
    void clearUniverses();
    std::vector<Universe *> *getUniverses();
    std::map<int, Universe *> getUniqueUniverses();
    std::map<int, Cell *> getAllCells(); // Lattice 中 Universe 相关的全部 Cell
    std::map<int, Universe *> getAllUniverses();
    Universe *getUniverse(int lat_x, int lat_y, int lat_z = 0);
    void updateUniverse(int lat_x, int lat_y, int lat_z, Universe *universe);
    void removeUniverse(Universe *universe);

    // The same layout along z-axis
    bool getNonUniform() const { return _non_uniform; }
    virtual int getNumX() const = 0;
    virtual int getNumY() const = 0;
    int getNumZ() const { return _num_z; }
    double getWidthZ() const { return _width_z; }
    const DoubleVec &getWidthsZ() const { return _widths_z; }
    const DoubleVec &getAccumulateZ() const { return _accumulate_z; }
    void setNumZ(int num_z) { _num_z = num_z; }
    void setWidthsZ(DoubleVec widthsz) { _widths_z = widthsz; }
    void setAccumulateZ(DoubleVec accumulatez) { _accumulate_z = accumulatez; }
    virtual void computeSizes() = 0;

    virtual void printString();

    //------------------------------------------------------
    // Interfaces to be implemented in each derivative class
    //------------------------------------------------------
    virtual int getNumLatticeCells() const = 0;
    int getMaxNumUniverses() const;

    virtual bool containsPoint(Point *point) = 0;                                            // 判断点是否在 Lattice 内
    virtual double minSurfaceDist(Point *point, double azim, double polar = M_PI / 2.0) = 0; // 计算局部点沿某方向到最近面的距离

    virtual double getMinX() const = 0;
    virtual double getMaxX() const = 0;
    virtual double getMinY() const = 0;
    virtual double getMaxY() const = 0;
    double getMinZ() const; // implemented
    double getMaxZ() const; // implemented

    virtual std::map<int, double>
    getUniqueRadius(std::map<int, Universe *> &unique_universes) = 0;

    virtual bool areValidIndices(int lat_x, int lat_y, int lat_z) const = 0;
    virtual bool isValidIndex(int index) const = 0;
    virtual std::array<int, 2> getLatXY(Point *point) = 0;
    virtual int getLatX(Point *point) = 0;
    virtual int getLatY(Point *point) = 0;
    int getLatZ(Point *point); // implemented

    int getLatXByIndex(int indx);                  // implemented
    int getLatYByIndex(int indx);                  // implemented
    int getLatZByIndex(int indx);                  // implemented
    std::array<int, 3> getLatXYZByIndex(int indx); // implemented

    virtual double getLocalPointX(const Point *p, int lat_x, int lat_y) = 0;
    virtual double getLocalPointY(const Point *p, int lat_x, int lat_y) = 0;
    double getLocalPointZ(const Point *p, int lat_z); // implemented
    // virtual double toLocalPoint(const Point *p, int lat_x, int lat_y, int lat_z) = 0;

    virtual int getLatticeCell(int lat_x, int lat_y, int lat_z) const = 0;
    virtual int getLatticeCell(Point *point) = 0;
    virtual int getLatticeSurface(int cell, Point *point, double azim, double polar) = 0;
    virtual int getLatticeSurfaceOTF(int cell, double z, int surface_2D) = 0;

    // Get the center point of a lattice cell
    virtual Point getLatticeCellCenter(int lat_x, int lat_y, int lat_z) = 0;
    Point getLatticeCellCenter(int index); // implemented

    // Get vertices of a lattice cell
    virtual std::vector<Point> getLatticeCellVertices(int index) = 0;

    virtual std::string toString() = 0;
    virtual std::string dimToString() = 0;

    // Virtual functions used with HexLattice and RecLattice
    virtual void setNumX(int num_x) {};
    virtual void setNumY(int num_y) {};
    virtual void setNumR(int num_r) {};
    // RecLattice
    virtual void setWidths(DoubleVec _cell_widths_x, DoubleVec _cell_widths_y, DoubleVec _cell_widths_z) {};
    virtual void setWidth(double _cell_width_x, double _cell_width_y, double _cell_width_z) {};
    // HexLattice
    virtual void setOrientation(std::string orientation) {};
    virtual void setOrientation(Orientation orientation) {};
    virtual void setWidth(double width_r, double width_z = std::numeric_limits<double>::infinity()) {}
    virtual void setWidths(double width_r, DoubleVec widths_z) {};
  };

  ///---------------------------------------------------------------------
  /// \class An iterator over lattice universes
  ///---------------------------------------------------------------------
  class LatticeIter
  {
  public:
    /// \brief Initialized from a lattice and the index of a lattice cell
    LatticeIter(Lattice &lat, size_t pos)
        : _lat(lat), _pos(pos)
    {
    }

    /// \brief Two LatticeIters equal iff they are at the same valid position
    bool operator==(const LatticeIter &rhs) { return (_pos == rhs._pos); }

    bool operator!=(const LatticeIter &rhs) { return !(*this == rhs); }

    /*
    重载解引用运算符（像指针的 *ptr）
    返回lattice中当前位置的一个指向Universe的指针
    &：引用，返回的是引用而非副本
    为什么返回引用：可以修改这个 Universe 指针，指向别的 Universe
    */
    Universe *&operator*() { return _lat._universes[_pos]; }

    /// \brief Move to the next valid lattice cell
    /// \details This overloaded operator iterates cells of the given
    ///          lattice and checks the validation
    /*
    返回 LatticeIter &：返回自身引用（标准做法）
    */
    LatticeIter &operator++()
    {
      // Find the next valid index
      while (_pos < _lat._universes.size())
      {
        ++_pos;
        if (_lat.isValidIndex(_pos)) // 调用 Lattice 的 isValidIndex() 检查索引是否有效
          return *this;              // 返回lattice中当前位置的一个指向Universe的指针
      }
      // The index past the last element
      _pos = _lat._universes.size();
      return *this;
    }

    /// \brief Return the current position of the iterator
    size_t getPos() const { return _pos; }

  protected:
    Lattice &_lat; /// 要迭代的 Lattice
    size_t _pos;   ///< Index of the current lattice cell
  };

  ///---------------------------------------------------------------------
  /// \class RecLattice Lattice.h "include/Lattice.h"
  /// \brief Represents a repeating 3D rectangular Lattice of Universes.
  ///---------------------------------------------------------------------
  class RecLattice : public Lattice
  {

  private:
    /** 沿x轴方向的Lattice单元格数量（整数）
     *  例如：_num_x = 5 表示x方向有5个格子排列 */
    int _num_x;

    /** 沿y轴方向的Lattice单元格数量（整数）
     *  例如：_num_y = 5 表示y方向有5个格子排列 */
    int _num_y;

    // --- 均匀网格（uniform lattice）：所有格子宽度相同 ---
    /** x方向每个格子的统一宽度（单位：cm）
     *  仅用于均匀晶格（所有格子宽度一样）
     *  例如：_width_x = 1.26 表示每个格子x方向都是1.26cm宽 */
    double _width_x;

    // --- 非均匀网格（non-uniform lattice）：每个格子宽度可以不同 ---
    /** x方向各个格子的宽度数组（DoubleVec是std::vector<double>的类型别名）
     *  长度为 _num_x，_widths_x[i] 存储第i个格子的宽度
     *  例如：_widths_x = {1.0, 1.5, 2.0} 表示3个格子分别宽1.0、1.5、2.0cm */
    DoubleVec _widths_x;

    /** x方向的累积宽度数组（前缀和数组，用于快速定位）
     *  长度为 _num_x+1，_accumulate_x[i] 是第 i 条纵向网格线相对于最小 x 的位置
     *  _accumulate_x[0] = -2（起点）
     *  _accumulate_x[i] = _widths_x[0] + _widths_x[1] + ... + _widths_x[i-1]（前i个格子的总宽）
     *  _accumulate_x[_num_x] = 整个lattice的总宽度
     *
     *  作用：给定一个x坐标，可通过二分查找快速确定它在哪个格子里
     *  例如：_widths_x = {1.0, 1.5, 2.0}
     *       则 _accumulate_x = {-2.0, -1.0, 0.5, 2.5}
     *       表示网格线位置在-2.0、-1.0、0.5、2.5cm处
     * */
    DoubleVec _accumulate_x;

    // --- Y方向的对应变量（含义与X方向完全相同）---
    /** y方向每个格子的统一宽度（cm），仅用于均匀晶格 */
    double _width_y;

    /** y方向各个格子的宽度数组，用于非均匀晶格 */
    DoubleVec _widths_y;

    /** y方向的累积宽度数组（前缀和），用于快速定位格子 */
    DoubleVec _accumulate_y;

  public:
    RecLattice(const int id = -1, const char *name = "");
    ~RecLattice() = default;

    //-----------------------------------------------------
    // Implemented interfaces
    //-----------------------------------------------------
    friend class LatticeIter;

    LatticeIter begin();

    /// \brief Compute the number of cells of the HexLattice
    int getNumLatticeCells() const
    {
      return _num_x * _num_y * _num_z;
    }

    bool containsPoint(Point *point);
    double minSurfaceDist(Point *point, double azim, double polar = M_PI / 2.0);

    std::map<int, double>
    getUniqueRadius(std::map<int, Universe *> &unique_universes);

    void setUniverses(int num_z, int num_y, int num_x, Universe **universes);

    double getMinX() const;
    double getMaxX() const;
    double getMinY() const;
    double getMaxY() const;

    bool areValidIndices(int lat_x, int lat_y, int lat_z) const;
    bool isValidIndex(int index) const;
    std::array<int, 2> getLatXY(Point *point);
    int getLatX(Point *point);
    int getLatY(Point *point);

    /// \brief Get the x-coord in the local system of a universe
    // 这两个函数把全局坐标 p 转成“当前 lattice 单元内部”的局部坐标：先用 getMinX() 得到 lattice 左边界，再加上 _accumulate_x[lat_x]（对应单元左边界距最左端的偏移）和该单元宽度的一半就得到该单元中心的全局 x 值；用 p->getX() 减去这个中心，就得到点相对于该单元中心的局部 x 坐标。y 方向完全对应。
    double getLocalPointX(const Point *p, int lat_x, int lat_y)
    {
      return p->getX() - (getMinX() + _widths_x[lat_x] / 2. + _accumulate_x[lat_x]);
    }

    /// \brief Get the y-coord in the local system of a universe
    double getLocalPointY(const Point *p, int lat_x, int lat_y)
    {
      return p->getY() - (getMinY() + _widths_y[lat_y] / 2. + _accumulate_y[lat_y]);
    }

    // double toLocalPoint(const Point *p, int lat_x, int lat_y, int lat_z);

    /// \brief Compute the uid of a lattice cell by its indices
    int getLatticeCell(int lat_x, int lat_y, int lat_z) const
    {
      return _num_y * _num_x * lat_z + _num_x * lat_y + lat_x;
      // 是把三维索引用成一维索引用的公式，遵循 k*_num_x*_num_y + j*_num_x + i 的行主序规则：lat_z 决定层偏移，lat_y 决定在层内的行偏移，lat_x 是列索引。
    }
    /// \brief Compute the uid of a lattice cell where a point resides
    int getLatticeCell(Point *point);

    int getLatticeSurface(int cell, Point *point, double azim = 0, double polar = M_PI / 2.0);
    int getLatticeSurfaceOTF(int cell, double z, int surface_2D);

    // Non-uniform lattice cells along z-axis
    void computeSizes();

    // Get the center point of a lattice cell
    Point getLatticeCellCenter(int lat_x, int lat_y, int lat_z);
    using Lattice::getLatticeCellCenter;

    // Get vertices of a lattice cell
    std::vector<Point> getLatticeCellVertices(int index);

    std::string toString();
    std::string dimToString();

    //-----------------------------------------------------
    // Member functions
    //-----------------------------------------------------
    int getNumX() const { return _num_x; }
    int getNumY() const { return _num_y; }
    double getWidthX() const { return _width_x; }
    double getWidthY() const { return _width_y; }
    const DoubleVec &getWidthsX() const { return _widths_x; }
    const DoubleVec &getWidthsY() const { return _widths_y; }
    const DoubleVec &getAccumulateX() const { return _accumulate_x; }
    const DoubleVec &getAccumulateY() const { return _accumulate_y; }

    void setNumX(int num_x) { _num_x = num_x; }
    void setNumY(int num_y) { _num_y = num_y; }
    void setNonUniform(bool non_uniform) { _non_uniform = non_uniform; }
    void setWidthsX(DoubleVec widthsx) { _widths_x = widthsx; }
    void setWidthsY(DoubleVec widthsy) { _widths_y = widthsy; }
    void setAccumulateX(DoubleVec accumulatex) { _accumulate_x = accumulatex; }
    void setAccumulateY(DoubleVec accumulatey) { _accumulate_y = accumulatey; }

    /* Set XYZ widths of uniform meshes */
    void setWidth(double width_x, double width_y,
                  double width_z = std::numeric_limits<double>::infinity());

    /* Set XYZ widths of non-uniform meshes */
    void setWidths(DoubleVec widths_x, DoubleVec widths_y,
                   DoubleVec widths_z);

    /* For debug use */
    void printLatticeSizes();
  };

  ///---------------------------------------------------------------------
  /// \file  Lattice.h "include/Lattice.h"
  /// \class HexLattice
  /// \brief Represents a repeating 3D hexagonal Lattice of Universes.
  ///---------------------------------------------------------------------
  class HexLattice : public Lattice
  {

  private:
    Orientation _orientation; ///< Orientation of the lattice
    int _num_r;               ///< The number of radial tiles
    double _width_r;          ///< The radial pitch of the lattice

  public:
    HexLattice(const int id = -1, const char *name = "");
    ~HexLattice() = default;

    //-----------------------------------------------------
    // Implemented interfaces
    //-----------------------------------------------------
    friend class LatticeIter;

    LatticeIter begin();

    /// \brief Compute the number of cells of the HexLattice
    int getNumLatticeCells() const
    {
      return (3 * _num_r * (_num_r - 1) + 1) * _num_z;
    }

    bool containsPoint(Point *point);
    double minSurfaceDist(Point *point, double azim, double polar = M_PI / 2.0);

    std::map<int, double>
    getUniqueRadius(std::map<int, Universe *> &unique_universes);

    void setUniverses(int num_z, int num_r, Universe **universes);
    void setUniversesY(int num_z, int num_r, Universe **universes);
    void setUniversesX(int num_z, int num_r, Universe **universes);

    double getMinX() const;
    double getMaxX() const;
    double getMinY() const;
    double getMaxY() const;

    bool areValidIndices(int lat_x, int lat_y, int lat_z) const;
    bool isValidIndex(int index) const;
    std::array<int, 2> getLatXY(Point *point);
    int getLatX(Point *point);
    int getLatY(Point *point);

    double getLocalPointX(const Point *p, int lat_x, int lat_y);
    double getLocalPointY(const Point *p, int lat_x, int lat_y);
    // double toLocalPoint(const Point *p, int lat_x, int lat_y, int lat_z);

    int getLatticeCell(int lat_x, int lat_y, int lat_z) const;
    int getLatticeCell(Point *point);
    int getLatticeSurface(int cell, Point *point, double azim, double polar = M_PI / 2.0);
    int getLatticeSurfaceOTF(int cell, double z, int surface_2D);

    // Non-uniform lattice cells along z-axis
    void computeSizes();

    // Get the center point of a lattice cell
    Point getLatticeCellCenterLocal(int lat_x, int lat_y, int lat_z);
    Point getLatticeCellCenter(int lat_x, int lat_y, int lat_z);
    using Lattice::getLatticeCellCenter;

    std::vector<Point> getLatticeCellVertices(int index);

    std::string toString();
    std::string dimToString();

    //-----------------------------------------------------
    // Member functions
    //-----------------------------------------------------
    /// \brief Returns the number of radial rings
    int getNumR() const { return _num_r; }

    /// \brief Returns the radial pitch
    double getWidthR() const { return _width_r; }

    /// \brief Returns the number of cells including invalid ones
    int getNumX() const { return 2 * _num_r - 1; }
    int getNumY() const { return 2 * _num_r - 1; }

    /// \brief Returns the orientation
    Orientation getOrientation() const { return _orientation; }

    bool isOrientationX() const { return _orientation == Orientation::x; }
    bool isOrientationY() const { return _orientation == Orientation::y; }

    void setOrientation(std::string);
    void setOrientation(Orientation);
    void setNumR(int num_r) { _num_r = num_r; }
    void setWidth(double width_r,
                  double width_z = std::numeric_limits<double>::infinity())
    {
      _width_r = width_r;
      _width_z = width_z;
    }
    // Set Z widths of non-uniform meshes
    void setWidths(double width_r, DoubleVec widths_z);
  };

} /* namespace antmoc */

#endif /* LATTICE_H_ */
