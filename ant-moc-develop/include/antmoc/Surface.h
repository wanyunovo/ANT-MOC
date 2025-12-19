/**
 * @文件Surface.h
 * @details Surface 类和子类。
 * @日期 2012 年 1 月 9 日
 * @作者 William Boyd，麻省理工学院，课程 22 (wboyd@mit.edu)
 */

#ifndef SURFACE_H_
#define SURFACE_H_

#include <algorithm>
#include <limits>
#include <map>
#include <string>
#include <vector>

#include "antmoc/constants.h"
#include "antmoc/enum_types.h"
#include "antmoc/math_utils.h"
#include "antmoc/Point.h"

namespace antmoc
{

  /** 转发声明 */
  class Cell;
  class HexLattice;
  class LocalCoords;

  int surface_id();
  void reset_surface_id();
  void maximize_surface_id(int surface_id);

  /**
   * @class Surface Surface.h“src/Surface.h”
   * @brief 表示 3D 中的一般表面。
   * @details Surface 类及其子类用于定义
   *使用构造实体进行 MOC 模拟的几何结构
   *几何（CSG）形式主义。光线追踪过程中使用表面
   *跨越几何体的特征轨迹。
   */
  class Surface
  {

  protected:
    /** 模拟中表面数量的静态计数器 */
    static int _n;

    /** 创建的每个 Surface 的单调递增唯一 ID */
    int _uid;

    /** 为每个创建的 Surface 定义一个用户定义的 id */
    int _id;

    /** 用户定义的 Surface 名称 */
    char *_name;

    /** 曲面的类型（平面、圆柱等） */
    surfaceType _surface_type;

    /** 用于此表面的边界条件类型
     *  （即真空或反射）*/
    boundaryType _boundary_type;

    /* 相邻细胞的向量 */
    // 键只有 +1 和 -1，分别表示 Surface 两侧的半空间；值是指向相应侧所有邻接 Cell 的动态数组。这样每个 Surface 记录“正侧有哪些 Cell、负侧有哪些 Cell”，便于射线追踪时判断穿越后会进入谁。
    // 理想情况下，一个封闭的曲面两侧各只有一个 Cell，但在真实几何里，多个 Cell 可能共用同一截面。例如把一个圆柱面作为多个环形 Cell 的公共边界时，圆柱面正侧可能有多个 Cell，负侧也可能有多个 Cell。每个半空间可以有任意多个 Cell
    // 看初始化函数，可以佐证
    std::map<int, std::vector<Cell *> *> _neighbors;

  public:
    Surface(const int id = 0, const char *name = "");
    virtual ~Surface();

    int getUid() const;
    int getId() const;
    char *getName() const;
    surfaceType getSurfaceType();
    boundaryType getBoundaryType();

    /**
     * @brief 返回轴方向的最小坐标
     *由 halfspace 和该表面定义的空间
     * @param axis 感兴趣的轴（0 = x，1 = y，2 = z）
     * @param halfspace 要考虑的半空间
     * @return轴方向最小坐标
     */
    double getMin(int axis, int halfspace);

    /**
     * @brief 返回轴方向的最大坐标
     *由 halfspace 和该表面定义的空间
     * @param axis 感兴趣的轴（0 = x，1 = y，2 = z）
     * @param halfspace 要考虑的半空间
     * @return轴方向最大坐标
     */
    double getMax(int axis, int halfspace);

    /**
     * @brief 返回此曲面半空间之一的最小 x 值。
     * @param halfspace 要考虑的 Surface 的 halfspace
     * @return最小x值
     */
    virtual double getMinX(int halfspace) = 0;

    /**
     * @brief 返回此曲面半空间之一的最大 x 值。
     * @param halfspace 要考虑的 Surface 的 halfspace
     * @返回最大x值
     */
    virtual double getMaxX(int halfspace) = 0;

    /**
     * @brief 返回此表面半空间之一的最小 y 值。
     * @param halfspace 要考虑的 Surface 的 halfspace
     * @return 最小y值
     */
    virtual double getMinY(int halfspace) = 0;

    /**
     * @brief 返回此表面半空间之一的最大 y 值。
     * @param halfspace 要考虑的 Surface 的 halfspace
     * @return 最大y值
     */
    virtual double getMaxY(int halfspace) = 0;

    /**
     * @brief 返回此曲面半空间之一的最小 z 值。
     * @param halfspace 要考虑的 Surface 的 halfspace
     * @return 最小z值
     */
    virtual double getMinZ(int halfspace) = 0;

    /**
     * @brief 返回此曲面半空间之一的最大 z 值。
     * @param halfspace 要考虑的 Surface 的 halfspace
     * @返回最大z值
     */
    virtual double getMaxZ(int halfspace) = 0;

    void setName(const char *name);
    void setBoundaryType(const boundaryType boundary_type);
    void addNeighborCell(int halfspace, Cell *cell);

    /**
     * @brief 使用表面的势方程评估点。
     * @details 此方法返回势值 f(x,y) f$
     *代表此 Surface 的函数 \f$f\f$。
     * @param 指向感兴趣的 Soint 的指针
     * @return 平面势方程中点的值。
     */
    virtual double evaluate(const Point *point) const = 0;

    /**
     * @brief 从给定的曲面中查找与该曲面的交点
     *由角度定义的点和轨迹。
     * @param point 指向兴趣点的指针
     * @param azim 方位角（以弧度为单位）
     * @param Polar 极角（以弧度为单位）
     * @parampoints 点数组来存储交叉点位置
     * @return 交点的数量（0或1）
     */
    virtual int intersection(Point *point, double azim, double polar,
                             Point *points) = 0;

    bool isPointOnSurface(Point *point);
    bool isCoordOnSurface(LocalCoords *coord);
    double getMinDistance(Point *point, double azim, double polar = M_PI_2);

    /**
     * @brief 将此 Surface 的属性转换为字符数组。
     * @details 返回的字符数组包含 Surface 的类型（即
     *PLANE) 和势方程中的系数。
     * @return 此 Surface 属性的字符数组
     */
    virtual std::string toString() = 0;

    void printString();
  };

  /**
   * @class 平面 Surface.h "src/Surface.h"
   * @brief 表示垂直于 xy 平面的平面。
   */
  class Plane : public Surface
  {

  protected:
    /** x 中线性项的系数 */
    double _A;

    /** y 中线性项的系数 */
    double _B;

    /** z 中的线性项的系数 */
    double _C;

    /** 常数偏移 */
    double _D;

    /** Plane 是 Surface 类的朋友 */
    friend class Surface;

    /** 飞机是 Z气缸类的朋友 */
    friend class ZCylinder;

  public:
    Plane(const double A, const double B, const double C, const double D,
          const int id = 0, const char *name = "");

    double getMinX(int halfspace);
    double getMaxX(int halfspace);
    double getMinY(int halfspace);
    double getMaxY(int halfspace);
    double getMinZ(int halfspace);
    double getMaxZ(int halfspace);
    double getA();
    double getB();
    double getC();
    double getD();

    double evaluate(const Point *point) const;
    int intersection(Point *point, double azim, double polar, Point *points);

    std::string toString();
  };

  /**
   * @class XPlane Surface.h“src/Surface.h”
   * @brief 表示垂直于 x 轴的平面。
   */
  class XPlane : public Plane
  {

  private:
    /** XPlane 沿 x 轴的位置 */
    double _x;

  public:
    XPlane(const double x, const int id = 0, const char *name = "");

    void setX(const double x);

    double getX();
    double getMinX(int halfspace);
    double getMaxX(int halfspace);

    std::string toString();
  };

  /**
   * @class YPlane Surface.h“src/Surface.h”
   * @brief 表示垂直于 y 轴的平面。
   */
  class YPlane : public Plane
  {

  private:
    /** YPlane 沿 y 轴的位置 */
    double _y;

  public:
    YPlane(const double y, const int id = 0, const char *name = "");

    void setY(const double y);

    double getY();
    double getMinY(int halfspace);
    double getMaxY(int halfspace);

    std::string toString();
  };

  /**
   * @class ZPlane Surface.h“src/Surface.h”
   * @brief 表示垂直于 z 轴的平面。
   */
  class ZPlane : public Plane
  {

  private:
    /** ZPlane 沿 z 轴的位置 */
    double _z;

  public:
    ZPlane(const double z, const int id = 0, const char *name = "");

    void setZ(const double z);

    double getZ();
    double getMinZ(int halfspace);
    double getMaxZ(int halfspace);

    std::string toString();
  };

  /**
   * @class ZCylinder Surface.h“src/Surface.h”
   * @brief 表示轴平行于 z 轴的圆柱体。
   */
  class ZCylinder : public Surface
  {

  private:
    /** ZCylinder 中心的一个点 */
    Point _center;

    /** ZCylinder 的半径 */
    double _radius;

    /** x 平方项的系数 */
    double _A;

    /** y 平方项的系数 */
    double _B;

    /** x 中的线性项的系数 */
    double _C;

    /** y 中的线性项的系数 */
    double _D;

    /** 常数偏移 */
    double _E;

    /** ZCylinder 是 Surface 类的朋友 */
    friend class Surface;

    /** ZCylinder 是 Plane 类的朋友 */
    friend class Plane;

  public:
    ZCylinder(const double x, const double y, const double radius,
              const int id = 0, const char *name = "");

    /// \brief 复制构造函数
    ZCylinder(const ZCylinder &zcylinder);

    double getX0();
    double getY0();
    double getRadius();
    double getMinX(int halfspace);
    double getMaxX(int halfspace);
    double getMinY(int halfspace);
    double getMaxY(int halfspace);
    double getMinZ(int halfspace);
    double getMaxZ(int halfspace);

    double evaluate(const Point *point) const;
    int intersection(Point *point, double azim, double polar, Point *points);

    std::string toString();
  };

  /// \class HexLatticePrism Surface.h "src/Surface.h"
  /// \brief 表示轴平行于 z 轴的六角晶格的表面。
  class HexLatticePrism : public Surface
  {

  private:
    ///< HexLatticePrism 中心的一个点 */
    Point _center;

    ///< 底层 HexLattice
    HexLattice *_lattice;

    ///< 底层晶格的包围柱面
    ZCylinder *_bounding_cylinder;

    ///< 包含包围圆柱体的包围格子
    HexLattice *_bounding_lattice;

    ///< HexLatticePrism 是 Surface 类的友元
    friend class Surface;

  public:
    HexLatticePrism(const double x, const double y,
                    const double width_r, const int num_r,
                    const char *orientation, const int id = 0,
                    const char *name = "");
    // 复制构造函数
    HexLatticePrism(const HexLatticePrism &);

    virtual ~HexLatticePrism();

    double getX0();
    double getY0();
    int getNumR();
    double getWidthR();
    double getHeight();
    double getMinX(int halfspace);
    double getMaxX(int halfspace);
    double getMinY(int halfspace);
    double getMaxY(int halfspace);
    double getMinZ(int halfspace);
    double getMaxZ(int halfspace);

    double evaluate(const Point *point) const;
    int intersection(Point *point, double azim, double polar, Point *points);

    std::string toString();
  };

} /* 命名空间 Antmoc */

#endif /* 表面高度 */
