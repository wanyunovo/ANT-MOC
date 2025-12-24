/**
 * @file Universe.h
 * @brief The Universe class.
 * @date January 9, 2012
 * @author William Boyd, MIT, Course 22 (wboyd@mit.edu)
 */

#ifndef UNIVERSE_H_
#define UNIVERSE_H_

#ifdef __cplusplus
#include <limits>
#include <map>
#include <set>
#include <vector>

#include "antmoc/constants.h"
#include "antmoc/enum_types.h"
#endif

namespace antmoc
{

  /** Forward declarations */
  class Cell;
  class LocalCoords;
  class Material;
  class Point;

  int universe_id();
  void reset_universe_id();
  void maximize_universe_id(int universe_id);

  /**
   * @class Universe Universe.h "include/Universe.h"
   * @brief A Universe represents an unbounded space in 3D.
   * @details A Universe contains cell which are bounded subspaces in 3D
   *          which together form the Universe. Universes allow
   *          for complex, repeating (i.e. lattices) geometries to be simply
   *          represented with as few data structures as possible.
   */
  class Universe
  {

  protected: // protected 访问权限：派生类（子类）可以访问，外部代码不能直接访问
    // ==================== ID 相关成员变量 ====================

    /** 静态成员变量：记录已创建的 Universe 总数 */
    // static 关键字：该变量被所有 Universe 对象共享，不属于某个具体对象
    // 所有 Universe 对象共用一个 _n，用于统计总数或生成唯一ID
    static int _n;

    /** 唯一ID (_uid)：每个 Universe 对象的唯一标识符 */
    // 单调递增：创建新对象时自动分配，保证全局唯一性
    int _uid;

    /** 为每个创建的 Universe 定义一个用户 ID */
    // 与 _uid 不同，这个ID可以由用户指定，可能重复
    int _id;

    /**  a univese id in the XML file or user-defined */
    // 从外部配置文件（如XML）读取的ID
    int _input_id;

    /** A user-defined name for the Surface */
    char *_name;

    /** The type of Universe (ie, SIMPLE or LATTICE) */
    universeType _type;

    /** A collection of Cell IDs and Cell pointers in this Universe */
    std::map<int, Cell *> _cells;

    /** 可裂变标志：该 Universe 是否包含可发生裂变的材料 */
    bool _fissionable;

    // ==================== 几何边界 ====================

    /** Universe 的空间范围极值：六个面的坐标 */
    double _min_x; // X轴最小坐标（左边界）
    double _max_x; // X轴最大坐标（右边界）
    double _min_y; // Y轴最小坐标（前边界）
    double _max_y; // Y轴最大坐标（后边界）
    double _min_z; // Z轴最小坐标（下边界）
    double _max_z; // Z轴最大坐标（上边界）

    // ==================== 状态标志 ====================

    /** 边界检查标志：边界是否已计算并保持最新 */
    bool _boundaries_inspected;

    /** The boundaryTypes of the universe */
    boundaryType _min_x_bound; // Universe 在最小 x 坐标处的边界条件类型（VACUUM 真空边界 或 REFLECTIVE 反射边界）
    boundaryType _max_x_bound; // Universe 在最大 x 坐标处的边界条件类型
    boundaryType _min_y_bound;
    boundaryType _max_y_bound;
    boundaryType _min_z_bound;
    boundaryType _max_z_bound;

  public:
    static constexpr int NOID = -1;
    bool anonymous() { return _input_id == NOID; }

    Universe(const int id = NOID, const char *name = "");
    virtual ~Universe();
    int getUid() const;
    int getId() const;
    char *getName() const;
    universeType getType();
    int getNumCells() const;
    double getMinX();
    double getMaxX();
    double getMinY();
    double getMaxY();
    double getMinZ();
    double getMaxZ();
    boundaryType getMinXBoundaryType();
    boundaryType getMaxXBoundaryType();
    boundaryType getMinYBoundaryType();
    boundaryType getMaxYBoundaryType();
    boundaryType getMinZBoundaryType();
    boundaryType getMaxZBoundaryType();

    Cell *getCell(int cell_id);
    std::map<int, Cell *> getCells() const;
    std::map<int, Cell *> getAllCells();
    std::map<int, Material *> getAllMaterials();
    std::map<int, Universe *> getAllUniverses();
    bool isFissionable();

    void resetBoundaries();
    void calculateBoundaries();
    void setName(const char *name);
    void setType(universeType type);
    void addCell(Cell *cell);
    void removeCell(Cell *cell);

    bool containsPoint(Point *point);
    Cell *findCell(LocalCoords *coords);
    void setFissionability(bool fissionable);
    void subdivideCells(double max_radius =
                            std::numeric_limits<double>::infinity());
    void buildNeighbors();
    int getInputId();

    virtual std::string toString();
    void printString();

    Universe *clone();
  };

  /**
   * @brief A helper struct for the Universe::findCell() method.
   * @details This is used to insert a Universe's Cells to the back of a vector
   *          of neighbor Cells in Universe::findCell() routine. This works in
   *          symbiosis with the pair_second method template defined below.
   * 把 map 的 value 提取出来
   */
  template <typename tPair>
  struct second_t
  {
    typename tPair::second_type operator()(const tPair &p) const
    {
      return p.second;
    }
  };

  /**
   * @brief A helper routine for the Universe::findCell() method.
   * @details This is used to insert a Universe's Cells to the back of a vector
   *         of neighbor Cells in Universe::findCell() routine. This works in
   *         symbiosis with the second_t struct template defined above.
   * @param map a std::map iterator
   * @return the second element in the iterator (e.g., map value)
   *
   *传入任意 map（或 map 迭代器意义上的容器）就返回一个 second_t 实例，调用处写法更简
   */
  template <typename tMap>
  second_t<typename tMap::value_type> pair_second(const tMap &map)
  {
    return second_t<typename tMap::value_type>();
  }

} /* namespace antmoc */

#endif /* UNIVERSE_H_ */
