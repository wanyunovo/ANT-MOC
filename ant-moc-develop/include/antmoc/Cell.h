/// \文件包含/Cell.h
/// \简介 Cell 类。
/// \日期 2012 年 1 月 18 日
/// 作者 William Boyd，麻省理工学院，课程 22 (wboyd@mit.edu)

#ifndef CELL_H_
#define CELL_H_

#ifdef __cplusplus
#include <limits>
#include <map>
#include <unordered_map>
#include <set>
#include <string>
#include <vector>

#include "antmoc/BoolParser.h"
#include "antmoc/enum_types.h"
#include "antmoc/math_utils.h"
#endif

namespace antmoc
{

  // 前向声明解决循环依赖
  class Material;
  class RecLattice;
  class LocalCoords;
  class Point;
  class Surface;
  class Universe;
  class CellIter;

  int cell_id();
  void reset_cell_id();
  void maximize_cell_id(int cell_id);

  using SurfaceMap = std::unordered_map<int, Surface *>;

  struct Halfspace
  {
    int halfspace; // +1 or -1
    Surface *surface;
  };

  /// \class Cell Cell.h "src/Cell.h"
  /// \brief 代表宇宙内部的一个细胞。
  class Cell
  {

  private:
    ///< 细胞数量的静态计数器
    static int _n;

    ///< 布尔解析器
    static BoolParser parser;

    ///< 创建的每个 Cell 的单调递增唯一 ID
    int _uid;

    ///< 为每个创建的 Cell 定义的 ID
    int _id;

    ///< 用户定义的 ID
    int _input_id;

    ///< 用户定义的 Surface 名称
    char *_name;

    ///< 单元类型（即 MATERIAL 或 FILL）
    cellType _cell_type;

    ///< 指向填充此单元格的材料或宇宙的指针
    /*
    void* 在 ant-moc-develop/include/antmoc/Cell.h (line 74) 中只是一个“无类型”指针，占位保存某个实体的地址，但编译器不知道它指向什么类型，不能直接解引用或调用成员；只能存/取再手动转型。
Cell 里 _fill 既可能指向 Material 也可能指向 Universe/Lattice，用 void* 就能灵活保存这几种情况，不过会失去类型安全，转换错了就会产生未定义行为。
    */
    void *_fill;

    ///< 根据重叠线段计算出的单元的体积/面积
    double _volume;

    ///< 几何中此单元的实例总数
    int _num_instances;

    ///< 指示单元格是否旋转的布尔值
    bool _rotated;

    ///< 具有以度为单位的角度的数组，用于围绕 x、y 和 z 旋转
    double _rotation[3];

    ///< 根据旋转角度定义的旋转矩阵
    double _rotation_matrix[9];

    ///< 指示是否翻译单元格的布尔值
    bool _translated;

    ///< 具有 x、y 和 z 方向平移的数组
    double _translation[3];

    ///< 细分此 Cell 的环数
    int _num_rings;

    ///< 细分该 Cell 的扇区数
    int _num_sectors;

    ///< 如果被另一个 Cell 克隆则为父 Cell
    Cell *_parent;

    ///< 指示这是一组半空间的交集的布尔值
    // Determine whether there is operator '|' or '~'
    bool _simple;

    ///< 区域的抽象语法树
    NodePtr _ast;

    ///< 可以从 _ast 生成的区域的逆波兰表示法
    PyVector<int> _rpn;

    ///< 带有指针的边界 Surface ID 的映射
    SurfaceMap _surfaces;

    ///< 相邻细胞的向量
    std::vector<Cell *> _neighbors;
    // 同一个 Cell 可以有多个邻居，来自不同的 Surface。每当 Surface 两侧有新组合，就通过Surface::addNeighborCell的双层循环把这些邻居关系追加进去
    /**
      为什么负侧每个 Cell可以把正侧所有 Cell 记为邻居，他们为什么相邻?
      因为“邻居”是从射线追踪/几何连通的角度来定义的：如果两个 Cell 共享同一个 Surface，射线穿过这个 Surface 时就会从一个 Cell 进入另一个 Cell，这就是邻接关系。Surface 的 +1 半空间包含“位于这面右侧/外侧/法向一端”的所有 Cell，-1 半空间则是左侧（反法向）所有 Cell。它们虽然各自可能不止一个，但只要这些 Cell 的边界里包含同一个 Surface，就说明“从这面穿过去能到对面的任意一个 Cell”。因此 addNeighborCell 把负侧的每个 Cell 与正侧的每个 Cell 互相登记为邻居，确保无论射线当前处在哪个 Cell，只要穿过该 Surface，就能直接查到对面的那些候选 Cell，继续追踪。
     */

    // The material cell has temperature
    float _temperature;
    // Cell表示的物理含义（棒、组件、堆芯等）
    cellPhy _cell_phy;

    void ringify(std::vector<Cell *> &subcells, double max_radius);
    void sectorize(std::vector<Cell *> &subcells);

  public:
    static constexpr int NOID = -1;
    bool anonymous() { return _input_id == NOID; }

    //-----------------------------------------------------
    // 迭代器
    //-----------------------------------------------------
    friend class CellIter;
    CellIter begin();
    CellIter end();

    Cell(int id = NOID, const char *name = "");
    virtual ~Cell();
    int getUid() const;
    int getId() const;
    int getInputId();
    char *getName() const;
    cellType getType() const;
    Material *getFillMaterial();
    Universe *getFillUniverse();
    double getVolume();
    int getNumInstances();
    bool isRotated();
    bool isTranslated();
    double getPhi(std::string units = "degrees");
    double getTheta(std::string units = "degrees");
    double getPsi(std::string units = "degrees");
    double *getRotationMatrix();
    double *getTranslation();
    void retrieveRotation(double *rotations, int num_axes,
                          std::string units = "degrees");
    void retrieveTranslation(double *translations, int num_axes);
    int getNumRings();
    int getNumSectors();
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
    int getNumSurfaces() const;
    SurfaceMap getSurfaces() const;
    std::vector<Cell *> getNeighbors() const;
    bool hasParent();
    Cell *getParent();
    Cell *getOldestAncestor();

    int getNumZCylinders();

    std::map<int, Cell *> getAllCells();
    std::map<int, Universe *> getAllUniverses();

    void setName(const char *name);
    void setFill(Material *fill);
    void setFill(Universe *fill);
    void setVolume(double volume);
    void incrementVolume(double volume);
    void setNumInstances(int num_instances);
    void incrementNumInstances();
    void setRotation(double *rotation, int num_axes, std::string units = "degrees");
    void setTranslation(double *translation, int num_axes);
    void setNumRings(int num_rings);
    void setNumSectors(int num_sectors);
    void setParent(Cell *parent);
    void addSurface(int halfspace, Surface *surface);
    void removeSurface(Surface *surface);
    void addNeighborCell(Cell *cell);

    std::string getRegion() const { return parser.toString(_rpn); }
    PyVector<int> getRPN() const { return _rpn; }
    void setSurfaces(const SurfaceMap surfaces) { _surfaces = surfaces; }
    void addRegionSurfaces(const std::string region, const SurfaceMap &surfaces);
    void generateRPN();

    /// \brief 指示标记是否代表半空间
    bool isValidHalfspace(int token)
    {
      return !parser.isOperator(token) && _surfaces.count(abs(token));
    }

    bool isFissionable();
    bool containsPoint(Point *point);
    bool containsSimple(Point *point);
    bool containsComplex(Point *point);
    bool containsCoords(LocalCoords *coords);
    double minSurfaceDist(Point *point, double azim, double polar);
    double minSurfaceDist(LocalCoords *coords);

    Cell *clone(bool clone_region = true);
    void subdivideCell(double max_radius);
    void buildNeighbors();

    // 额外参数和补充
    void setTemperature(float);
    float getTemperature();
    void setPhy(const char *);
    cellPhy getPhy();

    std::string toString();
    void printString();
  };

  /// \class 单元半空间上的迭代器
  /// \details 此迭代器仅迭代有效的标记。运营商都没有
  ///          删除的曲面也不会被迭代。
  ///          请注意，如果单元格的 RPN
  ///          在迭代内部发生改变。
  /**
   * 它遍历 _rpn（逆波兰式），只挑出那些合法的半空间 token，并返回一个 Halfspace 结构 {halfspace, surface}，其中 surface 通过 _surfaces[abs(token)] 取出。
   * 之所以不用直接遍历 SurfaceMap，是因为 _rpn 表示几何布尔表达式的顺序；Cell 可能有被删掉或组合成更复杂表达式的面，迭代器负责跳过无效项，只把真正定义边界的半空间提供给调用者。
   * 通过迭代器遍历cell，每次拿到的就是一个Halfspace
   */
  class CellIter
  {
  public:
    CellIter(Cell &cell, size_t pos = 0)
        : _cell(cell), _pos(pos)
    {
      // 从第一个有效的半空间开始
      while (_pos < _cell._rpn.size())
      {
        int token = _cell._rpn[_pos];
        if (_cell.isValidHalfspace(token))
          break;
        ++_pos;
      }
    }

    /// \brief 如果两个 CellIter 位于相同的有效位置，则它们相等
    bool operator==(const CellIter &rhs) { return (_pos == rhs._pos); }

    bool operator!=(const CellIter &rhs) { return !(*this == rhs); }

    /// \brief 取消引用迭代器
    /// \返回一个半空间
    Halfspace operator*()
    {
      int token = _cell._rpn[_pos];
      auto surface = _cell._surfaces[abs(token)];
      return {copysign(1, token), surface};
    }

    /// \brief 移至下一个有效标记
    CellIter &operator++()
    {
      // 查找下一个有效索引
      while (_pos < _cell._rpn.size() - 1)
      {
        ++_pos;
        int token = _cell._rpn[_pos];
        // 即在自增之后马上用 _pos 当索引去访问 _rpn。为了不越界，必须保证“自增之后仍落在有效下标内”，因此循环条件写成 < size() - 1：
        if (_cell.isValidHalfspace(token))
          return *this;
      }
      // 最后一个元素之后的索引
      _pos = _cell._rpn.size(); // 这里的 *this 只是返回对当前迭代器对象的引用，方便链式调用，并不会访问 _rpn[_pos]。真正的“解引用”是在 CellIter::operator*() 中才进行的。
      return *this;
    }

    /// \brief 返回迭代器的当前位置
    size_t getPos() const { return _pos; }

  protected:
    Cell &_cell; ///< 要迭代的单元格
    size_t _pos; ///< 当前半空间的索引
  };

} // 命名空间 Antmoc

#endif // 细胞 h
