#include "antmoc/Lattice.h"
#include "antmoc/log.h"
#include "antmoc/Cell.h"
#include "antmoc/LocalCoords.h"
#include "antmoc/string_utils.h"
#include "antmoc/Surface.h"

#define _USE_MATH_DEFINES
#include <algorithm>
#include <cmath>

namespace antmoc
{

  // Check two distances by coincidence tolerance
  inline bool coincident(double d1, double d2)
  {
    return std::abs(d1 - d2) < FLT_COINCIDENT;
  }

  //////////////////////////////////////////////////
  //                                              //
  //              Lattice: abstraction            //
  //                                              //
  //////////////////////////////////////////////////

  Lattice::Lattice(const int id, const char *name) : Universe(id, name)
  {

    _num_z = 0;
    _width_z = 0;
    _non_uniform = false;

    _type = LATTICE; // universeType
    _offset.setCoords(0.0, 0.0, 0.0);
  }

  /**
   * @brief Returns an iterator to the universe past the last one.
   */
  LatticeIter Lattice::end()
  {
    return LatticeIter(*this, _universes.size());
  }

  /// \brief Gets the maximum number of universes, including invalid ones
  int Lattice::getMaxNumUniverses() const
  {
    return getNumX() * getNumY() * getNumZ();
  }

  /// \brief Computes lat_x by index, including invalid cells
  /// \details This method is useful when we want to convert an index
  ///          back into the lattice cell position.
  int Lattice::getLatXByIndex(int indx)
  {
    return indx % getNumX();
  }

  /// \brief Computes lat_y by index, including invalid cells
  /// \details This method is useful when we want to convert an index
  ///          back into the lattice cell position.
  int Lattice::getLatYByIndex(int indx)
  {
    return (indx / getNumX()) % getNumY();
  }

  /// \brief Computes lat_z by index, including invalid cells
  /// \details This method is useful when we want to convert an index
  ///          back into the lattice cell position.
  int Lattice::getLatZByIndex(int indx)
  {
    return indx / (getNumX() * getNumY());
  }

  /// \brief Computes the position of a lattice cell, including invalid cells
  std::array<int, 3> Lattice::getLatXYZByIndex(int indx)
  {
    return {getLatXByIndex(indx),
            getLatYByIndex(indx),
            getLatZByIndex(indx)};
  }

  /// \brief Returns the center point of a lattice cell with respect to the center
  ///        of the center of parent universe
  Point Lattice::getLatticeCellCenter(int index)
  {

    // Converts the index into position
    auto lat_xyz = getLatXYZByIndex(index);
    return getLatticeCellCenter(lat_xyz[0], lat_xyz[1], lat_xyz[2]);
  }

  /**
   * @brief Clears memory for all of Universes pointers.
   */
  void Lattice::clearUniverses()
  {
    _universes.clear();
  }

  /**
   * @brief Set the offset in global coordinates for this Lattice.
   * @details A lattice is assumed to be a grid with the center/origin
   *          of the grid located in the center of the Lattice's parent universe.
   *          The offset represents the offset of the lattice center/origin with
   *          respect to the center of the parent universe. Therefore an offset of
   *          (-1,2,1) would move the center/origin of the lattice to the left
   *          1 cm, forward 2 cm, and up 1 cm.
   * @param x the offset in the x direction
   * @param y the offset in the y direction
   * @param z the offset in the z direction
   * 格子默认以父 universe 中心为原点，_offset 表示把格子的原点平移到哪里
   */
  void Lattice::setOffset(double x, double y, double z)
  {
    _offset.setX(x);
    _offset.setY(y);
    _offset.setZ(z);
  }

  /// \brief Set the offset in global coordinates for this Lattice
  void Lattice::setOffset(Point offset)
  {
    _offset = offset;
  }

  /**
   * @brief Returns the minimum reachable z-coordinate in the Lattice.
   * @return the minimum reachable z-coordinate
   * 所有把全局坐标转换成格子内部坐标、再转换回去的地方都要加减 _offset，所以必须存这个点。
   * 例如求几何范围就用 _offset ± 累计宽度一半 得到实际的 min/ma
   */
  double Lattice::getMinZ() const
  {
    return _offset.getZ() - _accumulate_z[_num_z] / 2.;
  }

  /**
   * @brief Returns the maximum reachable z-coordinate in the Lattice.
   * @return the maximum reachable z-coordinate
   */
  double Lattice::getMaxZ() const
  {
    return _offset.getZ() + _accumulate_z[_num_z] / 2.;
  }

  /**
   * @brief Get the z-coord of a point in the local system of a universe.
   * @param point a pointer to a point being evaluated.
   * @return the local z-coord
   */
  double Lattice::getLocalPointZ(const Point *p, int lat_z)
  {
    return p->getZ() - (getMinZ() + _widths_z[lat_z] / 2. + _accumulate_z[lat_z]);
  }

  /**
   * @brief Finds the Lattice cell z index that a point lies in.
   * @param point a pointer to a point being evaluated.
   * @return the Lattice cell z index.
   */
  int Lattice::getLatZ(Point *point)
  {

    /* Check to see if lattice is infinite in z direction */
    if (!is3D())
      return 0;

    int lat_z = -1;

    /* get the distance to the bottom surface */
    double dist_to_bottom = point->getZ() - getMinZ();

    /* Compute the z index for the Lattice cell this point is in */
    for (int i = 0; i < _num_z; i++)
    {
      if (dist_to_bottom >= _accumulate_z[i] &&
          dist_to_bottom < _accumulate_z[i + 1])
      {
        lat_z = i;
        break;
      }
    }

    /* Check if the Point is on the Lattice boundaries and if so adjust
     * z Lattice cell index */
    if (fabs(dist_to_bottom) < ON_SURFACE_THRESH)
      lat_z = 0;
    else if (fabs(dist_to_bottom - _accumulate_z[_num_z]) < ON_SURFACE_THRESH)
      lat_z = _num_z - 1;
    if (lat_z == -1)
      log::ferror("Trying to get lattice z index for point(z = %f) that is "
                  "outside lattice bounds. dist_to_bottom = %f is not within "
                  "[0.0, %f]",
                  point->getZ(), dist_to_bottom, _accumulate_z[_num_z]);

    return lat_z;
  }

  /**
   * @brief Finds the Cell within this Lattice that a LocalCoords is in.
   * @details This method first find the Lattice cell, then searches the
   *          Universe inside that Lattice cell. If LocalCoords is outside
   *          the bounds of the Lattice, this method will return NULL.
   * @param coords the LocalCoords of interest. Coordinates of coords and
   *        lattice._offset share the same origin.
   * @return a pointer to the Cell this LocalCoord is in or NULL
   * 链表视角：当前 coords：“我在 Lattice L 里，在第 (lat_x,lat_y,lat_z) 个小格子里，坐标是 (x,y,z)。”
   * 新的 next_coords：“在这个格子里填的 Universe U 的坐标系下，我的位置是 (next_x,next_y,next_z)。”
   * 又往下长了一节
   */
  Cell *Lattice::findCell(LocalCoords *coords)
  {

    if (coords->getX() != coords->getX() ||
        coords->getY() != coords->getY())
    {
      log::ferror("Nan is found in points when trying to find a lattice cell "
                  "(Lattice::findCell): x = %f, y = %f",
                  coords->getX(), coords->getY());
    }

    /* Set the LocalCoord to be a LAT type at this level */
    coords->setType(LAT); // 这一层是 Lattice 层

    /* Compute the x and y indices for the Lattice cell this coord is in */
    // 定位该点属于 lattice 网格中的哪一个单元
    auto lat_xy = getLatXY(coords->getPoint());
    int lat_x = lat_xy[0];
    int lat_y = lat_xy[1];
    int lat_z = getLatZ(coords->getPoint());

    /* Determin if we are referencing a valid lattice cell */
    if (!areValidIndices(lat_x, lat_y, lat_z))
      return nullptr;

    /* Compute local position of Point in the next level Universe. The offset of
       Lattice should be considered. */
    // 以该格子的局部坐标系为基准，算下一层 Universe 里的坐标
    double next_x = getLocalPointX(coords->getPoint(), lat_x, lat_y);
    double next_y = getLocalPointY(coords->getPoint(), lat_x, lat_y);
    double next_z = getLocalPointZ(coords->getPoint(), lat_z);

    /* Check for 2D problem or 2D lattice */
    if (!is3D())
      next_z = coords->getZ();

    /* Create a new LocalCoords object for the next level Universe */
    // 生成下一层 LocalCoords 结点
    LocalCoords *next_coords = coords->getNextCreate(next_x, next_y, next_z);
    Universe *univ = getUniverse(lat_x, lat_y, lat_z);
    next_coords->setUniverse(univ);

    /* Set Lattice indices */
    // 当前这一层记录自己是哪个 Lattice、哪个格子
    coords->setLattice(this);
    coords->setLatticeX(lat_x);
    coords->setLatticeY(lat_y);
    coords->setLatticeZ(lat_z);

    /* Search the next lowest level Universe for the Cell */
    // 最后递归调用子 Universe 的 findCell(next_coords)，继续向更低层搜索直到Cell被返回
    return univ->findCell(next_coords);
  }

  /**
   * @brief Return a 3D vector of the Universes in the Lattice.
   * @return 3D vector of Universes
   *
   * 获取填充 Lattice 的所有 Universe(三维 vector 容器)
   */
  std::vector<Universe *> *Lattice::getUniverses()
  {
    return &_universes;
  }

  /**
   * @brief Aggregates a list (vector) of the IDs of all Universes within
   *        the FILL type Cells filling this Universe.
   * @details Note that this method only searches the first level of Cells
   *          below this Universe within the nested Universe coordinate system.
   * @return a map of Universe keyed by the universe ID.
   *
   * 整理 Lattice 中包含的所有不重复 Universe 的标号以及指针,  并返回存储 Universe 标号与指针的 map 容器
   */
  std::map<int, Universe *> Lattice::getUniqueUniverses()
  {

    std::map<int, Universe *> unique_universes;

    // Iterate all universes
    for (auto u : *this)
      /* The lattice may reuse the same Universe pointer at different indices.
       * Using the Universe ID as the map key collapses duplicates. */
      if (u)
        unique_universes[u->getId()] = u;
    // 因为 map 的 key 唯一，重复的 ID 自动覆盖，也就是“去重”。

    return unique_universes;
  }

  /**
   * @brief Returns the std::map of Cell IDs and Cell pointers in this Lattice
   *        at all nested Universe levels.
   * @return std::map of Cell IDs and pointers
   */
  std::map<int, Cell *> Lattice::getAllCells()
  {

    std::map<int, Cell *> cells;
    std::map<int, Universe *> unique_universes = getUniqueUniverses();

    std::map<int, Cell *> nested_cells;

    for (auto &u : unique_universes)
    {
      /* Each Universe recursively reports its own Cells (possibly nested through
       * further fills) so we merge their maps into a single container. */
      nested_cells = u.second->getAllCells();
      cells.insert(nested_cells.begin(), nested_cells.end());
    }

    return cells;
  }

  /**
   * @brief Returns the std::map of all nested Universe IDs and Universe pointers
            filling this Lattice.
   * @return std::map of Universe IDs and pointers
   */
  std::map<int, Universe *> Lattice::getAllUniverses()
  {

    /* Initialize a map of all Universes contained by the Lattice in each
     * nested Universe level */
    std::map<int, Universe *> all_universes;

    /* Get all unique Universes contained in each of the Lattice cells */
    std::map<int, Universe *> unique_universes = getUniqueUniverses();

    /* Add the unique Universes filling each Lattice cell */
    all_universes.insert(unique_universes.begin(), unique_universes.end());

    /* Append all Universes containing each Cell to the map */
    std::map<int, Universe *> nested_universes;

    for (auto &u : unique_universes)
    {
      /* Drill down into each Universe to pick up any lower-level fills, then
       * merge them into the global set keyed by Universe ID. */
      nested_universes = u.second->getAllUniverses();
      all_universes.insert(nested_universes.begin(), nested_universes.end());
    }

    return all_universes;
  }

  /**
   * @brief Returns a pointer to the Universe within a specific Lattice cell.
   * @param lat_x the x index to the Lattice cell
   * @param lat_y the y index to the Lattice cell
   * @param lat_z the z index to the Lattice cell
   * @return pointer to a Universe filling the Lattice cell
   * 获取填充 Lattice 指定位置的 Universe
   */
  Universe *Lattice::getUniverse(int lat_x, int lat_y, int lat_z)
  {

    /* Checks that lattice indices are within the bounds of the lattice */
    if (!areValidIndices(lat_x, lat_y, lat_z))
      log::ferror("Cannot retrieve Universe from Lattice ID = %d: Index"
                  "out of bounds: Tried to access Cell x = %d, y = %d, z = %d but bounds"
                  "are %s",
                  _id, lat_x, lat_y, lat_z, dimToString().c_str());

    /* getLatticeCell flattens the (x,y,z) indices into the backing vector. */
    return _universes[getLatticeCell(lat_x, lat_y, lat_z)];
  }

  /**
   * @brief Update the Universe in a particular Lattice cell.
   * @details This method may only be used after an array of Universes
   *          has been assigned with the Lattice::setUniverses(...) method.
   * @param lat_x the Lattice cell index along x
   * @param lat_y the Lattice cell index along y
   * @param lat_z the Lattice cell index along z
   * @param universe the Universe to insert into the Lattice
   */
  void Lattice::updateUniverse(int lat_x, int lat_y, int lat_z,
                               Universe *universe)
  {

    if (_universes.size() == 0)
      log::ferror("Unable to set Universe %d in Lattice %d which "
                  "has not yet been assigned an array of Universes",
                  universe->getId(), _id);

    if (!areValidIndices(lat_x, lat_y, lat_z))
    {
      log::ferror("Unable to set Universe %d in Lattice %d with "
                  "Lattice cell index =(%d, %d, %d) which is outside the "
                  "array of Universes",
                  universe->getId(), _id,
                  lat_x, lat_y, lat_z);
    }

    /* Assign the Universe to the array */
    /* Because we store Universes in a linear vector, this assignment is the
     * single source of truth for what fills the specified lattice location. */
    _universes[getLatticeCell(lat_x, lat_y, lat_z)] = universe;
  }

  /**
   * @brief Removes all references to a Universe from the Lattice.
   * @param universe the Universe to remove
   */
  void Lattice::removeUniverse(Universe *universe)
  {

    /* Clear any Universes in the Lattice (from a previous run) */
    for (auto &u : *this)
      if (universe->getId() == u->getId())
        u = nullptr;
  }

  /**
   * @brief Subdivides all of the Material-filled Cells within this Lattice
   *        into rings and angular sectors aligned with the z-axis.
   * @param max_radius the maximum allowable radius used in the subdivisions
   *
   * 将 Lattice 中所有以材料填充的 Cell,按照垂直于 z 轴方向,进行环状和扇区分割
   */
  void Lattice::subdivideCells(double max_radius)
  {

    log::fdebug("Subdividing Cells for Lattice ID=%d "
                "with max radius %f",
                _id, max_radius);

    std::map<int, Universe *> universes = getUniqueUniverses();

    /* unique_radius is used as the maximum radius for the ringified Cells */
    std::map<int, double> unique_radius = getUniqueRadius(universes);

    /* Subdivide all Cells */
    for (auto &u : universes)
    {
      log::fdebug("univ_ID: %d, radius: %f, max_radius: %f",
                  u.first, unique_radius[u.first], max_radius);

      /* If the  local universe equivalent radius is smaller than max_radius
         parameter, over-ride it*/
      /* Each Universe decides how to subdivide its Cells, but we clamp the
       * target radius by both the lattice-specific radius and the user cap. */
      u.second->subdivideCells(std::min(unique_radius[u.first],
                                        max_radius));
    }
  }

  /**
   * @brief Builds collections of neighboring Cells for all Cells in each
   *        Universe in the Lattice for optimized ray tracing.
   */
  void Lattice::buildNeighbors()
  {

    /* Get list of unique Universes in this Lattice */
    std::map<int, Universe *> universes = getUniqueUniverses();

    /* Loop over each Universe and make recursive call */
    for (auto &u : universes)
      /* Delegate to each Universe so that Cell::buildNeighbors can walk down the
       * fill hierarchy and register adjacency information. */
      u.second->buildNeighbors();
  }

  /**
   * @brief Prints a string representation of all of the Lattice's attributes to
   *        the console.
   */
  void Lattice::printString()
  {
    log::fresult(toString().c_str());
  }

  //////////////////////////////////////////////////
  //                                              //
  //                  RecLattice                  //
  //                                              //
  //////////////////////////////////////////////////

  /**
   * @brief Constructor sets the user-specified and unique IDs for this Lattice.
   * @param id the user-specified optional Lattice (Universe) ID
   * @param name the user-specified optional Lattice (Universe) name
   */
  RecLattice::RecLattice(const int id, const char *name) : Lattice(id, name)
  {

    _lattice_type = latticeType::Rectangle;

    /* Default width and number of Lattice cells along each dimension */
    _num_x = 0;
    _num_y = 0;
    _width_x = 0;
    _width_y = 0;
  }

  /// \brief Get the begin iterator
  LatticeIter RecLattice::begin()
  {
    return LatticeIter(*this, 0);
  }

  /**
   * @brief Returns the minimum reachable x-coordinate in the Lattice.
   * @return the minimum reachable x-coordinate
   * _offset 存的是整个 lattice 中心在父级坐标系的位置；_accumulate_x[_num_x] 是 x 方向所有单元宽度的总和，即 lattice 的全宽
   * getMaxX()/getMinX() 的任务就是返回 lattice 在全局坐标下的左右边界。由于 _offset 位于 lattice 中心，把总宽度除以 2 再加减 _offset，就得到左右端点：
   * 最左：_offset.getX() - total_width/2
   * 最右：_offset.getX() + total_width/2。
   * y/z 方向同理，用 _accumulate_y[_num_y]/_accumulate_z[_num_z]。这样保证无论 lattice 被放在父几何哪个位置，只要知道中心和总宽，就能快速得到边界，后续判断点是否在 lattice 内、转换局部坐标等都依赖这些值。
   */
  double RecLattice::getMinX() const
  {
    return _offset.getX() - _accumulate_x[_num_x] / 2.; // offset.setX(min_x + (max_x - min_x)/2.0);相当于求min_x
  }

  /**
   * @brief Returns the maximum reachable x-coordinate in the Lattice.
   * @return the maximum reachable x-coordinate
   */
  double RecLattice::getMaxX() const
  {
    return _offset.getX() + _accumulate_x[_num_x] / 2.;
  }

  /**
   * @brief Returns the minimum reachable y-coordinate in the Lattice.
   * @return the minimum reachable y-coordinate
   */
  double RecLattice::getMinY() const
  {
    return _offset.getY() - _accumulate_y[_num_y] / 2.;
  }

  /**
   * @brief Returns the maximum reachable y-coordinate in the Lattice.
   * @return the maximum reachable y-coordinate
   */
  double RecLattice::getMaxY() const
  {
    return _offset.getY() + _accumulate_y[_num_y] / 2.;
  }

  /**
   * @brief Get the maximum equivalent radius of each unique universes. Equivalent
   *         radius are computed as the diagonal length to the cell boundary.
   * @param unique_universes The unique universes of this Lattice
   * @return a map of unique radius keyed by the universe ID.
   *  计算矩形 Lattice 内每种唯一 Universe 的“等效半径”，并返回一个 map<int,double>，键是 Universe ID，值是半径。
“等效半径”是指填充该 Universe 的矩形单元对角线的一半：sqrt( (width_x/2)^2 + (width_y/2)^2 )。它衡量了“从单元中心到边界最远点”的距离，后续在 subdivideCells 里会用它来控制环/扇形细分的最大半径。
   */
  std::map<int, double> RecLattice::getUniqueRadius(std::map<int, Universe *> &unique_universes)
  {

    std::map<int, double> unique_radius;

    /* Create and initialize the <universe ID, unique radius> map */
    for (auto &u : unique_universes)
      unique_radius[u.first] = 0.;

    /* Get the maximum equivalent radius of each unique universe */
    for (auto it = begin(); it != end(); ++it)
    {
      auto i = it.getPos() % _num_x;
      auto j = it.getPos() % (_num_x * _num_y) / _num_x;
      // 这里不懂
      /*
      pos % (_num_x * _num_y) 先把层数剥掉，只保留本层里的位置。比如 _num_x=4,_num_y=3，一层有 12 个格点；如果 pos=17，17 % 12 = 5，说明在当前层里是第 5 个位置（从 0 开始）。
再除以 _num_x（每一行 x 的数量）就得到 y 行号。继续上例：5 / 4 = 1，表示在 y=1 这一行，也就是第二行。
i 表示当前格点在 x 方向排第几个，用 “位置对 _num_x 取模”就能得到。
j 表示在 y 方向排第几个，先把位置按一层（_num_x * _num_y）取模，再除以 _num_x，相当于“去掉 x 的贡献后剩下的层内 y 索引”。
这样就可以用 _widths_x[i]、_widths_y[j] 查到当前格点在 x/y 方向的尺寸，从而计算该格点 (也即放在此处的 Universe) 的等效半径。*/
      auto u = *it;
      unique_radius[u->getId()] =
          std::max(unique_radius[u->getId()],
                   sqrt(_widths_x[i] * _widths_x[i] / 4.0 + _widths_y[j] * _widths_y[j] / 4.0));
    }
    return unique_radius;
  }

  /**
   * @brief Set the width of each Lattice cell.
   * @param width_x the width along the x-axis in centimeters
   * @param width_y the width along the y-axis in centimeters
   * @param width_z the width along the z-axis in centimeters
   */
  void RecLattice::setWidth(double width_x, double width_y, double width_z)
  {

    if (width_x <= 0 || width_y <= 0 || width_z <= 0)
      log::ferror("Unable to set the width of Lattice ID = %d "
                  "for x = %f, y = %f, and z = %f since they are not positive values",
                  _id, width_x, width_y, _width_z);

    _width_x = width_x;
    _width_y = width_y;
    _width_z = width_z;
  }

  /// \brief Sets the array of Universe pointers filling each Lattice cell.
  /// \details In OpenMOC, this is a helper method for SWIG to allow users to
  ///          assign Universes to a Lattice using a 2D Python list (list of
  ///          lists). In ANT-MOC, this method is used by geometry input class
  ///          to receive the universe layout of a lattice from an XML file.
  ///          The layout is provided as it looks like and converted to another
  ///          order in this method. That is, the input order is highest-z-
  ///          -first, which will be converted into lowest-z-section-first.
  ///          An example of that is C5G7 layout
  ///
  ///            __ the (0, 0, 0) universe in the input
  ///           /
  ///          u1 u2 u3   |
  ///          u2 u1 u3   | the highest section
  ///          u3 u3 u3   |
  ///
  ///          u1 u2 u3   |
  ///          u2 u1 u3   | the lowest section
  ///          u3 u3 u3   |
  ///           \__
  ///               the (0, 0, 0) universe after conversion
  ///
  /// \param num_z the number of Lattice cells along z
  /// \param num_y the number of Lattice cells along y
  /// \param num_x the number of Lattice cells along x
  /// \param universes the array of Universes for each Lattice cell
  /**
   * 函数的输入 universes 是一段已经排好序的一维数组（长度 num_x * num_y * num_z），来自 XML/脚本。它的排列规则是“先最高 z 层，再下一层…；在每一层里第一行表示该层的最上边，从左到右”。也就是说 universes[0] 对应 (x=0, y=0, z=最高层)，universes[_num_x-1] 是同一层最右列，接下来才是次上一行，行遍历完才换 z 层。
Lattice 内部维护的 _universes 才是后续查找 (findCell) 真正用的数据结构。它同样是一维数组，但假定索引公式 idx = k*_num_x*_num_y + j*_num_x + i 满足：k=0 = 最低 z 层，j=0 = 当前层最下那行，i=0 = 最左列。
setUniverses 要做的就是把输入顺序转换成内部顺序，于是三重循环：
k 是内部层索引（0=最低层）。要找对应输入层就用 _num_z-1-k。例如 _num_z=3，内部 k=0 属于第 0 层（最低层），输入里要取第 2 层（最高层），所以用 (_num_z-1-k)。
j 是内部行索引（0=最底行）。输入里 j=0 却代表最高行，所以同理用 _num_y-1-j 翻过来。
i 方向在输入和内部都是“从左到右”，无需翻转。
   */
  void RecLattice::setUniverses(int num_z, int num_y, int num_x,
                                Universe **universes)
  {

    // std::map<int, Universe*> unique_universes = getUniqueUniverses();

    /* Remove all Universes in the Lattice */
    // for (auto &u : unique_universes)
    //   removeUniverse(u.second);

    /* Clear any Universes in the Lattice (from a previous run) */
    clearUniverses();

    /* Set the Lattice dimensions */
    setNumX(num_x);
    setNumY(num_y);
    setNumZ(num_z);

    // The Lattice cells are assumed to be in row-major order starting from the
    // upper left corner at the highest z-section. This nested loop reorders the
    // Lattice cells to start from the lower left corner at the lowest z-section.
    // It actually represents the xyz coordinates
    // 把输入的 3D 网格（universes 指针数组）拷贝进 _universes，但顺序要翻转一下。
    _universes.resize(num_x * num_y * num_z, nullptr);
    for (int k = 0; k < _num_z; k++) // 对于层
    {
      for (int j = 0; j < _num_y; j++) // 对于行
      {
        for (int i = 0; i < _num_x; i++) // 对应列
        {
          auto idx_in = (_num_z - 1 - k) * _num_x * _num_y + (_num_y - 1 - j) * _num_x + i;
          auto idx = k * _num_x * _num_y + j * _num_x + i;
          // allow nullptr, delay the assignment
          _universes[idx] = universes[idx_in];
        }
      }
    }
    computeSizes();
  }

  /**
   * @brief Checks if a Point is within the bounds of a Lattice.
   * @param point a pointer to the Point of interest
   * @return true if the Point is in the bounds, false if not
   */
  bool RecLattice::containsPoint(Point *point)
  {

    /* Computes the Lattice bounds */
    double bound_x_max = getMaxX();
    double bound_x_min = getMinX();
    double bound_y_max = getMaxY();
    double bound_y_min = getMinY();
    double bound_z_max = getMaxZ();
    double bound_z_min = getMinZ();

    double x = point->getX();
    double y = point->getY();
    double z = point->getZ();

    /* If the Point is outside the x bounds */
    if (x > bound_x_max || x < bound_x_min)
      return false;

    /* If the Point is outside the y bounds */
    else if (y > bound_y_max || y < bound_y_min)
      return false;

    /* If the Point is outside the z bounds */
    else if (z > bound_z_max || z < bound_z_min)
      return false;

    /* If the Point is within the bounds */
    else
      return true;
  }

  /**
   * @brief Finds the distance to the nearest surface.
   * @details Knowing that a Lattice must be cartesian, this function computes
   *          the distance to the nearest boundary between lattice cells
   *          in the direction of the track.
   *          Returns distance to nearest Lattice cell boundary.
   * 该函数计算在轨迹方向上到晶格单元之间最近边界的距离。返回到最近晶格单元边界的距离
   * @param point a pointer to a starting point
   * @param azim the azimuthal angle of the track
   * @param polar the polar angle of the track
   * @return the distance to the nearest Lattice cell boundary
   */
  double RecLattice::minSurfaceDist(Point *point, double azim, double polar)
  {

    if (point->getX() != point->getX() ||
        point->getY() != point->getY())
    {
      log::ferror("Nan is found in points when try to compute the min distance "
                  "to surface: x = %f, y = %f, azim = %f, polar = %f",
                  point->getX(), point->getY(), azim, polar);
    }

    /* Compute the x, y, and z indices for the Lattice cell this point is in */
    int lat_x = getLatX(point);
    int lat_y = getLatY(point);
    int lat_z = getLatZ(point);

    /* Get unit vector components */
    double u_xy = sin(polar);
    double u_x = u_xy * cos(azim); // 射线在x方向上的单位分量
    double u_y = u_xy * sin(azim);
    double u_z = cos(polar);

    /* Get local coordinates of the point 当前晶格单元中的局部坐标*/
    double x = getLocalPointX(point, lat_x, lat_y);
    double y = getLocalPointY(point, lat_x, lat_y);
    double z = getLocalPointZ(point, lat_z);

    /* Get the min distance for X PLANE  */
    double dist_x;
    if (fabs(u_x) > FLT_EPSILON)
    { // copysign函数的功能是返回一个值，该值的大小与第一个参数相同，但符号与第二个参数相同。
      // 就是“取 0.5*width 的绝对值，符号取 u_x 的符号”
      double plane_x = std::copysign(0.5 * _widths_x[lat_x], u_x); // 确定当前点所在的晶格单元的下一个x方向的最近的平面位置
      dist_x = (plane_x - x) / u_x;                                // u_x 可以理解为射线在x方向上的速率，除以u_x  就是将平面距离转换为实际距离
    }
    else
    {
      dist_x = std::numeric_limits<double>::infinity();
    }

    /* Get the min distance for Y PLANE  */
    double dist_y;
    if (fabs(u_y) > FLT_EPSILON)
    {
      double plane_y = std::copysign(0.5 * _widths_y[lat_y], u_y);
      dist_y = (plane_y - y) / u_y;
    }
    else
    {
      dist_y = std::numeric_limits<double>::infinity();
    }

    /* Get the min distance for Z PLANE  */
    double dist_z;
    if (fabs(u_z) > FLT_EPSILON &&
        is3D())
    {
      double plane_z = std::copysign(0.5 * _widths_z[lat_z], u_z);
      dist_z = (plane_z - z) / u_z;
    }
    else
    {
      dist_z = std::numeric_limits<double>::infinity();
    }

    /* return shortest distance to next lattice cell */
    return std::min(dist_x, std::min(dist_y, dist_z));
  }

  /**
   * @brief Determine if the indices are outside the bound of the Lattice
   */
  bool RecLattice::areValidIndices(int lat_x, int lat_y, int lat_z) const
  {
    return lat_x >= 0 && lat_x < _num_x &&
           lat_y >= 0 && lat_y < _num_y &&
           lat_z >= 0 && lat_z < _num_z;
  }

  /**
   * @brief Determine if the index is outside the bound of the Lattice
   */
  bool RecLattice::isValidIndex(int index) const
  {
    return index >= 0 && index < getNumLatticeCells();
  }

  /// \brief Finds the Lattice cell x and y indices that a point lies in.
  /// \param point a pointer to a point being evaluated.
  /// \return the Lattice cell x and y indices.
  std::array<int, 2> RecLattice::getLatXY(Point *point)
  {
    return {getLatX(point), getLatY(point)};
  }

  /**
   * @brief Finds the Lattice cell x index that a point lies in.
   * @param point a pointer to a point being evaluated.
   * @return the Lattice cell x index.
   */
  int RecLattice::getLatX(Point *point)
  {

    int lat_x = -1;

    /* get the distance to the left surface */
    // 如果晶格最左面的绝对位置是 getMinX()，那么 dist_to_left = point->getX() - getMinX() 就代表该点沿 x 方向距离左边界多远
    double dist_to_left = point->getX() - getMinX();

    //_accumulate_x 中存的就是从左到右的累计宽度（第 0 格左边界=0，第 1 格右界=宽度1，第 2 格右界=宽度1+宽度2 ...）。for 循环 (Lattice.cpp (lines 808-815)) 找到满足 _accumulate_x[i] <= dist_to_left < _accumulate_x[i+1] 的区间，说明该点位于第 i 列。
    /* Compute the x index for the Lattice cell this point is in */
    for (int i = 0; i < _num_x; i++)
    {
      if (dist_to_left >= _accumulate_x[i] && dist_to_left < _accumulate_x[i + 1])
      {
        lat_x = i;
        break;
      }
    }

    /* Check if the Point is on the Lattice boundaries and if so adjust
     * x Lattice cell index */
    if (fabs(dist_to_left) < ON_SURFACE_THRESH)
      lat_x = 0;
    else if (fabs(dist_to_left - _accumulate_x[_num_x]) < ON_SURFACE_THRESH) // _accumulate_x 列表包含 _num_x+1 个值，代表每一条垂直网格线相对最左边界的累积距离（第 0 条线在 0，最后一条线在总宽度处），所以 _accumulate_x[_num_x] 是“最右边界线的位置”

      lat_x = _num_x - 1;
    if (lat_x == -1)
      log::ferror("Trying to get lattice x index for point(x = %f) that is "
                  "outside lattice bounds. dist_to_left = %f is not within "
                  "[0.0, %f]",
                  point->getX(), dist_to_left, _accumulate_x[_num_x]);

    return lat_x;
  }

  /**
   * @brief Finds the Lattice cell y index that a point lies in.
   * @param point a pointer to a point being evaluated.
   * @return the Lattice cell y index.
   */
  int RecLattice::getLatY(Point *point)
  {

    int lat_y = -1;

    /* get the distance to the bottom surface */
    double dist_to_bottom = point->getY() - getMinY();

    /* Compute the y index for the Lattice cell this point is in */
    for (int i = 0; i < _num_y; i++)
    {
      if (dist_to_bottom >= _accumulate_y[i] &&
          dist_to_bottom < _accumulate_y[i + 1])
      {
        lat_y = i;
        break;
      }
    }

    /* Check if the Point is on the Lattice boundaries and if so adjust
     * y Lattice cell index */
    if (fabs(dist_to_bottom) < ON_SURFACE_THRESH)
      lat_y = 0;
    else if (fabs(dist_to_bottom - _accumulate_y[_num_y]) < ON_SURFACE_THRESH)
      lat_y = _num_y - 1;
    if (lat_y == -1)
      log::ferror("Trying to get lattice y index for point(y = %f) that is "
                  "outside lattice bounds. dist_to_bottom = %f is not within "
                  "[0.0, %f]",
                  point->getY(), dist_to_bottom, _accumulate_y[_num_y]);

    return lat_y;
  }

  /**
   * @brief Converts a Lattice's attributes to a character array representation.
   * @return character array of this Lattice's attributes
   */
  std::string RecLattice::toString()
  {

    std::stringstream string;

    string << "Lattice ID = " << _id
           << ", name = " << _name
           << ", # cells along x = " << _num_x
           << ", # cells along y = " << _num_y
           << ", # cells along z = " << _num_z;

    if (_non_uniform)
    {
      string << "\nThis lattice is non-uniform.\nx widths: ";
      for (int i = 0; i < _num_x; i++)
        string << _widths_x[i] << "  ";
      string << "\ny widths: ";
      for (int i = 0; i < _num_y; i++)
        string << _widths_y[i] << "  ";
      string << "\nz widths: ";
      for (int i = 0; i < _num_z; i++)
        string << _widths_z[i] << "  ";
    }
    else
      string << ", x width = " << _width_x
             << ", y width = " << _width_y
             << ", z width = " << _width_z;

    string << ", type = Rectangle ";

    string << "\n\t\tUniverse IDs within this Lattice: ";

    for (auto it = begin(); it != end(); ++it)
    {
      auto u = *it;
      if (u)
        string << u->getId();
      else
        string << -1;
      string << ", ";

      if (it.getPos() % _num_x == 0)
        string << "\n\t\t";
    }

    return string.str();
  }

  /**
   * @brief Return a string containing the numbers of lattice cells
   */
  std::string RecLattice::dimToString()
  {

    std::stringstream string;

    string << "nx = " << _num_x
           << ", ny = " << _num_y
           << ", nz = " << _num_z;
    return string.str();
  }

  /**
   * @brief Finds the Lattice cell index that a point lies in.
   * @details Lattice cells are numbered starting with 0 in the lower left
   *          corner. Lattice cell IDs in all rows then increase monotonically
   *          from left to right. For example, the indices for a 4 x 4 lattice:
   *                  12  13  14  15
   *                  8    9  10  11
   *                  4    5   6   7
   *                  0    1   2   3
   * @param point a pointer to a point being evaluated.
   * @return the Lattice cell index.
   */
  int RecLattice::getLatticeCell(Point *point)
  {
    return getLatZ(point) * _num_x * _num_y + getLatY(point) * _num_x + getLatX(point);
  }

  /**
   * @brief Finds the Lattice cell surface that a point lies on.
   *        If the point is not on a surface, -1 is returned.
   * @details The surface indices for a lattice cell are 0 (left),
   *         1, (bottom), 2 (right), 3 (top), 4 (bottom-left corner),
   *         5 (bottom-right corner), 6 (top-right corner), and
   *         7 (top-left corner). The index returned takes into account
   *         the cell index and returns 8*cell_index + surface_index.
   * @param cell the cell index that the point is in.
   * @param point a pointer to a point being evaluated.
   * @return the Lattice surface index.
   *  查找2维点所在的Lattice单元格曲面。如果该点不在曲面上，则返回-1
   */
  int RecLattice::getLatticeSurface(int cell, Point *point, double azim, double polar)
  {

    int surface = -1;

    /* Get coordinates of point and cell boundaries */
    // unused
    // double x = point->getX();
    // double y = point->getY();
    // double z = point->getZ();
    int lat_x = (cell % (_num_x * _num_y)) % _num_x; // 获取三维全局索引
    int lat_y = (cell % (_num_x * _num_y)) / _num_x;
    int lat_z = cell / (_num_x * _num_y);

    /* Create planes representing the boundaries of the lattice cell 创建表示晶格单元边界的平面*/
    // 创建三个平面对象 XPlane/YPlane/ZPlane，初始放在原点；它们随后会被移动到对应单元面的真实坐标。
    XPlane xplane(0.0);
    YPlane yplane(0.0);
    ZPlane zplane(0.0);

    /* Bools indicating if point is on each surface 声明六个布尔 on_min_x 等，分别表示点是否落在该单元的六个面上。 */
    bool on_min_x, on_max_x, on_min_y, on_max_y, on_min_z, on_max_z;
    // 依次对 X/Y 方向的最小/最大面赋值：使用 _accumulate_x[lat_x] 这类累积距离 + 全局最小坐标 getMinX() 来得到当前单元左/右界位置，设置平面，再调用 isPointOnSurface(point) 检查点是否正好在那个面上
    /* Check if point is on X_MIN boundary */
    xplane.setX(_accumulate_x[lat_x] + getMinX());
    on_min_x = xplane.isPointOnSurface(point);

    /* Check if point is on X_MAX boundary */
    xplane.setX(_accumulate_x[lat_x + 1] + getMinX());
    on_max_x = xplane.isPointOnSurface(point);

    /* Check if point is on Y_MIN boundary */
    yplane.setY(_accumulate_y[lat_y] + getMinY());
    on_min_y = yplane.isPointOnSurface(point);

    /* Check if point is on Y_MAX boundary */
    yplane.setY(_accumulate_y[lat_y + 1] + getMinY());
    on_max_y = yplane.isPointOnSurface(point);

    /* Check if point is on Z_MIN boundary */
    on_min_z = false;
    if (is3D())
    {
      zplane.setZ(_accumulate_z[lat_z] + getMinZ());
      on_min_z = zplane.isPointOnSurface(point);
    }

    /* Check if point is on Z_MAX boundary */
    on_max_z = false;
    if (is3D())
    {
      zplane.setZ(_accumulate_z[lat_z + 1] + getMinZ());
      on_max_z = zplane.isPointOnSurface(point);
    }

    if (on_min_x)
    { // 通过嵌套的条件语句判断该点在哪些平面上
      /*
       通过嵌套 if/else 判断点落在的具体组合：优先看 X 方向，再看 Y，再看 Z。只要点在多个方向交线/角上，就返回不同枚举值，例如 SURFACE_X_MIN_Y_MIN_Z_MIN 代表该单元左-前-底角；这些枚举常量（注释如 // 18）对应局部面/棱/角编号。
      */
      if (on_min_y)
      {
        if (on_min_z)
          surface = SURFACE_X_MIN_Y_MIN_Z_MIN; // 18
        else if (on_max_z)
          surface = SURFACE_X_MIN_Y_MIN_Z_MAX; // 19
        else
          surface = SURFACE_X_MIN_Y_MIN; // 006
      }
      else if (on_max_y)
      {
        if (on_min_z)
          surface = SURFACE_X_MIN_Y_MAX_Z_MIN; // 20
        else if (on_max_z)
          surface = SURFACE_X_MIN_Y_MAX_Z_MAX; // 21
        else
          surface = SURFACE_X_MIN_Y_MAX; // 008
      }
      else
      {
        if (on_min_z)
          surface = SURFACE_X_MIN_Z_MIN; // 10
        else if (on_max_z)
          surface = SURFACE_X_MIN_Z_MAX; // 12
        else
          surface = SURFACE_X_MIN; // 000
      }
    }
    else if (on_max_x)
    {
      if (on_min_y)
      {
        if (on_min_z)
          surface = SURFACE_X_MAX_Y_MIN_Z_MIN; // 22
        else if (on_max_z)
          surface = SURFACE_X_MAX_Y_MIN_Z_MAX; // 23
        else
          surface = SURFACE_X_MAX_Y_MIN; // 007
      }
      else if (on_max_y)
      {
        if (on_min_z)
          surface = SURFACE_X_MAX_Y_MAX_Z_MIN; // 24
        else if (on_max_z)
          surface = SURFACE_X_MAX_Y_MAX_Z_MAX; // 25
        else
          surface = SURFACE_X_MAX_Y_MAX; // 009
      }
      else
      {
        if (on_min_z)
          surface = SURFACE_X_MAX_Z_MIN; // 11
        else if (on_max_z)
          surface = SURFACE_X_MAX_Z_MAX; // 13
        else
          surface = SURFACE_X_MAX; // 003
      }
    }
    else if (on_min_y)
    {
      if (on_min_z)
        surface = SURFACE_Y_MIN_Z_MIN; // 14
      else if (on_max_z)
        surface = SURFACE_Y_MIN_Z_MAX; // 16
      else
        surface = SURFACE_Y_MIN; // 001
    }
    else if (on_max_y)
    {
      if (on_min_z)
        surface = SURFACE_Y_MAX_Z_MIN; // 15
      else if (on_max_z)
        surface = SURFACE_Y_MAX_Z_MAX; // 17
      else
        surface = SURFACE_Y_MAX; // 004
    }
    else if (on_min_z)
      surface = SURFACE_Z_MIN; // 2
    else if (on_max_z)
      surface = SURFACE_Z_MAX; // 5

    if (surface != -1)
      surface = NUM_SURFACES * cell + surface; // 每个CMFD 网格都有26个面点边面，返回的是全局的CMFD surface数量

    return surface;
  }

  /**
   * @brief Finds the Lattice cell surface that a point lies on.
   *        If the point is not on a surface, -1 is returned.
   * 查找Z轴点所在的CMFD晶格单元表面。如果点不在表面上，则返回-1
   * @details The surface indices are defined in constants.h as they
   *          need to be consistent with the surface constant definitions
   *          used in Cmfd. The index returned takes into account
   *          the cell index and returns NUM_SURFACES*cell_index + surface_index.
   * @param cell the cell index that the point is in
   * @param z z coordinate of the point
   * @param surface_2D 2D surface considered
   * @return the Lattice surface index.
   * NUM_SURFACES=26，为：6个面，8个顶点，12个边
   */
  int RecLattice::getLatticeSurfaceOTF(int cell, double z, int surface_2D)
  {

    /* Determine min and max z boundaries of the cell 该CMFD所在的全局Z边界值*/
    int lat_z = cell / (_num_x * _num_y);
    double z_min = _accumulate_z[lat_z] + getMinZ();
    double z_max = z_min + _widths_z[lat_z];

    /* Check for z-surface crossing on 2D surface 检查surface_2D是否是2维的（2维面的索引只有10个，6面（包括了Z轴的两个面）和四个边）*/
    if (surface_2D % NUM_SURFACES > 9)
      log::ferror("Found a z-surface crossing on a 2D segment");

    /* Check min z boundary for crossing */
    if (fabs(z_min - z) < TINY_MOVE)
    { // 该点是否位于该CMFD的底面上
      int surface;
      switch (surface_2D % NUM_SURFACES)
      {
      case SURFACE_X_MIN: // 位于x-z左边 10
        surface = SURFACE_X_MIN_Z_MIN;
        break;
      case SURFACE_X_MAX: // 位于x-z右边 11
        surface = SURFACE_X_MAX_Z_MIN;
        break;
      case SURFACE_Y_MIN: // 位于y-z前边 14
        surface = SURFACE_Y_MIN_Z_MIN;
        break;
      case SURFACE_Y_MAX: // 位于y-z后边 14
        surface = SURFACE_Y_MAX_Z_MIN;
        break;
      case SURFACE_X_MIN_Y_MIN: // 18
        surface = SURFACE_X_MIN_Y_MIN_Z_MIN;
        break;
      case SURFACE_X_MIN_Y_MAX: // 20
        surface = SURFACE_X_MIN_Y_MAX_Z_MIN;
        break;
      case SURFACE_X_MAX_Y_MIN: // 22
        surface = SURFACE_X_MAX_Y_MIN_Z_MIN;
        break;
      case SURFACE_X_MAX_Y_MAX: // 24
        surface = SURFACE_X_MAX_Y_MAX_Z_MIN;
        break;
      default:
        surface = SURFACE_Z_MIN; // 2
      }
      return NUM_SURFACES * cell + surface; // 乘以26是给每个单元格一个基准值，对于单元格 0，它的基准值是 0，对于单元格 1，它的基准值是 26，后面的surface则为具体一个cell的哪一个面（包括点边面）
    }

    /* Check max z boundary for crossing */
    if (fabs(z_max - z) < TINY_MOVE)
    { // 该点是否位于该CMFD的顶面上
      int surface;
      switch (surface_2D % NUM_SURFACES)
      {
      case SURFACE_X_MIN:
        surface = SURFACE_X_MIN_Z_MAX;
        break;
      case SURFACE_X_MAX:
        surface = SURFACE_X_MAX_Z_MAX;
        break;
      case SURFACE_Y_MIN:
        surface = SURFACE_Y_MIN_Z_MAX;
        break;
      case SURFACE_Y_MAX:
        surface = SURFACE_Y_MAX_Z_MAX;
        break;
      case SURFACE_X_MIN_Y_MIN:
        surface = SURFACE_X_MIN_Y_MIN_Z_MAX;
        break;
      case SURFACE_X_MIN_Y_MAX:
        surface = SURFACE_X_MIN_Y_MAX_Z_MAX;
        break;
      case SURFACE_X_MAX_Y_MIN:
        surface = SURFACE_X_MAX_Y_MIN_Z_MAX;
        break;
      case SURFACE_X_MAX_Y_MAX:
        surface = SURFACE_X_MAX_Y_MAX_Z_MAX;
        break;
      default:
        surface = SURFACE_Z_MAX;
      }
      return NUM_SURFACES * cell + surface;
    }

    /* If no axial crossing, return the 2D surface */
    if (surface_2D == -1)
      return surface_2D;
    else                                                      // 如果没有轴向交叉，则返回二维曲面
      return NUM_SURFACES * cell + surface_2D % NUM_SURFACES; // 如果没有和z轴交叉，并且surface_2D有效，则直接这样计算（就没有加Z轴的数据）
  }

  /**
   * @brief Set widths of non-uniform meshes in x y z directions.
   * @details An example of how this may be called from Python illustrated below:
   *
   * @code
   *          RecLattice::setWidths([1.0,2.1,3.0], [4.0,5.1,6.0,7.0], [3.3,2.4])
   * @endcode
   *
   * @param widths_x x-direction widths of non-uniform meshes
   * @param widths_y y-direction widths of non-uniform meshes
   * @param widths_z z-direction widths of non-uniform meshes
   */
  void RecLattice::setWidths(DoubleVec widths_x,
                             DoubleVec widths_y, DoubleVec widths_z)
  {
    _non_uniform = true;
    _widths_x = widths_x;
    _widths_y = widths_y;
    _widths_z = widths_z;
  }

  /**
   * @brief Set _widths_x, _widths_y, _widths_z for uniform case, compute
   *        accumulate variables.
   */
  void RecLattice::computeSizes()
  {
    if (_non_uniform)
    {
      if ((_widths_x.size() - _num_x) || (_widths_y.size() - _num_y) ||
          (_widths_z.size() - _num_z)) // 检查 x、y 和 z 方向上非均匀网格尺寸数组的大小是否与对应方向上的Lattice网格的数量一致
        log::ferror("The sizes of non-uniform mesh widths are not consistent"
                    " with the sizes of filling Universes into Lattice");
    }
    else
    {
      _widths_x.resize(_num_x, _width_x);
      _widths_y.resize(_num_y, _width_y);
      _widths_z.resize(_num_z, _width_z);
    }

    /* Compute the accumulated lengths along each axis */
    _accumulate_x.resize(_num_x + 1, 0.0);
    _accumulate_y.resize(_num_y + 1, 0.0);
    _accumulate_z.resize(_num_z + 1, 0.0);

    // 长度为 _num_x+1 的前缀和，_accumulate_x[i] 是第 i 条纵向网格线相对于最小 x 的位置，_accumulate_x[_num_x] 则是整个晶格的总宽。对应的 _accumulate_y/_accumulate_z 处理 y/z 方向。
    for (int i = 0; i < _num_x; i++)
      _accumulate_x[i + 1] = _accumulate_x[i] + _widths_x[i];

    for (int i = 0; i < _num_y; i++)
      _accumulate_y[i + 1] = _accumulate_y[i] + _widths_y[i];

    for (int i = 0; i < _num_z; i++)
      _accumulate_z[i + 1] = _accumulate_z[i] + _widths_z[i];
  }

  /// \brief Returns the center point of a lattice cell with respect to the center
  ///        of the center of parent universe
  Point RecLattice::getLatticeCellCenter(int lat_x, int lat_y, int lat_z)
  {

    // Lattice cell center, local point
    // center 是在 lattice 自己的局部坐标系下计算的格子单元中心：左边界偏移 _accumulate_x[lat_x] 加上当前单元宽度的一半，y/z 一样
    Point center = {_accumulate_x[lat_x] + _widths_x[lat_x] / 2.,
                    _accumulate_y[lat_y] + _widths_y[lat_y] / 2.,
                    _accumulate_z[lat_z] + _widths_z[lat_z] / 2.};

    // Lattice center, local point
    // lat_center 是整个 lattice 的几何中心，同样用总累计宽的一半表示
    Point lat_center = {_accumulate_x[_num_x] / 2.,
                        _accumulate_y[_num_y] / 2.,
                        _accumulate_z[_num_z] / 2.};

    return _offset + center - lat_center;
    // 先把单元中心相对 lattice 中心的偏移量 center - lat_center 算出来，再加上 _offset（lattice 在父宇宙里的实际中心位置），就得到了该单元中心在父 Universe/global 坐标系中的点
  }

  /// \brief Returns coordinates of the 4 or 8 vertices of a lattice cell
  std::vector<Point> RecLattice::getLatticeCellVertices(int index)
  {

    std::vector<Point> vertices;

    // Gets the lattice cell center
    Point center = getLatticeCellCenter(index);

    // Defines a radius vector
    int lat_x = index % _num_x;
    double x = _widths_x[lat_x] / 2.;

    int lat_y = (index / _num_x) % _num_y;
    double y = _widths_y[lat_y] / 2.;

    // FIXME: support 2-D lattice
    int lat_z = index / (_num_x * _num_y);
    double z = _widths_z[lat_z] / 2.;

    for (auto z_sign : std::vector<int>{-1, 1})
    {
      vertices.push_back({x + center.getX(),
                          y + center.getY(),
                          z_sign * z + center.getZ()});

      vertices.push_back({-x + center.getX(),
                          y + center.getY(),
                          z_sign * z + center.getZ()});

      vertices.push_back({-x + center.getX(),
                          -y + center.getY(),
                          z_sign * z + center.getZ()});

      vertices.push_back({x + center.getX(),
                          -y + center.getY(),
                          z_sign * z + center.getZ()});
    }

    return vertices;
  }

  /// \brief For debug use.
  void RecLattice::printLatticeSizes()
  {
    int i;
    printf("non_uniform=%d, \nNum_XYZ: %2d, %2d, %2d\n", _non_uniform,
           _num_x, _num_y, _num_z);
    printf("offset: %f, %f, %f\n", _offset.getX(), _offset.getY(), _offset.getZ());
    printf("cell_width_XYZ: %f, %f, %f\n", _width_x, _width_y, _width_z);
    printf("cell_widths_XYZ:\n");
    for (i = 0; i < _num_x; i++)
      printf("i=%d, %f; ", i, _widths_x[i]);
    printf("\n");
    for (i = 0; i < _num_y; i++)
      printf("i=%d, %f; ", i, _widths_y[i]);
    printf("\n");
    for (i = 0; i < _num_z; i++)
      printf("i=%d, %f; ", i, _widths_z[i]);
    printf("\n");

    printf("accumulates_XYZ:\n");
    for (i = 0; i < _num_x + 1; i++)
      printf("i=%d, %f; ", i, _accumulate_x[i]);
    printf("\n");
    for (i = 0; i < _num_y + 1; i++)
      printf("i=%d, %f; ", i, _accumulate_y[i]);
    printf("\n");
    for (i = 0; i < _num_z + 1; i++)
      printf("i=%d, %f; ", i, _accumulate_z[i]);
    printf("\n");
  }

  //////////////////////////////////////////////////
  //                                              //
  //                  HexLattice                  //
  //                                              //
  //////////////////////////////////////////////////

  HexLattice::HexLattice(const int id, const char *name) : Lattice(id, name)
  {

    _lattice_type = latticeType::Hexagon;

    /* Default width and number of Lattice cells along each dimension */
    _orientation = Orientation::y;
    _num_r = 0;
    _width_r = 0;
  }

  /// \brief Get the begin iterator
  LatticeIter HexLattice::begin()
  {
    return LatticeIter(*this, std::min((size_t)_num_r - 1, _universes.size()));
  }

  /// \brief Returns the minimum reachable x-coordinate in the Lattice.
  double HexLattice::getMinX() const
  {
    double x = 0.;
    if (isOrientationY())
    {
      double cos30 = std::cos(M_PI / 6);
      x = _width_r * ((_num_r - 1) * cos30 + 0.5 / cos30);
    }
    else
    {
      x = _width_r * (_num_r - 0.5);
    }
    return _offset.getX() - x;
  }

  /// \brief Returns the maximum reachable x-coordinate in the Lattice.
  double HexLattice::getMaxX() const
  {
    double x = 0.;
    if (isOrientationY())
    {
      double cos30 = std::cos(M_PI / 6);
      x = _width_r * ((_num_r - 1) * cos30 + 0.5 / cos30);
    }
    else
    {
      x = _width_r * (_num_r - 0.5);
    }
    return _offset.getX() + x;
  }

  /// \brief Returns the minimum reachable y-coordinate in the Lattice.
  double HexLattice::getMinY() const
  {
    double x = 0.;
    if (isOrientationY())
    {
      x = _width_r * (_num_r - 0.5);
    }
    else
    {
      double cos30 = std::cos(M_PI / 6);
      x = _width_r * ((_num_r - 1) * cos30 + 0.5 / cos30);
    }
    return _offset.getY() - x;
  }

  /// \brief Returns the maximum reachable y-coordinate in the Lattice.
  double HexLattice::getMaxY() const
  {
    double x = 0.;
    if (isOrientationY())
    {
      x = _width_r * (_num_r - 0.5);
    }
    else
    {
      double cos30 = std::cos(M_PI / 6);
      x = _width_r * ((_num_r - 1) * cos30 + 0.5 / cos30);
    }
    return _offset.getY() + x;
  }

  std::map<int, double> HexLattice::getUniqueRadius(std::map<int, Universe *> &unique_universes)
  {

    std::map<int, double> unique_radius;

    /* Create and initialize the <universe ID, unique radius> map */
    for (auto &u : unique_universes)
      unique_radius[u.first] = 0.;

    /* Get the maximum equivalent radius of each unique universe */
    double radius = _width_r / std::sqrt(3);
    for (auto u : *this)
      unique_radius[u->getId()] = radius;

    return unique_radius;
  }

  /// \brief Checks if a Point is within the bounds of a Lattice.
  /// \param point a pointer to the Point of interest
  /// \return true if the Point is in the bounds, false if not
  bool HexLattice::containsPoint(Point *point)
  {

    // If the Point is outside the z bounds
    double z = point->getZ();
    if (z > getMaxZ() || z < getMinZ())
    {
      return false;
    }
    else
    {
      // Compute the lattice cell the point lies in
      auto lat_xy = getLatXY(point);
      int lat_x = lat_xy[0];
      int lat_y = lat_xy[1];
      int lat_z = getLatZ(point);

      return areValidIndices(lat_x, lat_y, lat_z);
    }
  }

  /// \brief Set the orientation of the HexLattice.
  /// \details The orientation is defined according to the direction
  ///          of the flat side.
  void HexLattice::setOrientation(std::string orientation)
  {
    stringutils::trim(orientation);
    stringutils::toUpper(orientation);

    if (orientation == "Y")
    {
      _orientation = Orientation::y;
    }
    else
    {
      _orientation = Orientation::x;
    }
  }

  /// \brief Set the orientation of the HexLattice.
  void HexLattice::setOrientation(Orientation orientation)
  {
    _orientation = orientation;
  }

  /// \brief Finds the distance to the nearest surface.
  /// \details See OpenMC. Different from the algorithm implemented in OpenMC,
  ///          the distance is computed based on the current lattice cell the
  ///          point lies in.
  double HexLattice::minSurfaceDist(Point *point, double azim, double polar)
  {

    if (point->getX() != point->getX() ||
        point->getY() != point->getY())
    {
      log::ferror("Nan is found in points when try to compute the min distance "
                  "to surface: x = %f, y = %f, azim = %f, polar = %f",
                  point->getX(), point->getY(), azim, polar);
    }

    // Short description of the direction vectors used here.  The beta, gamma, and
    // delta vectors point towards the flat sides of each hexagonal tile.
    // Y - orientation:
    //   basis0 = (1, 0)
    //   basis1 = (-1/sqrt(3), 1)   = +120 degrees from basis0
    //   beta   = (sqrt(3)/2, 1/2)  = +30 degrees from basis0
    //   gamma  = (sqrt(3)/2, -1/2) = -60 degrees from beta
    //   delta  = (0, 1)            = +60 degrees from beta
    // X - orientation:
    //   basis0 = (1/sqrt(3), -1)
    //   basis1 = (0, 1)            = +120 degrees from basis0
    //   beta   = (1, 0)            = +30 degrees from basis0
    //   gamma  = (1/2, -sqrt(3)/2) = -60 degrees from beta
    //   delta  = (1/2, sqrt(3)/2)  = +60 degrees from beta
    // The z-axis is considered the same as RecLattice.

    // Get unit vector components
    double u_xy = sin(polar);
    double u_x = u_xy * cos(azim);
    double u_y = u_xy * sin(azim);
    double u_z = cos(polar);

    // Compute the dot production
    double u_beta, u_gamma, u_delta;
    if (_orientation == Orientation::y)
    {
      u_beta = u_x * std::sqrt(3.0) / 2.0 + u_y / 2.0;
      u_gamma = u_x * std::sqrt(3.0) / 2.0 - u_y / 2.0;
      u_delta = u_y;
    }
    else
    {
      u_beta = u_x;
      u_gamma = u_x / 2.0 - u_y * std::sqrt(3.0) / 2.0;
      u_delta = u_x / 2.0 + u_y * std::sqrt(3.0) / 2.0;
    }

    // Compute the x, y, and z indices for the Lattice cell this point is in
    auto lat_xy = getLatXY(point);
    int lat_x = lat_xy[0];
    int lat_y = lat_xy[1];
    int lat_z = getLatZ(point);

    // Get local coordinates of the point
    double x = getLocalPointX(point, lat_x, lat_y);
    double y = getLocalPointY(point, lat_x, lat_y);
    double z = getLocalPointZ(point, lat_z);

    double dist = std::numeric_limits<double>::infinity();

    // Get the min distance in the beta direction
    if (fabs(u_beta) > FLT_EPSILON)
    {
      double plane = std::copysign(0.5 * _width_r, u_beta);

      // (x,y) \dot beta
      double l_beta;
      if (_orientation == Orientation::y)
      {
        l_beta = x * std::sqrt(3.0) / 2.0 + y / 2.0;
      }
      else
      {
        l_beta = x;
      }

      dist = (plane - l_beta) / u_beta;
    }

    /* Get the min distance in the gamma direction */
    if (fabs(u_gamma) > FLT_EPSILON)
    {
      double plane = std::copysign(0.5 * _width_r, u_gamma);

      // (x,y) \dot gamma
      double l_gamma;
      if (_orientation == Orientation::y)
      {
        l_gamma = x * std::sqrt(3.0) / 2.0 - y / 2.0;
      }
      else
      {
        l_gamma = x / 2.0 - y * std::sqrt(3.0) / 2.0;
      }

      double d_gamma = (plane - l_gamma) / u_gamma;
      dist = std::min(dist, d_gamma);
    }

    /* Get the min distance in the delta direction */
    if (fabs(u_delta) > FLT_EPSILON)
    {
      double plane = std::copysign(0.5 * _width_r, u_delta);

      // (x,y) \dot delta
      double l_delta;
      if (_orientation == Orientation::y)
      {
        l_delta = y;
      }
      else
      {
        l_delta = x / 2.0 + y * std::sqrt(3.0) / 2.0;
      }

      double d_delta = (plane - l_delta) / u_delta;
      dist = std::min(dist, d_delta);
    }

    /* Get the min distance along z-axis */
    if (fabs(u_z) > FLT_EPSILON && is3D())
    {
      double plane_z = std::copysign(0.5 * _widths_z[lat_z], u_z);
      double dist_z = (plane_z - z) / u_z;
      dist = std::min(dist, dist_z);
    }

    /* return shortest distance to next lattice cell */
    return dist;
  }

  /// \brief Determine if the indices are outside the bound of the Lattice
  /// \details This method assumes that a lattice cell should be within
  ///          the six sides of a hexgon formed by connecting the centers
  ///          of each cell on the boundary.
  bool HexLattice::areValidIndices(int lat_x, int lat_y, int lat_z) const
  {
    return ((lat_x >= 0) && (lat_y >= 0) && (lat_z >= 0) && (lat_x < 2 * _num_r - 1) && (lat_y < 2 * _num_r - 1) && (lat_x + lat_y > _num_r - 2) && (lat_x + lat_y < 3 * _num_r - 2) && (lat_z < _num_z));
  }

  bool HexLattice::isValidIndex(int index) const
  {
    int nx = 2 * _num_r - 1;
    int nxy = nx * nx; // nx == ny
    int lat_z = index / nxy;
    int lat_y = (index - nxy * lat_z) / nx;
    int lat_x = index - nxy * lat_z - nx * lat_y;
    return areValidIndices(lat_x, lat_y, lat_z);
  }

  /// \brief Finds the Lattice cell x index that a point lies in.
  /// \param point a pointer to a point being evaluated.
  /// \return the Lattice cell x index.
  int HexLattice::getLatX(Point *point)
  {
    auto i_xyz = getLatXY(point);
    return i_xyz[0];
  }

  /// \brief Finds the Lattice cell y index that a point lies in.
  /// \param point a pointer to a point being evaluated.
  /// \return the Lattice cell y index.
  int HexLattice::getLatY(Point *point)
  {
    auto i_xyz = getLatXY(point);
    return i_xyz[1];
  }

  /// \brief Get the x-coord of a point in the local system of a universe.
  ///        See OpenMC.
  /// \details First, move the point from universe system to lattice system.
  ///          Then, consider the contributions of lat_x, lat_y, and lat_z.
  ///          On x-y plane, (lat_x, lat_y, lat_z) = (0, 0, 0) is at the lower
  ///          left corner, which means the previous center is mapped to
  ///          (_num_r - 1, _num_r - 1, (_num_r - 1)/2)
  /// \param point a pointer to a point being evaluated.
  /// \return the local x-coord
  double HexLattice::getLocalPointX(const Point *p, int lat_x, int lat_y)
  {
    if (_orientation == Orientation::y)
    {
      // x_l = x_g - (center + pitch_x*cos(30)*index_x)
      return p->getX() - _offset.getX() - std::sqrt(3.) / 2. * _width_r * (lat_x - _num_r + 1);
    }
    else
    {
      // x_l = x_g - (center + pitch_x*index_x + pitch_y*sin(30)*index_y)
      return p->getX() - (_offset.getX() + (lat_x - _num_r + 1) * _width_r + (lat_y - _num_r + 1) * _width_r * 0.5);
    }
  }

  /// \brief Get the y-coord of a point in the local system of a universe.
  ///        See OpenMC.
  double HexLattice::getLocalPointY(const Point *p, int lat_x, int lat_y)
  {
    if (_orientation == Orientation::y)
    {
      // y_l = y_g - (center + pitch_x*index_x + pitch_y*sin(30)*index_y)
      return p->getY() - _offset.getY() - 0.5 * _width_r * (lat_x - _num_r + 1) - _width_r * (lat_y - _num_r + 1);
    }
    else
    {
      // y_l = y_g - (center + pitch_y*cos(30)*index_y)
      return p->getY() - _offset.getY() - std::sqrt(3.0) / 2.0 * (lat_y - _num_r + 1) * _width_r;
    }
  }

  /// \brief Return the index of a lattice cell
  /// \details Each of the lattice cells has an index but not all of them
  ///          are valid. For example, getLatticeCell(0, 0, 0) gives 0,
  ///          which is an invalid index.
  int HexLattice::getLatticeCell(int lat_x, int lat_y, int lat_z) const
  {
    return (2 * _num_r - 1) * (2 * _num_r - 1) * lat_z + (2 * _num_r - 1) * lat_y + lat_x;
  }

  int HexLattice::getLatticeCell(Point *point)
  {
    auto lat_xy = getLatXY(point);
    int lat_x = lat_xy[0];
    int lat_y = lat_xy[1];
    int lat_z = getLatZ(point);
    if (lat_x < 0 || lat_y < 0 || lat_z < 0 || lat_z > _num_z - 1 || lat_x > 2 * _num_r - 2 || lat_y > 2 * _num_r - 2)
    {
      return -1;
    }
    else
    {
      // log::finfo("Point is %.2f %.2f %.2f", point->getX(), point->getY(), point->getZ());
      return getLatticeCell(lat_x, lat_y, lat_z);
    }
  }

  // need fix, should choose the nearest surface, rather disthe < FLT_EPSILON
  /*
  int HexLattice::getLatticeSurface(int cell, Point* point, double azim, double polar) {

    int surface = -1;
    log::finfo("Azim:%lf", azim);
    // Get unit vector components
    double u_x = cos(azim);
    double u_y = sin(azim);

    // Compute the dot production
    double u_beta, u_gamma, u_delta;
    if (_orientation == Orientation::y) {
      u_beta = u_x * std::sqrt(3.0) / 2.0  + u_y / 2.0;
      u_gamma = u_x * std::sqrt(3.0) / 2.0  - u_y / 2.0;
      u_delta = u_y;
    }
    else {
      u_beta = u_x;
      u_gamma = u_x / 2.0  - u_y * std::sqrt(3.0) / 2.0;
      u_delta = u_x / 2.0  + u_y * std::sqrt(3.0) / 2.0;
    }

    // Compute the x, y, and z indices for the Lattice cell this point is in
    auto lat_xy = getLatXY(point);
    int lat_x = lat_xy[0];
    int lat_y = lat_xy[1];
    int lat_z = getLatZ(point);

    // Get local coordinates of the point
    double x = getLocalPointX(point, lat_x, lat_y);
    double y = getLocalPointY(point, lat_x, lat_y);
    double z = getLocalPointZ(point, lat_z);

    // HexLatticePrism hex_surface(point_x, point_y, _width_r, _num_r, _orientation);
    ZPlane zplane(0.0);

    // Bools indicating if point is on each surface
    bool on_min_beta, on_max_beta, on_min_gamma, on_max_gamma, on_min_delta, on_max_delta, on_min_z, on_max_z;

    double dist_beta = std::numeric_limits<double>::infinity();
    double dist_gamma = std::numeric_limits<double>::infinity();
    double dist_delta = std::numeric_limits<double>::infinity();

    // Get the min distance in the beta direction
    if (fabs(u_beta) > FLT_EPSILON) {
      double plane = std::copysign(0.5 * _width_r, u_beta);

      // (x,y) \dot beta
      double l_beta;
      if (_orientation == Orientation::y) {
        l_beta = x * std::sqrt(3.0) / 2.0  + y / 2.0;
      } else {
        l_beta = x;
      }

      dist_beta = (plane - l_beta) / u_beta;
      log::finfo("u_beta:%lf ,dist beta:%lf", u_beta, dist_beta);

      if(u_beta > FLT_EPSILON){
        if(fabs(dist_beta) < FLT_EPSILON)
          on_max_beta = true;
        else
          on_max_beta = false;
      }
      else{
        if(fabs(dist_beta) < FLT_EPSILON)
          on_min_beta = true;
        else
          on_min_beta = false;
      }
    }

    // Get the min distance in the gamma direction
    if (fabs(u_gamma) > FLT_EPSILON) {
      double plane = std::copysign(0.5 * _width_r, u_gamma);

      // (x,y) \dot gamma
      double l_gamma;
      if (_orientation == Orientation::y) {
        l_gamma = x * std::sqrt(3.0) / 2.0  - y / 2.0;
      } else {
        l_gamma = x / 2.0 - y * std::sqrt(3.0) / 2.0;
      }

      dist_gamma = (plane - l_gamma) / u_gamma;
      log::finfo("u_gamma:%lf ,dist gamma:%lf", u_gamma, dist_gamma);

      if(u_gamma > FLT_EPSILON){
        if(fabs(dist_gamma) < FLT_EPSILON)
          on_max_gamma = true;
        else
          on_max_gamma = false;
      }
      else{
        if(fabs(dist_gamma) < FLT_EPSILON)
          on_min_gamma = true;
        else
          on_min_gamma = false;
      }
    }

    // Get the min distance in the delta direction
    if (fabs(u_delta) > FLT_EPSILON) {
      double plane = std::copysign(0.5 * _width_r, u_delta);

      // (x,y) \dot delta
      double l_delta;
      if (_orientation == Orientation::y) {
        l_delta = y;
      } else {
        l_delta = x / 2.0 + y * std::sqrt(3.0) / 2.0;
      }

      dist_delta = (plane - l_delta) / u_delta;
      log::finfo("u_delta:%lf ,dist delta:%lf", u_delta, dist_delta);

      if(u_delta > FLT_EPSILON){
        if(fabs(dist_delta) < FLT_EPSILON)
          on_max_delta = true;
        else
          on_max_delta = false;
      }
      else{
        if(fabs(dist_delta) < FLT_EPSILON)
          on_min_delta = true;
        else
          on_min_delta = false;
      }
    }

    // Check if point is on Z_MIN boundary
    on_min_z = false;
    if (is3D()) {
      zplane.setZ(_accumulate_z[lat_z] + getMinZ());
      on_min_z = zplane.isPointOnSurface(point);
    }

    // Check if point is on Z_MAX boundary
    on_max_z = false;
    if (is3D()) {
      zplane.setZ(_accumulate_z[lat_z+1] + getMinZ());
      on_max_z = zplane.isPointOnSurface(point);
    }

    if (on_min_beta) {
      if (on_min_gamma) {
        if (on_min_z)
          surface = HEX_SURFACE_BETA_MIN_GAMMA_MIN_Z_MIN;
        else if (on_max_z)
          surface = HEX_SURFACE_BETA_MIN_GAMMA_MIN_Z_MAX;
        else
          surface = HEX_SURFACE_BETA_MIN_GAMMA_MIN;
      }
      else if (on_min_delta) {
        if (on_min_z)
          surface = HEX_SURFACE_BETA_MIN_DELTA_MIN_Z_MIN;
        else if (on_max_z)
          surface = HEX_SURFACE_BETA_MIN_DELTA_MIN_Z_MAX;
        else
          surface = HEX_SURFACE_BETA_MIN_DELTA_MIN;
      }
      else {
        if (on_min_z)
          surface = HEX_SURFACE_BETA_MIN_Z_MIN;
        else if (on_max_z)
          surface = HEX_SURFACE_BETA_MIN_Z_MAX;
        else
          surface = HEX_SURFACE_BETA_MIN;
      }
    }
    else if (on_max_beta) {
      if (on_max_gamma) {
        if (on_min_z)
          surface = HEX_SURFACE_BETA_MAX_GAMMA_MAX_Z_MIN;
        else if (on_max_z)
          surface = HEX_SURFACE_BETA_MAX_GAMMA_MAX_Z_MAX;
        else
          surface = HEX_SURFACE_BETA_MAX_GAMMA_MAX;
      }
      else if (on_max_delta) {
        if (on_min_z)
          surface = HEX_SURFACE_BETA_MAX_DELTA_MAX_Z_MIN;
        else if (on_max_z)
          surface = HEX_SURFACE_BETA_MAX_DELTA_MAX_Z_MAX;
        else
          surface = HEX_SURFACE_BETA_MAX_DELTA_MAX;
      }
      else {
        if (on_min_z)
          surface = HEX_SURFACE_BETA_MIN_Z_MIN;
        else if (on_max_z)
          surface = HEX_SURFACE_BETA_MIN_Z_MAX;
        else
          surface = HEX_SURFACE_BETA_MIN;
      }
    }
    else if (on_min_gamma) {
      if (on_max_delta) {
        if (on_min_z)
          surface = HEX_SURFACE_GAMMA_MIN_DELTA_MAX_Z_MIN;
        else if (on_max_z)
          surface = HEX_SURFACE_GAMMA_MIN_DELTA_MAX_Z_MAX;
        else
          surface = HEX_SURFACE_GAMMA_MIN_DELTA_MAX;
      }
      else {
        if (on_min_z)
          surface = HEX_SURFACE_GAMMA_MIN_Z_MIN;
        else if (on_max_z)
          surface = HEX_SURFACE_GAMMA_MIN_Z_MAX;
        else
          surface = HEX_SURFACE_GAMMA_MIN;
      }
    }
    else if (on_max_gamma) {
      if (on_min_delta) {
        if (on_min_z)
          surface = HEX_SURFACE_GAMMA_MAX_DELTA_MIN_Z_MIN;
        else if (on_max_z)
          surface = HEX_SURFACE_GAMMA_MAX_DELTA_MIN_Z_MAX;
        else
          surface = HEX_SURFACE_GAMMA_MAX_DELTA_MIN;
      }
      else if (on_min_z)
        surface = HEX_SURFACE_GAMMA_MAX_Z_MIN;
      else if (on_max_z)
        surface = HEX_SURFACE_GAMMA_MAX_Z_MAX;
      else
        surface = HEX_SURFACE_GAMMA_MAX;
    }
    else if (on_min_delta) {
      if (on_min_z)
        surface = HEX_SURFACE_DELTA_MIN_Z_MIN;
      else if (on_max_z)
        surface = HEX_SURFACE_DELTA_MIN_Z_MAX;
      else
        surface = HEX_SURFACE_DELTA_MIN;
    }
    else if (on_max_delta) {
      if (on_min_z)
        surface = HEX_SURFACE_DELTA_MAX_Z_MIN;
      else if (on_max_z)
        surface = HEX_SURFACE_DELTA_MAX_Z_MAX;
      else
        surface = HEX_SURFACE_DELTA_MAX;
    }
    else if (on_min_z)
      surface = HEX_SURFACE_Z_MIN;
    else if (on_max_z)
      surface = HEX_SURFACE_Z_MAX;

    if (surface != -1)
      surface = HEX_NUM_SURFACES * cell + surface;

    log::finfo("surface: %d", surface);
    return surface;
  }*/

  int HexLattice::getLatticeSurface(int cell, Point *point, double azim, double polar)
  {

    int surface = -1;

    // Get unit vector components
    double u_xy = sin(polar);
    double u_x = u_xy * cos(azim);
    double u_y = u_xy * sin(azim);
    double u_z = cos(polar);

    // Compute the dot production
    double u_beta, u_gamma, u_delta;
    if (_orientation == Orientation::y)
    {
      u_beta = u_x * std::sqrt(3.0) / 2.0 + u_y / 2.0;
      u_gamma = u_x * std::sqrt(3.0) / 2.0 - u_y / 2.0;
      u_delta = u_y;
    }
    else
    {
      u_beta = u_x;
      u_gamma = u_x / 2.0 - u_y * std::sqrt(3.0) / 2.0;
      u_delta = u_x / 2.0 + u_y * std::sqrt(3.0) / 2.0;
    }

    // Compute the x, y, and z indices for the Lattice cell this point is in
    auto lat_xy = getLatXY(point);
    int lat_x = lat_xy[0];
    int lat_y = lat_xy[1];
    int lat_z = getLatZ(point);

    // Get local coordinates of the point
    double x = getLocalPointX(point, lat_x, lat_y);
    double y = getLocalPointY(point, lat_x, lat_y);
    double z = getLocalPointZ(point, lat_z);

    // HexLatticePrism hex_surface(point_x, point_y, _width_r, _num_r, _orientation);
    ZPlane zplane(0.0);

    double dist_beta = std::numeric_limits<double>::infinity();
    double dist_gamma = std::numeric_limits<double>::infinity();
    double dist_delta = std::numeric_limits<double>::infinity();
    double dist = std::numeric_limits<double>::infinity();

    // Bools indicating if point is on each surface
    bool on_min_beta = false, on_max_beta = false,
         on_min_gamma = false, on_max_gamma = false,
         on_min_delta = false, on_max_delta = false,
         on_min_z = false, on_max_z = false;

    // Get the min distance in the beta direction
    if (fabs(u_beta) > FLT_EPSILON)
    {
      double plane = std::copysign(0.5 * _width_r, u_beta);

      // (x,y) \dot beta
      double l_beta;
      if (_orientation == Orientation::y)
      {
        l_beta = x * std::sqrt(3.0) / 2.0 + y / 2.0;
      }
      else
      {
        l_beta = x;
      }

      dist_beta = (plane - l_beta) / u_beta;
      dist = dist_beta;
    }

    // Get the min distance in the gamma direction
    if (fabs(u_gamma) > FLT_EPSILON)
    {
      double plane = std::copysign(0.5 * _width_r, u_gamma);

      // (x,y) \dot gamma
      double l_gamma;
      if (_orientation == Orientation::y)
      {
        l_gamma = x * std::sqrt(3.0) / 2.0 - y / 2.0;
      }
      else
      {
        l_gamma = x / 2.0 - y * std::sqrt(3.0) / 2.0;
      }

      dist_gamma = (plane - l_gamma) / u_gamma;
      dist = std::min(dist, dist_gamma);
    }

    // Get the min distance in the delta direction
    if (fabs(u_delta) > FLT_EPSILON)
    {
      double plane = std::copysign(0.5 * _width_r, u_delta);

      // (x,y) \dot delta
      double l_delta;
      if (_orientation == Orientation::y)
      {
        l_delta = y;
      }
      else
      {
        l_delta = x / 2.0 + y * std::sqrt(3.0) / 2.0;
      }

      dist_delta = (plane - l_delta) / u_delta;
      dist = std::min(dist, dist_delta);
    }

    // Check if point is on Z_MIN boundary
    if (is3D())
    {
      zplane.setZ(_accumulate_z[lat_z] + getMinZ());
      on_min_z = zplane.isPointOnSurface(point);
    }

    // Check if point is on Z_MAX boundary
    if (is3D())
    {
      zplane.setZ(_accumulate_z[lat_z + 1] + getMinZ());
      on_max_z = zplane.isPointOnSurface(point);
    }

    if (dist == dist_beta)
    {
      if (u_beta > 0)
        on_max_beta = true;
      else
        on_min_beta = true;
    }
    if (dist == dist_gamma)
    {
      if (u_gamma > 0)
        on_max_gamma = true;
      else
        on_min_gamma = true;
    }
    if (dist == dist_delta)
    {
      if (u_delta > 0)
        on_max_delta = true;
      else
        on_min_delta = true;
    }

    if (on_min_beta)
    {
      if (on_min_gamma)
      {
        if (on_min_z)
          surface = HEX_SURFACE_BETA_MIN_GAMMA_MIN_Z_MIN;
        else if (on_max_z)
          surface = HEX_SURFACE_BETA_MIN_GAMMA_MIN_Z_MAX;
        else
          surface = HEX_SURFACE_BETA_MIN_GAMMA_MIN;
      }
      else if (on_min_delta)
      {
        if (on_min_z)
          surface = HEX_SURFACE_BETA_MIN_DELTA_MIN_Z_MIN;
        else if (on_max_z)
          surface = HEX_SURFACE_BETA_MIN_DELTA_MIN_Z_MAX;
        else
          surface = HEX_SURFACE_BETA_MIN_DELTA_MIN;
      }
      else
      {
        if (on_min_z)
          surface = HEX_SURFACE_BETA_MIN_Z_MIN;
        else if (on_max_z)
          surface = HEX_SURFACE_BETA_MIN_Z_MAX;
        else
          surface = HEX_SURFACE_BETA_MIN;
      }
    }
    else if (on_max_beta)
    {
      if (on_max_gamma)
      {
        if (on_min_z)
          surface = HEX_SURFACE_BETA_MAX_GAMMA_MAX_Z_MIN;
        else if (on_max_z)
          surface = HEX_SURFACE_BETA_MAX_GAMMA_MAX_Z_MAX;
        else
          surface = HEX_SURFACE_BETA_MAX_GAMMA_MAX;
      }
      else if (on_max_delta)
      {
        if (on_min_z)
          surface = HEX_SURFACE_BETA_MAX_DELTA_MAX_Z_MIN;
        else if (on_max_z)
          surface = HEX_SURFACE_BETA_MAX_DELTA_MAX_Z_MAX;
        else
          surface = HEX_SURFACE_BETA_MAX_DELTA_MAX;
      }
      else
      {
        if (on_min_z)
          surface = HEX_SURFACE_BETA_MIN_Z_MIN;
        else if (on_max_z)
          surface = HEX_SURFACE_BETA_MIN_Z_MAX;
        else
          surface = HEX_SURFACE_BETA_MIN;
      }
    }
    else if (on_min_gamma)
    {
      if (on_max_delta)
      {
        if (on_min_z)
          surface = HEX_SURFACE_GAMMA_MIN_DELTA_MAX_Z_MIN;
        else if (on_max_z)
          surface = HEX_SURFACE_GAMMA_MIN_DELTA_MAX_Z_MAX;
        else
          surface = HEX_SURFACE_GAMMA_MIN_DELTA_MAX;
      }
      else
      {
        if (on_min_z)
          surface = HEX_SURFACE_GAMMA_MIN_Z_MIN;
        else if (on_max_z)
          surface = HEX_SURFACE_GAMMA_MIN_Z_MAX;
        else
          surface = HEX_SURFACE_GAMMA_MIN;
      }
    }
    else if (on_max_gamma)
    {
      if (on_min_delta)
      {
        if (on_min_z)
          surface = HEX_SURFACE_GAMMA_MAX_DELTA_MIN_Z_MIN;
        else if (on_max_z)
          surface = HEX_SURFACE_GAMMA_MAX_DELTA_MIN_Z_MAX;
        else
          surface = HEX_SURFACE_GAMMA_MAX_DELTA_MIN;
      }
      else if (on_min_z)
        surface = HEX_SURFACE_GAMMA_MAX_Z_MIN;
      else if (on_max_z)
        surface = HEX_SURFACE_GAMMA_MAX_Z_MAX;
      else
        surface = HEX_SURFACE_GAMMA_MAX;
    }
    else if (on_min_delta)
    {
      if (on_min_z)
        surface = HEX_SURFACE_DELTA_MIN_Z_MIN;
      else if (on_max_z)
        surface = HEX_SURFACE_DELTA_MIN_Z_MAX;
      else
        surface = HEX_SURFACE_DELTA_MIN;
    }
    else if (on_max_delta)
    {
      if (on_min_z)
        surface = HEX_SURFACE_DELTA_MAX_Z_MIN;
      else if (on_max_z)
        surface = HEX_SURFACE_DELTA_MAX_Z_MAX;
      else
        surface = HEX_SURFACE_DELTA_MAX;
    }
    else if (on_min_z)
      surface = HEX_SURFACE_Z_MIN;
    else if (on_max_z)
      surface = HEX_SURFACE_Z_MAX;

    if (surface != -1)
      surface = HEX_NUM_SURFACES * cell + surface;

    return surface;
  }

  int HexLattice::getLatticeSurfaceOTF(int cell, double z, int surface_2D)
  {

    /* Determine min and max z boundaries of the cell */
    // int lat_z = getLatZ(point);
    int num_x = _num_r * 2 - 1;
    int num_y = _num_r * 2 - 1;
    int lat_z = cell / (num_x * num_y);
    double z_min = _accumulate_z[lat_z] + getMinZ();
    double z_max = z_min + _widths_z[lat_z];

    /* Check for z-surface crossing on 2D surface */
    if (surface_2D % HEX_NUM_SURFACES > 13)
      log::ferror("Found a z-surface crossing on a 2D segment");

    /* Check min z boundary for crossing */
    if (fabs(z_min - z) < TINY_MOVE)
    {
      int surface;
      switch (surface_2D % HEX_NUM_SURFACES)
      {
      case HEX_SURFACE_BETA_MIN:
        surface = HEX_SURFACE_BETA_MIN_Z_MIN;
        break;
      case HEX_SURFACE_BETA_MAX:
        surface = HEX_SURFACE_BETA_MAX_Z_MIN;
        break;
      case HEX_SURFACE_GAMMA_MIN:
        surface = HEX_SURFACE_GAMMA_MIN_Z_MIN;
        break;
      case HEX_SURFACE_GAMMA_MAX:
        surface = HEX_SURFACE_GAMMA_MAX_Z_MIN;
        break;
      case HEX_SURFACE_DELTA_MIN:
        surface = HEX_SURFACE_DELTA_MIN_Z_MIN;
        break;
      case HEX_SURFACE_DELTA_MAX:
        surface = HEX_SURFACE_DELTA_MAX_Z_MIN;
        break;
      case HEX_SURFACE_BETA_MIN_GAMMA_MIN:
        surface = HEX_SURFACE_BETA_MIN_GAMMA_MIN_Z_MIN;
        break;
      case HEX_SURFACE_BETA_MIN_DELTA_MIN:
        surface = HEX_SURFACE_BETA_MIN_DELTA_MIN_Z_MIN;
        break;
      case HEX_SURFACE_BETA_MAX_GAMMA_MAX:
        surface = HEX_SURFACE_BETA_MAX_GAMMA_MAX_Z_MIN;
        break;
      case HEX_SURFACE_BETA_MAX_DELTA_MAX:
        surface = HEX_SURFACE_BETA_MAX_DELTA_MAX_Z_MIN;
        break;
      case HEX_SURFACE_GAMMA_MIN_DELTA_MAX:
        surface = HEX_SURFACE_GAMMA_MIN_DELTA_MAX_Z_MIN;
        break;
      case HEX_SURFACE_GAMMA_MAX_DELTA_MIN:
        surface = HEX_SURFACE_GAMMA_MAX_DELTA_MIN_Z_MIN;
        break;
      default:
        surface = HEX_SURFACE_Z_MIN;
      }
      return HEX_NUM_SURFACES * cell + surface;
    }

    /* Check max z boundary for crossing */
    if (fabs(z_max - z) < TINY_MOVE)
    {
      int surface;
      switch (surface_2D % HEX_NUM_SURFACES)
      {
      case HEX_SURFACE_BETA_MIN:
        surface = HEX_SURFACE_BETA_MIN_Z_MAX;
        break;
      case HEX_SURFACE_BETA_MAX:
        surface = HEX_SURFACE_BETA_MAX_Z_MAX;
        break;
      case HEX_SURFACE_GAMMA_MIN:
        surface = HEX_SURFACE_GAMMA_MIN_Z_MAX;
        break;
      case HEX_SURFACE_GAMMA_MAX:
        surface = HEX_SURFACE_GAMMA_MAX_Z_MAX;
        break;
      case HEX_SURFACE_DELTA_MIN:
        surface = HEX_SURFACE_DELTA_MIN_Z_MAX;
        break;
      case HEX_SURFACE_DELTA_MAX:
        surface = HEX_SURFACE_DELTA_MAX_Z_MAX;
        break;
      case HEX_SURFACE_BETA_MIN_GAMMA_MIN:
        surface = HEX_SURFACE_BETA_MIN_GAMMA_MIN_Z_MAX;
        break;
      case HEX_SURFACE_BETA_MIN_DELTA_MIN:
        surface = HEX_SURFACE_BETA_MIN_DELTA_MIN_Z_MAX;
        break;
      case HEX_SURFACE_BETA_MAX_GAMMA_MAX:
        surface = HEX_SURFACE_BETA_MAX_GAMMA_MAX_Z_MAX;
        break;
      case HEX_SURFACE_BETA_MAX_DELTA_MAX:
        surface = HEX_SURFACE_BETA_MAX_DELTA_MAX_Z_MAX;
        break;
      case HEX_SURFACE_GAMMA_MIN_DELTA_MAX:
        surface = HEX_SURFACE_GAMMA_MIN_DELTA_MAX_Z_MAX;
        break;
      case HEX_SURFACE_GAMMA_MAX_DELTA_MIN:
        surface = HEX_SURFACE_GAMMA_MAX_DELTA_MIN_Z_MAX;
        break;
      default:
        surface = HEX_SURFACE_Z_MAX;
      }
      return HEX_NUM_SURFACES * cell + surface;
    }

    /* If no axial crossing, return the 2D surface */
    if (surface_2D == -1)
      return surface_2D;
    else
      return HEX_NUM_SURFACES * cell + surface_2D % HEX_NUM_SURFACES;
  }

  /// \brief Finds the Lattice cell x and y indices that a point lies in.
  /// \param point a pointer to a point being evaluated.
  /// \return the Lattice cell x and y indices.
  std::array<int, 2> HexLattice::getLatXY(Point *point)
  {
    // Offset the xyz by the lattice center.
    double lx = point->getX() - _offset.getX();
    double ly = point->getY() - _offset.getY();

    int ix, iy;
    if (_orientation == Orientation::y)
    {
      // Convert coordinates into skewed bases.  The (x, alpha) basis is used to
      // find the index of the global coordinates to within 4 cells.
      double alpha = ly - lx / std::sqrt(3.0);
      ix = std::floor(lx / (0.5 * std::sqrt(3.0) * _width_r));
      iy = std::floor(alpha / _width_r);
    }
    else
    {
      // Convert coordinates into skewed bases.  The (alpha, y) basis is used to
      // find the index of the global coordinates to within 4 cells.
      double alpha = ly - lx * std::sqrt(3.0);
      ix = std::floor(-alpha / (std::sqrt(3.0) * _width_r));
      iy = std::floor(ly / (0.5 * std::sqrt(3.0) * _width_r));
    }

    // Add offset to indices (the center cell is (ix, iy) = (0, 0) but
    // the array is offset so that the indices never go below 0).
    ix += _num_r - 1;
    iy += _num_r - 1;

    // Calculate the (squared) distance between the point and the centers of
    // the four possible cells.  Regular hexagonal tiles form a Voronoi
    // tessellation so the xyz should be in the hexagonal cell that it is closest
    // to the center of.  This method is used over a method that uses the
    // remainders of the floor divisions above because it provides better finite
    // precision performance.  Squared distances are used because they are more
    // computationally efficient than normal distances.

    // COINCIDENCE CHECK
    // if a distance to center, d, is within the coincidence tolerance of the
    // current minimum distance, d_min, the point is on an edge or vertex.
    int dx = 0;
    int dy = 0;
    double d_min = std::numeric_limits<double>::infinity();
    for (int i = 0; i < 2; i++)
    {
      for (int j = 0; j < 2; j++)
      {
        // get local coordinates
        double loc_x = getLocalPointX(point, ix + j, iy + i);
        double loc_y = getLocalPointY(point, ix + j, iy + i);

        // calculate distance
        double d = loc_x * loc_x + loc_y * loc_y;

        // skip if the point is on an edge or vertex
        if (d < d_min && !coincident(d, d_min))
        {
          d_min = d;
          dx = j;
          dy = i;
        }
      }
    }

    // update outgoing indices
    ix += dx;
    iy += dy;

    return {ix, iy};
  }

  /// \brief Set universes for a HexLattice
  /// \param num_z the number of lattice cells along z-axis
  /// \param num_r the number of radial lattice cells
  /// \paran universes an array of pointers to universes
  void HexLattice::setUniverses(int num_z, int num_r, Universe **universes)
  {
    setNumR(num_r);
    setNumZ(num_z);

    clearUniverses();
    // Redundant spaces for invalid lattice cells
    _universes.resize(num_z * (2 * num_r - 1) * (2 * num_r - 1), nullptr);

    if (_orientation == Orientation::y)
    {
      setUniversesY(num_z, num_r, universes);
    }
    else
    {
      setUniversesX(num_z, num_r, universes);
    }

    // Compute accumulated arrays
    computeSizes();
  }

  /// \brief Set universes for a y-orientated lattice
  /// \details This algorithm comes from OpenMC.
  ///          Note that we support non-uniform z-layout for HexLattice.
  void HexLattice::setUniversesY(int num_z, int num_r, Universe **universes)
  {

    int in_idx = 0;
    for (int m = 0; m < num_z; m++)
    {
      // Initialize lattice indecies.
      int i_x = 1;
      int i_a = num_r - 1;

      // Map upper triangular region of hexagonal lattice which is found in the
      // first num_r-1 rows of the input.
      for (int k = 0; k < num_r - 1; k++)
      {
        // Walk the index to lower-left neighbor of last row start.
        i_x -= 1;

        // Iterate over the input columns.
        for (int j = 0; j < k + 1; j++)
        {
          int idx = getLatticeCell(i_x + num_r - 1, i_a + num_r - 1, num_z - 1 - m);
          _universes[idx] = universes[in_idx++];

          // Walk the index to the right neighbor (which is not adjacent).
          i_x += 2;
          i_a -= 1;
        }

        // Return the lattice index to the start of the current row.
        i_x -= 2 * (k + 1);
        i_a += (k + 1);
      }

      // Map the middle square region of the hexagonal lattice which is found in
      // the next 2*num_r-1 rows of the input.
      for (int k = 0; k < 2 * num_r - 1; k++)
      {
        if ((k % 2) == 0)
        {
          // Walk the index to the lower-left neighbor of the last row start.
          i_x -= 1;
        }
        else
        {
          // Walk the index to the lower-right neighbor of the last row start.
          i_x += 1;
          i_a -= 1;
        }

        // Iterate over the input columns.
        for (int j = 0; j < num_r - (k % 2); j++)
        {
          int idx = getLatticeCell(i_x + num_r - 1, i_a + num_r - 1, num_z - 1 - m);
          _universes[idx] = universes[in_idx++];

          // Walk the index to the right neighbor (which is not adjacent).
          i_x += 2;
          i_a -= 1;
        }

        // Return the lattice index to the start of the current row.
        i_x -= 2 * (num_r - (k % 2));
        i_a += num_r - (k % 2);
      }

      // Map the lower triangular region of the hexagonal lattice.
      for (int k = 0; k < num_r - 1; k++)
      {
        // Walk the index to the lower-right neighbor of the last row start.
        i_x += 1;
        i_a -= 1;

        // Iterate over the input columns.
        for (int j = 0; j < num_r - k - 1; j++)
        {
          int idx = getLatticeCell(i_x + num_r - 1, i_a + num_r - 1, num_z - 1 - m);
          _universes[idx] = universes[in_idx++];

          // Walk the index to the right neighbor (which is not adjacent).
          i_x += 2;
          i_a -= 1;
        }

        // Return lattice index to start of current row.
        i_x -= 2 * (num_r - k - 1);
        i_a += num_r - k - 1;
      }
    }
  }

  /// \brief Set universes for an x-orientated lattice
  /// \details This algorithm comes from OpenMC.
  void HexLattice::setUniversesX(int num_z, int num_r, Universe **universes)
  {

    int in_idx = 0;
    for (int m = 0; m < num_z; m++)
    {
      // Initialize lattice indecies.
      int i_a = -(num_r - 1);
      int i_y = num_r - 1;

      // Map upper region of hexagonal lattice which is found in the
      // first num_r-1 rows of the input.
      for (int k = 0; k < num_r - 1; k++)
      {

        // Iterate over the input columns.
        for (int j = 0; j < k + num_r; j++)
        {
          int idx = getLatticeCell(i_a + num_r - 1, i_y + num_r - 1, num_z - 1 - m);
          _universes[idx] = universes[in_idx++];

          // Move to the next right neighbour cell
          i_a += 1;
        }

        // Return the lattice index to the start of the current row.
        i_a = -(num_r - 1);
        i_y -= 1;
      }

      // Map the lower region from the centerline of cart to down side
      for (int k = 0; k < num_r; k++)
      {
        // Walk the index to the lower-right neighbor of the last row start.
        i_a = -(num_r - 1) + k;

        // Iterate over the input columns.
        for (int j = 0; j < 2 * num_r - k - 1; j++)
        {
          int idx = getLatticeCell(i_a + num_r - 1, i_y + num_r - 1, num_z - 1 - m);
          _universes[idx] = universes[in_idx++];

          // Move to the next right neighbour cell
          i_a += 1;
        }

        // Return lattice index to start of current row.
        i_y -= 1;
      }
    }
  }

  /**
   * @brief Set widths of non-uniform meshes in x y z directions.
   * @param width_r pitch on the x-y plane
   * @param widths_z z-direction widths of non-uniform meshes
   */
  void HexLattice::setWidths(double width_r, DoubleVec widths_z)
  {
    _non_uniform = true;
    _width_r = width_r;
    _widths_z = widths_z;
  }

  /// \brief Set _widths_z for uniform case, compute accumulate variables.
  void HexLattice::computeSizes()
  {
    if (_non_uniform)
    {
      if (_widths_z.size() - _num_z)
        log::ferror("The sizes of non-uniform mesh widths are not consistent"
                    " with the sizes of filling Universes into Lattice");
    }
    else
    {
      _widths_z.resize(_num_z, _width_z);
    }

    /* Compute the accumulated lengths along each axis */
    _accumulate_z.resize(_num_z + 1, 0.0);

    for (int i = 0; i < _num_z; i++)
      _accumulate_z[i + 1] = _accumulate_z[i] + _widths_z[i];
  }

  /// \brief Returns the center point of a lattice cell with respect to the center
  ///        of the lattice
  /// \details This method is designed for visualization. The most important part
  ///          is the change of basis. It actually consists of a shear and a scaling
  ///          transformations.
  ///
  ///          y-orientation:
  ///          A = ┍ r*cos(pi/6)  0 ┑ = ┍ r  0 ┑┍ cos(pi/6)  0 ┑
  ///              ┕ r*sin(pi/6)  r ┙   ┕ 0  r ┙┕ sin(pi/6)  1 ┙
  ///
  ///          x-orientation:
  ///          A = ┍ r  r*sin(pi/6) ┑ = ┍ r  0 ┑┍ 1  sin(pi/6) ┑
  ///              ┕ 0  r*cos(pi/6) ┙   ┕ 0  r ┙┕ 0  cos(pi/6) ┙
  Point HexLattice::getLatticeCellCenterLocal(int lat_x, int lat_y, int lat_z)
  {

    double r = _width_r;
    double cos30 = std::cos(M_PI / 6);
    double sin30 = 0.5;
    double x, y;

    // Apply the change of basis to a lattice cell
    if (_orientation == Orientation::y)
    {
      x = r * cos30 * lat_x;
      y = r * (sin30 * lat_x + lat_y);
    }
    else
    {
      x = r * (lat_x + sin30 * lat_y);
      y = r * cos30 * lat_y;
    }

    // Compute z coordinate
    double z = _accumulate_z[lat_z] + _widths_z[lat_z] / 2;

    return {x, y, z};
  }

  /// \brief Returns the center point of a lattice cell with respect to the center
  ///        of the center of parent universe
  Point HexLattice::getLatticeCellCenter(int lat_x, int lat_y, int lat_z)
  {

    // Lattice cell center, local point
    Point center = getLatticeCellCenterLocal(lat_x, lat_y, lat_z);
    // Lattice center, local point
    Point lat_center = getLatticeCellCenterLocal(_num_r - 1,
                                                 _num_r - 1,
                                                 lat_z);
    // Corrects the z-coordinate
    lat_center.setZ(_accumulate_z[_num_z] / 2.);

    return _offset + center - lat_center;
  }

  /// \brief Returns coordinates of the 6 or 12 vertices of a lattice cell
  /// \details First we got one vertex, and then we rotate it to get the rest
  ///          of the vertices.
  std::vector<Point> HexLattice::getLatticeCellVertices(int index)
  {

    double sin60 = std::sin(M_PI / 3);
    double cos60 = 1 / 2.;
    double tan60 = std::tan(M_PI / 3);

    std::vector<Point> vertices;

    // Gets the lattice cell center
    Point center = getLatticeCellCenter(index);

    // Defines a radius vector
    double x, y, z;
    if (_orientation == Orientation::y)
    {
      x = _width_r / (2. * sin60);
      y = 0;
    }
    else
    {
      x = _width_r / 2.;
      y = x / tan60;
    }

    // FIXME: support 2-D lattice
    int lat_z = index / ((2 * _num_r - 1) * (2 * _num_r - 1));
    z = _widths_z[lat_z] / 2.;

    for (auto z_sign : std::vector<int>{-1, 1})
    {
      for (int i = 0; i < 6; ++i)
      {
        vertices.push_back({x + center.getX(),
                            y + center.getY(),
                            z_sign * z + center.getZ()});
        // Rotation
        auto u = cos60 * x - sin60 * y;
        auto v = sin60 * x + cos60 * y;
        x = u;
        y = v;
      }
    }

    return vertices;
  }

  std::string HexLattice::toString()
  {

    std::stringstream string;

    string << "Lattice ID = " << _id
           << ", name = " << _name
           << ", # cells along r-axis = " << _num_r
           << ", # cells along z = " << _num_z;

    string << ", radial pitch = " << _width_r
           << ", axial  pitch = " << _width_z;

    string << ", type = Hexagon ";

    string << "\n\t\tUniverse IDs within this Lattice: ";

    for (auto it = begin(); it != end(); ++it)
    {
      auto u = *it;
      if (u)
        string << u->getId();
      else
        string << -1;
      string << ", ";

      // FIXME
      // if (it.getPos() % _num_x == 0)
      //  string << "\n\t\t";
    }

    return string.str();
  }

  std::string HexLattice::dimToString()
  {

    std::stringstream string;

    string << "nr = " << _num_r
           << ", nz = " << _num_z;
    return string.str();
  }

} /* namespace antmoc */
