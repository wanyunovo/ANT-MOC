#include "antmoc/LocalCoords.h"
#include "antmoc/Cell.h"
#include "antmoc/constants.h"
#include "antmoc/Lattice.h"
#include "antmoc/log.h"
#include "antmoc/Universe.h"

namespace antmoc
{

  /**
   * @brief 构造函数将 x、y 和 z 坐标和位置设置为坐标。
   * @param x x 坐标
   * @param y y 坐标
   * @param z z 坐标
   * @param first LocalCoords 是否是第一个，它将包含一个
   *所有下一个 LocalCoords 的数组
   */
  LocalCoords::LocalCoords(double x, double y, double z, bool first)
  {
    _coords.setCoords(x, y, z);
    _phi = 0.;
    _polar = M_PI_2;
    _universe = NULL;
    _lattice = NULL;
    _cell = NULL;
    _next = NULL;
    _prev = NULL;
    _version_num = 0;
    if (first)
    {
      // The first local coords will pre-allocate memory to save all subsequent local coords.
      _array_size = LOCAL_COORDS_LEN;
      _next_array = new LocalCoords[LOCAL_COORDS_LEN];
    }
    else
    {
      _array_size = 0;
      _next_array = NULL;
    }
    _position = -1;
  }

  /**
   * @briefDestructor。
   */
  LocalCoords::~LocalCoords()
  {
    prune();
    if (_position == -1)
      deleteArray();
  }

  /**
   * @brief 返回此 LocalCoords 的级别（UNIV 或 LAT）。
   * @return 嵌套的 Universe 级别（UNIV 或 LAT）
   */
  coordType LocalCoords::getType()
  {
    return _type;
  }

  /**
   * @brief 返回此 LocalCoords 所在的 Universe。
   * @回归宇宙
   */
  Universe *LocalCoords::getUniverse() const
  {
    return _universe;
  }

  /**
   * @brief 返回此 LocalCoords 所在的单元格。
   * @返回单元格
   */
  Cell *LocalCoords::getCell() const
  {
    return _cell;
  }

  /**
   * @brief 返回此 LocalCoords 所在的格子。
   * @返回格子
   */
  Lattice *LocalCoords::getLattice() const
  {
    return _lattice;
  }

  /**
   * @brief 返回此网格单元的第一个索引
   *LocalCoords 驻留。
   * @return第一个晶格单元索引
   */
  int LocalCoords::getLatticeX() const
  {
    return _lattice_x;
  }

  /**
   * @brief 返回此网格单元的第二个索引
   *LocalCoords 驻留。
   * @return第二个格子索引
   */
  int LocalCoords::getLatticeY() const
  {
    return _lattice_y;
  }

  /**
   * @brief 返回此网格单元的第三个索引
   *LocalCoords 驻留。
   * @return第三个格子索引
   */
  int LocalCoords::getLatticeZ() const
  {
    return _lattice_z;
  }

  /**
   * @brief 返回此 LocalCoords 位置的 x 坐标。
   * @return 此 LocalCoords 位置的 x 坐标
   */
  double LocalCoords::getX() const
  {
    return _coords.getX();
  }

  /**
   * @brief 返回此 LocalCoords 位置的 y 坐标。
   * @return 此 LocalCoords 位置的 y 坐标
   */
  double LocalCoords::getY() const
  {
    return _coords.getY();
  }

  /**
   * @brief 返回此 LocalCoords 位置的 z 坐标。
   * @return 此 LocalCoords 位置的 z 坐标
   */
  double LocalCoords::getZ() const
  {
    return _coords.getZ();
  }

  /**
   * @brief 返回相对于 x 轴的方向角（以弧度为单位）。
   * @返回以弧度表示的方向角
   */
  double LocalCoords::getPhi() const
  {
    return _phi;
  }

  /**
   * @brief 返回相对于 z 轴的方向角（以弧度为单位）。
   * @返回以弧度表示的方向角
   */
  double LocalCoords::getPolar() const
  {
    return _polar;
  }

  /**
   * @brief 返回指向包含此坐标的 Point 的指针
   *本地坐标。
   * @return 指向包含 x 和 y 坐标的 Point 的指针
   */
  Point *LocalCoords::getPoint()
  {
    return &_coords;
  }

  /**
   * @brief 返回指向下一个较低嵌套 Universe 的 LocalCoord 的指针
   *级别（如果存在）。
   * @return 指向下一个 LocalCoord 的指针
   */
  LocalCoords *LocalCoords::getNext() const
  {
    return _next;
  }

  /**
   * @brief 创建并返回指向下一个 LocalCoords 的指针（嵌套
   *更深）。
   * @param x x 坐标
   * @param y y 坐标
   * @param z z 坐标
   * @return 指向刚刚创建的下一个 LocalCoords 的指针
   * 调用者是“当前这一层”的 LocalCoords 结点 this。想要“往下一层走”，就调用 getNextCreate(下一层的 x, y, z)。
   *
   */
  LocalCoords *LocalCoords::getNextCreate(double x, double y, double z)
  {

    if (_next == NULL) //_next不为空就不会创建了
    {

      if (_position + 1 >= _array_size)
      {
        _next = new LocalCoords(x, y, z, true); // 新的结点自己再带一块新的数组。
        _next->setPrev(this);
      }
      else
      {
        _next = &_next_array[_position + 1];
        _next->setPrev(this);
        _next->setNext(NULL);
        _next->setArrayPosition(_next_array, _position + 1, _array_size);
        _next->getPoint()->setCoords(x, y, z);
        _next->setLattice(NULL);
        _next->setUniverse(NULL);
        _next->setCell(NULL);
        _next->setVersionNum(0);
        /*
        用 _next_array[_position+1] 这个现成的对象作为 _next。设置好它的 prev/next、_position 等信息。给它写入新的 (x,y,z)，并把 Universe/Lattice/Cell 先清空，让后面几何逻辑再填。
        */
      }
    }

    return _next;
  }

  /**
   * @brief 返回 _next_array 中的 LocalCoords 位置。
   * @return 该对象在底层 _next_array 中的位置
   */
  int LocalCoords::getPosition()
  {
    return _position;
  }

  /**
   * @brief 搜索 LocalCoords 对象以检测循环。
   * @details 如果 LocalCoords 表观长度更大，则假定存在循环
   *1000名会员
   */
  void LocalCoords::detectLoop()
  {
    int n = 0;

    LocalCoords *iter = _next;
    while (iter != NULL)
    {
      iter = iter->getNext();
      n++;
      if (n > 1000)
        log::ferror("Infinite loop of coords");
    }
    log::fdebug("The LocalCoords is: %s\n", toString().c_str());
    log::fdebug("The depth of the chain is %d \n", n);
  }

  /**
   * @brief 返回指向下一个更高嵌套 Universe 的 LocalCoord 的指针
   *级别（如果存在）。
   * @return 指向前一个 LocalCoord 的指针
   */
  LocalCoords *LocalCoords::getPrev() const
  {
    return _prev;
  }

  /**
   * @brief 返回 LocalCoords 对象的版本。
   * @details 版本号区分其他匹配的 FSR 密钥
   * @return 版本号
   */
  int LocalCoords::getVersionNum()
  {
    return _version_num;
  }

  /**
   * @brief 设置 LocalCoords 的类型（UNIV 或 LAT）。
   * @param type LocalCoords 的类型（UNIV 或 LAT）
   */
  void LocalCoords::setType(coordType type)
  {
    _type = type;
  }

  /**
   * @brief 设置此 LocalCoords 所在的 Universe。
   * @param宇宙宇宙
   */
  void LocalCoords::setUniverse(Universe *universe)
  {
    _universe = universe;
  }

  /**
   * @brief 设置此 LocalCoords 所在的单元格。
   * @param cell 单元格
   */
  void LocalCoords::setCell(Cell *cell)
  {
    _cell = cell;
  }

  /**
   * @brief 设置此 LocalCoords 所在的晶格。
   * @paramlattice 格子
   */
  void LocalCoords::setLattice(Lattice *lattice)
  {
    _lattice = lattice;
  }

  /**
   * @brief 设置该网格单元的行索引
   *LocalCoords 驻留。
   * @paramlattice_x 行晶格单元索引
   */
  void LocalCoords::setLatticeX(int lattice_x)
  {
    _lattice_x = lattice_x;
  }

  /**
   * @brief 设置该网格单元的列索引
   *LocalCoords 驻留。
   * @paramlattice_y 列晶格单元索引
   */
  void LocalCoords::setLatticeY(int lattice_y)
  {
    _lattice_y = lattice_y;
  }

  /**
   * @brief 设置该网格单元的 z 索引
   *LocalCoords 驻留。
   * @paramlattice_z z 晶格索引
   */
  void LocalCoords::setLatticeZ(int lattice_z)
  {
    _lattice_z = lattice_z;
  }

  /**
   * @brief 设置此 LocalCoords 的 x 坐标。
   * @param x x 坐标
   */
  void LocalCoords::setX(double x)
  {
    _coords.setX(x);
  }

  /**
   * @brief 设置此本地坐标的 y 坐标。
   * @param y y 坐标
   */
  void LocalCoords::setY(double y)
  {
    _coords.setY(y);
  }

  /**
   * @brief 设置此 LocalCoords 的 z 坐标。
   * @param z z 坐标
   */
  void LocalCoords::setZ(double z)
  {
    _coords.setZ(z);
  }

  /**
   * @brief 设置此 LocalCoords 的方位角。
   * @param phi 方位角
   */
  void LocalCoords::setPhi(double phi)
  {
    _phi = phi;
  }

  /**
   * @brief 设置此 LocalCoords 的极角。
   * @param Polar 极角
   */
  void LocalCoords::setPolar(double polar)
  {
    _polar = polar;
  }

  /**
   * @brief 设置指向下一个较低嵌套 Universe 上的 LocalCoords 的指针
   *级别。
   * @param next 指向下一个 LocalCoords 的指针
   */
  void LocalCoords::setNext(LocalCoords *next)
  {
    _next = next;
  }

  /**
   * @brief 将指针设置为下一个更高嵌套的 LocalCoords
   *宇宙级别。
   * @param prev 指向前一个 LocalCoords 的指针
   */
  void LocalCoords::setPrev(LocalCoords *prev)
  {
    _prev = prev;
  }

  /**
   * @brief 设置此 LocalCoords 在 LocalCoords 数组中的位置，
   *指向该数组的指针（下一个 LocalCoords）也是
   *已转移
   * @param array 指向下一个 LocalCoords 数组的指针
   * @param 该 LocalCoords 在所述数组中的位置索引
   * @param array_size 数组的大小
   */
  void LocalCoords::setArrayPosition(LocalCoords *array, int position,
                                     int array_size)
  {
    _next_array = array;
    _position = position;
    _array_size = array_size;
  }

  /**
   * @brief 设置 LocalCoords 对象的版本。
   * @details 版本号区分其他匹配的 FSR 密钥
   * @param version_num 版本号
   */
  void LocalCoords::setVersionNum(int version_num)
  {
    _version_num = version_num;
  }

  /**
   * @brief 查找并返回链表中的最后一个 LocalCoords
   *表示几何体最低层的局部坐标
   *嵌套宇宙。
   * @details 遍历 LocalCoords 的链接列表以查找位于
   *最低嵌套宇宙级别。
   * @return 指向列表中最后一个 LocalCoords 对象的指针
   */
  LocalCoords *LocalCoords::getLowestLevel()
  {
    LocalCoords *curr = this;

    /* 遍历链表 */
    while (curr->getNext() != NULL)
      curr = curr->getNext();

    return curr;
  }

  /**
   * @brief 查找并返回链表中的第一个 LocalCoords
   *表示几何体最高层的局部坐标
   *嵌套宇宙。
   * @details 遍历 LocalCoords 的链接列表以查找位于
   *最高嵌套宇宙级别。
   * @return 指向列表中第一个 LocalCoords 对象的指针
   */
  LocalCoords *LocalCoords::getHighestLevel()
  {
    LocalCoords *curr = this;

    /* 遍历链表 */
    while (curr->getPrev() != NULL)
      curr = curr->getPrev();

    return curr;
  }

  /**
   * @brief 翻译每个 LocalCoords 对象的所有 x,y,z 坐标
   *链接列表。
   * @details 这个方法会遍历整个链表并应用
   *对每个元素的翻译。
   * @param delta_x 我们希望将 x 移动的量
   * @param delta_y 我们希望将 y 移动的量
   * @param delta_z 我们希望将 z 移动的量
   */
  void LocalCoords::adjustCoords(double delta_x, double delta_y, double delta_z)
  {

    /* 沿着链表前进 */
    LocalCoords *curr = this;
    while (curr != NULL)
    {
      curr->setX(curr->getX() + delta_x);
      curr->setY(curr->getY() + delta_y);
      curr->setZ(curr->getZ() + delta_z);
      curr = curr->getNext();
    }

    /* 沿链表反转方向 */
    curr = _prev;
    while (curr != NULL)
    {
      curr->setX(curr->getX() + delta_x);
      curr->setY(curr->getY() + delta_y);
      curr->setZ(curr->getZ() + delta_z);
      curr = curr->getPrev();
    }
    return;
  }

  /**
   * @brief 更新链表中的最后一个元素（最低的元素）
   *嵌套宇宙的级别）与 a 具有相同的坐标
   *给定点。
   * @param point 指向兴趣点的指针
   */
  void LocalCoords::updateMostLocal(Point *point)
  {

    /* This function translates the lowest point in the linked list to a given position,
     *Correspondingly, other points in the linked list should also be translated by the same vector.
     *When using this function, ensure that the point is on the same layer as the lowest-level LocalCoords?
     */
    /* Get the lowest level coordinate */
    LocalCoords *curr = getLowestLevel();

    /* 将坐标平移适当的量 */
    double delta_x = point->getX() - curr->getX();
    double delta_y = point->getY() - curr->getY();
    double delta_z = point->getZ() - curr->getZ();
    adjustCoords(delta_x, delta_y, delta_z);

    return;
  }

  /**
   * @brief 删除当前节点以后(低层)的所有节点
   */
  void LocalCoords::prune()
  {

    LocalCoords *curr = getLowestLevel();
    LocalCoords *next = curr->getPrev();

    /* 迭代链接列表中该坐标下的 LocalCoords */
    while (curr != this)
    {
      next = curr->getPrev();
      if (curr->getPosition() == -1)
        // When curr is a starting local coords, delete the entire linked list headed by it at once
        curr->deleteArray();
      curr = next;
    }

    /* 将链表中的下一个 LocalCoord 设置为 null */
    setNext(NULL);
  }

  /**
   * @brief 删除下一个坐标的基础数组。
   */
  void LocalCoords::deleteArray()
  {
    if (_next_array != NULL)
    {
      delete[] _next_array;
      _next_array = NULL;
    }
  }

  /**
   * @brief
   * @details 将方法调用节点下面的所有节点(包含调用节点)拷贝到coords节点 。完成后，以当前节点和coords节点为首的两个链表中的节点值相同
   * @param coords 一个指向 LocalCoords 的指针
   *
   */
  void LocalCoords::copyCoords(LocalCoords *coords)
  {
    /* Copy the current node to the passed coords node,
     *Copy all nodes below the current node to the coords node,
     *After completion, the node values in the two linked lists headed by the current node and the coords node are the same.
     */

    LocalCoords *curr1 = this;
    LocalCoords *curr2 = coords;

    /* 修剪 LocalCoords 链表 */
    curr2->prune();

    /* 迭代此 LocalCoords 链表并创建一个
     * 输入 LocalCoords 的副本 */
    while (curr1 != NULL)
    {
      curr2->setX(curr1->getX());
      curr2->setY(curr1->getY());
      curr2->setZ(curr1->getZ());
      curr2->setUniverse(curr1->getUniverse());

      if (curr1->getType() == UNIV)
      {
        // There is only one cell under a universe
        curr2->setType(UNIV);
        curr2->setCell(curr1->getCell());
      }
      else
      {
        // There are many cells under a lattice, so you cannot just take one cell.
        curr2->setLattice(curr1->getLattice());
        curr2->setLatticeX(curr1->getLatticeX());
        curr2->setLatticeY(curr1->getLatticeY());
        curr2->setLatticeZ(curr1->getLatticeZ());
        curr2->setType(LAT);
      }

      curr1 = curr1->getNext();

      if (curr1 != NULL)
      {
        curr2 = curr2->getNextCreate(0, 0, 0);
      }
    }

    /* 删除旧坐标链表中的任何剩余部分 */
    if (curr2 != NULL)
      curr2->prune();
  }

  /**
   * @brief 将此 LocalCoords 的属性转换为字符数组
   *代表。
   * @return LocalCoord 属性的字符数组
   */
  std::string LocalCoords::toString()
  {

    std::stringstream string;
    LocalCoords *curr = this;

    /* 循环遍历列表中低于此的所有 LocalCoords */
    while (curr != NULL)
    {
      string << "LocalCoords: level = ";

      int univ_id = -1;
      if (curr->getUniverse() != NULL)
        univ_id = curr->getUniverse()->getId();
      int lat_id = -1;
      if (curr->getLattice() != NULL)
        lat_id = curr->getLattice()->getId();
      int cell_id = -1;
      if (curr->getCell() != NULL)
        cell_id = curr->getCell()->getId();

      if (curr->getType() == UNIV)
      {
        string << " UNIVERSE, x = " << curr->getX()
               << ", y = " << curr->getY()
               << ", z = " << curr->getZ();
        string << ", universe = " << univ_id;
        string << ", cell = " << cell_id;
      }
      else if (curr->getType() == LAT)
      {
        string << " LATTICE, x = " << curr->getX()
               << ", y = " << curr->getY()
               << ", z = " << curr->getZ()
               << ", universe = " << univ_id
               << ", lattice = " << lat_id
               << ", lattice_x = " << curr->getLatticeX()
               << ", lattice_y = " << curr->getLatticeY()
               << ", lattice_z = " << curr->getLatticeZ();
      }
      else
      {
        string << " NONE, x = " << curr->getX()
               << ", y = " << curr->getY()
               << ", z = " << curr->getZ()
               << ", universe = " << univ_id
               << ", lattice = " << lat_id
               << ", lattice_x = " << curr->getLatticeX()
               << ", lattice_y = " << curr->getLatticeY()
               << ", lattice_z = " << curr->getLatticeZ();
        string << ", cell = " << cell_id;
      }

      string << ", next:\n";
      curr = curr->getNext();
    }

    return string.str();
  }

} /* 命名空间 Antmoc */
