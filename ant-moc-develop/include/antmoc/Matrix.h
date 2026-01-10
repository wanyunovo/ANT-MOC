/**
 * @file Matrix.h
 * @brief A matrix object
 * @date May 5, 2015
 * @author Samuel Shaner, MIT, Course 22 (shaner@mit.edu)
 */

#ifndef MATRIX_H_
#define MATRIX_H_

#include <map>
#include <omp.h>
#include <vector>
#include "antmoc/constants.h"

namespace antmoc
{

  class Matrix
  {

  private:
    /** A vector of maps representing the matrix */
    std::vector<std::map<int, CMFD_PRECISION>> _LIL;
    // 列表-列表（List of Lists, LIL）表示法
    /*
    _LIL 结构示意：
  _LIL[0] = {(col_a, val_a), (col_b, val_b), ...}   // 第0行的非零元素
  _LIL[1] = {(col_c, val_c), ...}                   // 第1行的非零元素
  ...
  _LIL[n-1] = {...}                                  // 第n-1行的非零元素

  _LIL 是一个 vector<map<int, CMFD_PRECISION>>
每一行用一个 map 存储：键是列号，值是矩阵元素值
std::map 自动按键排序，且插入/查找效率为 O(log n)

    */

    /** The CSR（Compressed Sparse Row）是高效存储稀疏矩阵的标准格式，由三个数组组成：_A、_IA 和 _JA
     * - _A：存储矩阵中所有非零元素的值，按行优先顺序排列
     * - _IA：存储每行第一个非零元素在_A中的起始位置索引，注意是在_A数组中的位置，不是在原矩阵中
     * - _JA：存储_A中每个非零元素对应的列索引
     *
     * 举个栗子：有如下稀疏矩阵，CSR 格式的存储方式如下：
     * 1 0 0 0
     * 0 2 0 0
     * 3 0 4 0
     * 0 0 0 5
     * _A: {1, 2, 3, 4, 5}     //row = 2 IA[2] = 2 IA[3] = 4 JA[2] = 0/JA[3] = 2
     * _IA: {0, 1, 2, 4, 5}  最后一个元素值表示非零元素数量
     * _JA: {0, 1, 0, 2, 3}
     * _DIAG: {1, 2, 4, 5}
     *
     * 为什么需要两种格式？、
     * 格式	            优点	                       缺点	              适用阶段
     * LIL	    动态插入/修改 O(log n)	      遍历较慢，内存不连续	      构建阶段
     * CSR	    遍历快，内存连续，适合矩阵运算	      修改困难	           求解阶段
     */
    CMFD_PRECISION *_A;    // 存储矩阵中所有非零元素的值
    int *_IA;              // 存储每行第一个非零元素在_A中的起始位置索引,解释一下第四个为啥=4,因为第四行的非零元素为5,而5在_A中的索引位置是4
    int *_JA;              // 存储_A中每个非零元素对应的列索引
    CMFD_PRECISION *_DIAG; // 存储矩阵的对角线元素

    bool _modified; // 是否初始化了Matrix（Matrix构造方法是否被调用）
    int _num_x;
    int _num_y;
    int _num_z;
    int _num_groups;
    int _num_rows;

    /** OpenMP mutual exclusion locks for atomic cell updates
     * 为每个CMFD 单元分配一个锁，对同一 CMFD 单元的操作互斥
     */
    omp_lock_t *_cell_locks;

    void convertToCSR();
    void setNumX(int num_x);
    void setNumY(int num_y);
    void setNumZ(int num_z);
    void setNumGroups(int num_groups);

  public:
    Matrix(omp_lock_t *cell_locks, int num_x = 1, int num_y = 1, int num_z = 1,
           int num_groups = 1, int empty_cells = 0);
    virtual ~Matrix();

    /* Worker functions */
    void incrementValue(int cell_from, int group_from, int cell_to, int group_to,
                        CMFD_PRECISION val);
    void clear();
    void printString();
    void transpose();

    /* Getter functions */
    CMFD_PRECISION getValue(int cell_from, int group_from, int cell_to,
                            int group_to);
    CMFD_PRECISION *getA();
    int *getIA();
    int *getJA();
    CMFD_PRECISION *getDiag();
    int getNumX();
    int getNumY();
    int getNumZ();
    int getNumGroups();
    int getNumRows();
    int getNNZ();
    omp_lock_t *getCellLocks();

    /* Setter functions */
    void setValue(int cell_from, int group_from, int cell_to, int group_to,
                  CMFD_PRECISION val);
  };

} /* namespace antmoc */

#endif /* MATRIX_H_ */
