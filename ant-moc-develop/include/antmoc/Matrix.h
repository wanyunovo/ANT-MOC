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

class Matrix {

private:

  /** A vector of maps representing the matrix */
  std::vector< std::map<int, CMFD_PRECISION> > _LIL; //列表-列表（List of Lists, LIL）表示法

  /** The CSR（Compressed Sparse Row）压缩稀疏行 matrix variables ，举个栗子：有如下稀疏矩阵，CSR 格式的存储方式如下：
   * 1 0 0 0
   * 0 2 0 0
   * 3 0 4 0
   * 0 0 0 5
   * _A: {1, 2, 3, 4, 5}     //row = 2 IA[2] = 2 IA[3] = 4 JA[2] = 0/JA[3] = 2
   * _IA: {0, 1, 2, 4, 5}  最后一个元素值表示非零元素数量
   * _JA: {0, 1, 0, 2, 3}
   * _DIAG: {1, 2, 4, 5}
  */
  CMFD_PRECISION* _A;  //存储矩阵中所有非零元素的值
  int* _IA;   //存储每行第一个非零元素在_A中的起始位置索引,解释一下第四个为啥=4,因为第四行的非零元素为5,而5在_A中的索引位置是4
  int* _JA;   //存储_A中每个非零元素对应的列索引  
  CMFD_PRECISION* _DIAG;  //存储矩阵的对角线元素

  bool _modified;  //是否初始化了Matrix（Matrix构造方法是否被调用）
  int _num_x;
  int _num_y;
  int _num_z;
  int _num_groups;
  int _num_rows;

  /** OpenMP mutual exclusion locks for atomic cell updates */
  omp_lock_t* _cell_locks;

  void convertToCSR();
  void setNumX(int num_x);
  void setNumY(int num_y);
  void setNumZ(int num_z);
  void setNumGroups(int num_groups);

public:
  Matrix(omp_lock_t* cell_locks, int num_x=1, int num_y=1, int num_z=1,
         int num_groups=1, int empty_cells=0);
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
  CMFD_PRECISION* getA();
  int* getIA();
  int* getJA();
  CMFD_PRECISION* getDiag();
  int getNumX();
  int getNumY();
  int getNumZ();
  int getNumGroups();
  int getNumRows();
  int getNNZ();
  omp_lock_t* getCellLocks();

  /* Setter functions */
  void setValue(int cell_from, int group_from, int cell_to, int group_to,
                CMFD_PRECISION val);
};

} /* namespace antmoc */

#endif /* MATRIX_H_ */
