/**
 * @file Vector.h
 * @brief A vector object
 * @date May 5, 2015
 * @author Samuel Shaner, MIT, Course 22 (shaner@mit.edu)
 */

#ifndef SRC_VECTOR_H_
#define SRC_VECTOR_H_


#include <map>
#include <omp.h>
#include <vector>
#include "antmoc/pairwise_sum.h"

namespace antmoc
{

class Vector {
 private:
  /** A list of lists representing the vector */
  CMFD_PRECISION* _array;
  int _num_rows;  //Vector容器的大小=(_num_x*_num_y*_num_z-empty_cells)*_num_groups
  int _num_x;
  int _num_y;
  int _num_z;
  int _num_groups;

  /** OpenMP mutual exclusion locks for atomic cell updates */
  omp_lock_t* _cell_locks;

  void setNumX(int num_x);
  void setNumY(int num_y);
  void setNumZ(int num_z);
  void setNumGroups(int num_groups);

 public:
  Vector(omp_lock_t* cell_locks, int num_x = 1, int num_y = 1, int num_z = 1,
         int num_groups = 1, int empty_cells = 0);
  virtual ~Vector();

  /* Worker functions */
  void incrementValue(int cell, int group, CMFD_PRECISION val);
  void incrementValues(int cell, int group_start, int group_end,
                       CMFD_PRECISION* vals);
  void clear();
  void scaleByValue(CMFD_PRECISION val);
  void printString();
  void copyTo(Vector* vector);

  /* Getter functions */
  CMFD_PRECISION getValue(int cell, int group);
  CMFD_PRECISION* getArray();
  int getNumX();
  int getNumY();
  int getNumZ();
  int getNumGroups();
  int getNumRows();
  double getSum();
  omp_lock_t* getCellLocks();

  /* Setter functions */
  void setValue(int cell, int group, CMFD_PRECISION val);
  void setValues(int cell, int group_start, int group_end,
                 CMFD_PRECISION* vals);
  void setAll(CMFD_PRECISION val);
};

} /* namespace antmoc */

#endif  // SRC_VECTOR_H_
