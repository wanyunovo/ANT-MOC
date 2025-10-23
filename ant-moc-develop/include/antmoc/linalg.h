/**
 * @file linalg.h
 * @details This file contains a library of functions for performing linear
 *          algebra operations.
 * @date July 3, 2015
 * @author Samuel Shaner, MIT, Course 22 (shaner@mit.edu)
 */

/* File: linalg.h */

#ifndef LINALG_H_
#define LINALG_H_

#include "antmoc/constants.h"
#include "antmoc/mpi_utils.h"
#include "antmoc/Matrix.h"
#include "antmoc/Vector.h"

#include <map>
#include <vector>
#include <algorithm>

namespace antmoc
{

/**
 * @brief Verbose iteration information for the CMFD eigenvalue solver 
 */
struct ConvergenceData {

  /* The Max. prolongation factor of the CMFD */
  double pf;

  /* The initial residual of the CMFD eigenvalue problem */
  double cmfd_res_1;

  /* The final residual of the CMFD eigenvalue problem */
  double cmfd_res_end;

  /* The linear solver residual of the first CMFD eigenvalue iteration */
  double linear_res_1;

  /* The linear solver residual of the final CMFD eigenvalue iteration 最终CMFD特征值迭代的线性求解器残差*/
  double linear_res_end;

  /* The number of the CMFD eigenvalue iterations */
  int cmfd_iters;

  /* The number of linear iterations for the first CMFD eigenvalue iteration */
  int linear_iters_1;

  /* The number of linear iterations for the final CMFD eigenvalue iteration */
  int linear_iters_end;

  /* Constructor initializes convergence statistics to -1 */
  ConvergenceData() {
    pf = -1;
    cmfd_res_1 = -1;
    cmfd_res_end = -1;
    linear_res_1 = -1;
    linear_res_end = -1;
    cmfd_iters = -1;
    linear_iters_1 = -1;
    linear_iters_end = -1;
  }
};


/**
 * @brief Structure for communication of fluxes between neighbor domains
 * 相邻域之间通量通信的简要结构
 */
struct DomainCommunicator {

  int _num_domains_x;
  int _num_domains_y;
  int _num_domains_z;
  //当前域在x方向上的索引
  int _domain_idx_x;
  int _domain_idx_y;
  int _domain_idx_z;
  int _local_num_x;
  int _local_num_y;
  int _local_num_z;

  /* Sum of the starting CMFD global indexes of a domain, for the color 起始CMFD全局索引的总和*/
  int _offset;

  /* Number of connecting neighbors for each surface cell 每个表面单元的连接邻居数量 一个表示从当前表面出发的连接，一个表示到达当前表面的连接*/
  int** num_connections;

  /* Indexes of connecting neighbors for each surface cell 每个表面单元的连接邻居的索引*/
  int*** indexes;

  /* Surface numbers of connecting neighbors for each surface cell 每个表面单元的连接邻居的表面数*/
  int*** domains;

  /* Coupling coeffs between connecting neighbors and itself for each 
     surface cell */
  CMFD_PRECISION*** coupling_coeffs;

  /* Fluxes of connecting neighbors for each surface cell*/
  CMFD_PRECISION*** fluxes;

  /* Buffer for sending/receiving fluxes to/from connecting neighbors */
  CMFD_PRECISION** buffer;

  /* Map to the index of the boundary elements */
  std::map<int, int> mapLocalToSurface;

  int num_groups;
  bool stop;
#ifdef ENABLE_MPI_
  MPI_Comm _MPI_cart;
#endif
};


/**
 * @brief Get coupling fluxes and other information from neighbors. 
 * @details The information are transfered by reference.
 */
#ifdef ENABLE_MPI_
void getCouplingTerms(DomainCommunicator* comm, int color, int*& coupling_sizes,
                      int**& coupling_indexes, CMFD_PRECISION**& coupling_coeffs,
                      CMFD_PRECISION**& coupling_fluxes,
                      CMFD_PRECISION* curr_fluxes, int& offset);
#endif

double eigenvalueSolve(Matrix* A, Matrix* M, Vector* X, double k_eff,
                       double tol, double SOR_factor=1.5,
                       ConvergenceData* convergence_data = NULL,
                       DomainCommunicator* comm = NULL,
                       int empty_cells_num = 0);
bool linearSolve(Matrix* A, Matrix* M, Vector* X, Vector* B, double tol,
                 double SOR_factor=1.5,
                 ConvergenceData* convergence_data = NULL,
                 DomainCommunicator* comm = NULL,
                 int empty_cells_num = 0);
bool ddLinearSolve(Matrix* A, Matrix* M, Vector* X, Vector* B, double tol,
                   double SOR_factor, ConvergenceData* convergence_data,
                   DomainCommunicator* comm);
void matrixMultiplication(Matrix* A, Vector* X, Vector* B);
double computeRMSE(Vector* x, Vector* y, bool integrated,
                         DomainCommunicator* comm = NULL);

/* Hex Lattice Solve */
double hexeigenvalueSolve(Matrix* A, Matrix* M, Vector* X, double k_eff,
                       double tol, double SOR_factor=1.5,
                       ConvergenceData* convergence_data = NULL,
                       DomainCommunicator* comm = NULL,
                       int empty_cells_num = 0,
                       std::vector<bool> empty_fsrs_cells = std::vector<bool>(),
                       std::vector<int> logical_actual_map = std::vector<int>());
bool hexlinearSolve(Matrix* A, Matrix* M, Vector* X, Vector* B, double tol,
                    double SOR_factor=1.5,
                    ConvergenceData* convergence_data = NULL,
                    DomainCommunicator* comm = NULL,
                    int empty_cells_num = 0,
                    std::vector<bool> empty_fsrs_cells = std::vector<bool>(),
                    std::vector<int> logical_actual_map = std::vector<int>());
bool hexddLinearSolve(Matrix* A, Matrix* M, Vector* X, Vector* B, double tol,
                    double SOR_factor, ConvergenceData* convergence_data,
                    DomainCommunicator* comm,
                    int empty_cells_num = 0,
                    std::vector<bool> empty_fsrs_cells = std::vector<bool>(),
                    std::vector<int> logical_actual_map = std::vector<int>());
double hexComputeRMSE(Vector* x, Vector* y, bool integrated,
                      DomainCommunicator* comm = NULL,
                      int empty_cells_num = 0);


/**
 * @brief Transpose a 2D matrix.
 * @param matrix array to transpose
 * @param dim1 first dimension length
 * @param dim2 second dimension length
 */
template<typename T>
inline void matrix_transpose(T* matrix, int dim1, int dim2) {

  std::vector<T> temp(dim1 * dim2);

  for (int i=0; i < dim1; i++) {
    for (int j=0; j < dim2; j++)
      temp[i * dim1 + j] = matrix[j * dim1 + i];
  }

  std::copy(temp.begin(), temp.end(), matrix);
}

} /* namespace antmoc */

#endif /* LINALG_H_ */
