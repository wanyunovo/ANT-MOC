#include "antmoc/linalg.h"
#include "antmoc/log.h"

#include <fstream>
#include <fenv.h>
#include <cmath>
#include <omp.h>


namespace antmoc
{

/**
 * @brief Solves a generalized eigenvalue problem using a power iteration method  
 * 使用幂迭代方法（Power Iteration Method）来求解广义特征值问题的算法
 * @details This function takes in a loss + streaming Matrix (A),
 *          a fission gain Matrix (M), a flux Vector (X), a tolerance used
 *          for both the power method and linear solve convergence (tol), and
 *          a successive over-relaxation factor (SOR_factor) and computes the
 *          dominant eigenvalue and eigenvector using the Power method. The
 *          eigenvalue is returned and the input X Vector is modified in
 *          place to be the corresponding eigenvector.
 * @param A the loss + streaming Matrix object
 * @param M the fission gain Matrix object
 * @param X the flux Vector object    _new_flux  (其实老的通量复制过来的)
 * @param k_eff initial k_effective
 * @param tol the power method and linear solve source convergence threshold   _source_convergence_threshold
 * @param SOR_factor the successive over-relaxation factor
 * @param convergence_data a summary of how to solver converged
 * @param comm an MPI communicator for the domain-decomposed solver
 * @return k_eff the dominant eigenvalue
 */
double eigenvalueSolve(Matrix* A, Matrix* M, Vector* X, double k_eff,
                             double tol, double SOR_factor,
                             ConvergenceData* convergence_data,
                             DomainCommunicator* comm,
                             int empty_cells_num) {

  log::fverbose("Computing the Matrix-Vector eigenvalue...");
  tol = std::max(MIN_LINALG_TOLERANCE, tol);  ////确保容差 tol 不小于最小线性代数容差 MIN_LINALG_TOLERANCE

  /* Check for consistency of matrix and vector dimensions 
   维度一致性检查：依次检查矩阵 A 和 M 以及向量 X 的 X, Y, Z 和 Groups 维度是否一致。如果不一致，记录错误并停止执行。
  */
  if (A->getNumX() != M->getNumX() || A->getNumX() != X->getNumX())
    log::ferror("Cannot compute the Matrix-Vector eigenvalue with "
               "different x dimensions for the A matrix, M matrix, and X vector"
               ": (%d, %d, %d)", A->getNumX(), M->getNumX(), X->getNumX());
  else if (A->getNumY() != M->getNumY() || A->getNumY() != X->getNumY())
    log::ferror("Cannot compute the Matrix-Vector eigenvalue with "
               "different y dimensions for the A matrix, M matrix, and X vector"
               ": (%d, %d, %d)", A->getNumY(), M->getNumY(), X->getNumY());
  else if (A->getNumZ() != M->getNumZ() || A->getNumZ() != X->getNumZ())
    log::ferror("Cannot compute the Matrix-Vector eigenvalue with "
               "different z dimensions for the A matrix, M matrix, and X vector"
               ": (%d, %d, %d)", A->getNumZ(), M->getNumZ(), X->getNumZ());
  else if (A->getNumGroups() != M->getNumGroups() ||
           A->getNumGroups() != X->getNumGroups())
    log::ferror("Cannot compute the Matrix-Vector eigenvalue with "
               "different num groups  for the A matrix, M matrix, and X vector:"
               " (%d, %d, %d)", A->getNumGroups(), M->getNumGroups(),
               X->getNumGroups());

  /* Initialize variables 初始化变量*/
  omp_lock_t* cell_locks = X->getCellLocks();
  int num_x = X->getNumX();
  int num_y = X->getNumY();
  int num_z = X->getNumZ();
  int num_groups = X->getNumGroups();
  //分别初始化为旧的源向量和新的源向量，用于后续计算
  Vector old_source(cell_locks, num_x, num_y, num_z, num_groups, empty_cells_num);
  Vector new_source(cell_locks, num_x, num_y, num_z, num_groups, empty_cells_num);
  double residual;  //存储残差，用于判断收敛性
  int iter;

  /* Compute and normalize the initial source 初始源的计算与归一化
  //传过去的矩阵是M,M*X得到了新的old_source[里面的值发生了变化,本来全是0],计算初始源项，
  即矩阵 M 与向量 X 的乘积，结果保存在 old_source 中,对应B.30中M*粗网格通量
  */
  matrixMultiplication(M, X, &old_source); 
  double old_source_sum = old_source.getSum();  //old_source 中array的每一项求和
  int num_rows = X->getNumRows();
#ifdef ENABLE_MPI_
  if (comm != NULL) {
    //使用 MPI_Allreduce 进行全局求和操作，确保在所有进程间共享相同的源项和行数
    double temp_sum = old_source_sum;
    MPI_Allreduce(&temp_sum, &old_source_sum, 1, MPI_DOUBLE, MPI_SUM,
        comm->_MPI_cart);  
    int temp_rows = num_rows;
    MPI_Allreduce(&temp_rows, &num_rows, 1, MPI_INT, MPI_SUM,
        comm->_MPI_cart);
  }
#endif
  //归一化处理
  old_source.scaleByValue(num_rows / old_source_sum); //把old_source Vector中的array中所有值都乘以*num_rows / old_source_sum倍,理解为是一个归一化操作
  X->scaleByValue(num_rows * k_eff / old_source_sum); //把传过来的老通量的所有值也乘以*num_rows * k_eff / old_source_sum倍,注意老通量的值发生变化

  /* Power iteration Matrix-Vector solver 幂迭代的主循环*/
  double initial_residual = 0;
  bool solver_failure = false;  //线性求解器是否失败
  for (iter = 0; iter < MAX_LINALG_POWER_ITERATIONS; iter++) {  //主迭代循环，执行最多 MAX_LINALG_POWER_ITERATIONS 25000 次

    /* Solve X = A^-1 * old_source 首先使用常规的线性求解器 linearSolve 来解决 X = A^-1 * old_source 问题,对应B.30中的左侧项*/
    bool converged = false;
    if (!solver_failure)
      converged = linearSolve(A, M, X, &old_source, tol*1e-1, SOR_factor,
                              convergence_data, comm, empty_cells_num);

    /* If the solver failed, try the diagonally dominant solver 常规线性求解器未能收敛，请尝试对角占优求解器ddLinearSolve 进行冗余求解*/
    if (!converged) {
      solver_failure = true;
      converged = ddLinearSolve(A, M, X, &old_source, tol*1e-1, SOR_factor,
                                convergence_data, comm);
    }

    /* Check for divergence */
    if (!converged)
      return -1.0; //对角占优求解器也未能收敛，则认为算法发散，返回 -1.0 作为错误标志

    /* Compute the new source 使用更新后的 X 向量，计算新的源项 new_source*/
    matrixMultiplication(M, X, &new_source);

    /* Compute the sum of new sources 计算新源的总和*/
    double new_source_sum = new_source.getSum();
#ifdef ENABLE_MPI_
    if (comm != NULL) {
      //使用 MPI_Allreduce 确保 new_source_sum 是全局求和的结果
      double temp_sum = new_source_sum;
      MPI_Allreduce(&temp_sum, &new_source_sum, 1, MPI_DOUBLE, MPI_SUM,
          comm->_MPI_cart);
    }
#endif

    /* Compute and set keff */
    k_eff = new_source_sum / num_rows;   //计算新的特征值 k_eff

    /* Scale the new source by keff */
    new_source.scaleByValue(1.0 / k_eff);  //归一化 new_source，以确保其与 k_eff 一致

    /* Compute the residual */
    residual = computeRMSE(&new_source, &old_source, true, comm); //计算新旧源项之间的均方根误差（RMSE），作为当前迭代的残差
    if (iter == 0) {  //在第一次迭代时，记录初始残差，并将其用于后续的收敛性检查。如果初始残差过小，将其设为一个较小的值（1e-10）
      initial_residual = residual;
      if (initial_residual < 1e-14)
        initial_residual = 1e-10;
      if (convergence_data != NULL) {  //如果 convergence_data 非空，则更新收敛信息，记录初始残差及对应的线性求解迭代次数和残差
        convergence_data->cmfd_res_1 = residual;
        convergence_data->linear_iters_1 = convergence_data->linear_iters_end;
        convergence_data->linear_res_1 = convergence_data->linear_res_end;
      }
    }

    /* Copy the new source to the old source 将新的源项复制到旧的源项中，为下一次迭代做准备*/
    new_source.copyTo(&old_source);

    //输出当前迭代的 iter 值、特征值 k_eff 和残差 residual。
    log::fverbose_once("Matrix-Vector eigenvalue iter: %d, keff: %f, residual: "
               "%3.2e", iter, k_eff, residual);

    /* Check for convergence 收敛性检查
    如果当前残差相对于初始残差下降超过97%（residual / initial_residual < 0.03），
    或者当前残差小于最小线性代数容差 MIN_LINALG_TOLERANCE，
    并且迭代次数超过最小迭代次数 MIN_LINALG_POWER_ITERATIONS 25次，则认为算法收敛
    */
    if ((residual / initial_residual < 0.03 || residual < MIN_LINALG_TOLERANCE)
        && iter > MIN_LINALG_POWER_ITERATIONS) {
      if (convergence_data != NULL) { //记录最终的残差和迭代次数
        convergence_data->cmfd_res_end = residual;
        convergence_data->cmfd_iters = iter;
      }
      break;
    }
  }

  log::fverbose_once("Matrix-Vector eigenvalue solve iterations: %d", iter);
  if (iter == MAX_LINALG_POWER_ITERATIONS)  //如果达到最大迭代次数 MAX_LINALG_POWER_ITERATIONS 仍未收敛，记录错误日志并结束
    log::ferror("Eigenvalue solve failed to converge in %d iterations", iter);

  return k_eff;
}


double hexeigenvalueSolve(Matrix* A, Matrix* M, Vector* X, double k_eff,
                             double tol, double SOR_factor,
                             ConvergenceData* convergence_data,
                             DomainCommunicator* comm,
                             int empty_cells_num,
                             std::vector<bool> empty_fsrs_cells,
                             std::vector<int> logical_actual_map) {

  log::fverbose("Computing the Matrix-Vector eigenvalue...");
  tol = std::max(MIN_LINALG_TOLERANCE, tol);

  /* Check for consistency of matrix and vector dimensions */
  if (A->getNumX() != M->getNumX() || A->getNumX() != X->getNumX())
    log::ferror("Cannot compute the Matrix-Vector eigenvalue with "
               "different x dimensions for the A matrix, M matrix, and X vector"
               ": (%d, %d, %d)", A->getNumX(), M->getNumX(), X->getNumX());
  else if (A->getNumY() != M->getNumY() || A->getNumY() != X->getNumY())
    log::ferror("Cannot compute the Matrix-Vector eigenvalue with "
               "different y dimensions for the A matrix, M matrix, and X vector"
               ": (%d, %d, %d)", A->getNumY(), M->getNumY(), X->getNumY());
  else if (A->getNumZ() != M->getNumZ() || A->getNumZ() != X->getNumZ())
    log::ferror("Cannot compute the Matrix-Vector eigenvalue with "
               "different z dimensions for the A matrix, M matrix, and X vector"
               ": (%d, %d, %d)", A->getNumZ(), M->getNumZ(), X->getNumZ());
  else if (A->getNumGroups() != M->getNumGroups() ||
           A->getNumGroups() != X->getNumGroups())
    log::ferror("Cannot compute the Matrix-Vector eigenvalue with "
               "different num groups  for the A matrix, M matrix, and X vector:"
               " (%d, %d, %d)", A->getNumGroups(), M->getNumGroups(),
               X->getNumGroups());

  /* Initialize variables */
  omp_lock_t* cell_locks = X->getCellLocks();
  int num_x = X->getNumX();
  int num_y = X->getNumY();
  int num_z = X->getNumZ();
  int num_groups = X->getNumGroups();
  Vector old_source(cell_locks, num_x, num_y, num_z, num_groups, empty_cells_num);
  Vector new_source(cell_locks, num_x, num_y, num_z, num_groups, empty_cells_num);
  double residual;
  int iter;

  /*
  log::finfo("HexeigenvalueSolve Before scale X:");
  X->printString();*/

  /* Compute and normalize the initial source */
  matrixMultiplication(M, X, &old_source);

  /*
  log::finfo("HexeigenvalueSolve Before scale old_source:");
  old_source.printString();*/

  double old_source_sum = old_source.getSum();
  int num_rows = X->getNumRows();
#ifdef ENABLE_MPI_
  if (comm != NULL) {
    double temp_sum = old_source_sum;
    MPI_Allreduce(&temp_sum, &old_source_sum, 1, MPI_DOUBLE, MPI_SUM,
        comm->_MPI_cart);
    int temp_rows = num_rows;
    MPI_Allreduce(&temp_rows, &num_rows, 1, MPI_INT, MPI_SUM,
        comm->_MPI_cart);
  }
#endif
  old_source.scaleByValue(num_rows / old_source_sum);
  X->scaleByValue(num_rows * k_eff / old_source_sum);

  /*
  log::finfo("HexeigenvalueSolve After scale X:");
  X->printString();
  log::finfo("HexeigenvalueSolve After scale old_source:");
  old_source.printString();*/

  /* Power iteration Matrix-Vector solver */
  double initial_residual = 0;
  bool solver_failure = false;
  for (iter = 0; iter < MAX_LINALG_POWER_ITERATIONS; iter++) {

    /* Solve X = A^-1 * old_source */
    bool converged = false;
    if (!solver_failure)
      converged = hexlinearSolve(A, M, X, &old_source, tol*1e-1, SOR_factor,
                              convergence_data, comm, empty_cells_num,
                              empty_fsrs_cells, logical_actual_map);

    /* If the solver failed, try the diagonally dominant solver */
    if (!converged) {
      solver_failure = true;
      converged = hexddLinearSolve(A, M, X, &old_source, tol*1e-1, SOR_factor,
                              convergence_data, comm, empty_cells_num,
                              empty_fsrs_cells, logical_actual_map);
    }

    /* Check for divergence */
    if (!converged)
      return -1.0;

    /* Compute the new source */
    matrixMultiplication(M, X, &new_source);

    /*
    log::finfo("HexeigenvalueSolve Iteration %d:", iter);
    log::finfo("X is:");
    X->printString();
    log::finfo("old_source is:");
    new_source.printString();*/

    /* Compute the sum of new sources */
    double new_source_sum = new_source.getSum();
#ifdef ENABLE_MPI_
    if (comm != NULL) {
      double temp_sum = new_source_sum;
      MPI_Allreduce(&temp_sum, &new_source_sum, 1, MPI_DOUBLE, MPI_SUM,
          comm->_MPI_cart);
    }
#endif

    /* Compute and set keff */
    k_eff = new_source_sum / num_rows;

    /* Scale the new source by keff */
    new_source.scaleByValue(1.0 / k_eff);

    /*
    log::finfo("new_source is:");
    new_source.printString();*/

    /* Compute the residual */
    residual = hexComputeRMSE(&new_source, &old_source, true, comm, empty_cells_num);
    if (iter == 0) {
      initial_residual = residual;
      if (initial_residual < 1e-14)
        initial_residual = 1e-10;
      if (convergence_data != NULL) {
        convergence_data->cmfd_res_1 = residual;
        convergence_data->linear_iters_1 = convergence_data->linear_iters_end;
        convergence_data->linear_res_1 = convergence_data->linear_res_end;
      }
    }

    /* Copy the new source to the old source */
    new_source.copyTo(&old_source);

    log::fverbose_once("Matrix-Vector eigenvalue iter: %d, keff: %f, residual: "
               "%3.2e", iter, k_eff, residual);

    /* Check for convergence */
    if ((residual / initial_residual < 0.03 || residual < MIN_LINALG_TOLERANCE)
        && iter > MIN_LINALG_POWER_ITERATIONS) {
      if (convergence_data != NULL) {
        convergence_data->cmfd_res_end = residual;
        convergence_data->cmfd_iters = iter;
      }
      break;
    }
  }

  log::fverbose_once("Matrix-Vector eigenvalue solve iterations: %d", iter);
  if (iter == MAX_LINALG_POWER_ITERATIONS)
    log::ferror("Eigenvalue solve failed to converge in %d iterations", iter);

  return k_eff;
}


/**
 * @brief Solves a linear system using Red-Black Gauss Seidel with
 *        successive over-relaxation.
 * 求解线性方程组的函数，使用了红黑 Gauss-Seidel 方法，并结合了连续超松弛 (SOR) 技术
 * @details This function takes in a loss + streaming Matrix (A),
 *          a fission gain Matrix (M), a flux Vector (X), a source Vector (B),
 *          a source convergence tolerance (tol) and a successive
 *          over-relaxation factor (SOR_factor) and computes the
 *          solution to the linear system. The input X Vector is modified in
 *          place to be the solution vector.
 * @param A the loss + streaming Matrix object   A损失矩阵
 * @param M the fission gain Matrix object   M裂变矩阵
 * @param X the flux Vector object  // 老通量old_flux对象，输入时为初始猜测，输出时为解向量。已经进行归一化处理的
 * @param B the source Vector object  //old_source 为M*old_flux得来  也已经归一化了
 * @param tol the power method and linear solve source convergence threshold  源收敛阈值
 * @param SOR_factor the successive over-relaxation factor  //调整 Gauss-Seidel 更新步长,用于SOR算法
 * @param convergence_data a summary of how to solver converged  求解器收敛数据
 * @param comm an MPI communicator for the domain-decomposed solver
 */
bool linearSolve(Matrix* A, Matrix* M, Vector* X, Vector* B, double tol,
                 double SOR_factor, ConvergenceData* convergence_data,
                 DomainCommunicator* comm, int empty_cells_num) {

  bool success = true;
  tol = std::max(MIN_LINALG_TOLERANCE, tol);  //确保容差 tol 不小于系统设置的最小容差 MIN_LINALG_TOLERANCE

  /* Check for consistency of matrix and vector dimensions */
  if (A->getNumX() != B->getNumX() || A->getNumX() != X->getNumX() ||
      A->getNumX() != M->getNumX())
    log::ferror("Cannot perform linear solve with different x dimensions"
               " for the A matrix, M matrix, B vector, and X vector: "
               "(%d, %d, %d, %d)", A->getNumX(), M->getNumX(),
               B->getNumX(), X->getNumX());
  else if (A->getNumY() != B->getNumY() || A->getNumY() != X->getNumY() ||
           A->getNumY() != M->getNumY())
    log::ferror("Cannot perform linear solve with different y dimensions"
               " for the A matrix, M matrix, B vector, and X vector: "
               "(%d, %d, %d, %d)", A->getNumY(), M->getNumY(),
               B->getNumY(), X->getNumY());
  else if (A->getNumZ() != B->getNumZ() || A->getNumZ() != X->getNumZ() ||
           A->getNumZ() != M->getNumZ())
    log::ferror("Cannot perform linear solve with different z dimensions"
               " for the A matrix, M matrix, B vector, and X vector: "
               "(%d, %d, %d, %d)", A->getNumZ(), M->getNumZ(),
               B->getNumZ(), X->getNumZ());
  else if (A->getNumGroups() != B->getNumGroups() ||
           A->getNumGroups() != X->getNumGroups() ||
           A->getNumGroups() != M->getNumGroups())
    log::ferror("Cannot perform linear solve with different num groups"
               " for the A matrix, M matrix, B vector, and X vector: "
               "(%d, %d, %d, %d)", A->getNumGroups(), M->getNumGroups(),
               B->getNumGroups(), X->getNumGroups());

  /* Initialize variables 初始化变量*/
  double residual;    //当前的残差
  double min_residual = 1e6;  //最小残差初始化为一个大的数值
  int iter = 0;  //迭代次数初始为 0
  omp_lock_t* cell_locks = X->getCellLocks();
  //存储网格的 x、y、z 维度以及能群数量
  int num_x = X->getNumX();
  int num_y = X->getNumY();
  int num_z = X->getNumZ();
  int num_groups = X->getNumGroups();
  // unused
  // int num_rows = X->getNumRows();
  //Vector X_old(cell_locks, num_x, num_y, num_z, num_groups);  //FIXME delete
  //CMFD_PRECISION* x_old = X_old.getArray();
  // unused
  // CMFD_PRECISION* x_old = X_old.getArray();
  // 指向稀疏矩阵 A 的行起始位置和列索引
  int* IA = A->getIA(); 
  int* JA = A->getJA();
  CMFD_PRECISION* DIAG = A->getDiag();  //指向矩阵 A 的对角线元素数组
  CMFD_PRECISION* a = A->getA();   //指向矩阵 A 的非零元素数组
  //分别指向 X 和 B 向量的元素数组
  CMFD_PRECISION* x = X->getArray();  //old_flux的数组所有值
  CMFD_PRECISION* b = B->getArray();  //old_source的数组所有值
  //用于存储旧的和新的源项
  Vector old_source(cell_locks, num_x, num_y, num_z, num_groups, empty_cells_num);
  Vector new_source(cell_locks, num_x, num_y, num_z, num_groups, empty_cells_num);

  /* Compute initial source 计算初始源项*/
  matrixMultiplication(M, X, &old_source);  //将矩阵 M 与向量 X 相乘，结果存储在 old_source 中

#ifdef ENABLE_MPI_
  // Initialize communication buffers
  int* coupling_sizes = NULL;
  int** coupling_indexes = NULL;
  CMFD_PRECISION** coupling_coeffs = NULL;
  CMFD_PRECISION** coupling_fluxes = NULL;
#endif

  double initial_residual = 0;
  while (iter < MAX_LINEAR_SOLVE_ITERATIONS) {  //不超过10000

    /* Pass new flux to old flux */
    //X->copyTo(&X_old);

    // Iteration over red/black cells  红黑 Gauss-Seidel 迭代
    for (int color = 0; color < 2; color++) { //使用两种颜色（0 和 1）分别进行 Gauss-Seidel 迭代  
      int offset = 0;  //offset 用于调整单元格在不同线程中的计算顺序
#ifdef ENABLE_MPI_
      getCouplingTerms(comm, color, coupling_sizes, coupling_indexes,
                       coupling_coeffs, coupling_fluxes, x, offset);  //用于获取耦合项
#endif

#pragma omp parallel for collapse(2)
      for (int iz=0; iz < num_z; iz++) {
        for (int iy=0; iy < num_y; iy++) {
          for (int ix=(iy+iz+color+offset)%2; ix < num_x; ix+=2) {

            int cell = (iz*num_y + iy)*num_x + ix;  //计算当前单元格在整个网格中的索引
            int row_start = cell*num_groups;

            /* Find index into communicator buffers for cells on surfaces 
            判断当前单元格是否在域的表面，并根据需要获取该表面单元的索引 domain_surface_index*/
            bool on_surface = (iz==0) || (iz==num_z-1) || (iy==0) || (iy==num_y-1)
                 || (ix==0) || (ix==num_x-1);
            int domain_surface_index = -1;
            if (comm != NULL && on_surface)  //如果在域边界
              domain_surface_index = comm->mapLocalToSurface[cell];

            /* Contribution of off-diagonal terms, hard to SIMD vectorize 非对角项的贡献，难以进行SIMD矢量化*/
            for (int g=0; g < num_groups; g++) {  //对每一个能群 g，计算对应的行 row 的非对角线项的贡献

              int row = row_start + g;  
              x[row] = (1.0 - SOR_factor) * x[row] * (DIAG[row] / SOR_factor);
              log::fdebug("x[%d] is %.6f", row, x[row]);
              
              if (fabs(DIAG[row]) < FLT_EPSILON )
                log::ferror("A zero has been found on the diagonal of "
                            "the CMFD matrix cell [%d,%d,%d]=%d, group %d",
                            ix, iy, iz, cell, g);

              for (int i = IA[row]; i < IA[row+1]; i++) {  //IA存储每行第一个非零元素在a中的起始位置索引,循环次数为每行非零元素的个数

                int col = JA[i];  //JA存储a中每个非零元素对应的列索引 
                // Get the column index
                if (row != col) //!=代表不是位于对角线的元素
                  x[row] -= a[i] * x[col];  //a是损失矩阵A的非零值数组,x为old_flux中的第col个值,这句代码的作用是消除第row行方程中所有非对角线元素对当前行解x[row]的影响
                else
                  x[row] += b[row]; //b为old_source的array数组,对角线元素,增加x的对应值,x为old_flux的数组所有值
              }
#ifdef ENABLE_MPI_
              // Contribution of off node fluxes
              if (comm != NULL && on_surface) {
                int row_surf = domain_surface_index * num_groups + g;
                for (int i = 0; i < coupling_sizes[row_surf]; i++) {
                  int idx = coupling_indexes[row_surf][i] * num_groups + g;
                  int domain = comm->domains[color][row_surf][i];
                  CMFD_PRECISION flux = coupling_fluxes[domain][idx];
                  x[row] -= coupling_coeffs[row_surf][i] * flux;
                }
              }
#endif
              // Perform these operations separately, for performance 为了提高性能，请分别执行这些操作
              //(SOR_factor / DIAG[row]) 是一个调整因子，它将当前的 x[row] 值缩放至合适的大小，以确保 SOR 迭代按照预期进行
              x[row] *= (SOR_factor / DIAG[row]);
            }
          }
        }
      }
    }

    // Compute the new source
    matrixMultiplication(M, X, &new_source);  //计算新的源项,并赋值给new_source

    // Compute the residual 计算残差
    feclearexcept (FE_ALL_EXCEPT);  //清除程序中所有可能的浮点异常标志
    if (iter == 0) {   //如果是第一次迭代
      residual = computeRMSE(&new_source, &old_source, true, comm);
      if (convergence_data != NULL)
        convergence_data->linear_res_end = residual;
      initial_residual = residual;
    }

    // Increment the interations counter
    iter++;

    double ratio_residuals = 1;
    if (initial_residual > 0)
      ratio_residuals = residual / initial_residual;
    log::fdebug("SOR iter: %d, residual: %3.2e, initial residual: %3.2e"
               ", ratio = %3.2e, tolerance: %3.2e, end? %d", iter, residual,
               initial_residual, ratio_residuals, tol,
               (ratio_residuals < 0.1 || residual < tol) &&
               iter > MIN_LINEAR_SOLVE_ITERATIONS);

    // Compute residual only after minimum iteration number
    if (iter > MIN_LINEAR_SOLVE_ITERATIONS) {

      residual = computeRMSE(&new_source, &old_source, true, comm);

      // Record current minimum residual  记录较小的残差
      if (residual < min_residual)
        min_residual = residual;

      // Check for going off the rails 检查是否求解器出现发散（即求解未收敛）的情况
      int raised = fetestexcept (FE_INVALID); //如果发生了 FE_INVALID (浮点运算相关的异常)异常，则 raised 的值中相应的位会被设置为 1；如果没有发生，则相应位为 0
      if ((residual > 1e3 * min_residual && min_residual > 1e-10) || raised) {
        log::fwarn("linear solve divergent : res %e min_res %e NaN? %d",
                   residual, min_residual, raised);//输出的有这一句
        if (convergence_data != NULL)
          convergence_data->linear_iters_end = iter;
        success = false;
        break;   //代表求解失败,不收敛
      }

      // Check for convergence
      if (residual / initial_residual < 0.1 || residual < tol) {
        if (convergence_data != NULL)
          convergence_data->linear_iters_end = iter;
        break;  //代表已经收敛,求解成功
      }
    }

    // Copy the new source to the old source   迭代次数在25-9999,把新的源项赋值给旧的源项
    if (iter > MIN_LINEAR_SOLVE_ITERATIONS - 1 &&
        iter < MAX_LINEAR_SOLVE_ITERATIONS)
      new_source.copyTo(&old_source);  //将新的源项赋值为旧的源项
  }

  log::fdebug("linear solve iterations: %d", iter);

  // Check if the maximum iterations were reached
  if (iter == MAX_LINEAR_SOLVE_ITERATIONS) {  //如果到了最大迭代次数,还没有break收敛,则判定为求解失败
    matrixMultiplication(M, X, &new_source);
    double residual = computeRMSE(&new_source, &old_source, true, comm);
    log::fverbose("Ratio = %3.2e, tol = %3.2e", residual / initial_residual,
               tol);
    log::finfo("Linear solve failed to converge in %d iterations with "
               "initial residual %3.2e and final residual %3.2e", iter,
               initial_residual, residual);
    success = false;
  }
  return success;
}


bool hexlinearSolve(Matrix* A, Matrix* M, Vector* X, Vector* B, double tol,
                 double SOR_factor, ConvergenceData* convergence_data,
                 DomainCommunicator* comm, int empty_cells_num,
                 std::vector<bool> empty_fsrs_cells,
                 std::vector<int> logical_actual_map) {

  bool success = true;
  tol = std::max(MIN_LINALG_TOLERANCE, tol);

  /* Check for consistency of matrix and vector dimensions */
  if (A->getNumX() != B->getNumX() || A->getNumX() != X->getNumX() ||
      A->getNumX() != M->getNumX())
    log::ferror("Cannot perform linear solve with different x dimensions"
               " for the A matrix, M matrix, B vector, and X vector: "
               "(%d, %d, %d, %d)", A->getNumX(), M->getNumX(),
               B->getNumX(), X->getNumX());
  else if (A->getNumY() != B->getNumY() || A->getNumY() != X->getNumY() ||
           A->getNumY() != M->getNumY())
    log::ferror("Cannot perform linear solve with different y dimensions"
               " for the A matrix, M matrix, B vector, and X vector: "
               "(%d, %d, %d, %d)", A->getNumY(), M->getNumY(),
               B->getNumY(), X->getNumY());
  else if (A->getNumZ() != B->getNumZ() || A->getNumZ() != X->getNumZ() ||
           A->getNumZ() != M->getNumZ())
    log::ferror("Cannot perform linear solve with different z dimensions"
               " for the A matrix, M matrix, B vector, and X vector: "
               "(%d, %d, %d, %d)", A->getNumZ(), M->getNumZ(),
               B->getNumZ(), X->getNumZ());
  else if (A->getNumGroups() != B->getNumGroups() ||
           A->getNumGroups() != X->getNumGroups() ||
           A->getNumGroups() != M->getNumGroups())
    log::ferror("Cannot perform linear solve with different num groups"
               " for the A matrix, M matrix, B vector, and X vector: "
               "(%d, %d, %d, %d)", A->getNumGroups(), M->getNumGroups(),
               B->getNumGroups(), X->getNumGroups());

  /* Initialize variables */
  double residual;
  double min_residual = 1e6;
  int iter = 0;
  omp_lock_t* cell_locks = X->getCellLocks();
  int num_x = X->getNumX();
  int num_y = X->getNumY();
  int num_z = X->getNumZ();
  int num_groups = X->getNumGroups();
  
  // unused
  // int num_rows = X->getNumRows();
  //Vector X_old(cell_locks, num_x, num_y, num_z, num_groups);  //FIXME delete
  //CMFD_PRECISION* x_old = X_old.getArray();
  // unused
  // CMFD_PRECISION* x_old = X_old.getArray();
  int* IA = A->getIA();
  int* JA = A->getJA();
  CMFD_PRECISION* DIAG = A->getDiag();
  CMFD_PRECISION* a = A->getA();
  CMFD_PRECISION* x = X->getArray();
  CMFD_PRECISION* b = B->getArray();
  Vector old_source(cell_locks, num_x, num_y, num_z, num_groups, empty_cells_num);
  Vector new_source(cell_locks, num_x, num_y, num_z, num_groups, empty_cells_num);

  /* Compute initial source */
  matrixMultiplication(M, X, &old_source);

#ifdef ENABLE_MPI_
  // Initialize communication buffers
  int* coupling_sizes = NULL;
  int** coupling_indexes = NULL;
  CMFD_PRECISION** coupling_coeffs = NULL;
  CMFD_PRECISION** coupling_fluxes = NULL;
#endif

  double initial_residual = 0;
  while (iter < MAX_LINEAR_SOLVE_ITERATIONS) {

    /*
    log::finfo("HexlinearSolve Iterations %d:", iter);
    log::finfo("HexlinearSolve X befor iter:");
    X->printString();*/

    /* Pass new flux to old flux */
    //X->copyTo(&X_old);

    // Iteration over red/black cells
    for (int color = 0; color < 2; color++) {
      int offset = 0;
#ifdef ENABLE_MPI_
      getCouplingTerms(comm, color, coupling_sizes, coupling_indexes,
                       coupling_coeffs, coupling_fluxes, x, offset);
#endif

#pragma omp parallel for collapse(2)
      for (int iz=0; iz < num_z; iz++) {
        for (int iy=0; iy < num_y; iy++) {
          for (int ix=(iy+iz+color+offset)%2; ix < num_x; ix+=2) {

            int cell = (iz*num_y + iy)*num_x + ix;
            if (!empty_fsrs_cells[cell]) {
              cell = logical_actual_map[cell];
              int row_start = cell*num_groups;

              /* Find index into communicator buffers for cells on surfaces */
              bool on_surface = (iz==0) || (iz==num_z-1) || (iy==0) || (iy==num_y-1)
                  || (ix==0) || (ix==num_x-1);
              int domain_surface_index = -1;
              if (comm != NULL && on_surface)
                domain_surface_index = comm->mapLocalToSurface[cell];

              /* Contribution of off-diagonal terms, hard to SIMD vectorize */
              for (int g=0; g < num_groups; g++) {

                int row = row_start + g;
                x[row] = (1.0 - SOR_factor) * x[row] * (DIAG[row] / SOR_factor);
                //log::fdebug("x[%d] is %.2f, SOR_factor is %.2f, DIAG[%d] is %.2f", row, x[row], SOR_factor, row, DIAG[row]);
                
                if (fabs(DIAG[row]) < FLT_EPSILON )
                  log::ferror("A zero has been found on the diagonal of "
                              "the CMFD matrix cell %d, group %d", cell, g);

                for (int i = IA[row]; i < IA[row+1]; i++) {

                  // Get the column index
                  int col = JA[i];
                  if (row != col)
                    x[row] -= a[i] * x[col];
                  else
                    x[row] += b[row];
                }
#ifdef ENABLE_MPI_
                // Contribution of off node fluxes
                if (comm != NULL && on_surface) {
                  int row_surf = domain_surface_index * num_groups + g;
                  for (int i = 0; i < coupling_sizes[row_surf]; i++) {
                    int idx = coupling_indexes[row_surf][i] * num_groups + g;
                    int domain = comm->domains[color][row_surf][i];
                    CMFD_PRECISION flux = coupling_fluxes[domain][idx];
                    x[row] -= coupling_coeffs[row_surf][i] * flux;
                  }
                }
#endif
                // Perform these operations separately, for performance
                x[row] *= (SOR_factor / DIAG[row]);
              }
            }
          }
        }
      }
    }

    // Compute the new source
    matrixMultiplication(M, X, &new_source);

    /*
    log::finfo("HexlinearSolve X after iter:");
    X->printString();
    log::finfo("HexlinearSolve old_source");
    old_source.printString();
    log::finfo("HexlinearSolve new_source");
    new_source.printString();*/

    // Compute the residual
    feclearexcept (FE_ALL_EXCEPT);
    if (iter == 0) {
      residual = hexComputeRMSE(&new_source, &old_source, true, comm, empty_cells_num);
      if (convergence_data != NULL)
        convergence_data->linear_res_end = residual;
      initial_residual = residual;
    }

    // Increment the interations counter
    iter++;

    double ratio_residuals = 1;
    if (initial_residual > 0)
      ratio_residuals = residual / initial_residual;
    log::fdebug("SOR iter: %d, residual: %3.2e, initial residual: %3.2e"
               ", ratio = %3.2e, tolerance: %3.2e, end? %d", iter, residual,
               initial_residual, ratio_residuals, tol,
               (ratio_residuals < 0.1 || residual < tol) &&
               iter > MIN_LINEAR_SOLVE_ITERATIONS);

    // Compute residual only after minimum iteration number
    if (iter > MIN_LINEAR_SOLVE_ITERATIONS) {
      //log::fdebug("hexlinearSolve iterations:%d", iter);
      residual = hexComputeRMSE(&new_source, &old_source, true, comm, empty_cells_num);

      // Record current minimum residual
      if (residual < min_residual)
        min_residual = residual;
      
      log::fdebug("LinearSolve residual: %e, min_res: %e", residual, min_residual);

      // Check for going off the rails
      int raised = fetestexcept (FE_INVALID);
      if ((residual > 1e3 * min_residual && min_residual > 1e-10) || raised) {
        log::fwarn("linear solve divergent : res %e min_res %e NaN? %d",
                   residual, min_residual, raised);
        if (convergence_data != NULL)
          convergence_data->linear_iters_end = iter;
        success = false;
        break;
      }

      // Check for convergence
      if (residual / initial_residual < 0.1 || residual < tol) {
        if (convergence_data != NULL)
          convergence_data->linear_iters_end = iter;
        break;
      }
    }

    // Copy the new source to the old source
    if (iter > MIN_LINEAR_SOLVE_ITERATIONS - 1 &&
        iter < MAX_LINEAR_SOLVE_ITERATIONS)
      new_source.copyTo(&old_source);
  }

  log::fdebug("linear solve iterations: %d", iter);

  // Check if the maximum iterations were reached
  if (iter == MAX_LINEAR_SOLVE_ITERATIONS) {
    matrixMultiplication(M, X, &new_source);
    double residual = hexComputeRMSE(&new_source, &old_source, true, comm, empty_cells_num);
    log::fverbose("Ratio = %3.2e, tol = %3.2e", residual / initial_residual,
               tol);
    log::finfo("Linear solve failed to converge in %d iterations with "
               "initial residual %3.2e and final residual %3.2e", iter,
               initial_residual, residual);
    success = false;
  }
  return success;
}


#ifdef ENABLE_MPI_
/**
 * @brief Get coupling fluxes and other information from neighbors. 
 *        The information are transfered by reference.
 * @param comm Structure for communication of fluxes between neighbor domains
 * @param color red or black color
 * @param coupling_sizes Number of connecting neighbors for each surface cell
 * @param coupling_indexes Surface numbers of connecting neighbors for each 
 *        surface cell
 * @param coupling_coeffs Coupling coeffs between connecting neighbors and 
 *        itself for each surface cell
 * @param coupling_fluxes Fluxes of connecting neighbors for each surface cell
 * @param curr_fluxes CMFD cell fluxes of current iteration
 * @param offset Sum of the starting CMFD global indexes of a domain, for 
 *        calculation of the color
 */
void getCouplingTerms(DomainCommunicator* comm, int color, int*& coupling_sizes,
                      int**& coupling_indexes, CMFD_PRECISION**& coupling_coeffs,
                      CMFD_PRECISION**& coupling_fluxes,
                      CMFD_PRECISION* curr_fluxes, int& offset) {
  if (comm != NULL) {
    coupling_sizes = comm->num_connections[color];
    coupling_indexes = comm->indexes[color];
    coupling_coeffs = comm->coupling_coeffs[color];
    coupling_fluxes = comm->fluxes[color];

    offset = comm->_offset;

    MPI_Request requests[2*NUM_FACES];

    int nx = comm->_local_num_x;
    int ny = comm->_local_num_y;
    int nz = comm->_local_num_z;

    int ng = comm->num_groups;

    // Get numerical precision for communication
    MPI_Datatype flux_type;
    if (sizeof(CMFD_PRECISION) == 4)
      flux_type = MPI_FLOAT;
    else
      flux_type = MPI_DOUBLE;

    int sizes[NUM_FACES];
    for (int coord=0; coord < 3; coord++) {
      for (int d=0; d<2; d++) {

        int dir = 2*d-1;
        int surf = coord + 3*d;
        int op_surf = surf - 3*dir;
        int source, dest;

        MPI_Cart_shift(comm->_MPI_cart, coord, dir, &source, &dest);

        // Pack MPI buffer
        int size = 0;
        if (surf == SURFACE_X_MIN) {
          size = ny * nz * ng;
          for (int i=0; i < nz; i++)
            for (int j=0; j < ny; j++)
              for (int g=0; g < ng; g++)
                comm->buffer[surf][ng*(i*ny+j)+g] =
                  curr_fluxes[ng*((i*ny + j)*nx) + g];
        }
        else if (surf == SURFACE_X_MAX) {
          size = ny * nz * ng;
          for (int i=0; i < nz; i++)
            for (int j=0; j < ny; j++)
              for (int g=0; g < ng; g++)
                comm->buffer[surf][ng*(i*ny+j)+g] =
                  curr_fluxes[ng*((i*ny + j)*nx + nx-1) + g];
        }
        else if (surf == SURFACE_Y_MIN) {
          size = nx * nz * ng;
          for (int i=0; i < nz; i++)
            for (int j=0; j < nx; j++)
              for (int g=0; g < ng; g++)
                comm->buffer[surf][ng*(i*nx+j)+g] =
                  curr_fluxes[ng*(i*nx*ny + j) + g];
        }
        else if (surf == SURFACE_Y_MAX) {
          size = nx * nz * ng;
          for (int i=0; i < nz; i++)
            for (int j=0; j < nx; j++)
              for (int g=0; g < ng; g++)
                comm->buffer[surf][ng*(i*nx+j)+g] =
                  curr_fluxes[ng*(i*nx*ny + j + nx*(ny-1)) + g];
        }
        else if (surf == SURFACE_Z_MIN) {
          size = nx * ny * ng;
          for (int i=0; i < ny; i++)
            for (int j=0; j < nx; j++)
              for (int g=0; g < ng; g++)
                comm->buffer[surf][ng*(i*nx+j)+g] =
                  curr_fluxes[ng*(i*nx + j)+g];
        }
        else if (surf == SURFACE_Z_MAX) {
          size = nx * ny * ng;
          for (int i=0; i < ny; i++)
            for (int j=0; j < nx; j++)
              for (int g=0; g < ng; g++)
                comm->buffer[surf][ng*(i*nx+j)+g] =
                  curr_fluxes[ng*(i*nx + j + nx*ny*(nz-1)) + g];
        }

        sizes[surf] = size;

        // Post send
        MPI_Isend(comm->buffer[surf], size, flux_type,
                  dest, 0, comm->_MPI_cart, &requests[2*surf]);

        // Post receive
        MPI_Irecv(&comm->buffer[op_surf][size], size, flux_type,
                  source, 0, comm->_MPI_cart, &requests[2*op_surf+1]);
      }
    }

    // Block for communication round to complete
    bool round_complete = false;
    while (!round_complete) {

      round_complete = true;
      int flag;
      MPI_Status send_stat;
      MPI_Status recv_stat;

      for (int coord=0; coord < 3; coord++) {
        for (int d=0; d<2; d++) {
          int surf = coord + 3*d;

          MPI_Test(&requests[2*surf], &flag, &send_stat);
          if (flag == 0)
            round_complete = false;

          MPI_Test(&requests[2*surf+1], &flag, &recv_stat);
          if (flag == 0)
            round_complete = false;
          else
            // Copy received data into coupling_fluxes
            for (int i=0; i < sizes[surf]; i++)
              coupling_fluxes[surf][i] = comm->buffer[surf][sizes[surf]+i];
        }
      }
    }
  }
}
#endif


/**
 * @brief Performs a matrix vector multiplication.
 * @details This function takes in a Matrix (A), a variable Vector (X),
 *          and a solution Vector (B) and computes the matrix vector product.
 *          The solution Vector is modified in place.
 * @param A a Matrix object
 * @param X the variable Vector object
 * @param B the solution Vector object
 */
void matrixMultiplication(Matrix* A, Vector* X, Vector* B) {

  /* Check for consistency of matrix and vector dimensions */
  if (A->getNumX() != B->getNumX() || A->getNumX() != X->getNumX())
    log::ferror("Cannot perform matrix multiplication  with different x "
               "dimensions for the A matrix, B vector, and X vector: "
               "(%d, %d, %d)", A->getNumX(), B->getNumX(), X->getNumX());
  else if (A->getNumY() != B->getNumY() || A->getNumY() != X->getNumY())
    log::ferror("Cannot perform matrix multiplication with different y "
               "dimensions for the A matrix, B vector, and X vector: "
               "(%d, %d, %d)", A->getNumY(), B->getNumY(), X->getNumY());
  else if (A->getNumZ() != B->getNumZ() || A->getNumZ() != X->getNumZ())
    log::ferror("Cannot perform matrix multiplication with different z "
               "dimensions for the A matrix, B vector, and X vector: "
               "(%d, %d, %d)", A->getNumZ(), B->getNumZ(), X->getNumZ());
  else if (A->getNumGroups() != B->getNumGroups() ||
           A->getNumGroups() != X->getNumGroups())
    log::ferror("Cannot perform matrix multiplication with different "
               "num groups for the A matrix, M matrix, and X vector: "
               "(%d, %d, %d)", A->getNumGroups(), B->getNumGroups(),
               X->getNumGroups());

  B->setAll(0.0);
  int* IA = A->getIA();
  int* JA = A->getJA();
  CMFD_PRECISION* a = A->getA();
  CMFD_PRECISION* x = X->getArray();
  CMFD_PRECISION* b = B->getArray();
  int num_rows = X->getNumRows();

#pragma omp parallel for
  for (int row = 0; row < num_rows; row++) {
    for (int i = IA[row]; i < IA[row+1]; i++) {
      b[row] += a[i] * x[JA[i]];
    }
  }
}


/**
 * @brief Computes the Root Mean Square Error of two Vectors.
 * 计算两个向量的均方根误差。
 * @details This function takes in two vectors (X and Y) and computes the
 *          Root Mean Square Error(均方根误差) of the Vector Y with respect to Vector X.
 *          The boolean integrated must also be given to indicate whether the
 *          operation on the vector should be group-wise integrated before
 *          performing the RMSE operation.
 * @param X a Vector object
 * @param Y a second Vector object
 * @param integrated a boolean indicating whether to group-wise integrate.
 * @param comm an MPI communicator to exchange residuals between domains
 */
double computeRMSE(Vector* X, Vector* Y, bool integrated,
                   DomainCommunicator* comm) {

  /* Check for consistency of vector dimensions */
  if (X->getNumX() != Y->getNumX() || X->getNumY() != Y->getNumY() ||
      X->getNumZ() != Y->getNumZ() || X->getNumGroups() != Y->getNumGroups())
    log::ferror("Cannot compute RMSE with different vector dimensions: "
               "(%d, %d, %d, %d) and (%d, %d, %d, %d)",
               X->getNumX(), X->getNumY(), X->getNumZ(), X->getNumGroups(),
               Y->getNumX(), Y->getNumY(), Y->getNumZ(), Y->getNumGroups());

  double rmse;
  double sum_residuals = 0.;
  int norm;
  int num_x = X->getNumX();
  int num_y = X->getNumY();
  int num_z = X->getNumZ();
  int num_groups = X->getNumGroups();
  omp_lock_t* cell_locks = X->getCellLocks();

  if (integrated) {

    double new_source, old_source;
    CMFD_PRECISION residual[num_x * num_y * num_z]
         __attribute__ ((aligned(VEC_ALIGNMENT)));//指定数组的内存对齐方式,VEC_ALIGNMENT 是一个预定义的值，表示要对齐的字节数
    memset(residual, 0, num_x * num_y * num_z * sizeof(CMFD_PRECISION)); //将 residual 数组中的所有元素初始化为零

    /* Compute the RMSE */
#pragma omp parallel for private(new_source, old_source)
    for (int i = 0; i < num_x*num_y*num_z; i++) {
      new_source = 0.0;
      old_source = 0.0;
#pragma omp simd reduction(+:new_source,old_source)
      for (int g = 0; g < num_groups; g++) {
        new_source += X->getValue(i, g);
        old_source += Y->getValue(i, g);
      }
      if (fabs(old_source) > FLUX_EPSILON)
        residual[i] = std::pow((new_source - old_source) / old_source, 2);//这里应该/old_source,还是new_source?
    }

    // Sum residuals
#pragma omp simd reduction(+:sum_residuals) aligned(residual)
    for (int i = 0; i < num_x*num_y*num_z; i++)
      sum_residuals += residual[i];

    norm = num_x * num_y * num_z;
  }
  else {

    //FIXME Incurs a memory allocation, uses unnecessary locks
    Vector residual(cell_locks, num_x, num_y, num_z, num_groups);

    /* Compute the RMSE */
#pragma omp parallel for
    for (int i = 0; i < num_x*num_y*num_z; i++) {
      for (int g = 0; g < num_groups; g++) {
        if (fabs(X->getValue(i, g)) > FLUX_EPSILON)
          residual.setValue(i, g, std::pow((X->getValue(i, g) - Y->getValue(i, g)) /
                                           X->getValue(i, g), 2));
      }
    }
    sum_residuals = residual.getSum();
    norm = num_x * num_y * num_z * num_groups;
  }

#ifdef ENABLE_MPI_
  if (comm != NULL) {
    double temp_residual = sum_residuals;
    MPI_Allreduce(&temp_residual, &sum_residuals, 1, MPI_DOUBLE, MPI_SUM,
        comm->_MPI_cart);
    int temp_norm = norm;
    MPI_Allreduce(&temp_norm, &norm, 1, MPI_INT, MPI_SUM, comm->_MPI_cart);
  }
#endif

  /* Error check residual components */
  if (sum_residuals < 0.0) {
    log::fwarn("CMFD Residual mean square error %6.4f less than zero",
               sum_residuals);
    sum_residuals = 0.0;
  }
  if (norm <= 0) {
    log::fwarn("CMFD residual norm %d less than one", norm);
    norm = 1;
  }

  /* Compute RMS residual error */
  rmse = sqrt(sum_residuals / norm);  //对应公式HINTS-6.37,不应该先开平方根,再除以norm吗?

#ifdef ENABLE_MPI_
  if (comm != NULL) {
    double temp_rmse = rmse;
    MPI_Allreduce(&temp_rmse, &rmse, 1, MPI_DOUBLE, MPI_MAX, comm->_MPI_cart);
  }
#endif

  return rmse;
}


/**
 * @brief Computes the Root Mean Square Error of two Vectors.
 * @details This function takes in two vectors (X and Y) and computes the
 *          Root Mean Square Error of the Vector Y with respect to Vector X.
 *          The boolean integrated must also be given to indicate whether the
 *          operation on the vector should be group-wise integrated before
 *          performing the RMSE operation.
 * @param X a Vector object
 * @param Y a second Vector object
 * @param integrated a boolean indicating whether to group-wise integrate.
 * @param comm an MPI communicator to exchange residuals between domains
 */
double hexComputeRMSE(Vector* X, Vector* Y, bool integrated,
                      DomainCommunicator* comm, int empty_cells_num) {

  /* Check for consistency of vector dimensions */
  if (X->getNumX() != Y->getNumX() || X->getNumY() != Y->getNumY() ||
      X->getNumZ() != Y->getNumZ() || X->getNumGroups() != Y->getNumGroups())
    log::ferror("Cannot compute RMSE with different vector dimensions: "
               "(%d, %d, %d, %d) and (%d, %d, %d, %d)",
               X->getNumX(), X->getNumY(), X->getNumZ(), X->getNumGroups(),
               Y->getNumX(), Y->getNumY(), Y->getNumZ(), Y->getNumGroups());

  double rmse;
  double sum_residuals = 0.;
  int norm;
  int num_x = X->getNumX();
  int num_y = X->getNumY();
  int num_z = X->getNumZ();
  int num_groups = X->getNumGroups();
  omp_lock_t* cell_locks = X->getCellLocks();

  if (integrated) {

    double new_source, old_source;
    CMFD_PRECISION residual[num_x * num_y * num_z - empty_cells_num]
         __attribute__ ((aligned(VEC_ALIGNMENT)));  //FIXME Overflow for large cases?
    memset(residual, 0, (num_x * num_y * num_z - empty_cells_num) * sizeof(CMFD_PRECISION));

    /* Compute the RMSE */
#pragma omp parallel for private(new_source, old_source)
    for (int i = 0; i < num_x*num_y*num_z-empty_cells_num; i++) {
      new_source = 0.0;
      old_source = 0.0;
#pragma omp simd reduction(+:new_source,old_source)
      for (int g = 0; g < num_groups; g++) {
        new_source += X->getValue(i, g);
        old_source += Y->getValue(i, g);
      }
      if (fabs(old_source) > FLUX_EPSILON) {
        residual[i] = std::pow((new_source - old_source) / old_source, 2);
        //log::fdebug("cell(%d) new_source:%3.2e, old_source:%3.2e, residual[%d] = %3.2e", i, new_source, old_source, i, residual[i]);
      }
    }

    // Sum residuals
#pragma omp simd reduction(+:sum_residuals) aligned(residual)
    for (int i = 0; i < num_x*num_y*num_z-empty_cells_num; i++)
      sum_residuals += residual[i];

    norm = num_x * num_y * num_z - empty_cells_num;
  }
  else {

    //FIXME Incurs a memory allocation, uses unnecessary locks
    Vector residual(cell_locks, num_x, num_y, num_z, num_groups, empty_cells_num);

    /* Compute the RMSE */
#pragma omp parallel for
    for (int i = 0; i < num_x*num_y*num_z-empty_cells_num; i++) {
      for (int g = 0; g < num_groups; g++) {
        if (fabs(X->getValue(i, g)) > FLUX_EPSILON)
          residual.setValue(i, g, std::pow((X->getValue(i, g) - Y->getValue(i, g)) /
                                           X->getValue(i, g), 2));
      }
    }
    sum_residuals = residual.getSum();
    norm = (num_x * num_y * num_z - empty_cells_num) * num_groups;
  }

#ifdef ENABLE_MPI_
  if (comm != NULL) {
    double temp_residual = sum_residuals;
    MPI_Allreduce(&temp_residual, &sum_residuals, 1, MPI_DOUBLE, MPI_SUM,
        comm->_MPI_cart);
    int temp_norm = norm;
    MPI_Allreduce(&temp_norm, &norm, 1, MPI_INT, MPI_SUM, comm->_MPI_cart);
  }
#endif

  /* Error check residual components */
  if (sum_residuals < 0.0) {
    log::fwarn("CMFD Residual mean square error %6.4f less than zero",
               sum_residuals);
    sum_residuals = 0.0;
  }
  if (norm <= 0) {
    log::fwarn("CMFD residual norm %d less than one", norm);
    norm = 1;
  }

  /* Compute RMS residual error */
  rmse = sqrt(sum_residuals / norm);

#ifdef ENABLE_MPI_
  if (comm != NULL) {
    double temp_rmse = rmse;
    MPI_Allreduce(&temp_rmse, &rmse, 1, MPI_DOUBLE, MPI_MAX, comm->_MPI_cart);
  }
#endif
  //log::fdebug("RMSE is %.2f", rmse);
  return rmse;
}


/**
 * @brief Solves a linear system using the linear solver above, but makes the
 *        loss and streaming matrix diagonally dominant first, to increase
 *        likelihood of convergence.
 * 使用上述线性求解器求解线性系统，但首先使损失和流矩阵对角占优，以增加收敛的可能性
 * @details This function takes in a loss + streaming Matrix (A),
 *          a fission gain Matrix (M), a flux Vector (X), a source Vector (B),
 *          a source convergence tolerance (tol) and a successive
 *          over-relaxation factor (SOR_factor) and makes (A) diagonally
 *          dominant before calling the linear solve routine to compute the
 *          solution to the linear system. The input X Vector is modified in
 *          place to be the solution vector. The transformation to make (A)
 *          diagonally dominant is compensated by another matrix multiplication.
 * @param A the loss + streaming Matrix object
 * @param M the fission gain Matrix object
 * @param X the flux Vector object
 * @param B the source Vector object
 * @param tol the power method and linear solve source convergence threshold
 * @param SOR_factor the successive over-relaxation factor
 * @param convergence_data a summary of the convergence performance
 * @param comm a communicator for exchanging data through MPI
 */
bool ddLinearSolve(Matrix* A, Matrix* M, Vector* X, Vector* B, double tol,
                   double SOR_factor, ConvergenceData* convergence_data,
                   DomainCommunicator* comm) {

  /* Create vector for stabilizing flux */
  omp_lock_t* cell_locks = X->getCellLocks();
  int num_x = X->getNumX();
  int num_y = X->getNumY();
  int num_z = X->getNumZ();
  int num_groups = X->getNumGroups();
  int num_rows = X->getNumRows();

  Vector dd(cell_locks, num_x, num_y, num_z, num_groups);
  dd.setAll(0.0);
  
  CMFD_PRECISION* dd_array = dd.getArray();
  CMFD_PRECISION* x = X->getArray();

  /* Stabilize matrix A to be diagonally dominant */
  CMFD_PRECISION* a = A->getA(); //指向矩阵 A 的非零元素数组
  // unused
  // CMFD_PRECISION* a_diag = A->getDiag();
  // 指向稀疏矩阵 A 的行起始位置和列索引
  int* IA = A->getIA();
  int* JA = A->getJA();

#ifdef ENABLE_MPI_
  // Initialize communication buffers
  int* coupling_sizes = NULL;
  int** coupling_indexes = NULL;
  CMFD_PRECISION** coupling_coeffs = NULL;
  CMFD_PRECISION** coupling_fluxes = NULL;
#endif

  /* Loop over cells */
  for (int color = 0; color < 2; color++) {
    int offset = 0;
#ifdef ENABLE_MPI_
    getCouplingTerms(comm, color, coupling_sizes, coupling_indexes,
                     coupling_coeffs, coupling_fluxes, x, offset);
#endif
#pragma omp parallel for collapse(2)
    for (int iz=0; iz < num_z; iz++) {
      for (int iy=0; iy < num_y; iy++) {
        for (int ix=(iy+iz+color+offset)%2; ix < num_x; ix+=2) {

          int cell = (iz*num_y + iy)*num_x + ix; //计算当前单元格在整个网格中的索引

          /* Find index into communicator buffers for cells on surfaces 
          判断当前单元格是否在域的表面，并根据需要获取该表面单元的索引 domain_surface_index*/
          bool on_surface = (iz==0) || (iz==num_z-1) || (iy==0) || (iy==num_y-1)
               || (ix==0) || (ix==num_x-1);
          int domain_surface_index = -1;
          if (comm != NULL && on_surface) //如果在域边界
            domain_surface_index = comm->mapLocalToSurface[cell];

          /* Determine whether each group's row is diagonally dominant 
          确定每组的行是否对角占优*/
          for (int e = 0; e < num_groups; e++) {

            int row = cell * num_groups + e;

            /* Add local off-diagonal elements */
            double row_sum = 0.0;
            int diag_ind = -1;
            for (int idx = IA[row]; idx < IA[row+1]; idx++) { //IA存储每行第一个非零元素在a中的起始位置索引,循环次数为每行非零元素的个数
              if (JA[idx] != row)  //代表不是位于对角线的元素
                row_sum += fabs(a[idx]);
              else
                diag_ind = idx;
            }

            /* Add off-node off-diagonal elements 添加对角线元素外的节点*/
#ifdef ENABLE_MPI_
            if (comm != NULL && on_surface) {
              int row_surf = domain_surface_index * num_groups + e;
              int* coupling_sizes = comm->num_connections[color];
              CMFD_PRECISION** coupling_coeffs = comm->coupling_coeffs[color];
              for (int idx=0; idx < coupling_sizes[row_surf]; idx++)
                row_sum += fabs(coupling_coeffs[row_surf][idx]);
            }
#endif

            /* Check for diagonal dominance 检查对角线优势
            这里的 row_sum 是对非对角线元素的绝对值求和。
            如果这些非对角线元素的和超过了当前对角线元素的值，
            就需要增加对角线元素，使其对角占优*/
            if (row_sum > a[diag_ind])  //如果不对角占优，将差值存储在 dd 向量中
              dd.incrementValue(cell, e, row_sum - a[diag_ind]);
          }
        }
      }
    }
  }

  /* Adjust matrix A to be diagonally dominant 将矩阵A调整为对角占优*/
#pragma omp parallel for
  for (int i=0; i < num_x*num_y*num_z; i++) {
    for (int e=0; e < num_groups; e++) {
      A->incrementValue(i, e, i, e, dd.getValue(i,e));
    }
  }

  /* Create a vector for the remainder right hand side 
  右侧向量 RHS 来存储调整后的右端项*/
  Vector RHS(cell_locks, num_x, num_y, num_z, num_groups);
  CMFD_PRECISION* rhs_array = RHS.getArray();
  CMFD_PRECISION* b = B->getArray();

  /* Keep track of sources */
  Vector old_source(cell_locks, num_x, num_y, num_z, num_groups);
  Vector new_source(cell_locks, num_x, num_y, num_z, num_groups);

  /* Calculate source */
  matrixMultiplication(M, X, &new_source);

  /* Iterate to get the total solution */
  double initial_residual = 0.0;
  double residual = 0.0;
  double min_residual = 1e6;
  for (int iter=0; iter < MAX_LINEAR_SOLVE_ITERATIONS; iter++) {

    // Copy the new source to the old source
    new_source.copyTo(&old_source);

    for (int row=0; row < num_rows; row++)
      rhs_array[row] = b[row] + dd_array[row] * x[row];

    bool converged = linearSolve(A, M, X, &RHS, tol, SOR_factor,
                                 convergence_data, comm);
    if (!converged)
      log::ferror("Stabilized linear solver inner iteration failed"
                 " to converge");

    // Compute the new source
    matrixMultiplication(M, X, &new_source);

    // Compute the residual
    residual = computeRMSE(&new_source, &old_source, true, comm);
    if (iter == 0){
      initial_residual = residual;
      if (initial_residual < 1e-14)
        initial_residual = 1e-10;
    }

    // Record current minimum residual
    if (residual < min_residual)
      min_residual = residual;

    // Check for going off the rails
    if (residual > 1e3 * min_residual && min_residual > 1e-10) {
      log::fwarn("Inner linear solve divergent.");
      log::finfo("Residual = %6.4e, Min Res = %6.4e", residual, min_residual);
      return false;
    }

    // Check for convergence
    if (residual / initial_residual < 0.1 || residual < tol)
      break;
  }

  /* Reset matrix A 在完成求解后，将之前增加的对角占优部分从矩阵 A 中减去，恢复其原始状态*/
#pragma omp parallel for
  for (int i=0; i < num_x*num_y*num_z; i++) {
    for (int e=0; e < num_groups; e++) {
      A->incrementValue(i, e, i, e, -dd.getValue(i,e));
    }
  }

  return true;
}


bool hexddLinearSolve(Matrix* A, Matrix* M, Vector* X, Vector* B, double tol,
                   double SOR_factor, ConvergenceData* convergence_data,
                   DomainCommunicator* comm, int empty_cells_num,
                   std::vector<bool> empty_fsrs_cells,
                   std::vector<int> logical_actual_map) {

  /* Create vector for stabilizing flux */
  omp_lock_t* cell_locks = X->getCellLocks();
  int num_x = X->getNumX();
  int num_y = X->getNumY();
  int num_z = X->getNumZ();
  int num_groups = X->getNumGroups();
  int num_rows = X->getNumRows();

  Vector dd(cell_locks, num_x, num_y, num_z, num_groups, empty_cells_num);
  dd.setAll(0.0);

  CMFD_PRECISION* dd_array = dd.getArray();
  CMFD_PRECISION* x = X->getArray();

  /* Stabilize matrix A to be diagonally dominant */
  CMFD_PRECISION* a = A->getA();
  // unused
  // CMFD_PRECISION* a_diag = A->getDiag();
  int* IA = A->getIA();
  int* JA = A->getJA();

#ifdef ENABLE_MPI_
  // Initialize communication buffers
  int* coupling_sizes = NULL;
  int** coupling_indexes = NULL;
  CMFD_PRECISION** coupling_coeffs = NULL;
  CMFD_PRECISION** coupling_fluxes = NULL;
#endif

  /* Loop over cells */
  for (int color = 0; color < 2; color++) {
    int offset = 0;
#ifdef ENABLE_MPI_
    getCouplingTerms(comm, color, coupling_sizes, coupling_indexes,
                     coupling_coeffs, coupling_fluxes, x, offset);
#endif
#pragma omp parallel for collapse(2)
    for (int iz=0; iz < num_z; iz++) {
      for (int iy=0; iy < num_y; iy++) {
        for (int ix=(iy+iz+color+offset)%2; ix < num_x; ix+=2) {

          int cell = (iz*num_y + iy)*num_x + ix;
          if (!empty_fsrs_cells[cell]) {
            cell = logical_actual_map[cell];

            /* Find index into communicator buffers for cells on surfaces */
            bool on_surface = (iz==0) || (iz==num_z-1) || (iy==0) || (iy==num_y-1)
                || (ix==0) || (ix==num_x-1);
            int domain_surface_index = -1;
            if (comm != NULL && on_surface)
              domain_surface_index = comm->mapLocalToSurface[cell];

            /* Determine whether each group's row is diagonally dominant */
            for (int e = 0; e < num_groups; e++) {

              int row = cell * num_groups + e;

              /* Add local off-diagonal elements */
              double row_sum = 0.0;
              int diag_ind = -1;
              for (int idx = IA[row]; idx < IA[row+1]; idx++) {
                if (JA[idx] != row)
                  row_sum += fabs(a[idx]);
                else
                  diag_ind = idx;
              }

              /* Add off-node off-diagonal elements */
#ifdef ENABLE_MPI_
              if (comm != NULL && on_surface) {
                int row_surf = domain_surface_index * num_groups + e;
                int* coupling_sizes = comm->num_connections[color];
                CMFD_PRECISION** coupling_coeffs = comm->coupling_coeffs[color];
                for (int idx=0; idx < coupling_sizes[row_surf]; idx++)
                  row_sum += fabs(coupling_coeffs[row_surf][idx]);
              }
#endif
              /* Check for diagonal dominance */
              if (row_sum > a[diag_ind])
                dd.incrementValue(cell, e, row_sum - a[diag_ind]);
            }
          }
        }
      }
    }
  }

  /* Adjust matrix A to be diagonally dominant */
#pragma omp parallel for
  for (int i=0; i < num_x*num_y*num_z-empty_cells_num; i++) {
    for (int e=0; e < num_groups; e++) {
      A->incrementValue(i, e, i, e, dd.getValue(i,e));
    }
  }

  /* Create a vector for the remainder right hand side */
  Vector RHS(cell_locks, num_x, num_y, num_z, num_groups, empty_cells_num);
  CMFD_PRECISION* rhs_array = RHS.getArray();
  CMFD_PRECISION* b = B->getArray();

  /* Keep track of sources */
  Vector old_source(cell_locks, num_x, num_y, num_z, num_groups, empty_cells_num);
  Vector new_source(cell_locks, num_x, num_y, num_z, num_groups, empty_cells_num);

  /* Calculate source */
  matrixMultiplication(M, X, &new_source);

  /* Iterate to get the total solution */
  double initial_residual = 0.0;
  double residual = 0.0;
  double min_residual = 1e6;
  for (int iter=0; iter < MAX_LINEAR_SOLVE_ITERATIONS; iter++) {

    // Copy the new source to the old source
    new_source.copyTo(&old_source);

    for (int row=0; row < num_rows; row++)
      rhs_array[row] = b[row] + dd_array[row] * x[row];

    bool converged = hexlinearSolve(A, M, X, &RHS, tol, SOR_factor,
                                  convergence_data, comm, empty_cells_num,
                                  empty_fsrs_cells, logical_actual_map);
    if (!converged)
      log::ferror("Stabilized linear solver inner iteration failed"
                 " to converge");

    // Compute the new source
    matrixMultiplication(M, X, &new_source);

    // Compute the residual
    residual = hexComputeRMSE(&new_source, &old_source, true, comm, empty_cells_num);
    if (iter == 0){
      initial_residual = residual;
      if (initial_residual < 1e-14)
        initial_residual = 1e-10;
    }

    // Record current minimum residual
    if (residual < min_residual)
      min_residual = residual;
    log::fdebug("ddLinearSolve residual: %e, min_res: %e", residual, min_residual);

    // Check for going off the rails
    if (residual > 1e3 * min_residual && min_residual > 1e-10) {
      log::fwarn("Inner linear solve divergent.");
      log::finfo("Residual = %6.4e, Min Res = %6.4e", residual, min_residual);
      return false;
    }

    // Check for convergence
    if (residual / initial_residual < 0.1 || residual < tol)
      break;
  }

  /* Reset matrix A */
#pragma omp parallel for
  for (int i=0; i < num_x*num_y*num_z-empty_cells_num; i++) {
    for (int e=0; e < num_groups; e++) {
      A->incrementValue(i, e, i, e, -dd.getValue(i,e));
    }
  }

  return true;
}


} /* namespace antmoc */
