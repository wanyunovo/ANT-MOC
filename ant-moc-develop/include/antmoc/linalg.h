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
   * @brief CMFD 特征值求解器的详细迭代信息
   */
  struct ConvergenceData
  {

    /* CMFD 的最大延拓因子（Prolongation Factor）
     * 【物理含义】：
     * - CMFD 求解后得到粗网通量，需要"延长"回 MOC 细网
     * - pf = max(新MOC通量 / 旧MOC通量) 在所有 FSR 和能群上
     * - 衡量 CMFD 对 MOC 通量的修正幅度
     */
    double pf;

    /* CMFD 特征值问题的初始残差（外层第 1 次迭代）
     */
    double cmfd_res_1;

    /* The final residual of the CMFD eigenvalue problem (residual at the last outer iteration)CMFD 特征值问题的最终残差（外层最后一次迭代）
     */
    double cmfd_res_end;

    /* 内层线性求解器第 1 次迭代时的残差
     */
    double linear_res_1;

    /* 内层线性求解器在外层最后一次迭代时的残差
     */
    double linear_res_end;

    /* 外层 CMFD 特征值迭代的总次数
     * 【典型值】：
     * - 正常收敛：5 ~ 20 次
     * - 慢收敛：20 ~ 50 次（可能需要更好的初值或预处理）
     * - 不收敛：达到最大限制（如 100 次）仍未满足容差
     */
    int cmfd_iters;

    /* 外层第 1 次迭代时，内层线性求解器的迭代次数
     * 【参考值】：
     * - 使用 SOR 方法：10 ~ 50 次
     * - 使用 Krylov 方法（如 GMRES）：5 ~ 20 次
     */
    int linear_iters_1;

    /* 外层最后一次迭代时，内层线性求解器的迭代次数
     * 【预期行为】：
     * - 通常 linear_iters_end < linear_iters_1（因为外层接近收敛，初值更好）
     * - 如果反而增加，说明矩阵性质在变差（需检查物理合理性）
     */
    int linear_iters_end;

    /* 构造函数：初始化所有字段为 -1
     * 【-1 的含义】：
     * - 作为"未赋值"或"未计算"的标记
     * - 在打印日志时，如果看到 -1，说明该数据未被求解器填充
     * - 用于调试：快速发现代码路径中是否正确调用了求解器
     */
    ConvergenceData()
    {
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
   * [MPI并行] 域通信器结构体
   *
   * 【背景】域分解并行计算原理：
   * - 将整个反应堆几何区域划分为多个子区域（domain）
   * - 每个MPI进程负责计算一个子区域
   * - 相邻子区域的边界网格需要交换通量信息
   * - 这个结构体管理所有与边界通信相关的数据
   */
  struct DomainCommunicator
  {

    /* ==================== 域分解拓扑信息 ==================== */

    /* 全局域分解配置：整个几何在X/Y/Z方向上被切分成多少个域 */
    int _num_domains_x; // X方向域的总数
    int _num_domains_y; // Y方向域的总数
    int _num_domains_z; // Z方向域的总数

    /* 当前进程负责的域的索引（从0开始编号）
     * 例如：8个进程负责一个2×2×2的域分解，则 idx 范围是 [0,1]×[0,1]×[0,1]
     */
    int _domain_idx_x; // 当前域在X方向上的索引
    int _domain_idx_y; // 当前域在Y方向上的索引
    int _domain_idx_z; // 当前域在Z方向上的索引

    /* 当前域的负责的CMFD网格数量（不包括幽灵网格） */
    int _local_num_x; // 本域X方向的CMFD网格数
    int _local_num_y; // 本域Y方向的CMFD网格数
    int _local_num_z; // 本域Z方向的CMFD网格数

    /* ==================== 全局索引偏移量 ==================== */

    /* 全局索引起始偏移
     * - 作用：将本地索引转换为全局索引
     * - 例如：如果前面的域有100个网格，则当前域的 _offset=100
     * - 全局索引 = 本地索引 + _offset
     */
    int _offset;

    /* ==================== 边界连接信息（处理域间耦合）==================== */

    /* Number of connecting neighbors for each surface cell
     * 二维数组：记录每个边界单元连接的邻居域数量
     * - 第1维：边界单元索引
     * - 第2维：连接类型（如：入/出方向）
     * 注：边界单元可能与多个邻居域相连（如角点/边缘处）
     */
    int **num_connections;

    /* Indexes of connecting neighbors for each surface cell
     * 三维数组：存储与边界单元相连的邻居单元的全局索引
     * - 第1维：边界单元索引
     * - 第2维：连接类型
     * - 第3维：具体哪个邻居的索引
     */
    int ***indexes;

    /* Surface numbers of connecting neighbors for each surface cell
     * 三维数组：存储连接发生在邻居的哪个表面
     * - 用于确定耦合发生在邻居域的哪个面（如X+面、Y-面等）
     */
    int ***domains;

    /* Coupling coeffs between connecting neighbors and itself for each surface cell
     * 三维数组：耦合系数矩阵
     * - 描述边界单元与邻居单元之间的物理耦合强度
     * - 来自扩散方程离散化的非对角项
     */
    CMFD_PRECISION ***coupling_coeffs;

    /* Fluxes of connecting neighbors for each surface cell
     * 三维数组：存储从邻居域接收到的边界通量值
     * - 用于构建本域线性方程组的右端项
     */
    CMFD_PRECISION ***fluxes;

    /* Buffer for sending/receiving fluxes to/from connecting neighbors
     * 二维数组：MPI通信缓冲区
     * - 发送缓冲区：打包本域边界的通量数据
     * - 接收缓冲区：解包从邻居域收到的通量数据
     */
    CMFD_PRECISION **buffer;

    /* Map to the index of the boundary elements
     * 映射表：本地单元ID -> 边界数组中的索引
     * - 用于快速查找某个单元在边界通信数组中的位置
     * - std::map<本地网格ID, 边界数组索引>
     */
    std::map<int, int> mapLocalToSurface;

    /* 能群数量（用于分配通信缓冲区大小）*/
    int num_groups;

    /* 停止标志（用于控制迭代终止）*/
    bool stop;

#ifdef ENABLE_MPI_
    /* MPI笛卡尔拓扑通信子
     * - 将MPI进程组织成3D网格结构，便于相邻域通信
     * - 支持自动确定邻居进程的rank（进程编号）
     */
    MPI_Comm _MPI_cart;
#endif
  };

/**
 * @brief Get coupling fluxes and other information from neighbors.
 * @details The information are transfered by reference.
 *
 * @param comm 域通信器指针
 * @param color 颜色标记（用于红黑排序等算法）
 * @param coupling_sizes 输出：耦合数据的大小
 * @param coupling_indexes 输出：耦合单元的索引
 * @param coupling_coeffs 输出：耦合系数
 * @param coupling_fluxes 输出：耦合通量
 * @param curr_fluxes 输入：当前域的通量
 * @param offset 输出：索引偏移量
 */
#ifdef ENABLE_MPI_
  void getCouplingTerms(DomainCommunicator *comm, int color, int *&coupling_sizes,
                        int **&coupling_indexes, CMFD_PRECISION **&coupling_coeffs,
                        CMFD_PRECISION **&coupling_fluxes,
                        CMFD_PRECISION *curr_fluxes, int &offset);
#endif

  double eigenvalueSolve(Matrix *A, Matrix *M, Vector *X, double k_eff,
                         double tol, double SOR_factor = 1.5,
                         ConvergenceData *convergence_data = NULL,
                         DomainCommunicator *comm = NULL,
                         int empty_cells_num = 0);
  bool linearSolve(Matrix *A, Matrix *M, Vector *X, Vector *B, double tol,
                   double SOR_factor = 1.5,
                   ConvergenceData *convergence_data = NULL,
                   DomainCommunicator *comm = NULL,
                   int empty_cells_num = 0);
  bool ddLinearSolve(Matrix *A, Matrix *M, Vector *X, Vector *B, double tol,
                     double SOR_factor, ConvergenceData *convergence_data,
                     DomainCommunicator *comm);
  void matrixMultiplication(Matrix *A, Vector *X, Vector *B);
  double computeRMSE(Vector *x, Vector *y, bool integrated,
                     DomainCommunicator *comm = NULL);

  /* Hex Lattice Solve */
  double hexeigenvalueSolve(Matrix *A, Matrix *M, Vector *X, double k_eff,
                            double tol, double SOR_factor = 1.5,
                            ConvergenceData *convergence_data = NULL,
                            DomainCommunicator *comm = NULL,
                            int empty_cells_num = 0,
                            std::vector<bool> empty_fsrs_cells = std::vector<bool>(),
                            std::vector<int> logical_actual_map = std::vector<int>());
  bool hexlinearSolve(Matrix *A, Matrix *M, Vector *X, Vector *B, double tol,
                      double SOR_factor = 1.5,
                      ConvergenceData *convergence_data = NULL,
                      DomainCommunicator *comm = NULL,
                      int empty_cells_num = 0,
                      std::vector<bool> empty_fsrs_cells = std::vector<bool>(),
                      std::vector<int> logical_actual_map = std::vector<int>());
  bool hexddLinearSolve(Matrix *A, Matrix *M, Vector *X, Vector *B, double tol,
                        double SOR_factor, ConvergenceData *convergence_data,
                        DomainCommunicator *comm,
                        int empty_cells_num = 0,
                        std::vector<bool> empty_fsrs_cells = std::vector<bool>(),
                        std::vector<int> logical_actual_map = std::vector<int>());
  double hexComputeRMSE(Vector *x, Vector *y, bool integrated,
                        DomainCommunicator *comm = NULL,
                        int empty_cells_num = 0);

  /**
   * @brief Transpose a 2D matrix.
   * @param matrix array to transpose
   * @param dim1 first dimension length
   * @param dim2 second dimension length
   */
  template <typename T>
  inline void matrix_transpose(T *matrix, int dim1, int dim2)
  {

    std::vector<T> temp(dim1 * dim2);

    for (int i = 0; i < dim1; i++)
    {
      for (int j = 0; j < dim2; j++)
        temp[i * dim1 + j] = matrix[j * dim1 + i];
    }

    std::copy(temp.begin(), temp.end(), matrix);
  }

} /* namespace antmoc */

#endif /* LINALG_H_ */
