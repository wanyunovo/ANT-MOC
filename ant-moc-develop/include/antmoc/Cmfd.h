/**
 * @file Cmfd.h
 * @brief The Cmfd class.
 * @date October 14, 2013
 * @author Sam Shaner, MIT, Course 22 (shaner@mit.edu)
 */

#ifndef CMFD_H_
#define CMFD_H_
/*
作用：这是 C/C++ 的标准惯例。防止这个头文件被同一个源文件（.cpp）多次包含（include）。
原理：如果没有这几行，当 A 文件包含了 Cmfd.h，B 文件也包含了 Cmfd.h，而 C 文件同时包含了 A 和 B 时，编译器会报错说 Cmfd 类被重复定义了。加上这个保护符后，第一次包含时定义 CMFD_H_，第二次包含时因为 CMFD_H_ 已经存在，就会跳过整个文件内容。
*/

#ifdef __cplusplus // C++ 环境检查与数学常量
#define _USE_MATH_DEFINES

#include <memory> //引入 C++ 智能指针支持（如 std::shared_ptr），用于自动管理内存，防止内存泄漏。
#include <omp.h>  //引入 OpenMP 库。这是用于多线程并行计算的标准库（利用单台机器上的多个 CPU 核心）。

#include "constants.h"
#include "enum_types.h"
#include "linalg.h"
#include "log.h"
#include "container_utils.h"
#include "math_utils.h"
#include "Point.h"
#include "Universe.h"

#ifdef ENABLE_MPI_
#include <mpi.h> //如果在编译命令中加了 -DENABLE_MPI_（通常在超级计算机或集群上编译时），编译器就会引入 <mpi.h>，开启多进程通信功能。
#endif

#endif

/** Optimization macro for 3D calculations to avoid branch statements
 * 如果在编译时定义了 THREED，则定义宏 _SOLVE_3D 为 true。
 * 目的：注释提到“avoid branch statements”（避免分支语句）。在高性能计算的核心循环中，if (is_3d) { ... } else { ... } 这样的判断非常耗时。通过宏定义，编译器可以在编译阶段就把不需要的代码（比如 2D 的逻辑）直接剔除，生成更高效的机器码，而不需要在运行时反复判断。
 */
#ifdef THREED
#define _SOLVE_3D (true)
#endif

namespace antmoc
{

  /** Forward declarations */
  class Geometry;
  class RecLattice;
  class HexLattice;
  class Lattice;
  class LocalCoords;
  class Material;
  class Matrix;
  class Point;
  class Quadrature;
  class Timer;
  class Vector;
  struct segment;

  using QuadraturePtr = std::shared_ptr<Quadrature>; // 含义：给 std::shared_ptr<Quadrature> 起个短名字叫 QuadraturePtr。
                                                     // 作用：简化代码书写，提高可读性。std::shared_ptr 是 C++ 11 引入的智能指针，用于自动管理内存，当没有任何地方引用该对象时，自动释放内存。

  /** Comparator for sorting k-nearest stencil std::pair objects 升序排序
   * CMFD 算法中有一个步骤是寻找某个网格周围“最近的 K 个邻居”（K-Nearest Neighbors），用于插值更新通量。这些邻居通常存储为 (ID, 距离) 的对子。
   * 作用：这是一个自定义的比较函数，用于 std::sort 等排序算法。它告诉排序算法：比较两个邻居时，只看它们的距离（second 元素），距离小的排在前面（升序）
   */
  inline bool
  stencilCompare(const std::pair<int, double> &firstElem,
                 const std::pair<int, double> &secondElem)
  {
    return firstElem.second < secondElem.second;
  }

#undef track_flux

/** Indexing macro for the angular fluxes for each polar angle and energy
 *  group for either the forward or reverse direction for a given Track
 * 为给定轨道的每个极角和正向或反向能量组的角通量建立索引宏
 *
 * */
#define track_flux(p, e) (track_flux[(p) * _num_moc_groups + (e)])

  /**
   * @class Cmfd Cmfd.h "src/Cmfd.h"
   * @brief A class for Coarse Mesh Finite Difference (CMFD) acceleration.
   */
  class Cmfd
  {

  private:
    /** Pointer to polar quadrature object 求积组不是一组包含方位角和极角的组合以及他们的权值，这里为什么说是极角求积组？*/
    QuadraturePtr _quadrature;

    /** Pointer to geometry object
     *
     */
    Geometry *_geometry;

    /** The keff eigenvalue */
    double _k_eff;

    /** The A (destruction) matrix 对应CMFD文档中的B.29公式
     * 损失矩阵
     */
    Matrix *_A;

    /** The M (production) matrix */
    Matrix *_M;

    /** The old source vector B.32，CMFD求解前的通量*/
    Vector *_old_source;

    /** The new source vector */
    Vector *_new_source;

    /* ==================== 并行计算：域边界通信缓冲区（MPI多进程通信用）==================== */
    /*
     * 【背景知识】MPI域分解并行计算：
     * - MPI将整个计算区域切分成多个子区域（domain），每个进程负责一个子区域
     * - 相邻子区域的边界处需要交换数据（称为"边界通信"）
     * 发送端进程               接收端进程
     *   ↓                       ↓
     * 准备发送缓冲区 → MPI_Send → MPI_Recv → 存入接收缓冲区
     * (_send_xxx)            (_boundary_xxx)
     *
     * 【三层指针含义】CMFD_PRECISION ***变量 等价于 三维数组
     * - 第1维：面的索引（如X方向左面、右面，Y方向前面、后面，Z方向上面、下面，共6个面）
     * - 第2维：该面上的边界CMFD网格索引（如X面的网格数 = _local_num_y × _local_num_z）
     * - 第3维：能群索引（ncg = num_cmfd_groups，粗网有限差分的能量组数）
     */

    // -------------------- 接收缓冲区：存储从相邻进程收到的边界数据 --------------------
    CMFD_PRECISION ***_boundary_volumes;          // 边界单元的体积（或2D情况下的面积）
    CMFD_PRECISION ***_boundary_reaction;         // 边界单元的反应率（吸收、散射等核反应截面）
    CMFD_PRECISION ***_boundary_diffusion;        // 边界单元的扩散系数（描述中子扩散快慢）
    CMFD_PRECISION ***_old_boundary_flux;         // 边界单元的上一次迭代通量值（用于迭代收敛判断）
    CMFD_PRECISION ***_boundary_surface_currents; // 边界表面的中子流（净流入/流出量）

    // -------------------- 发送缓冲区：准备发送给相邻进程的本域边界数据 --------------------
    CMFD_PRECISION ***_send_volumes;   // 待发送：本域边界单元的体积
    CMFD_PRECISION ***_send_reaction;  // 待发送：本域边界单元的反应率
    CMFD_PRECISION ***_send_diffusion; // 待发送：本域边界单元的扩散系数
    CMFD_PRECISION ***_send_currents;  // 待发送：本域边界表面的中子流

    // -------------------- 分裂中子流特殊处理缓冲区（处理网格角点和边的中子流分配）--------------------
    /*
     *split current = 分裂中子流（角点/边上的中子流量按比例分配给相邻网格）
     *surface current = 表面中子流（通过网格表面的中子净流量）
     *net current = 净中子流（流入 - 流出）
     * 【分裂中子流（split current）概念】：
     * - 在直角网格中：
     *   · 角点：被4个单元共享（2D）或8个单元共享（3D）
     *   · 边：被2个单元共享（2D）或4个单元共享（3D）
     *   · 面：被2个单元共享（仅3D有面的概念）
     * - 角点/边上的中子流需要"分裂"按比例分配给相邻单元，避免重复计算或遗漏
     * - 例如2D角点：中子流 = 25%给左上 + 25%给右上 + 25%给左下 + 25%给右下
     * - 例如3D边：中子流 = 25%分配给围绕该边的4个单元
     */
    CMFD_PRECISION *_send_split_current_data;       // 一维数组：打包所有待发送的分裂中子流（连续内存便于MPI传输）
    CMFD_PRECISION *_receive_split_current_data;    // 一维数组：接收所有分裂中子流的连续缓冲区
    CMFD_PRECISION **_send_split_currents_array;    // 二维数组：将上述一维数组重新组织成按面索引的格式
    CMFD_PRECISION **_receive_split_currents_array; // 二维数组：将接收数据按面索引组织
    CMFD_PRECISION ***_off_domain_split_currents;   // 三维数组：存储本域外（相邻域）传来的分裂中子流
    CMFD_PRECISION ***_received_split_currents;     // 三维数组：最终整理好的接收到的分裂中子流数据

    /* A tally buffer for threads to tally temporary surface currents
     * [OpenMP并行计算] 线程私有缓冲区
     * - 为了避免多个线程同时写入同一个变量导致冲突（数据竞争），给每个线程分配独立的内存空间
     * - 二维数组：第一维是线程ID，第二维是该线程计算的数据
     * - 计算完成后，再将所有线程的数据汇总
     *  每个线程提供一个独立的缓冲区，用于记录每个线程计算出的临时表面电流值
     */
    CMFD_PRECISION **_temporary_currents;

    /** Vector representing the flux for each cmfd cell and cmfd enegy group at
     * the end of a CMFD solve 
     * 表示cmfd求解结束时每个cmfd单元和cmfd能量组的通量的向量 新通量
     * */
    Vector *_new_flux;

    /** Vector representing the flux for each cmfd cell and cmfd enegy group at
     * the beginning of a CMFD solve CMFD求
     * 解前的通量向量(由moc数据累加平均得到) 旧通量
     * */
    Vector *_old_flux;

    /** The corrected diffusion coefficients from the previous iteration 上一次迭代的表面修正扩散系数(B.26)*/
    Vector *_old_dif_surf_corr;

    /** Whether the old diffusion coefficient has been set
     * 标志位：是否已经存在表面修正扩散系数（第一步迭代时通常为false）
     */
    bool _old_dif_surf_valid;

    /** Gauss-Seidel SOR relaxation factor
     * [数值计算] SOR（逐次超松弛）迭代法的松弛因子
     * - 值通常在 (1, 2) 之间
     * - 作用：加快线性方程组求解的收敛速度
     */
    double _SOR_factor;

    /** cmfd source convergence threshold
     * 源收敛阈值：当两次迭代的源项差异小于此值时，认为计算收敛，停止迭代
     */
    double _source_convergence_threshold;

    /** Number of cells in x direction
     * 整个几何X方向上的CMFD粗网格数量
     */
    int _num_x;

    /** Number of cells in y direction */
    int _num_y;

    /** Number of cells in z direction */
    int _num_z;

    /** Number of radial tiles 径向六边形数量层数，一层1个，2层6个，三层12个*/
    int _num_r;

    /** Number of energy groups
     * MOC使用的能群总数
     */
    int _num_moc_groups;

    /** Number of polar angles
     * 角度离散参数：极角数量
     */
    int _num_polar;

    /** Number of azimuthal angles
     * 角度离散参数：方位角数量
     */
    int _num_azim;

    /** Number of energy groups used in cmfd solver. Note that cmfd supports
     * energy condensation from the MOC
     * CMFD（粗网计算）使用的能群数量
     * 通常 _num_cmfd_groups < _num_moc_groups，即进行了"能群归并"（Energy Condensation）以加速计算
     */
    int _num_cmfd_groups;

    /** Number of energy groups used in hex cmfd solver. Note that cmfd supports
     * energy condensation from the MOC
     * 六边形CMFD求解器的能群数量
     */
    int _hex_num_groups;

    /** Coarse energy indices for fine energy groups
     * 能群映射索引数组
     * - 数组大小 = CMFD能群数 + 1
     * - 定义了MOC细能群如何归并到CMFD粗能群的边界
     * - 例如：CMFD第1群包含MOC第0~5群，则 indices[0]=0, indices[1]=6
     */
    int *_group_indices;

    /** Map of MOC groups to CMFD groups
     * 详细映射表：数组下标是MOC能群ID，值是对应的CMFD能群ID
     * map[moc_group_id] = cmfd_group_id
     * 给定一个 MOC 细群 ID，快速查到它属于哪个 CMFD 粗群
     */
    int *_group_indices_map;

    /** Number of energy groups in the backup CMFD solver
     * 备用的 CMFD 求解器中使用的能群数
     */
    int _num_backup_groups;

    /** Map of MOC groups to backup CMFD group structure
     * 二维向量：定义MOC能群到备用求解器能群的映射结构
     */
    std::vector<std::vector<int>> _backup_group_structure;

    /** Map of CMFD groups to backup CMFD group structure
     * 用于将主 CMFD 求解器中的 CMFD 能群映射到备用 CMFD 求解器中的能群的映射数组
     */
    int *_cmfd_group_to_backup_group;

    /** If the user specified fine-to-coarse group indices
     * 标志位：用户是否在输入文件中显式指定了能群归并关系
     * 即细化能量组到粗化能量组的映射关系
     */
    bool _user_group_indices;

    /** If a linear source approximation is used
     * 标志位：是否使用线性源近似
     * - true: 线性源近似（源在网格内线性变化，精度高）
     */
    bool _linear_source;

    /** If diffusion coefficients are limited by the flux
     * 是否使用通量限制来限制扩散系数
     * - 防止在通量梯度极大的地方计算出非物理的扩散系数
     * - 增强数值稳定性
     */
    bool _flux_limiting;

    /** Whether to rebalance the computed sigma-t to be consistent with the MOC
     *  solution on every sweep
     * 是否在每次进行扫描计算时，程序会根据最新的 MOC 解调整总截面(Sigma-t)值，
     * 确保总截面与中子通量分布保持一致性。这可以提高计算的准确性和稳定性
     *  */
    bool _balance_sigma_t; //

    /** Number of FSRs 是整个几何的FSR还是一个CMFD中所有的FSR，是整个几何的FSR数量*/
    long _num_FSRs;

    /** The volumes (areas) for each FSR
     * 数组：存储每个FSR的体积（3D）或面积（2D）
     * FP_PRECISION 是浮点精度宏（通常是 float 或 double）
     */
    FP_PRECISION *_FSR_volumes;

    /** Pointers to Materials for each FSR
     * 指针数组：存储每个FSR对应的材料对象指针
     * _FSR_materials[fsr_id] 返回对应FSR的Material指针
     */
    Material **_FSR_materials;

    /** The FSR scalar flux in each energy group
     * 数组：存储每个FSR在每个能群的标通量
     * 大小通常是 _num_FSRs * _num_moc_groups
     */
    FP_PRECISION *_FSR_fluxes;

    /** The FSR source in each energy group
     * 数组：存储每个FSR的源项（裂变源 + 散射源）
     */
    FP_PRECISION *_FSR_sources;

    /** The source region flux moments (x, y, and z) for each energy group 存储每个能量组在源区域中的通量矩*/
    FP_PRECISION *_flux_moments;

    /** Array of CMFD cell volumes
     * 向量：存储每个CMFD粗网格的体积
     */
    Vector *_volumes;

    /** Array of material pointers for CMFD cell materials CMFD材料数组（二维指针数组）*/
    Material **_materials;

    /** Physical dimensions of the geometry and each CMFD cell
     * 几何体的物理尺寸（总长、总宽、总高）
     */
    double _width_x;
    double _width_y;
    double _width_z;
    double _width_r;
    DoubleVec _widths_z; // Z方向各层的宽度向量

    /** 一个CMFD网格的长宽高
     * 均匀网格情况下的单个网格尺寸
     */
    double _cell_width_x;
    double _cell_width_y;
    double _cell_width_z;

    /** Distance of each mesh from the left-lower-bottom most point
     * _accumulate_x[i] 表示第i个网格右边界距离晶格最左边界的距离
     */
    std::vector<double> _accumulate_x;
    std::vector<double> _accumulate_y;
    std::vector<double> _accumulate_z;

    /** Physical dimensions of non-uniform CMFD meshes (for whole geometry)
     * 非均匀网格尺寸数组：存储每个网格的具体尺寸
     */
    std::vector<double> _cell_widths_x;
    std::vector<double> _cell_widths_y;
    std::vector<double> _cell_widths_z;

    /** True if the cmfd meshes are non-uniform
     * 标志位：是否使用非均匀CMFD网格
     */
    bool _non_uniform;

    /** True if the cmfd mesh has been adjusted to fit the domain decomposition
     * 标志位：网格是否为了适应MPI域分解而进行了调整
     */
    bool _widths_adjusted_for_domains;

    /** Array of geometry boundaries
     * 边界条件数组：存储几何体各个面的边界类型（如反射、真空等）
     */
    boundaryType *_boundaries;

    /** Array of surface currents for each CMFD cell
     * Vector数组：存储每个CMFD单元各个表面的中子流
     */
    Vector *_surface_currents;

    /** Array of total current from starting boundary fluxes
     * Vector数组：存储从边界入射通量产生的总中子流
     */
    Vector *_starting_currents;

    /** Array of net currents of all CMFD cells
     * Vector数组：存储每个CMFD单元的净中子流（Net Current = 流出 - 流入）
     */
    Vector *_net_currents;

    /** Array of surface currents on all faces + edges and corners used in
        debugging
        调试用：存储所有面、边、角上的中子流数据
     */
    Vector *_full_surface_currents;

    /** Array of surface currents on edges and corners for each CMFD cell
     * 映射表：存储边和角上的中子流（用于更精细的分析）
     * std::map<int, ...> 键是边/角的唯一ID
     */
    std::map<int, CMFD_PRECISION> _edge_corner_currents;

    /** Vector of vectors of FSRs containing in each cell
     * 嵌套向量：建立 CMFD粗网格 -> FSR的索引关系
     * _cell_fsrs[cmfd_cell_id] 返回该粗网格内包含的所有FSR ID列表
     */
    std::vector<std::vector<long>> _cell_fsrs;

    /** vector for storing whether the CMFD cells are empty or not.
     * 布尔向量：标记CMFD网格是否为空（即不包含任何FSR，通常位于几何体外的填充区域）
     */
    std::vector<bool> _empty_fsrs_cells;
    int _empty_cells_num; // 空网格计数

    /** logical cmfd index to actual cmfd cell index从逻辑cmfd索引到实际cmfd单元索引*/
    std::vector<int> _logical_actual_map;

    /** Pointer to Lattice object representing the CMFD mesh
     * 指向代表整个 CMFD 网格系统的 Lattice 对象的指针
     */
    Lattice *_lattice;

    /** Orientation of the lattice
     * CMFD lattice的方向
     */
    Orientation _orientation;

    /** Bool type, which is used to determine the shape of the Lattice
     * 标志位：是否启用六边形Lattice（Hexagonal Lattice）
     */
    bool _hexlattice_enable;

    /** Flag indicating whether to update the MOC flux
     * 标志位：是否需要 CMFD
     如果没有定义 CMFD，或者用户显式关闭了通量更新（加速），就直接退出，不初始化了。
     */
    bool _flux_update_on;

    /** Flag indicating whether to us centroid updating
     * 标志位：是否使用质心更新策略（一种改进的通量更新方法）
     */
    bool _centroid_update_on;

    /** Flag indicating whether to check neutron balance on every CMFD solve
     * 标志位：是否在每次CMFD求解后检查中子守恒（调试和验证用）
     */
    bool _check_neutron_balance;

    /** Number of cells to used in updating MOC flux
     * 参数 K：用于K近邻插值更新MOC通量
     * 用于更新MOC通量的CMFD数量
     */
    int _k_nearest;

    /** Relaxation factor to use for corrected diffusion coefficients 用于修正扩散系数的松弛因子*/
    double _relaxation_factor;

    /** Map storing the k-nearest stencil for each fsr
     * 映射表：存储每个FSR的K个最近邻CMFD网格及其权重
     * - 用于将CMFD粗网格结果插值回FSR细网格
     * - pair<int, double> 表示 <网格ID, 权重>
     */
    std::map<int, std::vector<std::pair<int, double>>>
        _k_nearest_stencils;

    /** OpenMP mutual exclusion locks for atomic CMFD cell operations
     * [OpenMP并行] 互斥锁数组
     * - 每个CMFD网格一个锁
     * - 防止多个线程同时修改同一个网格的数据
     */
    omp_lock_t *_cell_locks;

    /** OpenMP mutual exclusion lock for edge/corner current tallies
     * [OpenMP并行] 边/角中子流计数的OpenMP互斥锁
     */
    omp_lock_t _edge_corner_lock;

#ifndef THREED
    /** Flag indicating whether the problem is 2D or 3D
     * 标志位：指示当前是2D还是3D计算（当未定义THREED宏时使用）
     */
    bool _SOLVE_3D;
#endif

    /** Array of azimuthal track spacings
     * 几何参数：方位角特征线的间距
     */
    double *_azim_spacings;

    /** 2D array of polar track spacings
     * 几何参数：极角特征线的间距（二维数组）
     */
    double **_polar_spacings;

    /** Whether to use axial interpolation for flux update ratios
     * 是否在更新通量比率时使用轴向插值
     * 0：不插值
     * 1：使用 FSR 在轴向上的平均值进行插值
     * 2：使用质心（centroid）的 z 坐标值进行插值
     */
    int _use_axial_interpolation;

    /** Axial interpolation constants
     * 预计算的轴向插值常数
     */
    std::vector<double *> _axial_interpolants;

    /* Structure to contain information about the convergence of the CMFD
     * 结构体指针：存储CMFD收敛过程数据（如残差历史、迭代次数）
     */
    ConvergenceData *_convergence_data;

    /* MPI communicator to transfer buffers, mainly currents at interfaces
     * [MPI并行] 域通信器对象
     * - 负责管理跨进程的数据传输
     * - 处理域分解边界上的中子流交换
     */
    DomainCommunicator *_domain_communicator;

    /* Buffer to contain received data
     * [MPI并行] 接收缓冲区：存储从其他进程接收的数据
     */
    CMFD_PRECISION *_inter_domain_data;

    /* Buffer to contain sent data from domain
     * [MPI并行] 发送缓冲区：存储准备发送给其他进程的数据
     */
    CMFD_PRECISION *_send_domain_data;

    /* For each face (1st dimension of the array), will contain data received
     * [MPI并行] 二维数组，对于每个面（数组的第一个维度），将包含接收到的数据
     */
    CMFD_PRECISION **_domain_data_by_surface;

    /* For each face (1st dimension of the array), will contain data to send
     * [MPI并行] 二维数组，对于每个面（数组的第一个维度），将包含准备发送的数据
     */
    CMFD_PRECISION **_send_data_by_surface;

    /* Map of the indexes to each boundary in the tally arrays
     * 映射表：边界索引映射，用于快速查找边界数据在数组中的位置
     */
    std::vector<std::map<int, int>> _boundary_index_map;

    /* The number of on-domain cells in the x-direction
     * 本地（当前进程负责的）子区域在X方向的网格数
     一个domain域中x方向上CMFD单元数
     */
    int _local_num_x;

    /* The number of on-domain cells in the y-direction */
    int _local_num_y;

    /* The number of on-domain cells in the z-direction */
    int _local_num_z;

    /* Size of _tally_memory array
     * 统计内存总大小
     */
    long _total_tally_size;

    /* 1D array that contains all tallies (diffusion, reaction and volume)
     */
    CMFD_PRECISION *_tally_memory;

    /* 2D array that contains reaction rates in each cell and group
     * 指针数组：指向 _tally_memory 中的反应率数据段
     */
    CMFD_PRECISION **_reaction_tally;

    /* 2D array that contains volume tallies of each cell
     * 指针数组：指向 _tally_memory 中的体积数据段
     */
    CMFD_PRECISION **_volume_tally;

    /* 2D array that contains diffusion tallies for each cell and groups
     * 指针数组：指向 _tally_memory 中的扩散系数数据段
     */
    CMFD_PRECISION **_diffusion_tally;

    /* Boolean to check if tallies are allocated
     * 标志位：内存是否已分配
     */
    bool _tallies_allocated;

    /* Boolean to check if the domain communicator (for domain decomposed CMFD)
     * has been allocated
     * 标志位：MPI通信器是否已初始化
     */
    bool _domain_communicator_allocated;

    /** A timer to record timing data for a simulation
     * 计时器对象：用于性能分析
     */
    Timer *_timer;

    /** A one-group backup CMFD solver
     * 备用CMFD求解器指针（通常用于加速主求解器的收敛，或作为预处理）
     */
    Cmfd *_backup_cmfd;

    /* Private worker functions */

    /* 计算Larsen有效扩散系数因子
     * - 用于修正扩散系数，使其更符合输运理论的结果
     * - dif_coef: 原始扩散系数
     * - delta: 网格尺寸
     */
    CMFD_PRECISION computeLarsensEDCFactor(CMFD_PRECISION dif_coef,
                                           CMFD_PRECISION delta);

    /* 构建CMFD线性方程组的矩阵 (Ax = b)
     * - moc_iteration: 当前MOC迭代步数
     * - 计算矩阵系数 A 和 M
     */
    void constructMatrices(int moc_iteration);
    void hexConstructMatrices(int moc_iteration); // 六边形几何版本

    /* 归并截面 (Fine -> Coarse)
     * - 将MOC细网格的截面加权平均，得到CMFD粗网格的截面
     * - 权重通常是通量和体积
     */
    void collapseXS();
    void hexCollapseXS(); // 六边形版本

    /* 更新MOC通量 (Coarse -> Fine)
     * - 将CMFD计算得到的粗网格通量变化，调制回MOC细网格通量
     * - 实现加速收敛的关键步骤
     */
    void updateMOCFlux();
    void updateHexMOCFlux(); // 六边形版本

    /* 重新缩放通量
     * - 确保全堆总中子产生率（或总功率）守恒
     */
    void rescaleFlux();

    /* 分裂顶点中子流
     * - 处理角点处的中子流分配问题
     */
    void splitVertexCurrents();
    void splitVertexCurrentsHex();

    /* 分裂边中子流
     * - 处理边处的中子流分配问题
     */
    void splitEdgeCurrents();
    void splitEdgeCurrentsHex();

    /* 获取顶点分裂关联的表面列表
     * - 辅助函数：确定一个顶点连接了哪些表面
     */
    void getVertexSplitSurfaces(int cell, int vertex, std::vector<int> *surfaces);
    void getVertexSplitSurfacesHex(int cell, int vertex, std::vector<int> *surfaces);

    /* 获取边分裂关联的表面列表 */
    void getEdgeSplitSurfaces(int cell, int edge, std::vector<int> *surfaces);
    void getEdgeSplitSurfacesHex(int cell, int edge, std::vector<int> *surfaces);

    /* 初始化材料属性 */
    void initializeMaterials();
    void initializeHexMaterials();

    /* 初始化中子流数组 */
    void initializeCurrents();
    void initializeCurrentsHex();

    /* 生成K近邻模板
     * - 预计算每个FSR周围最近的K个CMFD网格，用于插值
     */
    void generateKNearestStencils();
    void generateKNearestStencilsHex();

    /* 方向与表面索引转换工具函数 */
    int convertDirectionToSurface(int *direction);
    void convertSurfaceToDirection(int surface, int *direction);
    void convertSurfaceToDirectionHex(int surface, int *direction);
    std::string getSurfaceNameFromDirection(int *direction);
    std::string getSurfaceNameFromSurface(int surface);

    /* Private getter functions */
    /* 获取相邻网格ID
     * - global=true: 返回全局ID
     * - neighbor=true: 即使在不同域也返回ID
     */
    int getCellNext(int cell_id, int surface_id, bool global = true,
                    bool neighbor = false);

    /* 根据模板获取网格ID */
    int getCellByStencil(int cell_id, int stencil_id);
    int getHexCellByStencil(int cell_id, int stencil_id);
    int getHexCellByStencilPre(int cell_id, int stencil_id);

    /* 获取通量比率 (新通量/旧通量) */
    CMFD_PRECISION getFluxRatio(int cell_id, int group, int fsr);
    CMFD_PRECISION getHexFluxRatio(int cell_id, int group, int fsr);

    /* 获取更新比率 (用于MOC更新) */
    CMFD_PRECISION getUpdateRatio(int cell_id, int moc_group, int fsr);
    CMFD_PRECISION getHexUpdateRatio(int cell_id, int moc_group, int fsr);

    /* 计算到质心的距离 (用于插值权重计算) */
    double getDistanceToCentroid(Point *centroid, int cell_id,
                                 int stencil_index);
    double getDistanceToCentroidHex(Point *centroid, int cell_id,
                                    int stencil_index);
    double getDistanceToCentroidHexPre(Point *centroid, int cell_id,
                                       int stencil_index);

    /* 获取表面扩散系数
     * - correction=true: 返回经过CMFD修正后的系数
     */
    CMFD_PRECISION getSurfaceDiffusionCoefficient(int cmfd_cell, int surface,
                                                  int group, int moc_iteration,
                                                  bool correction);
    CMFD_PRECISION getHexSurfaceDiffusionCoefficient(int cmfd_cell, int surface,
                                                     int group, int moc_iteration,
                                                     bool correction);

    /* 获取网格中心扩散系数 */
    CMFD_PRECISION getDiffusionCoefficient(int cmfd_cell, int group);

    /* 获取表面几何宽度 */
    CMFD_PRECISION getSurfaceWidth(int surface);
    CMFD_PRECISION getHexSurfaceWidth(int surface);
    CMFD_PRECISION getPerpendicularSurfaceWidth(int surface);
    CMFD_PRECISION getHexPerpendicularSurfaceWidth(int surface);

    /* 获取表面方向 (正向/负向) */
    int getSense(int surface);

    /* ID转换函数 */
    int getLocalCMFDCell(int cmfd_cell);  // 全局ID -> 本地ID
    int getGlobalCMFDCell(int cmfd_cell); // 本地ID -> 全局ID
    int getCellColor(int cmfd_cell);      // 获取网格颜色(用于红黑迭代等并行算法)

    /* 打包MPI通信缓冲区 */
    void packBuffers();

#ifdef ENABLE_MPI_
    /* 幽灵网格交换 (Ghost Cell Exchange)
     * - MPI通信核心函数
     * - 交换相邻域边界上的网格数据
     */
    void ghostCellExchange();

    /* 通信分裂中子流
     * - faces=true: 交换面上的数据
     * - faces=false: 交换边/角上的数据
     */
    void communicateSplits(bool faces);
#endif

    /* 解包分裂中子流数据 */
    void unpackSplitCurrents(bool faces);

    /* 复制完整表面中子流 (用于调试或输出) */
    void copyFullSurfaceCurrents();
    void copyHexFullSurfaceCurrents();

    /* 检查中子守恒 (Neutron Balance Check)
     * - 验证流入+产生 = 流出+吸收
     * - pre_split: 是否在分裂电流处理前检查
     */
    void checkNeutronBalance(bool pre_split = true);
    void checkNeutronBalanceHex(bool pre_split = true);

    /* 打印延长因子 (Prolongation Factors)
     * - 调试用：查看粗网到细网的插值系数
     */
    void printProlongationFactors(int iteration);

  public:
    Cmfd();
    virtual ~Cmfd();

    /* Worker functions */
    double computeKeff(int moc_iteration);
    void initialize();
    void initializeHex();
    void initializeCellMap();
    void initializeGroupMap();
    void allocateTallies();
    void allocateHexTallies();
    void initializeLattice(Point *offset);
    void initializeBackupCmfdSolver();
    void copyCurrentsToBackup();
    int findCmfdCell(LocalCoords *coords);
    int findCmfdSurface(int cell_id, LocalCoords *coords, double azim, double polar);
    int findCmfdSurfaceOTF(int cell_id, double z, int surface_2D);
    void addFSRToCell(int cell_id, long fsr_id);
    void printCellFSRs();
    void zeroCurrents();
    void zeroHexCurrents();
    void tallyCurrent(segment *curr_segment, float *track_flux,
                      int azim_index, int polar_index, bool fwd);
    void tallyHexCurrent(segment *curr_segment, float *track_flux,
                         int azim_index, int polar_index, bool fwd);
    void tallyStartingCurrent(Point *point, double delta_x, double delta_y,
                              double delta_z, float *track_flux, double weight);
    void recordNetCurrents();
    void printInputParamsSummary();
    void printTimerReport();
    void checkBalance();

    /* Get parameters */
    int getNumCmfdGroups();
    int getNumMOCGroups();
    int getNumCells();
    int getCmfdGroup(int group);
    int getBoundary(int side);
    Lattice *getLattice();
    int getNumX();
    int getNumY();
    int getNumZ();
    Vector *getLocalCurrents();
    CMFD_PRECISION ***getBoundarySurfaceCurrents();
    int convertFSRIdToCmfdCell(long fsr_id);
    int convertGlobalFSRIdToCmfdCell(long global_fsr_id);
    std::vector<std::vector<long>> *getCellFSRs();
    bool isFluxUpdateOn();
    bool isCentroidUpdateOn();
    bool isSigmaTRebalanceOn();

    /* Set parameters */
    void setSORRelaxationFactor(double SOR_factor);
    void setCMFDRelaxationFactor(double relaxation_factor);
    void setGeometry(Geometry *geometry);
    void setWidthX(double width);
    void setWidthY(double width);
    void setWidthZ(double width);
    void setWidthsZ(DoubleVec widths);
    void setWidthR(double width);
    void setNumX(int num_x);
    void setNumY(int num_y);
    void setNumZ(int num_z);
    void setNumR(int num_r);
    void setNumFSRs(long num_fsrs);
    void setNumMOCGroups(int num_moc_groups);
    void setBoundary(int side, boundaryType boundary);
    void setLatticeStructure(int num_x, int num_y, int num_z = 1);
    void setFluxUpdateOn(bool flux_update_on);
    void setCentroidUpdateOn(bool centroid_update_on);
    void setGroupStructure(std::vector<std::vector<int>> group_indices);
    void setSourceConvergenceThreshold(double source_thresh);
    void setQuadrature(QuadraturePtr quadrature);
    void setKNearest(int k_nearest);
    void setSolve3D(bool solve_3d);
    void setAzimSpacings(const std::vector<double> &azim_spacings,
                         int num_azim);
    void setPolarSpacings(const std::vector<std::vector<double>> &
                              polar_spacings,
                          int num_azim, int num_polar);
    void setKeff(double k_eff);
    void setBackupGroupStructure(std::vector<std::vector<int>> group_indices);
    void setOrientation(Orientation orientation);
    void setOrientation(std::string);

#ifdef ENABLE_MPI_
    void setNumDomains(int num_x, int num_y, int num_z);
    void setDomainIndexes(int idx_x, int idx_y, int idx_z);
#endif
    void setConvergenceData(ConvergenceData *convergence_data);
    void useAxialInterpolation(int interpolate);

    /* Methods to try to fix stability issues */
    void useFluxLimiting(bool flux_limiting);
    void enforceBalanceOnDiagonal(int cmfd_cell, int group);
    void rebalanceSigmaT(bool balance_sigma_t);

    /* Set FSR parameters */
    void setFSRMaterials(Material **FSR_materials);
    void setFSRVolumes(FP_PRECISION *FSR_volumes);
    void setFSRFluxes(FP_PRECISION *scalar_flux);
    void setFSRSources(FP_PRECISION *sources);
    void setCellFSRs(std::vector<std::vector<long>> *cell_fsrs);
    void setFluxMoments(FP_PRECISION *flux_moments);

    /* Set XYZ widths of non-uniform meshes */
    void setWidths(std::vector<std::vector<double>> widths);

    /* Set the bool _hexlattice_enable */
    void setHexLatticeEnable(bool hexlattice_enable);

    /* Set _hex_num_groups */
    void setHexGroups(int hex_num_groups)
    {
      _hex_num_groups = hex_num_groups;
    }
    int getHexGroups()
    {
      return _hex_num_groups;
    }

    /* Return _hexlattice_enable */
    bool GetHexLatticeEnable();

    /// \brief Returns the orientation
    Orientation getOrientation() const
    {
      return _orientation;
    }

    void convertVertexToSurfaces(int surface, int *remainder_surfaces, int *partial_surfaces);

    void convertEdgeToSurfaces(int surface, int *partial_surfaces);

    void findEmptyCmfdCells();

    void printHex();

    bool CellinHexLattice(int cell);
    bool CellinXYZBoundary(int x, int y, int z, int surface);
    bool CellinXYZBoundary(int cell, int surface);
    bool CellinXYBoundary(int cell);
    bool CellinXYBoundaryWithStencil(int cell, int stencil_index);

    void GetXYSurfaces(int *cmfd_surfaces);
  };

} /* namespace antmoc */

#endif /* CMFD_H_ */
