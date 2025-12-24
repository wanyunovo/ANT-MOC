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

    /** Pointer to geometry object */
    Geometry *_geometry;

    /** The keff eigenvalue */
    double _k_eff;

    /** The A (destruction) matrix 对应CMFD文档中的B.29公式*/
    Matrix *_A;

    /** The M (production) matrix */
    Matrix *_M;

    /** The old source vector B.32，CMFD求解前的通量*/
    Vector *_old_source;

    /** The new source vector */
    Vector *_new_source;

    /* Domain boundary communication buffers    15个变量*/
    CMFD_PRECISION ***_boundary_volumes;
    CMFD_PRECISION ***_boundary_reaction;
    CMFD_PRECISION ***_boundary_diffusion;
    // 三维大小:NUM_FACES(面数),某一个面对应的边界CMFD网格数(比如x面的CMFD网格数即为:_local_num_y * _local_num_z),ncg
    // storage_per_cell * num_boundary_cells + ncg * num_boundary_cells;其中storage_per_cell=((2 + NUM_FACES) * ncg + 1)
    // num_boundary_cells=所有面的CMFD网格数,
    CMFD_PRECISION ***_old_boundary_flux;
    CMFD_PRECISION ***_boundary_surface_currents;

    CMFD_PRECISION ***_send_volumes;
    CMFD_PRECISION ***_send_reaction;
    CMFD_PRECISION ***_send_diffusion;
    CMFD_PRECISION ***_send_currents;

    CMFD_PRECISION *_send_split_current_data;
    CMFD_PRECISION *_receive_split_current_data;
    CMFD_PRECISION **_send_split_currents_array;
    CMFD_PRECISION **_receive_split_currents_array;
    CMFD_PRECISION ***_off_domain_split_currents;
    CMFD_PRECISION ***_received_split_currents;

    /* A tally buffer for threads to tally temporary surface currents 每个线程提供一个独立的缓冲区，用于记录每个线程计算出的临时表面电流值*/
    CMFD_PRECISION **_temporary_currents;

    /** Vector representing the flux for each cmfd cell and cmfd enegy group at
     * the end of a CMFD solve 表示cmfd求解结束时每个cmfd单元和cmfd能量组的通量的向量*/
    Vector *_new_flux;

    /** Vector representing the flux for each cmfd cell and cmfd enegy group at
     * the beginning of a CMFD solve CMFD求解前的网格通量(由moc数据累加平均得到)*/
    Vector *_old_flux;

    /** The corrected diffusion coefficients from the previous iteration 上一次迭代的表面修正扩散系数*/
    Vector *_old_dif_surf_corr;

    /** Whether the old diffusion coefficient has been set 是否已经有了表面修正扩散系数*/
    bool _old_dif_surf_valid;

    /** Gauss-Seidel SOR relaxation factor 高斯-赛德尔 SOR 松弛因子*/
    double _SOR_factor;

    /** cmfd source convergence threshold cmfd源收敛阈值*/
    double _source_convergence_threshold;

    /** Number of cells in x direction 整个几何上x方向上CMFD网格的数量*/
    int _num_x;

    /** Number of cells in y direction */
    int _num_y;

    /** Number of cells in z direction */
    int _num_z;

    /** Number of radial tiles 径向六边形数量层数，一层1个，2层6个，三层12个*/
    int _num_r;

    /** Number of energy groups */
    int _num_moc_groups;

    /** Number of polar angles */
    int _num_polar;

    /** Number of azimuthal angles */
    int _num_azim;

    /** Number of energy groups used in cmfd solver. Note that cmfd supports
     * energy condensation from the MOC */
    int _num_cmfd_groups;

    /** Number of energy groups used in hex cmfd solver. Note that cmfd supports
     * energy condensation from the MOC */
    int _hex_num_groups;

    /** Coarse energy indices for fine energy groups 用于存储CMFD的索引对应全部MOC能群数的第几个(下标从0开始,左闭右开),数组大小为CMFD能群数+1*/
    int *_group_indices;

    /** Map of MOC groups to CMFD groups 用于将细化MOC能量组映射到对应的粗化CMFD能量组，这个映射可以是一对一的，也可以是一对多的，取决于细化和粗化能量组之间的关系。*/
    int *_group_indices_map;

    /** Number of energy groups in the backup CMFD solver 备用的 CMFD 求解器中使用的能量组的数量*/
    int _num_backup_groups;

    /** Map of MOC groups to backup CMFD group structure 是一个二维 std::vector，其元素类型为 int，用于将 MOC 组映射到备用 CMFD组的结构或布局。*/
    std::vector<std::vector<int>> _backup_group_structure;

    /** Map of CMFD groups to backup CMFD group structure 用于将主 CMFD 求解器中的 CMFD 能群映射到备用 CMFD 求解器中的能群结构*/
    int *_cmfd_group_to_backup_group;

    /** If the user specified fine-to-coarse group indices 用户是否已经指定了细化能量组到粗化能量组的映射关系*/
    bool _user_group_indices;

    /** If a linear source approximation is used 是否使用线性源近似*/
    bool _linear_source;

    /** If diffusion coefficients are limited by the flux 是否使用通量限制来限制扩散系数*/
    bool _flux_limiting;

    /** Whether to rebalance the computed sigma-t to be consistent with the MOC
     *  solution on every sweep
     * 是否在每次进行扫描计算时，程序会根据最新的 MOC 解调整总截面值，
     * 确保总截面与中子通量分布保持一致性。这可以提高计算的准确性和稳定性
     *  */
    bool _balance_sigma_t; //

    /** Number of FSRs 是整个几何的FSR还是一个CMFD中所有的FSR，是整个几何的*/
    long _num_FSRs;

    /** The volumes (areas) for each FSR */
    FP_PRECISION *_FSR_volumes;

    /** Pointers to Materials for each FSR */
    Material **_FSR_materials;

    /** The FSR scalar flux in each energy group */
    FP_PRECISION *_FSR_fluxes;

    /** The FSR source in each energy group */
    FP_PRECISION *_FSR_sources;

    /** The source region flux moments (x, y, and z) for each energy group 存储每个能量组在源区域中的通量矩*/
    FP_PRECISION *_flux_moments;

    /** Array of CMFD cell volumes */
    Vector *_volumes;

    /** Array of material pointers for CMFD cell materials CMFD材料数组（二维的）*/
    Material **_materials;

    /** Physical dimensions of the geometry and each CMFD cell 整个几何上X方向的宽度*/
    double _width_x;
    double _width_y;
    double _width_z;
    double _width_r;
    DoubleVec _widths_z;

    /** 一个CMFD网格的长宽高*/
    double _cell_width_x;
    double _cell_width_y;
    double _cell_width_z;

    /** Distance of each mesh from the left-lower-bottom most point */
    // 从晶格最左边界开始到每个网格边界的距离
    std::vector<double> _accumulate_x;
    std::vector<double> _accumulate_y;
    std::vector<double> _accumulate_z;

    /** Physical dimensions of non-uniform CMFD meshes (for whole geometry) */
    std::vector<double> _cell_widths_x;
    std::vector<double> _cell_widths_y;
    std::vector<double> _cell_widths_z;

    /** True if the cmfd meshes are non-uniform 判断是CMFD网格是均匀的还是非均匀的*/
    bool _non_uniform;

    /** True if the cmfd mesh has been adjusted to fit the domain decomposition */
    bool _widths_adjusted_for_domains;

    /** Array of geometry boundaries */
    boundaryType *_boundaries;

    /** Array of surface currents for each CMFD cell 存储每个CMFD单元格表面的中子流信息的Vector数组*/
    Vector *_surface_currents;

    /** Array of total current from starting boundary fluxes 存储从开始边界通量得到的总电流的Vector数组*/
    Vector *_starting_currents;

    /** Array of net currents of all CMFD cells 存储所有CMFD单元格的净电流的Vector数组*/
    Vector *_net_currents;

    /** Array of surface currents on all faces + edges and corners used in
        debugging 用于存储所有面、边和角上的表面电流的Vector数组*/
    Vector *_full_surface_currents;

    /** Array of surface currents on edges and corners for each CMFD cell */
    std::map<int, CMFD_PRECISION> _edge_corner_currents;

    /** Vector of vectors of FSRs containing in each cell 这个是二维vector数组，第一维确定是具体的哪个CMFD单元，第二个是单个CMFD中所关联的所有FSR，*/
    std::vector<std::vector<long>> _cell_fsrs;

    /** vector for storing whether the CMFD cells are empty or not. */
    std::vector<bool> _empty_fsrs_cells; // 判断CMFD中的FSR是否为空（在边界外）
    int _empty_cells_num;                // CMFD中空单元格数量

    /** logical cmfd index to actual cmfd cell index从逻辑cmfd索引到实际cmfd单元索引*/
    std::vector<int> _logical_actual_map;

    /** Pointer to Lattice object representing the CMFD mesh 指向Lattice对象的CMFD*/
    Lattice *_lattice;

    /** Orientation of the lattice CMFD晶格的方向*/
    Orientation _orientation;

    /** Bool type, which is used to determine the shape of the Lattice */
    bool _hexlattice_enable;

    /** Flag indicating whether to update the MOC flux
     * 在求解过程中，通常需要多次迭代来达到收敛状态。每次迭代中，
     * 可能需要更新MOC通量以反映新的物理状态或边界条件的变化。
     * 如果当前迭代需要更新MOC通量，则 _flux_update_on 设为 true。
     *
     */
    bool _flux_update_on;

    /** Flag indicating whether to us centroid updating 是否进行质心更新操作*/
    bool _centroid_update_on;

    /** Flag indicating whether to check neutron balance on every CMFD solve */
    bool _check_neutron_balance;

    /** Number of cells to used in updating MOC flux 用于更新MOC通量的CMFD数*/
    int _k_nearest;

    /** Relaxation factor to use for corrected diffusion coefficients 用于修正扩散系数的松弛因子*/
    double _relaxation_factor;

    /** Map storing the k-nearest stencil for each fsr 存储每个 FSR 的最近邻单元及其权重信息，用于更新 MOC 通量*/
    std::map<int, std::vector<std::pair<int, double>>>
        _k_nearest_stencils;

    /** OpenMP mutual exclusion locks for atomic CMFD cell operations  CMFD互斥锁*/
    omp_lock_t *_cell_locks;

    /** OpenMP mutual exclusion lock for edge/corner current tallies 边/角中子流计数的OpenMP互斥锁*/
    omp_lock_t _edge_corner_lock;

#ifndef THREED
    /** Flag indicating whether the problem is 2D or 3D 指示问题是二维还是三维的标志*/
    bool _SOLVE_3D;
#endif

    /** Array of azimuthal track spacings */
    double *_azim_spacings;

    /** 2D array of polar track spacings 极角轨迹间距的2D阵列*/
    double **_polar_spacings;

    /** Whether to use axial interpolation for flux update ratios
     * 是否在更新通量比率时使用轴向插值
     * 0：不使用插值。1：使用 FSR 在轴向上的平均值进行插值。2：使用质心（centroid）的 z 坐标值进行插值
     * */
    int _use_axial_interpolation;

    /** Axial interpolation constants */
    std::vector<double *> _axial_interpolants;

    /* Structure to contain information about the convergence of the CMFD 有关CMFD收敛性信息的（ConvergenceData类型的）结构*/
    ConvergenceData *_convergence_data;

    /* MPI communicator to transfer buffers, mainly currents at interfaces MPI通讯器（DomainCommunicator）结构，用于传输缓冲区，主要是CMFD界面处的中子流 */
    DomainCommunicator *_domain_communicator;

    /* Buffer to contain received data 包含接收数据的缓冲区*/
    CMFD_PRECISION *_inter_domain_data;

    /* Buffer to contain sent data from domain 包含域中发送数据的缓冲区*/
    CMFD_PRECISION *_send_domain_data;

    /* For each face (1st dimension of the array), will contain data received 对于每个面（阵列的第一个维度），将包含接收到的数据*/
    CMFD_PRECISION **_domain_data_by_surface;

    /* For each face (1st dimension of the array), will contain data to send //对于每个面（数组的第一个维度），将包含要发送的数据*/
    CMFD_PRECISION **_send_data_by_surface;

    /* Map of the indexes to each boundary in the tally arrays 界面的索引映射信息*/
    std::vector<std::map<int, int>> _boundary_index_map;

    /* The number of on-domain cells in the x-direction 一个domain域中x方向上CMFD单元数*/
    int _local_num_x;

    /* The number of on-domain cells in the y-direction */
    int _local_num_y;

    /* The number of on-domain cells in the z-direction */
    int _local_num_z;

    /* Size of _tally_memory array */
    long _total_tally_size;

    /* 1D array that contains all tallies (diffusion, reaction and volume) */
    CMFD_PRECISION *_tally_memory;

    /* 2D array that contains reaction rates in each cell and group */
    CMFD_PRECISION **_reaction_tally;

    /* 2D array that contains volume tallies of each cell */
    CMFD_PRECISION **_volume_tally;

    /* 2D array that contains diffusion tallies for each cell and groups */
    CMFD_PRECISION **_diffusion_tally;

    /* Boolean to check if tallies are allocated */
    bool _tallies_allocated;

    /* Boolean to check if the domain communicator (for domain decomposed CMFD)
     * has been allocated */
    bool _domain_communicator_allocated;

    /** A timer to record timing data for a simulation */
    Timer *_timer;

    /** A one-group backup CMFD solver 一个单能群的备用 CMFD求解器*/
    Cmfd *_backup_cmfd;

    /* Private worker functions */
    CMFD_PRECISION computeLarsensEDCFactor(CMFD_PRECISION dif_coef,
                                           CMFD_PRECISION delta);
    void constructMatrices(int moc_iteration);
    void hexConstructMatrices(int moc_iteration);
    void collapseXS();
    void hexCollapseXS();
    void updateMOCFlux();
    void updateHexMOCFlux();
    void rescaleFlux();
    void splitVertexCurrents();
    void splitVertexCurrentsHex();
    void splitEdgeCurrents();
    void splitEdgeCurrentsHex();
    void getVertexSplitSurfaces(int cell, int vertex, std::vector<int> *surfaces);
    void getVertexSplitSurfacesHex(int cell, int vertex, std::vector<int> *surfaces);
    void getEdgeSplitSurfaces(int cell, int edge, std::vector<int> *surfaces);
    void getEdgeSplitSurfacesHex(int cell, int edge, std::vector<int> *surfaces);
    void initializeMaterials();
    void initializeHexMaterials();
    void initializeCurrents();
    void initializeCurrentsHex();
    void generateKNearestStencils();
    void generateKNearestStencilsHex();
    int convertDirectionToSurface(int *direction);
    void convertSurfaceToDirection(int surface, int *direction);
    void convertSurfaceToDirectionHex(int surface, int *direction);
    std::string getSurfaceNameFromDirection(int *direction);
    std::string getSurfaceNameFromSurface(int surface);

    /* Private getter functions */
    int getCellNext(int cell_id, int surface_id, bool global = true,
                    bool neighbor = false);
    int getCellByStencil(int cell_id, int stencil_id);
    int getHexCellByStencil(int cell_id, int stencil_id);
    int getHexCellByStencilPre(int cell_id, int stencil_id);
    CMFD_PRECISION getFluxRatio(int cell_id, int group, int fsr);
    CMFD_PRECISION getHexFluxRatio(int cell_id, int group, int fsr);
    CMFD_PRECISION getUpdateRatio(int cell_id, int moc_group, int fsr);
    CMFD_PRECISION getHexUpdateRatio(int cell_id, int moc_group, int fsr);
    double getDistanceToCentroid(Point *centroid, int cell_id,
                                 int stencil_index);
    double getDistanceToCentroidHex(Point *centroid, int cell_id,
                                    int stencil_index);
    double getDistanceToCentroidHexPre(Point *centroid, int cell_id,
                                       int stencil_index);
    CMFD_PRECISION getSurfaceDiffusionCoefficient(int cmfd_cell, int surface,
                                                  int group, int moc_iteration,
                                                  bool correction);
    CMFD_PRECISION getHexSurfaceDiffusionCoefficient(int cmfd_cell, int surface,
                                                     int group, int moc_iteration,
                                                     bool correction);
    CMFD_PRECISION getDiffusionCoefficient(int cmfd_cell, int group);
    CMFD_PRECISION getSurfaceWidth(int surface);
    CMFD_PRECISION getHexSurfaceWidth(int surface);
    CMFD_PRECISION getPerpendicularSurfaceWidth(int surface);
    CMFD_PRECISION getHexPerpendicularSurfaceWidth(int surface);
    int getSense(int surface);
    int getLocalCMFDCell(int cmfd_cell);  // TODO: optimize, document
    int getGlobalCMFDCell(int cmfd_cell); // TODO: optimize, document
    int getCellColor(int cmfd_cell);      // TODO: optimize, document
    void packBuffers();
#ifdef ENABLE_MPI_
    void ghostCellExchange();
    void communicateSplits(bool faces);
#endif
    void unpackSplitCurrents(bool faces);
    void copyFullSurfaceCurrents();
    void copyHexFullSurfaceCurrents();
    void checkNeutronBalance(bool pre_split = true);
    void checkNeutronBalanceHex(bool pre_split = true);
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
