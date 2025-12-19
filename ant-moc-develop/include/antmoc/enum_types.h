/**
 * @文件 enum_types.h
 * @details 类型枚举。
 *在HPC&DE、科技大学重构
 * @日期 2019 年 4 月 18 日
 * @作者 William Boyd，麻省理工学院，课程 22 (wboyd@mit.edu)
 * @作者王安，科技大学 (wangan@xs.ustb.edu.cn)
 */

#ifndef ENUM_TYPES_H_
#define ENUM_TYPES_H_

#if defined(__INTEL_COMPILER)
#define BETTER_ENUMS_NO_CONSTEXPR
#endif

#define BETTER_ENUMS_DEFAULT_CONSTRUCTOR(Enum) \
  public:                                      \
    Enum() = default;

#include <enum.h>
#include "antmoc/string_utils.h"

namespace antmoc {


/// \brief 将 better-enum 对象作为字符串返回。
template <typename T>
std::string enumToString(T e) {
  return stringutils::underscoreToSpace(e._to_string());
}


/**
 * @enum 边界类型
 * @brief 表面边界条件的类型。
 */
enum boundaryType {
  /** 真空边界条件 */
  VACUUM,

  /** 反射边界条件 */
  REFLECTIVE,

  /** 周期性边界条件 */
  PERIODIC,

  /* Boundary between two domains (only in domain-decomposed geometry) 两个域之间的界面*/
  INTERFACE,

  /** No boundary type (typically an interface between flat source regions) 在FSR之间的表面是无类型的表面*/
  BOUNDARY_NONE
};


/// \enum 分段类型
/// \brief ANT-MOC 支持的轨道分段类型。
BETTER_ENUM(segmentationType, char,

  /**< 显式 2D 段（2D 模拟所需） */
  EXPLICIT_2D,

  /**显式 3D 追踪，所有 3D 轨迹分段都显式存下来； */
  EXPLICIT_3D,

  /**< 通过 3D 轨迹形成轴向动态 3D 分段 */
  OTF_TRACKS,

  /**< 通过 z 堆栈形成轴向动态 3D 段 */
  OTF_STACKS
)


/// \enum surface type
/// \brief ANT-MOC supported surface types.
enum surfaceType {
  PLANE,            ///< Ordinary aircraft
  XPLANE,           ///< plane perpendicular to the x-axis
  YPLANE,           ///< plane perpendicular to the y-axis
  ZPLANE,           ///< plane perpendicular to the z-axis
  ZCYLINDER,        ///< Cylinder with axis parallel to z-axis
  QUADRATIC,        ///< Generalized Quadratic Surface
  LATTICEPRISM,     ///< lattice surface
  RECLATTICEPRISM,  ///< RecLattice 的表面
  HEXLATTICEPRISM   ///< Surface of hexagonal lattice
};


/**
 * @enum UniverseType
 * @brief 宇宙的类型
 */
enum universeType{

  /** 一个简单的非重复宇宙 */
  SIMPLE,

  /** 矩形格子中的宇宙集合 */
  LATTICE
};


/**
 * @enumlatticeType
 * @brief 格子的类型
 */
enum class latticeType{

  /** 矩形布局格子 */
  Rectangle,

  /** 六边形布局格子 */
  Hexagon,

  /** 没有任何布局的空格子 */
  EMPTY
};


/**
 * @enum cellType
 * @brief 单元格的类型。
*/
enum cellType {

  /** 一个由材质填充的单元格 */
  MATERIAL,

  /** 一个细胞充满了宇宙 */
  FILL,

  /** 细胞还没有被任何东西填满 */
  UNFILLED
};


/**
 * @enum cellPhy
 * @brief 细胞的物理意义
*/
enum cellPhy {

  /** 一个单元用于一个引脚 */
  PIN, 

  /** 装配单元 */
  ASSEMBLY, 

  /** 反应堆堆芯的单元 */
  CORE,

  /** 控制棒的单元 */
  CROD,

  /** 为他人准备的牢房 */
  ELSE
};


/**
 * @enum 求解器模式
 * @brief MOC 求解器使用的求解模式。
 */
enum solverMode {

  /** 正向通量分布 */
  FORWARD,

  /** 伴随通量分布 */
  ADJOINT,
};


/**
 * @enum残差类型
 * @brief 用于收敛标准的残差类型。
 */
enum residualType {

  /** 标量通量分布的残差 */
  SCALAR_FLUX,

  /** 裂变源分布的残差 */
  FISSION_SOURCE,

  /** 总源分布的残差 */
  TOTAL_SOURCE,
};


/**
 * @enum 稳定类型
 * @brief 在源迭代中使用的稳定类型
 */
BETTER_ENUM(stabilizationType, char,

  /** 一般对角稳定 */
  DIAGONAL,

  /** 山本的群体稳定 */
  YAMAMOTO,

  /** 标量通量更新的全局阻尼 */
  GLOBAL,

  NONE
)


/// \enum 正交类型
/// \brief ANT-MOC 支持的正交集类型。
BETTER_ENUM(quadratureType, char,
  TABUCHI_YAMAMOTO,
  LEONARD,
  GAUSS_LEGENDRE,
  EQUAL_WEIGHT,
  EQUAL_ANGLE
)


/// \enum 求解器类型
/// \brief ANT-MOC 支持的求解器类型。
BETTER_ENUM(solverType, char,
  CPU_SOLVER,
  CPU_LS_SOLVER,
  GPU_SOLVER
)


/// \参见 OpenMC
enum class Orientation {
  x,  ///< 晶格平行于 y 轴的平面
  y   ///< 晶格平行于 x 轴的平面
};


/// \enum xsFileLayout
/// \brief 横截面文件的布局
BETTER_ENUM(XSFileLayout, char,
  NAMED,
  COMPRESSED
)


/// \enumtallyMeshType
/// \brief 计数网格的类型
BETTER_ENUM(tallyMeshType, char,
  RECTANGLE, /**< 矩形布局计数网格 */
  HEXAGON    /**< 六边形布局计数网格 */
)


/// \enum trackMappingType
/// \简介 轨迹映射算法的类型
BETTER_ENUM(trackMappingType, char,
  BLOCK,          /**< 区块分布的映射算法 */
  CYCLIC_TRACK,   /**< 循环轨迹分布的映射算法 */
  ANGLE,          /**< 角度分解的映射算法 */
  AUTO            /**< 自动选择最佳算法 */
)


} /* 命名空间 Antmoc */

#endif /* 枚举类型 h */
