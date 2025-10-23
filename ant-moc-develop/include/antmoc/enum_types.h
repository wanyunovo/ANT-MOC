/**
 * @file enum_types.h
 * @details The type enums.
 *          Refactored in HPC&DE, USTB
 * @date April 18, 2019
 * @author William Boyd, MIT, Course 22 (wboyd@mit.edu)
 * @author An Wang, USTB (wangan@xs.ustb.edu.cn)
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


/// \brief Return better-enum object as string.
template <typename T>
std::string enumToString(T e) {
  return stringutils::underscoreToSpace(e._to_string());
}


/**
 * @enum boundaryType
 * @brief The types of boundary conditions for Surfaces.
 */
enum boundaryType {
  /** A vacuum boundary condition */
  VACUUM,

  /** A reflective boundary condition */
  REFLECTIVE,

  /** A periodic boundary condition */
  PERIODIC,

  /* Boundary between two domains (only in domain-decomposed geometry) 两个域之间的界面*/
  INTERFACE,

  /** No boundary type (typically an interface between flat source regions) 在FSR之间的表面是无类型的表面*/
  BOUNDARY_NONE
};


/// \enum segmentationType
/// \brief The types of Track segmentation supported by ANT-MOC.
BETTER_ENUM(segmentationType, char,

  /**< Explicit 2D segments (required for 2D simulations) */
  EXPLICIT_2D,

  /**< Explicit 3D segments */
  EXPLICIT_3D,

  /**< Axial on-the-fly 3D segment formation by 3D track */
  OTF_TRACKS,

  /**< Axial on-the-fly 3D segment formation by z-stack */
  OTF_STACKS
)


/// \enum surfaceType
/// \brief The types of surfaces supported by ANT-MOC.
enum surfaceType {
  PLANE,            ///< A general plane
  XPLANE,           ///< A plane perpendicular to the x-axis
  YPLANE,           ///< A plane perpendicular to the y-axis
  ZPLANE,           ///< A plane perpendicular to the z-axis
  ZCYLINDER,        ///< A cylinder with axis parallel to the z-axis
  QUADRATIC,        ///< A generalized quadratic surface
  LATTICEPRISM,     ///< The surface of a lattice
  RECLATTICEPRISM,  ///< The surface of a RecLattice
  HEXLATTICEPRISM   ///< The surface of a HexLattice
};


/**
 * @enum universeType
 * @brief The type of universe
 */
enum universeType{

  /** A simple non-repeating Universe */
  SIMPLE,

  /** A collection of Universes in a rectangular Lattice */
  LATTICE
};


/**
 * @enum latticeType
 * @brief The type of Lattice
 */
enum class latticeType{

  /** A rectangular layout Lattice */
  Rectangle,

  /** A hexagonal layout Lattice */
  Hexagon,

  /** An empty Lattice without any layout */
  EMPTY
};


/**
 * @enum cellType
 * @brief The type of cell.
*/
enum cellType {

  /** A cell filled by a Material */
  MATERIAL,

  /** A cell filled by a Universe */
  FILL,

  /** A cell not yet filled by anything */
  UNFILLED
};


/**
 * @enum cellPhy
 * @brief The physical meaning of cell
*/
enum cellPhy {

  /** A cell for a pin */
  PIN, 

  /** A cell for an assembly */
  ASSEMBLY, 

  /** A cell for a reactor core */
  CORE,

  /** A cell for a control rod */
  CROD,

  /** A cell for others */
  ELSE
};


/**
 * @enum solverMode
 * @brief The solution mode used by the MOC solver.
 */
enum solverMode {

  /** The forward flux distribution */
  FORWARD,

  /** The adjoint flux distribution */
  ADJOINT,
};


/**
 * @enum residualType
 * @brief The type of residual used for the convergence criterion.
 */
enum residualType {

  /** A residual on the scalar flux distribution */
  SCALAR_FLUX,

  /** A residual on the fission source distribution */
  FISSION_SOURCE,

  /** A residual on the total source distribution */
  TOTAL_SOURCE,
};


/**
 * @enum stabilizationType
 * @brief The type of stabilization to use on source iteration
 */
BETTER_ENUM(stabilizationType, char,

  /** General diagonal stabilization */
  DIAGONAL,

  /** Yamamoto's groupwise stabilization */
  YAMAMOTO,

  /** Global damping on the scalar flux update */
  GLOBAL,

  NONE
)


/// \enum quadratureType
/// \brief The types of quadrature sets supported by ANT-MOC.
BETTER_ENUM(quadratureType, char,
  TABUCHI_YAMAMOTO,
  LEONARD,
  GAUSS_LEGENDRE,
  EQUAL_WEIGHT,
  EQUAL_ANGLE
)


/// \enum solverType
/// \brief The types of solvers supported by ANT-MOC.
BETTER_ENUM(solverType, char,
  CPU_SOLVER,
  CPU_LS_SOLVER,
  GPU_SOLVER
)


/// \see OpenMC
enum class Orientation {
  x,  ///< Flat side of lattice parallel to y-axis
  y   ///< Flat side of lattice parallel to x-axis
};


/// \enum xsFileLayout
/// \brief The layout of cross-section file
BETTER_ENUM(XSFileLayout, char,
  NAMED,
  COMPRESSED
)


/// \enum tallyMeshType
/// \brief The type of tally mesh
BETTER_ENUM(tallyMeshType, char,
  RECTANGLE, /**< A rectangular layout tally mesh */
  HEXAGON    /**< A hexagonal layout tally mesh */
)


/// \enum trackMappingType
/// \brief The type of track mapping algorithms
BETTER_ENUM(trackMappingType, char,
  BLOCK,          /**< Mapping algorithm for block distribution */
  CYCLIC_TRACK,   /**< Mapping algorithm for cyclic track distribution */
  ANGLE,          /**< Mapping algorithm for angular decomposition */
  AUTO            /**< Automatically choose the best algorithm */
)


} /* namespace antmoc */

#endif /* ENUM_TYPES_H_ */
