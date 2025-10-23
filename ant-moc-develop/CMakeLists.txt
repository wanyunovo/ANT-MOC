cmake_minimum_required(VERSION 3.16 FATAL_ERROR)

# Project name
project(antmoc VERSION 0.1.14)

# ==============================================================================
# Set up options
# ==============================================================================
option(USE_INTERNAL_ALL "Build all dependencies from sources" OFF)
option(USE_EXTERNAL_BETTER_ENUMS "Use external better-enums" OFF)
option(USE_EXTERNAL_CXXOPTS "Use external cxxopts" OFF)

option(BUILD_SHARED_LIBS "Builds shared libraries" ON)
option(ENABLE_ALL_WARNINGS "Enables all warnings" ON)
option(ENABLE_DEBUG "Enables debugging mode" ON)
option(ENABLE_IPO "Enables interprocedural optimization" OFF)
option(ENABLE_PROFILING "Enables profiling with gprof" OFF)
option(ENABLE_TESTS "Builds tests" OFF)
option(ENABLE_TESTS_UPDATE     "Build tests for updating oracles" OFF)
option(ENABLE_TESTS_VISUALIZE  "Build tests for visualization"    OFF)
option(ENABLE_COVERAGE "Enables code coverage" OFF)

option(ENABLE_OPENMP "Enables OpenMP" ON)
option(ENABLE_MPI "Enables MPI" ON)
option(ENABLE_HIP "Enables HIP" OFF)
option(ENABLE_DOUBLE_PRECISION "Enables double precision" ON)
option(ENABLE_CYCLIC "Enables cyclic track generation" OFF)
option(ENABLE_CMFD "Enables CMFD acceleration" ON)
option(ENABLE_LINEAR_SOURCE "Enables linear source approximation" OFF)

if (ENABLE_COVERAGE)
  set(ENABLE_DEBUG ON CACHE BOOL "" FORCE)
  set(ENABLE_TESTS ON CACHE BOOL "" FORCE)
endif ()

if (USE_INTERNAL_ALL)
  set(USE_EXTERNAL_BETTER_ENUMS OFF CACHE BOOL "" FORCE)
  set(USE_EXTERNAL_CXXOPTS OFF CACHE BOOL "" FORCE)
endif ()

# Define global variables
set(antmoc_depends)
set(antmoc_defines)
set(antmoc_includes)
set(antmoc_compile_flags)
set(antmoc_link_flags)
set(antmoc_cxx_std cxx_std_11)
set(antmoc_exports)

# Include CMake modules
include(CMakeDependentOption)
include(cmake/ANTMOCMacros.cmake)
include(cmake/SetupCompilers.cmake)
include(cmake/SetupUserOptions.cmake)

# ==============================================================================
# Set up dependencies
# ==============================================================================
# Header directories
list(APPEND antmoc_includes
            $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
            $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
            $<INSTALL_INTERFACE:include>)

# Set up all the dependencies
include(cmake/SetupDependencies.cmake)

# Required third-party packages
list(APPEND antmoc_depends
            better-enums
            cxxopts
            hdf5
            fmt::fmt
            toml11::toml11
            tinyxml2::tinyxml2)

# MPI
if (ENABLE_MPI)
  list(APPEND antmoc_depends mpi)
  list(APPEND antmoc_defines ENABLE_MPI_ USE_MPI_IN_PLACE)
endif ()

# OpenMP
if (ENABLE_OPENMP)
  list(APPEND antmoc_depends openmp)
endif ()

# HIP
if (ENABLE_HIP)
  set(antmoc_cxx_std cxx_std_17)
  list(APPEND antmoc_depends hiplibs)
  list(APPEND antmoc_defines ENABLE_HIP_)
endif ()

# ==============================================================================
# Add antmoc libraries
# ==============================================================================
set(antmoc_sources
    src/fileio/BoolParser.cpp
    src/fileio/file_utils.cpp
    src/fileio/FSRDataHandler.cpp
    src/fileio/FSRDataHandlerHDF5.cpp
    src/fileio/GeoInput.cpp
    src/fileio/GeoInputXml.cpp
    src/fileio/GeoInputXMLUtils.cpp
    src/fileio/HDF5Handler.cpp
    src/fileio/MaterialHandler.cpp
    src/fileio/MaterialHandlerHDF5.cpp
    src/fileio/Mesh.cpp
    src/fileio/tally_utils.cpp
    src/fileio/VecParser.cpp
    src/fileio/xml_utils.cpp
    src/geometry/Cell.cpp
    src/geometry/FSRDataTypes.cpp
    src/geometry/Geometry.cpp
    src/geometry/Lattice.cpp
    src/geometry/LocalCoords.cpp
    src/geometry/Point.cpp
    src/geometry/Surface.cpp
    src/geometry/Universe.cpp
    src/logger/log.cpp
    src/logger/Progress.cpp
    src/material/Material.cpp
    src/options/ConfigInputCLI.cpp
    src/options/ConfigInput.cpp
    src/options/ConfigInputFile.cpp
    src/options/Option.cpp
    src/options/ParentNodeStack.cpp
    src/options/toml_utils.cpp
    src/solver/Cmfd.cpp
    src/solver/CPULSSolver.cpp
    src/solver/CPUSolver.cpp
    src/solver/ExpEvaluator.cpp
    src/solver/linalg.cpp
    src/solver/Matrix.cpp
    src/solver/MOCKernel.cpp
    src/solver/openmp_utils.cpp
    src/solver/Solver.cpp
    src/solver/TrackTraversingAlgorithms.cpp
    src/solver/TraverseSegments.cpp
    src/solver/Vector.cpp
    src/timer/Timer.cpp
    src/track/Quadrature.cpp
    src/track/Track3D.cpp
    src/track/Track.cpp
    src/track/TrackGenerator3D.cpp
    src/track/TrackGenerator.cpp
    src/track/TrackLoadBalance.cpp
    src/utils/math_utils.cpp
    src/utils/mpi_utils.cpp
    src/utils/pairwise_sum.cpp
    src/utils/ParallelHashMap.cpp
    src/utils/string_utils.cpp
)

set(antmoc_hip_sources
    src/geometry/hip/DeviceFSRDataTypes.cpp
    src/geometry/hip/DevicePoint.cpp
    src/material/hip/DeviceMaterial.cpp
    src/solver/hip/DeviceMOCKernel.cpp
    #src/solver/hip/GPUSolver.cpp
    src/track/hip/DeviceQuadrature.cpp
    src/track/hip/DeviceRayTracer.cpp
    src/track/hip/DeviceTrack.cpp
    src/track/hip/DeviceTrack3D.cpp
    src/utils/hip/hip_info.cpp
    src/utils/hip/hip_utils.cpp
    src/utils/hip/hip_utils_instantiator.cpp
)

if (ENABLE_HIP)
  list(APPEND antmoc_sources ${antmoc_hip_sources})
endif ()

antmoc_add_library(NAME          libantmoc
                   SOURCES       ${antmoc_sources}
                   INCLUDES      ${antmoc_includes}
                   DEPENDS_ON    ${antmoc_depends}
                   COMPILE_FLAGS ${antmoc_compile_flags}
                   LINK_FLAGS    ${antmoc_link_flags}
                   DEFINES       ${antmoc_defines}
                   CXX_STD       ${antmoc_cxx_std})

set_target_properties(libantmoc PROPERTIES PREFIX "")

# ==============================================================================
# Add tests
# ==============================================================================
if (ENABLE_COVERAGE)
  include(cmake/CodeCoverage.cmake)
endif ()
if (ENABLE_TESTS)
  enable_testing()
  antmoc_add_library(NAME        antmoc_test_main
                     SOURCES     tests/run_tests.cpp
                     INCLUDES    tests/include/
                     DEPENDS_ON  libantmoc googletest)
  add_subdirectory(tests)
endif ()

# ==============================================================================
# Add antmoc applications 
# ==============================================================================
set(antmoc_app_sources ${PROJECT_SOURCE_DIR}/apps/antmoc.cpp)
antmoc_add_executable(NAME       antmoc
                      SOURCES    ${antmoc_app_sources}
                      DEPENDS_ON libantmoc)

# Debugging
if (ENABLE_DEBUG)
  print_targets_properties(${antmoc_depends} libantmoc antmoc)
endif ()
# ==============================================================================
# Installation
# ==============================================================================
include(CMakePackageConfigHelpers)

install(TARGETS               libantmoc antmoc
        EXPORT                antmocTargets
        DESTINATION           cmake/antmoc/
        LIBRARY  DESTINATION  lib
        ARCHIVE  DESTINATION  lib
        RUNTIME  DESTINATION  bin
        INCLUDES DESTINATION  include)

list(REMOVE_DUPLICATES antmoc_exports)
foreach(dep ${antmoc_exports})
  # If the target is EXPORTABLE, add it to the export set
  get_target_property(_is_imported ${dep} IMPORTED)
  if(NOT ${_is_imported})
    install(TARGETS     ${dep}
            EXPORT      antmocTargets
            DESTINATION lib)
    # Namespace target to avoid conflicts
    set_target_properties(${dep} PROPERTIES EXPORT_NAME antmoc::${dep})
  endif()
endforeach()

configure_package_config_file(antmocConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/antmocConfig.cmake
    INSTALL_DESTINATION cmake/antmoc/)

install(DIRECTORY tools/ DESTINATION tools)
install(DIRECTORY include/ DESTINATION include FILES_MATCHING PATTERN *.h)
install(FILES       ${CMAKE_CURRENT_BINARY_DIR}/antmocConfig.cmake
        DESTINATION cmake/antmoc/)

install(EXPORT      antmocTargets
        FILE        antmocTargets.cmake
        NAMESPACE   antmoc::
        DESTINATION cmake/antmoc/)

# ==============================================================================
# Print summary
# ==============================================================================
message("

======================================================
               Configuration Summary
------------------------------------------------------

    External better-enums     : ${USE_EXTERNAL_BETTER_ENUMS}
    External cxxopts          : ${USE_EXTERNAL_CXXOPTS}
    HDF5                      : ${HAS_HDF5} (Parallel: ${HAS_PHDF5})

    Shared libraries          : ${BUILD_SHARED_LIBS}
    IPO                       : ${ENABLE_IPO}
    Warnings                  : ${ENABLE_ALL_WARNINGS}
    Debugging                 : ${ENABLE_DEBUG}
    Profiling                 : ${ENABLE_PROFILING}
    Code coverage             : ${ENABLE_COVERAGE}
    Tests                     : ${ENABLE_TESTS}
    Tests updating            : ${ENABLE_TESTS_UPDATE}
    Tests visualizing         : ${ENABLE_TESTS_VISUALIZE}

    OpenMP                    : ${ENABLE_OPENMP}
    MPI                       : ${ENABLE_MPI}
    HIP                       : ${ENABLE_HIP}
    Double precision          : ${ENABLE_DOUBLE_PRECISION}
    Cyclic tracks             : ${ENABLE_CYCLIC}
    CMFD                      : ${ENABLE_CMFD}
    Linear source             : ${ENABLE_LINEAR_SOURCE}

======================================================
")
