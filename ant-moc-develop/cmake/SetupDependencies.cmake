# ==============================================================================
# Set up dependencies
# ==============================================================================
set(antmoc_depends_dir ${PROJECT_SOURCE_DIR}/dependencies)

# OpenMP
if (ENABLE_OPENMP)
  find_package(OpenMP 3.1 REQUIRED)
  list(APPEND antmoc_exports openmp)
  antmoc_add_library(NAME       openmp
                     DEPENDS_ON OpenMP::OpenMP_CXX)
endif ()

# MPI
if (ENABLE_MPI)
  find_package(MPI 3.0 REQUIRED)
  list(APPEND antmoc_exports mpi)
  antmoc_add_library(NAME       mpi
                     DEPENDS_ON MPI::MPI_C)
endif ()

# HIP
if (ENABLE_HIP)
  include(cmake/SetupHIP.cmake)
  find_package(hip REQUIRED)
  find_package(rocthrust REQUIRED)
  list(APPEND antmoc_exports hiplibs)
  antmoc_add_library(NAME          hiplibs
                     DEPENDS_ON    hip::device roc::rocthrust
                     COMPILE_FLAGS -fgpu-rdc
                     LINK_FLAGS    -fgpu-rdc)
endif ()

# better-enums, required
if (USE_EXTERNAL_BETTER_ENUMS)
    set(BETTER_ENUMS_DIR $ENV{BETTER_ENUMS_DIR}
        CACHE STRING "better-enums directory")
else ()
    antmoc_get_submodule(
        NAME            better-enums
        HEADER_ONLY
        DEPS_DIR        ${antmoc_depends_dir}
        GIT_REPOSITORY  https://gitee.com/antmoc/better-enums.git
        GIT_TAG         0.11.3)

    set(BETTER_ENUMS_DIR ${antmoc_depends_dir}/better-enums/)
endif ()
list(APPEND antmoc_exports better-enums)
antmoc_add_library(NAME     better-enums
                   INCLUDES $<BUILD_INTERFACE:${BETTER_ENUMS_DIR}>)

# cxxopts, required
if (USE_EXTERNAL_CXXOPTS)
    find_package(cxxopts 3.0.0 REQUIRED)
else ()
    antmoc_get_submodule(
        NAME            cxxopts
        DEPS_DIR        ${antmoc_depends_dir}
        GIT_REPOSITORY  https://gitee.com/antmoc/cxxopts.git
        GIT_TAG         v3.0.0)

    set(CXXOPTS_BUILD_TESTS OFF CACHE BOOL "")
    set(CXXOPTS_BUILD_EXAMPLES OFF CACHE BOOL "")
    add_subdirectory(${antmoc_depends_dir}/cxxopts)
endif ()

# fmt (fmt::fmt), required
find_package(fmt 6.0.0 REQUIRED)

# toml11 (toml11::toml11), required
find_package(toml11 3.5.0 REQUIRED)

# TinyXML2 (tinyxml2::tinyxml2), required
find_package(tinyxml2 REQUIRED CONFIG)

# HDF5 (hdf5-shared or hdf5-static), required
find_package(HDF5 1.10.0 REQUIRED)
set(HAS_HDF5 ON CACHE BOOL "HDF5 found")

list(APPEND antmoc_exports hdf5)
antmoc_import_library(NAME    hdf5
                      PACKAGE HDF5
                      MAYBE   hdf5 hdf5-shared hdf5-static)

# Find parallel HDF5
find_program(EXTERNAL_PHDF5 h5pcc)
cmake_dependent_option(HAS_PHDF5 "" ON "EXTERNAL_PHDF5" OFF)
if (ENABLE_MPI AND NOT HAS_PHDF5)
    message(FATAL_ERROR "Could not find parallel HDF5 (MPI=${ENABLE_MPI})")
elseif (NOT ENABLE_MPI AND HAS_PHDF5)
    message(FATAL_ERROR "Could not find serial HDF5 (MPI=${ENABLE_MPI})")
endif ()

# Googletest
if (ENABLE_TESTS)
    find_package(GTest 1.10.0 REQUIRED CONFIG)
    list(APPEND antmoc_exports googletest)
    antmoc_add_library(NAME       googletest
                       DEPENDS_ON GTest::gtest GTest::gmock)
endif ()
