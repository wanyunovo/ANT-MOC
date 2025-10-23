set(ENABLE_ALL_WARNINGS ON CACHE BOOL "")
set(ENABLE_DEBUG ON CACHE BOOL "")
set(ENABLE_TESTS ON CACHE BOOL "")

set(ENABLE_MPI ON CACHE BOOL "")
set(ENABLE_HIP ON CACHE BOOL "")

# Build shared libraries to avoid linking errors
set(BUILD_SHARED_LIBS ON CACHE BOOL "")

# CMake targets, e.g. hip::device, require compiler hipcc
set(CMAKE_C_COMPILER hipcc CACHE STRING "")
set(CMAKE_CXX_COMPILER hipcc CACHE STRING "")

set(_clang_compile_flags -Wno-unused-command-line-argument
                         -Wno-unused-parameter
                         -Wno-unused-variable)
set(_clang_link_flags)
set(_clang_common_flags)
set(ANTMOC_CLANG_COMPILE_FLAGS ${_clang_common_flags} ${_clang_compile_flags} CACHE STRING "")
set(ANTMOC_CLANG_LINK_FLAGS ${_clang_common_flags} ${_clang_link_flags} CACHE STRING "")

# Use local config files
set(ANTMOC_HIP_CONFIG_DIR ${CMAKE_CURRENT_LIST_DIR}/dtk CACHE STRING "")

# The only available GPU architecture is gfx906
set(GPU_TARGETS gfx906 CACHE STRING "")
if (NOT ${CMAKE_VERSION} VERSION_LESS "3.21")
  set(CMAKE_HIP_ARCHITECTURES gfx906 CACHE STRING "")
  set(HIP_SEPARABLE_COMPILATION ON CACHE BOOL "" )
endif ()

# Flags controlling builds
set(ANTMOC_HOST_CONFIG_LOADED ON CACHE BOOL "")
