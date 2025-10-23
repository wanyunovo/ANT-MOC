# Compiler identification
set(gcc_like_cxx "$<COMPILE_LANG_AND_ID:CXX,AppleClang,Clang,GNU>")
set(gcc_like_c   "$<COMPILE_LANG_AND_ID:C,AppleClang,Clang,GNU>")
set(intel_cxx    "$<COMPILE_LANG_AND_ID:CXX,Intel>")
set(intel_c      "$<COMPILE_LANG_AND_ID:C,Intel>")
set(gcc_cxx      "$<COMPILE_LANG_AND_ID:CXX,GNU>")
set(clang_cxx    "$<COMPILE_LANG_AND_ID:CXX,AppleClang,Clang>")

# Extensions
set(CMAKE_CXX_EXTENSIONS OFF CACHE BOOL "")

if (CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
  #list(APPEND antmoc_depends mkl)
  list(APPEND antmoc_defines INTEL)
endif ()

if (ENABLE_DEBUG)
  message(STATUS "Debug: adding flags for debugging")
  list(APPEND antmoc_defines ENABLE_DEBUG_)
  list(APPEND antmoc_compile_flags -O1 -g)
else ()
  list(APPEND antmoc_compile_flags
       -O3
       "$<${intel_cxx}:$<BUILD_INTERFACE:-xHost -ansi-alias -no-prec-div>>")

  # Interprocedural optimization
  if (ENABLE_IPO)
    message(STATUS "IPO: adding flags for IPO")
    list(APPEND antmoc_compile_flags
         "$<${gcc_like_cxx}:$<BUILD_INTERFACE:-ipo>>"
         "$<${intel_cxx}:$<BUILD_INTERFACE:-flto>>")
  endif (ENABLE_IPO)

    # Optimization reporting flags
  if (ENABLE_OPT_REPORT)
    message(STATUS "Optimization report: adding flags for reporting compile optimization")
    list(APPEND antmoc_compile_flags
         "$<${gcc_like_cxx}:$<BUILD_INTERFACE:-fopt-info-all;-openmp-report>>"
         "$<${intel_cxx}:$<BUILD_INTERFACE:-qopt-report=5,-openmp-report>>")
  endif (ENABLE_OPT_REPORT)
endif (ENABLE_DEBUG)

if (ENABLE_ALL_WARNINGS)
  list(APPEND antmoc_compile_flags -Wall -Wextra)
endif (ENABLE_ALL_WARNINGS)

# Code coverage with gcov
if (ENABLE_COVERAGE)
    list(APPEND antmoc_compile_flags --coverage)
    list(APPEND antmoc_link_flags --coverage)
endif(ENABLE_COVERAGE)

#===============================================================================
# Extra flags
#===============================================================================
# Extra flags for GCC
if (ANTMOC_GCC_COMPILE_FLAGS)
  list(APPEND antmoc_compile_flags
       "$<${gcc_cxx}:$<BUILD_INTERFACE:${ANTMOC_GCC_COMPILE_FLAGS}>>")
endif ()
if (ANTMOC_GCC_LINK_FLAGS)
  list(APPEND antmoc_link_flags
       "$<${gcc_cxx}:$<BUILD_INTERFACE:${ANTMOC_GCC_LINK_FLAGS}>>")
endif ()

# Extra flags for Clang
if (ANTMOC_CLANG_COMPILE_FLAGS)
  list(APPEND antmoc_compile_flags
       "$<${clang_cxx}:$<BUILD_INTERFACE:${ANTMOC_CLANG_COMPILE_FLAGS}>>")
endif ()
if (ANTMOC_CLANG_LINK_FLAGS)
  list(APPEND antmoc_link_flags
       "$<${clang_cxx}:$<BUILD_INTERFACE:${ANTMOC_CLANG_LINK_FLAGS}>>")
endif ()
