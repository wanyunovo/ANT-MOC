# ------------------------------------------------------------------
# Cyclic track generation
# This is the track initialization algorithm for original antmoc
# ------------------------------------------------------------------
if (ENABLE_CYCLIC)
    message(STATUS "Track generation: cyclic")
    list(APPEND antmoc_defines USTB_)
else ()
    list(APPEND antmoc_defines MIT_)
endif ()

# ------------------------------------------------------------------
# Profiling
# ------------------------------------------------------------------
if (ENABLE_PROFILING STREQUAL "ON")
    message(STATUS "Profiling: adding flags for gprof")
    list(APPEND antmoc_defines ENABLE_PROFILE_)
    list(APPEND antmoc_compile_flags -pg)
    list(APPEND antmoc_link_flags -pg)
endif ()

# ------------------------------------------------------------------
# CMFD
# ------------------------------------------------------------------
# FIXME
if (ENABLE_CMFD)
    list(APPEND antmoc_defines ENABLECMFD)
endif ()

# ------------------------------------------------------------------
# Floating-point precision
# ------------------------------------------------------------------
if (ENABLE_DOUBLE_PRECISION)
    list(APPEND antmoc_defines
                FP_PRECISION=double
                CMFD_PRECISION=double
                LINALG_TOL=1E-15
                VEC_LENGTH=4)
else ()
    list(APPEND antmoc_defines
                FP_PRECISION=float
                CMFD_PRECISION=float
                LINALG_TOL=1E-7
                VEC_LENGTH=8)
endif ()
# The vector alignment used in solvers when allocating aligned data
# structures using MM_MALLOC and MM_FREE
list(APPEND antmoc_defines VEC_ALIGNMENT=64)

# ------------------------------------------------------------------
# Source type
# ------------------------------------------------------------------
if (ENABLE_LINEAR_SOURCE)
    list(APPEND antmoc_defines LINEARSOURCE)
endif ()
