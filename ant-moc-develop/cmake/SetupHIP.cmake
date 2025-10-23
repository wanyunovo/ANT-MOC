# Set HIP paths for antmoc
if (NOT DEFINED ROCM_PATH)
    if (DEFINED ENV{ROCM_PATH})
        set(ROCM_PATH $ENV{ROCM_PATH} CACHE PATH "ROCm path")
    elseif (DEFINED ENV{HIP_PATH})
        set(ROCM_PATH $ENV{HIP_PATH}/.. CACHE PATH "ROCm path")
    else ()
        set(ROCM_PATH /opt/rocm CACHE PATH "ROCm path")
    endif ()
endif ()

# Set ROCm and HIP paths for CMake
if (ANTMOC_HIP_CONFIG_DIR)
  list(APPEND CMAKE_PREFIX_PATH ${ANTMOC_HIP_CONFIG_DIR})
else ()
  list(APPEND CMAKE_PREFIX_PATH ${ROCM_PATH} ${ROCM_PATH}/hip)
endif ()
