# Linux - set a default path for ROCM_PATH
if(NOT DEFINED ROCM_PATH)
  set(ROCM_PATH /opt/rocm)
endif()

#If HIP is not installed under ROCm, need this to find HSA assuming HSA is under ROCm
if(DEFINED ENV{ROCM_PATH})
  set(ROCM_PATH "$ENV{ROCM_PATH}")
endif()

macro(set_and_check _var _file)
  set(${_var} "${_file}")
  if(NOT EXISTS "${_file}")
    message(FATAL_ERROR "File or directory ${_file} referenced by variable ${_var} does not exist !")
  endif()
endmacro()

####################################################################################
include(CheckCXXCompilerFlag)
include(CMakeFindDependencyMacro)

#Number of parallel jobs by default is 1
if(NOT DEFINED HIP_CLANG_NUM_PARALLEL_JOBS)
  set(HIP_CLANG_NUM_PARALLEL_JOBS 1)
endif()
set(HIP_COMPILER "clang")
set(HIP_RUNTIME "rocclr")

set_and_check( hip_INCLUDE_DIR "${ROCM_PATH}/include" )
set_and_check( hip_INCLUDE_DIRS "${hip_INCLUDE_DIR}" )
set_and_check( hip_LIB_INSTALL_DIR "${ROCM_PATH}/lib" )
set_and_check( hip_BIN_INSTALL_DIR "${ROCM_PATH}/bin" )

set_and_check(hip_HIPCC_EXECUTABLE "${hip_BIN_INSTALL_DIR}/hipcc")
set_and_check(hip_HIPCONFIG_EXECUTABLE "${hip_BIN_INSTALL_DIR}/hipconfig")

if(HIP_COMPILER STREQUAL "clang")
    set(HIP_CLANG_ROOT "${ROCM_PATH}/llvm")
  if(NOT HIP_CXX_COMPILER)
    set(HIP_CXX_COMPILER ${CMAKE_CXX_COMPILER})
  endif()
  if(HIP_CXX_COMPILER MATCHES ".*hipcc")
    execute_process(COMMAND ${HIP_CXX_COMPILER} --version
                    OUTPUT_STRIP_TRAILING_WHITESPACE
                    OUTPUT_VARIABLE HIP_CLANG_CXX_COMPILER_VERSION_OUTPUT)
    if(HIP_CLANG_CXX_COMPILER_VERSION_OUTPUT MATCHES "InstalledDir:[ \t]*([^\n]*)")
      get_filename_component(HIP_CLANG_ROOT "${CMAKE_MATCH_1}" DIRECTORY)
    endif()
  elseif (HIP_CXX_COMPILER MATCHES ".*clang\\+\\+")
    get_filename_component(HIP_CLANG_ROOT "${HIP_CXX_COMPILER}" DIRECTORY)
    get_filename_component(HIP_CLANG_ROOT "${HIP_CLANG_ROOT}" DIRECTORY)
  endif()
  file(GLOB HIP_CLANG_INCLUDE_SEARCH_PATHS ${HIP_CLANG_ROOT}/lib/clang/*/include)
  find_path(HIP_CLANG_INCLUDE_PATH stddef.h
      HINTS
          ${HIP_CLANG_INCLUDE_SEARCH_PATHS}
      NO_DEFAULT_PATH)
  find_dependency(AMDDeviceLibs)
  set(AMDGPU_TARGETS "gfx900;gfx906;gfx908" CACHE STRING "AMD GPU targets to compile for")
  set(GPU_TARGETS "${AMDGPU_TARGETS}" CACHE STRING "GPU targets to compile for")
endif()

find_dependency(amd_comgr)

include( "${CMAKE_CURRENT_LIST_DIR}/hip-targets.cmake" )

#Using find_dependency to locate the dependency for the packages
#This makes the cmake generated file xxxx-targets to supply the linker libraries
# without worrying other transitive dependencies
find_dependency(hsa-runtime64)
find_dependency(Threads)

set(_IMPORT_PREFIX "${ROCM_PATH}/hip")

#if HSA is not under ROCm then provide CMAKE_PREFIX_PATH=<HSA_PATH>
find_path(HSA_HEADER hsa/hsa.h
  PATHS
    "${_IMPORT_PREFIX}/../include"
    /opt/rocm/include
)

if (HSA_HEADER-NOTFOUND)
  message (FATAL_ERROR "HSA header not found! ROCM_PATH environment not set")
endif()

# Right now this is only supported for amd platforms
set_target_properties(hip::host PROPERTIES
  INTERFACE_COMPILE_DEFINITIONS "__HIP_PLATFORM_HCC__=1;__HIP_PLATFORM_AMD__=1"
)

if(HIP_RUNTIME MATCHES "rocclr")
  set_target_properties(hip::galaxyhip PROPERTIES
    INTERFACE_COMPILE_DEFINITIONS "__HIP_ROCclr__=1"
  #  INTERFACE_INCLUDE_DIRECTORIES "${_IMPORT_PREFIX}/include;${HSA_HEADER}"
  #  INTERFACE_SYSTEM_INCLUDE_DIRECTORIES "${_IMPORT_PREFIX}/include;${HSA_HEADER}"
    INTERFACE_INCLUDE_DIRECTORIES "${_IMPORT_PREFIX}/include"
    INTERFACE_SYSTEM_INCLUDE_DIRECTORIES "${_IMPORT_PREFIX}/include"
  )
  set_target_properties(hip::device PROPERTIES
    INTERFACE_COMPILE_DEFINITIONS "__HIP_ROCclr__=1"
    INTERFACE_INCLUDE_DIRECTORIES "${_IMPORT_PREFIX}/../include"
    INTERFACE_SYSTEM_INCLUDE_DIRECTORIES "${_IMPORT_PREFIX}/../include"
  )
else()
  set_target_properties(hip::hip_hcc_static PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${_IMPORT_PREFIX}/include;${HSA_HEADER}"
    INTERFACE_SYSTEM_INCLUDE_DIRECTORIES "${_IMPORT_PREFIX}/include;${HSA_HEADER}")

  set_target_properties(hip::hip_hcc PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${_IMPORT_PREFIX}/include;${HSA_HEADER}"
    INTERFACE_SYSTEM_INCLUDE_DIRECTORIES "${_IMPORT_PREFIX}/include;${HSA_HEADER}"
  )
  set_target_properties(hip::device PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${_IMPORT_PREFIX}/../include"
    INTERFACE_SYSTEM_INCLUDE_DIRECTORIES "${_IMPORT_PREFIX}/../include"
  )
endif()

if(HIP_COMPILER STREQUAL "clang")
  get_property(compilePropIsSet TARGET hip::device PROPERTY INTERFACE_COMPILE_OPTIONS SET)

  if (NOT compilePropIsSet AND HIP_CXX_COMPILER MATCHES ".*clang\\+\\+")
    set_property(TARGET hip::device APPEND PROPERTY
      INTERFACE_COMPILE_OPTIONS -mllvm -amdgpu-early-inline-all=true -mllvm -amdgpu-function-calls=false
    )
  endif()

  if (NOT compilePropIsSet)
    if (EXISTS ${AMD_DEVICE_LIBS_PREFIX}/amdgcn/bitcode)
      set_property(TARGET hip::device APPEND PROPERTY
        INTERFACE_COMPILE_OPTIONS -xhip
      )
    else()
      set_property(TARGET hip::device APPEND PROPERTY
        INTERFACE_COMPILE_OPTIONS -xhip --hip-device-lib-path=${AMD_DEVICE_LIBS_PREFIX}/lib
      )
    endif()
  endif()

  set_property(TARGET hip::device APPEND PROPERTY
     INTERFACE_LINK_LIBRARIES --hip-link
  )

  set_property(TARGET hip::device APPEND PROPERTY
    INTERFACE_INCLUDE_DIRECTORIES "${HIP_CLANG_INCLUDE_PATH}/.."
  )

  set_property(TARGET hip::device APPEND PROPERTY
    INTERFACE_SYSTEM_INCLUDE_DIRECTORIES "${HIP_CLANG_INCLUDE_PATH}/.."
  )

  foreach(GPU_TARGET ${GPU_TARGETS})
      if (NOT compilePropIsSet)
      set_property(TARGET hip::device APPEND PROPERTY
        INTERFACE_COMPILE_OPTIONS "--cuda-gpu-arch=${GPU_TARGET}"
      )
      endif()
      set_property(TARGET hip::device APPEND PROPERTY
        INTERFACE_LINK_LIBRARIES "--cuda-gpu-arch=${GPU_TARGET}"
      )
  endforeach()
  #Add support for parallel build and link
  if(${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang")
    check_cxx_compiler_flag("-parallel-jobs=1" HIP_CLANG_SUPPORTS_PARALLEL_JOBS)
  endif()
  if(HIP_CLANG_NUM_PARALLEL_JOBS GREATER 1)
    if(${HIP_CLANG_SUPPORTS_PARALLEL_JOBS} )
      if (NOT compilePropIsSet)
      set_property(TARGET hip::device APPEND PROPERTY
        INTERFACE_COMPILE_OPTIONS -parallel-jobs=${HIP_CLANG_NUM_PARALLEL_JOBS} -Wno-format-nonliteral
      )
      endif()
      set_property(TARGET hip::device APPEND PROPERTY
        INTERFACE_LINK_LIBRARIES -parallel-jobs=${HIP_CLANG_NUM_PARALLEL_JOBS}
      )
    else()
      message("clang compiler doesn't support parallel jobs")
    endif()
  endif()

  # Add support for __fp16 and _Float16, explicitly link with compiler-rt
  set_property(TARGET hip::host APPEND PROPERTY
    INTERFACE_LINK_LIBRARIES "-L\"${HIP_CLANG_INCLUDE_PATH}/../lib/linux\" -lclang_rt.builtins-x86_64"
  )
  set_property(TARGET hip::device APPEND PROPERTY
    INTERFACE_LINK_LIBRARIES "-L\"${HIP_CLANG_INCLUDE_PATH}/../lib/linux\" -lclang_rt.builtins-x86_64"
  )
endif()

set(hip_LIBRARIES hip::host hip::device)
set(hip_LIBRARY ${hip_LIBRARIES})

set(HIP_INCLUDE_DIR ${hip_INCLUDE_DIR})
set(HIP_INCLUDE_DIRS ${hip_INCLUDE_DIRS})
set(HIP_LIB_INSTALL_DIR ${hip_LIB_INSTALL_DIR})
set(HIP_BIN_INSTALL_DIR ${hip_BIN_INSTALL_DIR})
set(HIP_LIBRARIES ${hip_LIBRARIES})
set(HIP_LIBRARY ${hip_LIBRARY})
set(HIP_HIPCC_EXECUTABLE ${hip_HIPCC_EXECUTABLE})
set(HIP_HIPCONFIG_EXECUTABLE ${hip_HIPCONFIG_EXECUTABLE})
