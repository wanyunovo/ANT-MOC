# ==============================================================================
#
# Print properties of targets for debugging
#
#     print_targets_properties(target1 target2 ...)
#
# Targets are processed one by one to provide more readable output.
#
# ==============================================================================
include(CMakePrintHelpers)

function(print_targets_properties)
  list(APPEND _props_interface
       TYPE
       INTERFACE_INCLUDE_DIRECTORIES
       INTERFACE_COMPILE_DEFINITIONS
       INTERFACE_COMPILE_OPTIONS
       INTERFACE_LINK_OPTIONS
       INTERFACE_LINK_LIBRARIES)

  # properties of other types of libraries
  list(APPEND _props_concrete
       VERSION
       ${_props_interface}
       INCLUDE_DIRECTORIES
       COMPILE_DEFINITIONS
       COMPILE_OPTIONS
       LINK_OPTIONS
       LINK_LIBRARIES)

  foreach(_parsed_target IN LISTS ARGN)
    # Check target type
    get_target_property(_target_type ${_parsed_target} TYPE)
    if (_target_type STREQUAL "INTERFACE_LIBRARY")
      set(_all_props ${_props_interface})
    else ()
      set(_all_props ${_props_concrete})
    endif ()

    # Check existed properties
    set(_existed_props)

    foreach(_prop IN LISTS _all_props)
      # Generate a property list
      get_target_property(_prop_value ${_parsed_target} ${_prop})
      if (_prop_value)
        list(APPEND _existed_props ${_prop})
      endif ()
    endforeach()

    cmake_print_properties(TARGETS    ${_parsed_target}
                           PROPERTIES ${_existed_props})
  endforeach()
endfunction()

# ==============================================================================
# Function
#     antmoc_get_submodule(
#         NAME           modulename
#         DEPS_DIR       dir    # Directory for dependencies
#         HEADER_ONLY           # Header-only projects
#         GIT_REPOSITORY url    # URL of the git repository
#         GIT_TAG        tag    # Git branch name, tag or commit hash
#     )
#
# Fetch a submodule.
# ==============================================================================
include(FetchContent)

function(antmoc_get_submodule)
  set(options HEADER_ONLY)
  set(one_value_args NAME DEPS_DIR TYPE GIT_REPOSITORY GIT_TAG)
  set(multi_value_args NONE_MULTI)
  cmake_parse_arguments(module "${options}"
                               "${one_value_args}"
                               "${multi_value_args}"
                               ${ARGN})
  if (NOT module_NAME)
    message(FATAL_ERROR "No module name specified for antmoc_get_submodule")
  endif ()

  if (NOT module_GIT_REPOSITORY)
    message(FATAL_ERROR "No repo specified for antmoc_get_submodule")
  endif ()

  if (NOT module_GIT_TAG)
    set(module_GIT_TAG master)
  endif ()

  if (NOT module_DEPS_DIR)
    set(module_DEPS_DIR ${PROJECT_SOURCE_DIR}/dependencies/)
  endif ()

  set(module_dir "${module_DEPS_DIR}/${module_NAME}")

  if ((NOT EXISTS "${module_dir}" AND module_HEADER_ONLY) OR
      (NOT EXISTS "${module_dir}/CMakeLists.txt" AND NOT module_HEADER_ONLY))
    FetchContent_Populate(
        ${module_NAME}
        GIT_REPOSITORY ${module_GIT_REPOSITORY}
        GIT_TAG        ${module_GIT_TAG}
        SOURCE_DIR     ${module_dir})
  endif ()
endfunction()

# ==============================================================================
# Set properties for a target
# ==============================================================================
macro(antmoc_patch_target)
  set(singleValueArgs NAME SCOPE CXX_STD)
  set(multiValueArgs DEPENDS_ON
                     INCLUDES
                     COMPILE_FLAGS
                     LINK_FLAGS
                     DEFINES)
  cmake_parse_arguments(arg "${options}"
                            "${singleValueArgs}"
                            "${multiValueArgs}"
                            ${ARGN})

  set(_scope PUBLIC)
  if (arg_SCOPE)
    set(_scope ${arg_SCOPE})
  endif ()

  if (arg_DEPENDS_ON)
    target_link_libraries(${arg_NAME}
        ${_scope} ${arg_DEPENDS_ON})
  endif ()

  if (arg_INCLUDES)
    target_include_directories(${arg_NAME}
        ${_scope} ${arg_INCLUDES})
  endif ()

  if (arg_COMPILE_FLAGS)
    target_compile_options(${arg_NAME}
        ${_scope} ${arg_COMPILE_FLAGS})
  endif ()

  if (arg_LINK_FLAGS)
    target_link_options(${arg_NAME}
        ${_scope} ${arg_LINK_FLAGS})
  endif ()

  if (arg_DEFINES)
    target_compile_definitions(${arg_NAME}
        ${_scope} ${arg_DEFINES})
  endif ()

  if (arg_CXX_STD)
    target_compile_features(${arg_NAME}
        ${_scope} ${arg_CXX_STD})
  endif ()
endmacro()

# ==============================================================================
# Add an executable
# ==============================================================================
macro(antmoc_add_executable)
  set(singleValueArgs NAME)
  set(multiValueArgs SOURCES
                     INCLUDES
                     DEPENDS_ON
                     COMPILE_FLAGS
                     LINK_FLAGS
                     DEFINES)
  cmake_parse_arguments(arg "${options}"
                            "${singleValueArgs}"
                            "${multiValueArgs}"
                            ${ARGN} )

  add_executable(${arg_NAME} ${arg_SOURCES})

  antmoc_patch_target(NAME          ${arg_NAME}
                      DEPENDS_ON    ${arg_DEPENDS_ON}
                      INCLUDES      ${arg_INCLUDES}
                      COMPILE_FLAGS ${arg_COMPILE_FLAGS}
                      LINK_FLAGS    ${arg_LINK_FLAGS}
                      DEFINES       ${arg_DEFINES}
                      SCOPE         PRIVATE)
endmacro()

# ==============================================================================
# Add a library
# ==============================================================================
macro(antmoc_add_library)
  set(singleValueArgs NAME SCOPE CXX_STD)
  set(multiValueArgs SOURCES
                     DEPENDS_ON
                     INCLUDES
                     COMPILE_FLAGS
                     LINK_FLAGS
                     DEFINES)
  cmake_parse_arguments(arg "${options}"
                            "${singleValueArgs}"
                            "${multiValueArgs}"
                            ${ARGN} )

  set(_scope PUBLIC)
  if (arg_SCOPE)
    set(_scope PRIVATE)
  endif ()

  if (arg_SOURCES)
    add_library(${arg_NAME} ${arg_SOURCES})
  else ()
    add_library(${arg_NAME} INTERFACE)
    set(_scope INTERFACE)
  endif ()

  antmoc_patch_target(NAME          ${arg_NAME}
                      DEPENDS_ON    ${arg_DEPENDS_ON}
                      INCLUDES      ${arg_INCLUDES}
                      COMPILE_FLAGS ${arg_COMPILE_FLAGS}
                      LINK_FLAGS    ${arg_LINK_FLAGS}
                      DEFINES       ${arg_DEFINES}
                      SCOPE         ${_scope}
                      CXX_STD       ${arg_CXX_STD})

  if (arg_DEPENDS_ON)
    add_dependencies(${arg_NAME} ${arg_DEPENDS_ON})
  endif ()
endmacro()

# ==============================================================================
# Import a library
# ==============================================================================
macro(antmoc_import_library)
  set(singleValueArgs NAME EXPORTABLE)
  set(multiValueArgs PACKAGE MAYBE)
  cmake_parse_arguments(arg "${options}"
                            "${singleValueArgs}"
                            "${multiValueArgs}"
                            ${ARGN} )
  if (NOT TARGET ${arg_NAME})
    # create a target
    if (${arg_EXPORTABLE})
      add_library(${arg_NAME} INTERFACE)
    else ()
      add_library(${arg_NAME} INTERFACE IMPORTED)
    endif ()

    # set package name for searching
    set(_package ${arg_NAME})
    if (arg_PACKAGE)
      set(_package ${arg_PACKAGE})
    endif ()

    set(_use_target OFF)

    foreach(_target IN LISTS ${arg_MAYBE})
      if (TARGET _target)
        message(STATUS "${arg_NAME}: using target ${_target}")
        set(_use_target ON)
        antmoc_patch_target(NAME      ${arg_NAME}
                            DEPENDS_ON ${_target})
        break()
      endif ()
    endforeach()

    if (NOT _use_target)
      set(_defines "${_package}_DEFINITIONS")
      set(_includes "${_package}_INCLUDE_DIRS")
      set(_libraries "${_package}_LIBRARIES")
      message(STATUS "${_defines}: ${${_defines}}")
      message(STATUS "${_includes}: ${${_includes}}")
      message(STATUS "${_libraries}: ${${_libraries}}")
      antmoc_patch_target(NAME       ${arg_NAME}
                          SCOPE      INTERFACE
                          DEFINES    ${${_defines}}
                          INCLUDES   ${${_includes}}
                          DEPENDS_ON ${${_libraries}})
    endif ()
  endif ()
endmacro()

# ==============================================================================
#
# Discover tests and register them to CTest.
# Names of tests are determined by gtest_discover_tests, whereas the
# root name is the name of the source file without its suffix.
#
#     antmoc_discover_tests(
#         SOURCES        test1.cpp test2.cpp ...
#         PREFIX         test_prefix
#         MAIN           test_main      # Default to gtest_main
#         EXTRA_ARGS     arg1 arg2 ...  # Extra arguments for googletest
#         EXTRA_INCLUDES arg1 arg2 ...  # Provide extra includes
#         MATERIAL_DIR   dir            # Define macro for tests, default to cases
#         GEOMETRY_DIR   dir            # Define macro for tests, default to tests/geometry
#         WORKING_DIRECTORY dir         # Default to current directory
#         TARGET_LIST    var            # Return the list of targets to var
#     )
#
# ==============================================================================
include(GoogleTest)


function(antmoc_discover_tests)

    set(options NONE)
    set(one_value_args   PREFIX MAIN TARGET_LIST MATERIAL_DIR GEOMETRY_DIR WORKING_DIRECTORY)
    set(multi_value_args SOURCES EXTRA_ARGS EXTRA_INCLUDES)
    cmake_parse_arguments(test "${options}"
                               "${one_value_args}"
                               "${multi_value_args}"
                               ${ARGN})

    if (NOT test_SOURCES)
        message(FATAL_ERROR "Testing: no source files")
    endif()

    if (NOT test_MAIN)
        set(test_MAIN gtest_main)
    endif()

    if (NOT test_MATERIAL_DIR)
        set(test_MATERIAL_DIR ${PROJECT_SOURCE_DIR}/tests/materials)
    endif()

    if (NOT test_GEOMETRY_DIR)
        set(test_GEOMETRY_DIR ${PROJECT_SOURCE_DIR}/tests/geometry)
    endif()

    if (NOT test_WORKING_DIRECTORY)
        set(test_WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
    endif()

    foreach (testfile ${test_SOURCES})

        # creates an testing executable for CTest
        get_filename_component(testname ${testfile} NAME_WE)
        set(testname "${test_PREFIX}${testname}")
        add_executable(${testname} ${testfile})

        # link the entry to it
        target_link_libraries(${testname} PRIVATE ${test_MAIN})

        # add extra includes in addition to what test_MAIN have
        target_include_directories(${testname} PRIVATE ${test_EXTRA_INCLUDES})

        # add extra definitions
        target_compile_definitions(${testname}
                                   PRIVATE MATERIAL_DIR=${test_MATERIAL_DIR}
                                           GEOMETRY_DIR=${test_GEOMETRY_DIR}
                                           WORKING_DIR=${test_WORKING_DIRECTORY})

        list(REMOVE_DUPLICATES test_EXTRA_ARGS)

        # register the test to CTest and set the working directory to root
        gtest_discover_tests(${testname}
                             TEST_PREFIX       ${test_PREFIX}
                             WORKING_DIRECTORY ${test_WORKING_DIRECTORY}
                             EXTRA_ARGS        ${test_EXTRA_ARGS})

        set_target_properties(${testname} PROPERTIES FOLDER tests)

        # append the executable to a list
        set(test_targets ${test_targets} ${testname})

    endforeach()

    if (test_TARGET_LIST)
        set(${test_TARGET_LIST} ${test_targets} PARENT_SCOPE)
    endif()

endfunction()


# ==============================================================================
#
# Add tests rather than discovering and register them to CTest.
# Names of tests are determined by omitting file suffix.
# This function allows single source file which contains only a main function.
#
#     antmoc_add_tests(
#         SOURCES        test1.cpp test2.cpp ...
#         PREFIX         test_prefix
#         MAIN           test_main      # Defaults to none
#         MPI_RANKS      4              # Defaults to none
#         EXTRA_ARGS     arg1 arg2 ...  # Extra arguments for test command
#         EXTRA_LIBS     arg1 arg2 ...  # Provice extra libraries
#         EXTRA_INCLUDES arg1 arg2 ...  # Provide extra includes
#         MATERIAL_DIR   dir            # Define macro for tests, default to cases
#         GEOMETRY_DIR   dir            # Define macro for tests, default to tests/geometry
#         WORKING_DIRECTORY dir         # Default to current directory
#         TARGET_LIST    var            # Return the list of targets to var
#     )
#
# ==============================================================================

function(antmoc_add_tests)

    set(options NONE)
    set(one_value_args   PREFIX MAIN MPI_RANKS TARGET_LIST MATERIAL_DIR GEOMETRY_DIR WORKING_DIRECTORY)
    set(multi_value_args SOURCES EXTRA_ARGS EXTRA_LIBS EXTRA_INCLUDES)
    cmake_parse_arguments(test "${options}"
                               "${one_value_args}"
                               "${multi_value_args}"
                               ${ARGN})

    if (NOT test_SOURCES)
        message(FATAL_ERROR "Testing: no source files")
    endif()

    if (NOT test_MATERIAL_DIR)
        set(test_MATERIAL_DIR ${PROJECT_SOURCE_DIR}/tests/materials)
    endif()

    if (NOT test_GEOMETRY_DIR)
        set(test_GEOMETRY_DIR ${PROJECT_SOURCE_DIR}/tests/geometry)
    endif()

    if (NOT test_WORKING_DIRECTORY)
        set(test_WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
    endif()

    if (test_MPI_RANKS)
        set(test_driver_wrapper "mpirun")
        set(test_driver_args    "-np" "${test_MPI_RANKS}")
        if (test_MPI_RANKS GREATER 1)
            execute_process(COMMAND "${test_driver_wrapper} ${test_driver_args}"
                            RESULT_VARIABLE test_mpi_oversubscribe)
            # oversubscribe slots for OpenMPI
            if (test_no_oversubscribe)
                message(STATUS "Testing: using --oversubscribe for MPI tests")
                set(test_driver_args ${test_driver_args} "--oversubscribe")
            endif()
        endif()
    endif()

    foreach (testfile ${test_SOURCES})

        # creates an testing executable for CTest
        get_filename_component(testname ${testfile} NAME_WE)
        set(testname "${test_PREFIX}${testname}")
        add_executable(${testname} ${testfile})

        # link the entry to it
        target_link_libraries(${testname} PRIVATE ${test_MAIN} ${test_EXTRA_LIBS})

        # add extra includes in addition to what test_MAIN have
        target_include_directories(${testname} PRIVATE ${test_EXTRA_INCLUDES})

        # add extra definitions
        target_compile_definitions(${testname}
                                   PRIVATE MATERIAL_DIR=${test_MATERIAL_DIR}
                                           GEOMETRY_DIR=${test_GEOMETRY_DIR}
                                           WORKING_DIR=${test_WORKING_DIRECTORY})

        set_target_properties(${testname} PROPERTIES FOLDER tests)

        list(REMOVE_DUPLICATES test_EXTRA_ARGS)

        # register the test to CTest and set the working directory to root
        add_test(NAME
                    ${testname}
                 COMMAND
                    ${test_driver_wrapper}
                            ${test_driver_args}
                            $<TARGET_FILE:${testname}> ${test_EXTRA_ARGS}
                 WORKING_DIRECTORY
                    ${test_WORKING_DIRECTORY}
        )

        # append the executable to a list
        set(test_targets ${test_targets} ${testname})

    endforeach()

    if (test_TARGET_LIST)
        set(${test_TARGET_LIST} ${test_targets} PARENT_SCOPE)
    endif()

endfunction()
