# ==============================================================================
# Register tests to CTest
# ==============================================================================
# antmoc application
if (NOT ENABLE_MPI)
    antmoc_add_tests(SOURCES           ${PROJECT_SOURCE_DIR}/apps/antmoc.cpp
                     PREFIX            test_
                     EXTRA_ARGS        -c cases/c5g7/simple-lattice/settings.xml -i 1 -o " "
                     EXTRA_LIBS        libantmoc
                     WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
                     TARGET_LIST       antmoc_apptest_targets)
else()
    antmoc_add_tests(SOURCES           ${PROJECT_SOURCE_DIR}/apps/antmoc.cpp
                     PREFIX            test_mpi_
                     MPI_RANKS         2
                     EXTRA_ARGS        -c cases/c5g7/simple-lattice/settings.xml -i 1 -o " "
                     EXTRA_LIBS        libantmoc
                     WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
                     TARGET_LIST       antmoc_apptest_targets)
endif()

set(antmoc_test_targets antmoc_test_main ${antmoc_apptest_targets})

# unit and regression tests
add_subdirectory(unittests)
add_subdirectory(regression)

# ==============================================================================
# Test coverage
# ==============================================================================
if (ENABLE_COVERAGE)
    if (ENABLE_MPI)
      message(FATAL_ERROR "Test coverage is not available when MPI is enabled")
    endif ()
    message(STATUS "Test coverage: enabled")
    message(STATUS "Test coverage: dependencies are\n${antmoc_test_targets}")

    include(ProcessorCount)
    ProcessorCount(n_proc)
    if (n_proc EQUAL 0)
        set(n_proc 1)
    endif()

    setup_target_for_coverage_lcov(
        NAME            coverage
        EXECUTABLE      ctest --progress --output-on-failure -j ${n_proc}
        DEPENDENCIES    ${antmoc_test_targets}
        EXTRACT         "apps/*" "include/*" "src/*"
        TITLE           "ANTMOC"
    )
endif()
