
message(STATUS "Testing: building regression tests for antmoc")

#===============================================================================
# List sources
#===============================================================================
# tests under the same directory share the same input/output files
set(antmoc_regressiontests
    test_boundary_mixed_3x3/test_boundary_mixed_3x3.cpp
    test_boundary_reflective_3x3/test_boundary_reflective_3x3.cpp
    test_boundary_vacuum_3x3/test_boundary_vacuum_3x3.cpp
    test_config_toml/test_config_toml.cpp
    test_config_xml/test_config_xml.cpp
    test_fluxes_symmetry_3x3/test_fluxes_symmetry_3x3.cpp
    test_HexLattice_vacuum/test_HexLattice_vacuum.cpp
    test_Lattice_control_rod/test_Lattice_control_rod.cpp
    test_Lattice_tally_mesh/test_Lattice_tally_mesh.cpp
    test_reflector/test_reflector.cpp
)

# list tests associated with Periodic Track Decomposition
set(antmoc_regressiontests_mpi_track
    test_boundary_mixed_3x3/mpiTrack_test_boundary_mixed_3x3.cpp
    test_boundary_mixed_3x3/mpiTrack_test_boundary_mixed_3x3_cyclicTrackMapping.cpp
    test_boundary_reflective_3x3/mpiTrack_test_boundary_reflective_3x3.cpp
    test_boundary_reflective_3x3/mpiTrack_test_boundary_reflective_3x3_cyclicTrackMapping.cpp
    test_boundary_vacuum_3x3/mpiTrack_test_boundary_vacuum_3x3.cpp
    test_boundary_vacuum_3x3/mpiTrack_test_boundary_vacuum_3x3_cyclicTrackMapping.cpp
)


#===============================================================================
# Extra options for all tests
#===============================================================================
if (ENABLE_TESTS_UPDATE)
  list(APPEND extra_options "--update-results")
endif()
if (ENABLE_TESTS_VISUALIZE)
  list(APPEND extra_options "--visualize")
endif()

#===============================================================================
# Discover serial tests when MPI is disabled,
# Register parallel tests when MPI is enabled.
#===============================================================================
if (NOT ENABLE_MPI)
    antmoc_discover_tests(SOURCES           ${antmoc_regressiontests}
                          PREFIX            regression_
                          MAIN              antmoc_test_main
                          EXTRA_ARGS        ${extra_options}
                          MATERIAL_DIR      ${PROJECT_SOURCE_DIR}/tests/materials/
                          GEOMETRY_DIR      ${PROJECT_SOURCE_DIR}/tests/geometry/
                          WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                          TARGET_LIST       antmoc_regressiontest_targets)
else()
    antmoc_add_tests(SOURCES           ${antmoc_regressiontests_mpi_track}
                     PREFIX            regression_np2_
                     MAIN              antmoc_test_main
                     MPI_RANKS         2
                     MATERIAL_DIR      ${PROJECT_SOURCE_DIR}/tests/materials/
                     GEOMETRY_DIR      ${PROJECT_SOURCE_DIR}/tests/geometry/
                     WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                     TARGET_LIST       antmoc_regressiontest_targets)

    antmoc_add_tests(SOURCES           ${antmoc_regressiontests_mpi_track}
                     PREFIX            regression_np3_
                     MAIN              antmoc_test_main
                     MPI_RANKS         3
                     MATERIAL_DIR      ${PROJECT_SOURCE_DIR}/tests/materials/
                     GEOMETRY_DIR      ${PROJECT_SOURCE_DIR}/tests/geometry/
                     WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                     TARGET_LIST       antmoc_regressiontest_targets)
endif()

set(antmoc_test_targets ${antmoc_test_targets} ${antmoc_regressiontest_targets} PARENT_SCOPE)
