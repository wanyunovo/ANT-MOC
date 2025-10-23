
message(STATUS "Testing: building unittests for antmoc")

# list sources
set(antmoc_unittests
    test_BoolParser.cpp
    test_Cell_containsPoints.cpp
    test_Cell_distance.cpp
    test_ConfigInputCLI_parser.cpp
    test_FileUtils.cpp
    test_GeoInputXml_buildGlobal.cpp
    test_GeoInputXml_buildHexLattice.cpp
    test_GeoInputXml_buildHexLatticePrism.cpp
    test_GeoInputXml_buildLattice2D.cpp
    test_GeoInputXml_buildLattice3D.cpp
    test_GeoInputXml_buildPrimitives.cpp
    test_GeoInputXml_buildQuadraticSurfaces.cpp
    test_GeoInputXml_helpers.cpp
    test_GeoInputXml_layoutUtils.cpp
    test_HexLatticePrism.cpp
    test_initializePeriodicTracks.cpp
    test_initializeReflectiveTracks.cpp
    test_Lattice_cellCenters.cpp
    test_Lattice_cellVertices.cpp
    test_Lattice_variableWidths.cpp
    test_MaterialIdConversion.cpp
    test_MathUtils.cpp
    test_Point_arithmetic.cpp
    test_PyVector.cpp
    test_StringUtils.cpp
    test_TallyUtils.cpp
    test_VecParser.cpp
    test_XMLUtils.cpp
    test_ZCylinder.cpp
)

if (ENABLE_HIP)
    list(APPEND antmoc_unittests
        test_hip_DeviceFSRDataTypes.cpp
        test_hip_DeviceMaterial.cpp
        test_hip_DevicePoint.cpp
        test_hip_DeviceQuadrature.cpp
        test_hip_DeviceRayTracer.cpp
        test_hip_DeviceTracks.cpp
        test_hip_TrackIndexTypes.cpp
    )
endif ()

# discover serial tests when MPI is disabled
if (NOT ENABLE_MPI)
    antmoc_discover_tests(SOURCES     ${antmoc_unittests}
                          PREFIX      unit_
                          MAIN        antmoc_test_main
                          TARGET_LIST antmoc_unittest_targets)
else()
    antmoc_add_tests(SOURCES     ${antmoc_unittests}
                     PREFIX      unit_mpi_
                     MAIN        antmoc_test_main
                     MPI_RANKS   2
                     TARGET_LIST antmoc_unittest_targets)
endif()

set(antmoc_test_targets ${antmoc_test_targets} ${antmoc_unittest_targets} PARENT_SCOPE)
