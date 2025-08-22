# This is the top level CMake file for the SCOREC build
cmake_minimum_required(VERSION 3.8)
project(simLandIceMeshGen VERSION 0.0.1 LANGUAGES CXX)

enable_testing()
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

option(USE_SIMMODSUITE "Enable SimModSuite support" ON)

if(USE_SIMMODSUITE)
  find_package(SimModSuite MODULE REQUIRED)
  message(STATUS "Building with SimModSuite support")
else()
  message(STATUS "Building without SimModSuite support")
endif()

find_package(BLAS REQUIRED)


add_library(splineInterp splineInterpolation.cc BSpline.cc Expression.cc PolyNomial.cc)
target_link_libraries(splineInterp PUBLIC BLAS::BLAS)

if(USE_SIMMODSUITE)
  add_library(landIceMeshGen landIceMeshGen.cc)
  target_include_directories(landIceMeshGen PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
  target_include_directories(landIceMeshGen PUBLIC ${SIMMODSUITE_INCLUDE_DIR})
  target_include_directories(landIceMeshGen PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/Quadtree/include)
  target_link_libraries(landIceMeshGen PRIVATE ${SIMMODSUITE_LIBS})
  target_link_libraries(landIceMeshGen PRIVATE splineInterp)
  target_compile_definitions(landIceMeshGen PRIVATE USE_SIMMODSUITE)
endif()

set(testDataDir ${CMAKE_CURRENT_SOURCE_DIR}/contour_test_data)

if(USE_SIMMODSUITE)
  function(addExe name)
    add_executable(${name} ${name}.cc)
    target_link_libraries(${name} PRIVATE landIceMeshGen)
  endfunction()

  addExe(modelGen)

  add_test(NAME landIceMeshGen_gis
  COMMAND ./modelGen ${testDataDir}/gis.msh gis 1 100 100 1 "m")
set_tests_properties(landIceMeshGen_gis PROPERTIES
  PASS_REGULAR_EXPRESSION "Number of vertices in model: 925")
set_tests_properties(landIceMeshGen_gis PROPERTIES
  PASS_REGULAR_EXPRESSION "mesh faces in surface: 3326")

add_test(NAME landIceMeshGen_ais
  COMMAND ./modelGen ${testDataDir}/ais.vtk ais 1 100 100 1 "m")

add_test(NAME landIceMeshGen_oneSquare_jig
  COMMAND ./modelGen ${testDataDir}/oneSquare.msh one_jig 1e-6 100 100 1 "km")
set_tests_properties(landIceMeshGen_oneSquare_jig PROPERTIES
  PASS_REGULAR_EXPRESSION "mesh faces in surface: 4")

add_test(NAME landIceMeshGen_twoSquares_jig
  COMMAND ./modelGen ${testDataDir}/twoSquares.msh two_jig 1e-6 100 100 1 "km")
set_tests_properties(landIceMeshGen_twoSquares_jig PROPERTIES
  PASS_REGULAR_EXPRESSION "mesh faces in surface: 18")

add_test(NAME landIceMeshGen_oneSquare_vtk
  COMMAND ./modelGen ${testDataDir}/oneSquare.vtk one_vtk 1e-6 120 30 1 "km")
set_tests_properties(landIceMeshGen_oneSquare_vtk PROPERTIES
  PASS_REGULAR_EXPRESSION "mesh faces in surface: 4")

add_test(NAME landIceMeshGen_twoSquares_vtk
  COMMAND ./modelGen ${testDataDir}/twoSquares.vtk two_vtk 1e-6 120 30 1 "km")
set_tests_properties(landIceMeshGen_twoSquares_vtk PROPERTIES
  PASS_REGULAR_EXPRESSION "mesh faces in surface: 18")

add_test(NAME landIceMeshGen_oval_vtk
  COMMAND ./modelGen ${testDataDir}/oval_wBbox.vtk oval 1e-8 120 120 1 "km")
set_tests_properties(landIceMeshGen_oval_vtk PROPERTIES
  PASS_REGULAR_EXPRESSION "Number of vertices in model: 6")
endif()

add_library(curveReader curveReader.cc)
target_include_directories(curveReader PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_library(modelGen2d modelGen2d.cc)
target_include_directories(modelGen2d PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(modelGen2d PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/Quadtree/include)
target_link_libraries(modelGen2d PRIVATE splineInterp)

if(USE_SIMMODSUITE)
  addExe(testCreateEdges)
  target_link_libraries(testCreateEdges PRIVATE curveReader)
endif()

if(USE_SIMMODSUITE)
  add_test(NAME curve4277
    COMMAND ./testCreateEdges ${testDataDir}/testCurve4277.csv curve4277)
  set_tests_properties(curve4277 PROPERTIES
    PASS_REGULAR_EXPRESSION "Number of vertices in model: 3")
  set_tests_properties(curve4277 PROPERTIES
    PASS_REGULAR_EXPRESSION "Number of edges in model: 3")

  add_test(NAME curve1293
    COMMAND ./testCreateEdges ${testDataDir}/testCurve1293.csv curve1293)
  set_tests_properties(curve1293 PROPERTIES
    PASS_REGULAR_EXPRESSION "Number of vertices in model: 2")
  set_tests_properties(curve1293 PROPERTIES
    PASS_REGULAR_EXPRESSION "Number of edges in model: 2")

  add_test(NAME curve1655
    COMMAND ./testCreateEdges ${testDataDir}/testCurve1655.csv curve1655)
  set_tests_properties(curve1655 PROPERTIES
    PASS_REGULAR_EXPRESSION "Number of vertices in model: 2")
  set_tests_properties(curve1655 PROPERTIES
    PASS_REGULAR_EXPRESSION "Number of edges in model: 2")

  add_test(NAME intersection2339_2344
    COMMAND ./testCreateEdges
    ${testDataDir}/testIntersectionEdge2339Edge2344.csv
    intersection2339_2344)

  add_test(NAME edge522
    COMMAND ./testCreateEdges ${testDataDir}/testEdge522.csv edge522)
  set_tests_properties(edge522 PROPERTIES
    PASS_REGULAR_EXPRESSION "Number of vertices in model: 4")
  set_tests_properties(edge522 PROPERTIES
    PASS_REGULAR_EXPRESSION "Number of edges in model: 4")

  add_test(NAME dblMdlVtx
    COMMAND ./testCreateEdges ${testDataDir}/testDoubleMdlVtx.csv dblMdlVtx)
  set_tests_properties(dblMdlVtx PROPERTIES
    PASS_REGULAR_EXPRESSION "Number of vertices in model: 3")
  set_tests_properties(dblMdlVtx PROPERTIES
    PASS_REGULAR_EXPRESSION "Number of edges in model: 3")

  add_test(NAME fjord
    COMMAND ./testCreateEdges ${testDataDir}/fjord.csv fjord)
endif()

# Create function for executables that only use curveReader and modelGen2d (no SimModSuite)
function(addCurveExe name)
  add_executable(${name} ${name}.cc)
  target_link_libraries(${name} PRIVATE curveReader modelGen2d)
endfunction()

# Build curve-only executables when SimModSuite is disabled
if(NOT USE_SIMMODSUITE)
  addCurveExe(testOnCurve)
  add_test(NAME onCurve_edge522
    COMMAND ./testOnCurve ${testDataDir}/testEdge522.csv 4 0 150)
  
  addCurveExe(testDiscoverTopo)
  add_test(NAME testDiscoverTopo_dblMdlVtx
    COMMAND ./testDiscoverTopo ${testDataDir}/testDoubleMdlVtxWBbox.csv)
  add_test(NAME testDiscoverTopo_L11pts
    COMMAND ./testDiscoverTopo ${testDataDir}/L_11pts.csv)
  set_tests_properties(testDiscoverTopo_L11pts PROPERTIES
    PASS_REGULAR_EXPRESSION "number of model vertices: 3")
endif()

if(USE_SIMMODSUITE)
  addExe(testOnCurve)
  target_link_libraries(testOnCurve PRIVATE curveReader)
  add_test(NAME onCurve_edge522
    COMMAND ./testOnCurve ${testDataDir}/testEdge522.csv 4 0 150)

  addExe(testDiscoverTopo)
  target_link_libraries(testDiscoverTopo PRIVATE curveReader)
  add_test(NAME testDiscoverTopo_dblMdlVtx
    COMMAND ./testDiscoverTopo ${testDataDir}/testDoubleMdlVtxWBbox.csv)
  add_test(NAME testDiscoverTopo_L11pts
    COMMAND ./testDiscoverTopo ${testDataDir}/L_11pts.csv)
  set_tests_properties(testDiscoverTopo_L11pts PROPERTIES
    PASS_REGULAR_EXPRESSION "number of model vertices: 3")

  addExe(testSplineFitting)
  target_link_libraries(testSplineFitting PRIVATE curveReader splineInterp)
  
  add_test(NAME testSplineFitting_straight_edge
           COMMAND ./testSplineFitting ${testDataDir}/testCurveStraightLine.csv 9)
  add_test(NAME testSplineFitting_single_edge
           COMMAND ./testSplineFitting ${testDataDir}/testCurveSplineFitting.csv 12.2685)
  add_test(NAME testSplineFitting_sine_edge
           COMMAND ./testSplineFitting ${testDataDir}/testCurveSineSegment.csv 10.24227846)
endif()

