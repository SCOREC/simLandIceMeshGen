# This is the top level CMake file for the SCOREC build
cmake_minimum_required(VERSION 3.8)
project(simLandIceMeshGen VERSION 0.0.1 LANGUAGES CXX)

enable_testing()
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

find_package(SimModSuite MODULE REQUIRED)
find_package(BLAS REQUIRED)


add_library(splineInterp splineInterpolation.cc BSpline.cc Expression.cc PolyNomial.cc)
target_link_libraries(splineInterp PRIVATE BLAS::BLAS)

add_library(landIceMeshGen landIceMeshGen.cc)
target_include_directories(landIceMeshGen PUBLIC ${CURRENT_CMAKE_SOURCE_DIR})
target_include_directories(landIceMeshGen PUBLIC ${SIMMODSUITE_INCLUDE_DIR})
target_link_libraries(landIceMeshGen PRIVATE ${SIMMODSUITE_LIBS})
target_link_libraries(landIceMeshGen PRIVATE splineInterp)

function(addExe name)
  add_executable(${name} ${name}.cc)
  target_link_libraries(${name} PRIVATE landIceMeshGen)
endfunction()

addExe(modelGen)
set(testDataDir ${CMAKE_CURRENT_SOURCE_DIR}/contour_test_data)
add_test(NAME landIceMeshGen_gis
  COMMAND ./modelGen ${testDataDir}/gis.msh gis 1 120 40 1)
set_tests_properties(landIceMeshGen_gis PROPERTIES
  PASS_REGULAR_EXPRESSION "mesh faces in surface: 172874")

add_test(NAME landIceMeshGen_oneSquare_jig
  COMMAND ./modelGen ${testDataDir}/oneSquare.msh one_jig 1e-6 120 30 1)
set_tests_properties(landIceMeshGen_oneSquare_jig PROPERTIES
  PASS_REGULAR_EXPRESSION "mesh faces in surface: 4")

add_test(NAME landIceMeshGen_twoSquares_jig
  COMMAND ./modelGen ${testDataDir}/twoSquares.msh two_jig 1e-6 120 30 1)
set_tests_properties(landIceMeshGen_twoSquares_jig PROPERTIES
  PASS_REGULAR_EXPRESSION "mesh faces in surface: 18")

add_test(NAME landIceMeshGen_oneSquare_vtk
  COMMAND ./modelGen ${testDataDir}/oneSquare.vtk one_vtk 1e-6 120 30 1)
set_tests_properties(landIceMeshGen_oneSquare_vtk PROPERTIES
  PASS_REGULAR_EXPRESSION "mesh faces in surface: 4")

add_test(NAME landIceMeshGen_twoSquares_vtk
  COMMAND ./modelGen ${testDataDir}/twoSquares.vtk two_vtk 1e-6 120 30 1)
set_tests_properties(landIceMeshGen_twoSquares_vtk PROPERTIES
  PASS_REGULAR_EXPRESSION "mesh faces in surface: 18")

add_test(NAME landIceMeshGen_oval_vtk
  COMMAND ./modelGen ${testDataDir}/oval_wBbox.vtk oval 1e-8 120 120 1)
set_tests_properties(landIceMeshGen_oval_vtk PROPERTIES
  PASS_REGULAR_EXPRESSION "Number of vertices in model: 14")

add_library(curveReader curveReader.cc)
target_include_directories(curveReader PUBLIC ${CURRENT_CMAKE_SOURCE_DIR})

addExe(testCreateEdges)
target_link_libraries(testCreateEdges PRIVATE curveReader)

add_test(NAME curve4277
  COMMAND ./testCreateEdges ${testDataDir}/testCurve4277.csv curve4277)
set_tests_properties(curve4277 PROPERTIES
  PASS_REGULAR_EXPRESSION "Number of vertices in model: 6")
set_tests_properties(curve4277 PROPERTIES
  PASS_REGULAR_EXPRESSION "Number of edges in model: 6")

add_test(NAME curve1293
  COMMAND ./testCreateEdges ${testDataDir}/testCurve1293.csv curve1293)
set_tests_properties(curve1293 PROPERTIES
  PASS_REGULAR_EXPRESSION "Number of vertices in model: 5")
set_tests_properties(curve1293 PROPERTIES
  PASS_REGULAR_EXPRESSION "Number of edges in model: 5")

add_test(NAME curve1655
  COMMAND ./testCreateEdges ${testDataDir}/testCurve1655.csv curve1655)
set_tests_properties(curve1655 PROPERTIES
  PASS_REGULAR_EXPRESSION "Number of vertices in model: 4")
set_tests_properties(curve1655 PROPERTIES
  PASS_REGULAR_EXPRESSION "Number of edges in model: 4")

add_test(NAME intersection2339_2344
  COMMAND ./testCreateEdges
  ${testDataDir}/testIntersectionEdge2339Edge2344.csv
  intersection2339_2344)

add_test(NAME edge522
  COMMAND ./testCreateEdges ${testDataDir}/testEdge522.csv edge522)
set_tests_properties(edge522 PROPERTIES
  PASS_REGULAR_EXPRESSION "Number of vertices in model: 8")
set_tests_properties(edge522 PROPERTIES
  PASS_REGULAR_EXPRESSION "Number of edges in model: 8")

addExe(testOnCurve)
target_link_libraries(testOnCurve PRIVATE curveReader)
add_test(NAME onCurve_edge522
  COMMAND ./testOnCurve ${testDataDir}/testEdge522.csv 4 0 150)

addExe(testSplineFitting)
add_test(NAME testSplineFitting_single_edge
         COMMAND ./testSplineFitting)

