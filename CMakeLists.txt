# This is the top level CMake file for the SCOREC build
cmake_minimum_required(VERSION 3.8)
project(simLandIceMeshGen VERSION 0.0.1 LANGUAGES CXX)

enable_testing()
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

find_package(SimModSuite MODULE REQUIRED)
find_package(BLAS REQUIRED)

add_library(splineInterp splineInterpolation.cc BSpline.cc Expression.cc PolyNomial.cc)
target_link_libraries(splineInterp PRIVATE BLAS::BLAS)

function(addExe name)
  add_executable(${name} ${name}.cc)
  target_include_directories(${name} PRIVATE ${SIMMODSUITE_INCLUDE_DIR})
  target_link_libraries(${name} PRIVATE ${SIMMODSUITE_LIBS})
  target_link_libraries(${name} PRIVATE splineInterp)
endfunction()

addExe(landIceMeshGen)
set(testDataDir ${CMAKE_CURRENT_SOURCE_DIR}/contour_test_data)
add_test(NAME landIceMeshGen_gis
  COMMAND ./landIceMeshGen ${testDataDir}/gis.msh gis 0.5 1000)
set_tests_properties(landIceMeshGen_gis PROPERTIES
  PASS_REGULAR_EXPRESSION "mesh faces in surface: 172874")

add_test(NAME landIceMeshGen_oneSquare_jig
  COMMAND ./landIceMeshGen ${testDataDir}/oneSquare.msh one_jig 0.5 1)
set_tests_properties(landIceMeshGen_oneSquare_jig PROPERTIES
  PASS_REGULAR_EXPRESSION "mesh faces in surface: 4")

add_test(NAME landIceMeshGen_twoSquares_jig
  COMMAND ./landIceMeshGen ${testDataDir}/twoSquares.msh two_jig 0.5 1)
set_tests_properties(landIceMeshGen_twoSquares_jig PROPERTIES
  PASS_REGULAR_EXPRESSION "mesh faces in surface: 18")

add_test(NAME landIceMeshGen_oneSquare_vtk
  COMMAND ./landIceMeshGen ${testDataDir}/oneSquare.vtk one_vtk 0.5 1)
set_tests_properties(landIceMeshGen_oneSquare_jig PROPERTIES
  PASS_REGULAR_EXPRESSION "mesh faces in surface: 4")

add_test(NAME landIceMeshGen_twoSquares_vtk
  COMMAND ./landIceMeshGen ${testDataDir}/twoSquares.vtk two_vtk 0.5 1)
set_tests_properties(landIceMeshGen_twoSquares_jig PROPERTIES
  PASS_REGULAR_EXPRESSION "mesh faces in surface: 18")

add_test(NAME landIceMeshGen_oval_vtk
  COMMAND ./landIceMeshGen ${testDataDir}/oval_wBbox.vtk oval 0.0001 100)
set_tests_properties(landIceMeshGen_oval_vtk PROPERTIES
  PASS_REGULAR_EXPRESSION "Number of vertices in model: 5")


addExe(testCurveCreation)

addExe(testSplineFitting)
add_test(NAME testSplineFitting_single_edge
         COMMAND ./testSplineFitting)

