# This is the top level CMake file for the SCOREC build
cmake_minimum_required(VERSION 3.8)
project(simLandIceMeshGen VERSION 0.0.1 LANGUAGES CXX)

enable_testing()
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

find_package(SimModSuite MODULE REQUIRED)
find_package(BLAS REQUIRED)

function(addExe name)
  add_executable(${name} ${name}.cc)
  target_include_directories(${name} PRIVATE ${SIMMODSUITE_INCLUDE_DIR})
  target_link_libraries(${name} PRIVATE ${SIMMODSUITE_LIBS})
endfunction()

addExe(landIceMeshGen)
add_test(NAME landIceMeshGen_oneSquare_jig
  COMMAND ./landIceMeshGen ${CMAKE_CURRENT_SOURCE_DIR}/oneSquare.msh one_jig 0.5)
add_test(NAME landIceMeshGen_twoSquares_jig
  COMMAND ./landIceMeshGen ${CMAKE_CURRENT_SOURCE_DIR}/twoSquares.msh two_jig 0.5)
add_test(NAME landIceMeshGen_oneSquare_vtk
  COMMAND ./landIceMeshGen ${CMAKE_CURRENT_SOURCE_DIR}/oneSquare.vtk one_vtk 0.5)
add_test(NAME landIceMeshGen_twoSquares_vtk
  COMMAND ./landIceMeshGen ${CMAKE_CURRENT_SOURCE_DIR}/twoSquares.vtk two_vtk 0.5)

addExe(testCurveCreation)

add_library(splineInterp splineInterpolation.cc BSpline.cc Expression.cc PolyNomial.cc)
target_link_libraries(splineInterp PRIVATE BLAS::BLAS)

addExe(testSplineFitting)
target_link_libraries(testSplineFitting PRIVATE splineInterp)
add_test(NAME testSplineFitting_single_edge
         COMMAND ./testSplineFitting)

