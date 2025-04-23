# This is the top level CMake file for the SCOREC build
cmake_minimum_required(VERSION 3.8)
project(simLandIceMeshGen VERSION 0.0.1 LANGUAGES CXX)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

find_package(SimModSuite MODULE REQUIRED)
find_package(BLAS REQUIRED)

function(addExe name)
  add_executable(${name} ${name}.cc)
  target_include_directories(${name} PRIVATE ${SIMMODSUITE_INCLUDE_DIR})
  target_link_libraries(${name} PRIVATE ${SIMMODSUITE_LIBS})
endfunction()

addExe(landIceMeshGen)

addExe(testCurveCreation)

add_library(splineInterp splineInterpolation.cc BSpline.cc Expression.cc PolyNomial.cc)
target_link_libraries(splineInterp PRIVATE BLAS::BLAS)

addExe(testSplineFitting)
target_link_libraries(testSplineFitting PRIVATE splineInterp)
enable_testing()
add_test(NAME testSplineFitting_single_edge
         COMMAND ./testSplineFitting)

