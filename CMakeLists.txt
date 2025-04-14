# This is the top level CMake file for the SCOREC build
cmake_minimum_required(VERSION 3.8)
project(simLandIceMeshGen VERSION 0.0.1 LANGUAGES CXX)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

find_package(SimModSuite MODULE REQUIRED)

function(addTest name)
  add_executable(${name} ${name}.cc)
  target_include_directories(${name} PRIVATE ${SIMMODSUITE_INCLUDE_DIR})
  target_link_libraries(${name} PRIVATE ${SIMMODSUITE_LIBS})
endfunction()

addTest(landIceMeshGen)
addTest(testCurveCreation)
addTest(testSplineFitting)

