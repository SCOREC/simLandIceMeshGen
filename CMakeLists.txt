# This is the top level CMake file for the SCOREC build
cmake_minimum_required(VERSION 3.8)
project(simLandIceMeshGen VERSION 0.0.1 LANGUAGES CXX)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

find_package(SimModSuite MODULE REQUIRED)

add_executable(landIceMeshGen landIceMeshGen.cc)
target_include_directories(landIceMeshGen PRIVATE ${SIMMODSUITE_INCLUDE_DIR})
target_link_libraries(landIceMeshGen PRIVATE ${SIMMODSUITE_LIBS})

