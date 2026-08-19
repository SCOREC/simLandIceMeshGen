#!/bin/bash

set -e

usage="Usage: ./run_BSplineKokkos2D_scaling_test.sh /path/to/build_dir"

#Check commandline arguments
echo $#
echo $1
if [[ $# -lt 1 ]]; then
  echo $usage && exit 1
fi

BUILD_DIR=$1

#Check if the build dir provided exists
[[ ! -d $BUILD_DIR ]] && "BUILD_DIR ${BUILD_DIR} does not exist... exiting" && exit 1

#Compiling with the build provided
echo "Compiling with the build provided"
cmake --build $BUILD_DIR -j 10

#Running tests
RUN_DIR="${BUILD_DIR}/testKokkos2DScalingTest"
echo "Running Scaling Test with uniform number of pts per spline"

#1000 splines, 10 pts per spline, 50000 paracoords
${RUN_DIR} 1000 10 50000 uniform "scalingTestResult.csv"

#2000 splines, 10 pts per spline, 50000 paracoords
${RUN_DIR} 2000 10 50000 uniform "scalingTestResult.csv"

#6000 splines, 10 pts per spline, 50000 paracoords
${RUN_DIR} 6000 10 50000 uniform "scalingTestResult.csv"

echo "-------- END OF UNIFORM TEST --------"

echo "Running Scaling Test with variable sized pts per spline, sampled from a gaussian distribution"

#1000 splines, 10 pts per spline, 50000 paracoords
${RUN_DIR} 1000 10 50000 gaussian "scalingTestResult.csv"

#2000 splines, 10 pts per spline, 50000 paracoords
${RUN_DIR} 2000 10 50000 gaussian "scalingTestResult.csv"

#6000 splines, 10 pts per spline, 50000 paracoords
${RUN_DIR} 6000 10 50000 gaussian "scalingTestResult.csv"
echo "-------- END OF GAUSSIAN TEST --------"


