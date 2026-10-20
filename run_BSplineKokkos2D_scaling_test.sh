#!/bin/bash

set -e

usage="Usage: ./run_BSplineKokkos2D_scaling_test.sh /path/to/build_dir <disableSerial: yes -or- no> <number of splines>"

#Check commandline arguments
if [[ $# -lt 3 ]]; then
  echo $usage && exit 1
fi

BUILD_DIR=$1
DISABLE_SERIAL=$2
NUM_SPLINES=$3
MODE=$4

#Check if the build dir provided exists
[[ ! -d $BUILD_DIR ]] && "BUILD_DIR ${BUILD_DIR} does not exist... exiting" && exit 1

#Compiling with the build provided
echo "Compiling with the build provided"
cmake --build $BUILD_DIR -j 10

#Running tests
RUN_DIR="${BUILD_DIR}/testKokkos2DScalingTest"
echo "Running Scaling Test with uniform number of pts per spline"

#10 pts per spline average, 100 paracoords
${RUN_DIR} ${NUM_SPLINES} 10 100 ${MODE} "scalingTestResult.csv" ${DISABLE_SERIAL}


