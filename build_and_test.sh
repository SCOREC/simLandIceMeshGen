#!/bin/bash

# Script to build and test simLandIceMeshGen with and without SimModSuite

set -e  # Exit on any error

usage="Usage: ./build_and_test.sh /path/to/simLandIceMeshGen /path/to/simmodsuite/env"
[[ $# -ne 2 ]] && echo $usage && exit 1

PROJECT_DIR=$PWD
SOURCE_DIR=$1
SIMMODSUITE_ENV=$2
[[ ! -d $SOURCE_DIR ]] && "SOURCE_DIR ${SOURCE_DIR} does not exist... exiting" && exit 1
[[ ! -d $SIMMODSUITE_ENV ]] && "SIMMODSUITE_ENV ${SIMMODSUITE_ENV} does not exist... exiting" && exit 1

echo "=== Building and Testing simLandIceMeshGen ==="
echo "Project directory: $PROJECT_DIR"
echo "Source directory: $SOURCE_DIR"
echo

# Function to clean build directory
clean_build() {
    local build_dir=$1
    if [ -d "$build_dir" ]; then
        echo "Cleaning $build_dir..."
        rm -rf "$build_dir"
    fi
}

# Function to build and test
build_and_test() {
    local config=$1
    local build_dir=$2
    local cmake_options=$3
    local env_setup=$4
    
    echo "=== Building with $config ==="
    
    # Clean previous build
    clean_build "$build_dir"
    
    # Configure
    echo "Configuring..."
    if [ -n "$env_setup" ]; then
        eval "$env_setup && cmake -S \"$SOURCE_DIR\" -B \"$build_dir\" $cmake_options"
    else
        cmake -S "$SOURCE_DIR" -B "$build_dir" $cmake_options
    fi
    
    # Build
    echo "Building..."
    cmake --build "$build_dir"
    
    # Test
    echo "Running tests..."
    if [ -n "$env_setup" ]; then
        eval "$env_setup && ctest --test-dir \"$build_dir\" --output-on-failure"
    else
        ctest --test-dir "$build_dir" --output-on-failure
    fi
    
    echo "✅ $config build and test completed successfully"
    echo
}

# Build without SimModSuite
echo "🔧 Testing WITHOUT SimModSuite support..."
build_and_test "SimModSuite DISABLED" "$PROJECT_DIR/buildSimLandIceMeshGen" "-DUSE_SIMMODSUITE=OFF"

# Build with SimModSuite
echo "🔧 Testing WITH SimModSuite support..."
build_and_test "SimModSuite ENABLED" "$PROJECT_DIR/buildSimLandIceMeshGenWithSim" "-DUSE_SIMMODSUITE=ON" "source ${SIMMODSUITE_ENV}"

echo "🎉 All builds and tests completed successfully!"
echo
echo "Summary:"
echo "- SimModSuite DISABLED: buildSimLandIceMeshGen/"
echo "- SimModSuite ENABLED:  buildSimLandIceMeshGenWithSim/"
