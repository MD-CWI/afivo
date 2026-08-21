#!/usr/bin/env bash

set -e

# Set names/directories
TARGET_DIR=`pwd`"/silo"
SILO_TARNAME="Silo-4.12.1.tar.gz"
SILO_DIRNAME="Silo-4.12.1"
BUILD_DIR="build"

# Do compilation etc. in build directory
mkdir -p ${BUILD_DIR}
cd ${BUILD_DIR}

# Extract
if [ ! -d ${SILO_DIRNAME} ]; then
    tar -xzf ${SILO_TARNAME}
fi

# Configure with CMake (out-of-source build)
cmake -S ${SILO_DIRNAME} -B ${SILO_DIRNAME}/cmake-build \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX=${TARGET_DIR} \
    -DSILO_ENABLE_SHARED=OFF \
    -DSILO_ENABLE_HDF5=ON \
    -DSILO_ENABLE_BROWSER=OFF \
    -DSILO_ENABLE_SILEX=OFF \
    -DSILO_ENABLE_FORTRAN=ON \
    -DBUILD_TESTING=OFF

# Build and install
cmake --build ${SILO_DIRNAME}/cmake-build --parallel
cmake --install ${SILO_DIRNAME}/cmake-build
