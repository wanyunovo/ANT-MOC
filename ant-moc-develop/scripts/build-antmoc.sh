#!/bin/bash
set -e

BUILD_DIR=build

cmake -S . -B $BUILD_DIR \
  -DCMAKE_BUILD_TYPE=Release \
  -DENABLE_ALL_WARNINGS=ON \
  -DENABLE_MPI=ON \
  -DENABLE_TESTS=ON

cmake --build $BUILD_DIR -j$(nproc)

