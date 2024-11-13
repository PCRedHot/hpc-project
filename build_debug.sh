#!/bin/bash

# Create build directory if it doesn't exist
mkdir -p build_debug

# Navigate to build directory
cd build_debug

# Run cmake with Debug type
cmake -DCMAKE_BUILD_TYPE=Debug ..

# Build the project
cmake --build .