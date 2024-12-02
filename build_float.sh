#!/bin/bash

# Create build directory if it doesn't exist
mkdir -p build_float

# Navigate to build directory
cd build_float

# Run cmake with Release type
cmake -DCMAKE_BUILD_TYPE=Release -DPRECISION=float ..

# Build the project
cmake --build .