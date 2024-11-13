#!/bin/bash

# Create build directory if it doesn't exist
mkdir -p build

# Navigate to build directory
cd build

# Run cmake with Release type
cmake -DCMAKE_BUILD_TYPE=Release ..

# Build the project
cmake --build .