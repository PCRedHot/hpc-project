# Minimum CMake version
cmake_minimum_required(VERSION 3.10)

# Project name and C++ standard
project(hpc_project)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Add the executable
add_executable(program main.cpp)

# Add src subdirectory
add_subdirectory(src)

# Link the library to the program
target_link_libraries(program PRIVATE fin_diff_lib)

# Catch2 Testing
add_subdirectory(tests)
