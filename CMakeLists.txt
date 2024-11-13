# Minimum CMake version
cmake_minimum_required(VERSION 3.10)

# Project name and C++ standard
project(hpc_project)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Debug mode
if (CMAKE_BUILD_TYPE MATCHES Debug)
    add_definitions(-D__DEBUG__)
endif()
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

add_subdirectory(src)

# Add the executable
add_executable(program main.cpp)
# Link the library to the program
target_link_libraries(program PRIVATE fin_diff_lib)

# Option to compile tests
option(COMPILE_TESTS "Compile the tests subdirectory" ON)

if(COMPILE_TESTS)
  add_subdirectory(tests)
endif()

