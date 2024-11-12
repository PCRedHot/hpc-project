# Minimum CMake version
cmake_minimum_required(VERSION 3.10)

# Project name and C++ standard
project(hpc_project)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Set the default build type to Debug if not specified
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Debug CACHE STRING "Choose the type of build (Debug or Release)" FORCE)
endif()

# Add the executable
add_executable(program main.cpp)

# Add src subdirectory
add_subdirectory(src)

# Link the library to the program
target_link_libraries(program PRIVATE fin_diff_lib)

# Catch2 Testing
add_subdirectory(tests)

# Set output directories based on build type
set_target_properties(program PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY_DEBUG ${CMAKE_BINARY_DIR}/Debug/bin
    RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_BINARY_DIR}/Release/bin
    ARCHIVE_OUTPUT_DIRECTORY_DEBUG ${CMAKE_BINARY_DIR}/Debug/lib
    ARCHIVE_OUTPUT_DIRECTORY_RELEASE ${CMAKE_BINARY_DIR}/Release/lib
    LIBRARY_OUTPUT_DIRECTORY_DEBUG ${CMAKE_BINARY_DIR}/Debug/lib
    LIBRARY_OUTPUT_DIRECTORY_RELEASE ${CMAKE_BINARY_DIR}/Release/lib
)

# Add debug-specific compiler options and definitions
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  target_compile_options(program PRIVATE -g)
  target_compile_definitions(program PRIVATE DEBUG)
endif()