Include(FetchContent)

FetchContent_Declare(
  Catch2
  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
  GIT_TAG        v3.4.0 # or a later release
)

FetchContent_MakeAvailable(Catch2)


# Add the test executable
add_executable(tests
    test_mesh2d.cpp
    test_matrix.cpp
    test_matrix_legacy.cpp
    test_discrete_2d.cpp
    test_exprtk.cpp
    test_solver.cpp
    test_param_parse.cpp
)

# Link the Catch2 library and your static library to the test executable
target_link_libraries(tests PRIVATE Catch2::Catch2WithMain fin_diff_lib)

# Include directories for the tests
target_include_directories(tests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../src)

set_target_properties(tests PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/tests
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests
)

message(STATUS "tests CMAKE_CURRENT_SOURCE_DIR: ${CMAKE_CURRENT_SOURCE_DIR}")
message(STATUS "tests WORKING_DIRECTORY: ${CMAKE_BINARY_DIR}/tests")

# Copy the data folder to the build directory
add_custom_command(TARGET tests POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
    ${CMAKE_CURRENT_SOURCE_DIR}/data ${CMAKE_BINARY_DIR}/tests/data
)

# Define platform-specific macros
if (CMAKE_SYSTEM_NAME STREQUAL "Windows")
    target_compile_definitions(tests PRIVATE PLATFORM_WINDOWS)
elseif (CMAKE_SYSTEM_NAME STREQUAL "Linux" OR CMAKE_SYSTEM_NAME STREQUAL "Darwin")
    target_compile_definitions(tests PRIVATE PLATFORM_UNIX)
endif()


# Enable testing
enable_testing()

# Discover tests
include(CTest)
include(Catch)
catch_discover_tests(tests)