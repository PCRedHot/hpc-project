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
    # test_utils.cpp
)

# Link the Catch2 library and your static library to the test executable
target_link_libraries(tests PRIVATE Catch2::Catch2WithMain fin_diff_lib)

# Include directories for the tests
target_include_directories(tests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../src)

# Enable testing
enable_testing()

# Discover tests
include(CTest)
include(Catch)
catch_discover_tests(tests)