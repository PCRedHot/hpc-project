# Add static library for package
add_library(fin_diff_lib STATIC 
    mesh2d.cpp
    lhs2d.cpp
    matrix.cpp
    )

# Include directories for this target
target_include_directories(fin_diff_lib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include/internal)

# Include directories for the project
include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}/include/internal
  ${CMAKE_CURRENT_SOURCE_DIR}/include/external
)