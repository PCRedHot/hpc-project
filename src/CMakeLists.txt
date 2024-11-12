# Add static library for package
add_library(fin_diff_lib STATIC 
    mesh2d.cpp
    )

# Include directories for this target
target_include_directories(fin_diff_lib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include/internal)
target_include_directories(fin_diff_lib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include/external)

# Include directories for the project
include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}/include/internal
  ${CMAKE_CURRENT_SOURCE_DIR}/include/external
)