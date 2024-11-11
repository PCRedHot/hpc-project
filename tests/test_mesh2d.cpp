#include "catch2/catch_all.hpp"

#include "mesh2d.hpp"
#include <iostream>

using namespace fin_diff;

TEST_CASE("Mesh2D Creation test", "[mesh2d]") {
    std::vector<std::pair<double, double>> coords = {
        {0.0, 0.0}, {1.0, 0.0}, {0.0, 1.0}, {1.0, 1.0}};
    std::vector<std::pair<size_t, size_t>> lines = {
        {0, 1}, {1, 3}, {3, 2}, {2, 0}  // Line elements forming a quad
    };

    Mesh2D mesh(coords, lines);
    REQUIRE(mesh.get_num_points() == 4);
    REQUIRE(mesh.get_num_elements() == 4);

    const auto &mesh_coords = mesh.get_coordinates();
    REQUIRE(mesh_coords.size() == 4);
    REQUIRE(mesh_coords[0] == std::make_pair(0.0, 0.0));
    REQUIRE(mesh_coords[3] == std::make_pair(1.0, 1.0));

    const auto &mesh_connectivity = mesh.get_lines();
    REQUIRE(mesh_connectivity.size() == 4);
    REQUIRE(mesh_connectivity[0] == std::make_pair<size_t, size_t>(0, 1));
    REQUIRE(mesh_connectivity[3] == std::make_pair<size_t, size_t>(2, 0));
}

TEST_CASE("RectMesh2D Creation test 1", "[rect_mesh2d_1]") {
    int Nx = 10;
    int Ny = 10;

    RectMesh2D mesh(10, 10);
    REQUIRE(mesh.get_Nx() == 10);
    REQUIRE(mesh.get_Ny() == 10);
    REQUIRE(mesh.get_num_points() == 100);
    REQUIRE(mesh.get_num_elements() == (Nx - 1) * (Ny - 1) * 2 + (Nx - 1) + (Ny - 1));

    const auto &coords = mesh.get_coordinates();
    REQUIRE(coords.size() == 100);
    REQUIRE(coords[0] == std::make_pair(0.0, 0.0));
    REQUIRE(coords[99] == std::make_pair(1.0, 1.0));

    const auto &lines = mesh.get_lines();
    REQUIRE(lines.size() == (Nx - 1) * (Ny - 1) * 2 + (Nx - 1) + (Ny - 1));

    mesh.write_to_vtk("test_mesh[rect_mesh2d_1].vtk");
}

TEST_CASE("RectMesh2D Creation test 2", "[rect_mesh2d_2]") {
    int Nx = 2;
    int Ny = 4;

    RectMesh2D mesh(Nx, Ny);
    REQUIRE(mesh.get_Nx() == 2);
    REQUIRE(mesh.get_Ny() == 4);
    REQUIRE(mesh.get_num_points() == 8);

    double dx = 1.0 / (Nx - 1);
    double dy = 1.0 / (Ny - 1);

    const auto &coords = mesh.get_coordinates();

    REQUIRE(coords.size() == 8);
    for (int i = 0; i < 8; i++) {
        REQUIRE(coords[i] == std::make_pair(i % Nx * dx, i / Nx * dy));
    }

    const auto &lines = mesh.get_lines();
    REQUIRE(lines.size() == 10);    

    mesh.write_to_vtk("test_mesh[rect_mesh2d_2].vtk");

    // UNSCOPED_INFO("Test case start");
    // for (const auto& conn : lines) {
    //     UNSCOPED_INFO(conn.first << ", " << conn.second);
    // }
    // CHECK(false);
}