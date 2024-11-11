#include "mesh2d.hpp"

#include <cmath>
#include <fstream>
#include <iostream>

namespace fin_diff {

Mesh2D::Mesh2D(const std::vector<std::pair<double, double>>& coordinates,
               const std::vector<std::pair<size_t, size_t>>& lines)
    : N_pt(coordinates.size()),
      N_el(lines.size()),
      coordinates(coordinates),
      lines(lines) {
    for (const auto& line : lines) {
        double dx =
            coordinates[line.second].first - coordinates[line.first].first;
        double dy =
            coordinates[line.second].second - coordinates[line.first].second;
        line_distance.push_back(std::sqrt(dx * dx + dy * dy));
    }
    std::cout << "Mesh2D Created" << std::endl;
}

Mesh2D::~Mesh2D() { std::cout << "Mesh2D Deleted" << std::endl; }

size_t Mesh2D::get_num_points() const { return N_pt; }

size_t Mesh2D::get_num_elements() const { return N_el; }

const std::vector<std::pair<double, double>>& Mesh2D::get_coordinates() const {
    return coordinates;
}

const std::vector<std::pair<size_t, size_t>>& Mesh2D::get_lines() const {
    return lines;
}

// const std::vector<std::pair<size_t, size_t>>& Mesh2D::get_faces() const {
//     return faces;
// }

const std::vector<double>& Mesh2D::get_line_distance() const {
    return line_distance;
}

void Mesh2D::write_to_vtk(const std::string& filename) const {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    file << "# vtk DataFile Version 3.0\n";
    file << "Mesh2D Data\n";
    file << "ASCII\n";
    file << "DATASET UNSTRUCTURED_GRID\n";
    file << "POINTS " << N_pt << " float\n";
    for (const auto& coord : coordinates) {
        file << coord.first << " " << coord.second
             << " 0.0\n";  // 2D points with z=0
    }

    size_t num_cells = lines.size();
    file << "CELLS " << num_cells << " " << 3 * num_cells << "\n";
    for (const auto& conn : lines) {
        file << "2 " << conn.first << " " << conn.second << "\n";
    }

    file << "CELL_TYPES " << num_cells << "\n";
    for (size_t i = 0; i < num_cells; ++i) {
        file << "3\n";  // VTK_LINE
    }

    file.close();
    std::cout << "Data written to " << filename << std::endl;
}

RectMesh2D::RectMesh2D(size_t Nx, size_t Ny)
    : Mesh2D(generate_coordinates(Nx, Ny), generate_lines(Nx, Ny)),
      Nx(Nx),
      Ny(Ny) {
    std::cout << "RectMesh2D Created" << std::endl;
}

RectMesh2D::~RectMesh2D() { std::cout << "RectMesh2D Deleted" << std::endl; }

size_t RectMesh2D::get_Nx() const { return Nx; }

size_t RectMesh2D::get_Ny() const { return Ny; }

const std::vector<std::pair<double, double>>& RectMesh2D::get_coordinates()
    const {
    return coordinates;
}

std::vector<std::pair<double, double>> RectMesh2D::generate_coordinates(
    size_t Nx, size_t Ny) {
    // Throw error if Nx or Ny is less than 2
    if (Nx < 2 || Ny < 2) {
        throw std::invalid_argument("Nx and Ny must be greater than 1");
    }

    std::vector<std::pair<double, double>> coords;
    coords.reserve(Nx * Ny);

    double dx = 1.0 / (Nx - 1);
    double dy = 1.0 / (Ny - 1);

    for (size_t j = 0; j < Ny; j++) {
        for (size_t i = 0; i < Nx; i++) {
            coords.emplace_back(i * dx, j * dy);
        }
    }
    return coords;
}

std::vector<std::pair<size_t, size_t>> RectMesh2D::generate_lines(size_t Nx,
                                                                  size_t Ny) {
    std::vector<std::pair<size_t, size_t>> conn;
    conn.reserve((Nx - 1) * (Ny - 1) * 2 + (Nx - 1) + (Ny - 1));

    for (size_t j = 0; j < Ny - 1; j++) {
        for (size_t i = 0; i < Nx - 1; i++) {
            size_t p0 = j * Nx + i;
            size_t p1 = p0 + 1;
            size_t p2 = p0 + Nx;

            conn.push_back({p0, p1});  // Line element
            conn.push_back({p0, p2});  // Line element
        }

        // right edge
        conn.push_back({Nx * j + Nx - 1, Nx * j + 2 * Nx - 1});
    }

    // Add boundary lines on the top edge
    for (size_t i = 0; i < Nx - 1; i++) {
        size_t p0 = Nx * (Ny - 1) + i;
        size_t p1 = p0 + 1;
        conn.push_back({p0, p1});  // Line element
    }

    return conn;
}

}  // namespace fin_diff