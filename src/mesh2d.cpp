#include "mesh2d.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <unordered_map>

#ifdef PRECISION_FLOAT
using PrecisionType = float;
#else
using PrecisionType = double;
#endif

namespace fin_diff {

    Mesh2D::Mesh2D(const std::vector<std::pair<PrecisionType, PrecisionType>>& c,
                   const std::vector<std::pair<size_t, size_t>>& l,
                   const std::vector<std::array<size_t, 4>>& f)
        : N_pt(c.size()), N_el(l.size()), coordinates(c), lines(l), faces(f) {
        // Count number of lines connected to each point
        std::unordered_map<size_t, size_t> line_count;

        for (const auto& line : lines) {
            size_t point1 = line.first;
            size_t point2 = line.second;

            line_count[point1]++;
            line_count[point2]++;
        }

        // Set boundary if point is connected to 3 or less lines
        boundary.reserve(N_pt);
        for (size_t i = 0; i < N_pt; i++) {
            bool is_boundary_pt = line_count[i] <= 3;
            boundary.push_back(is_boundary_pt);

            if (is_boundary_pt) {
                num_boundary_points++;
            }
        }
#ifdef __DEBUG__
        std::cout << "Mesh2D Created" << std::endl;
#endif
    }

    Mesh2D::Mesh2D(const std::vector<std::pair<PrecisionType, PrecisionType>>& c,
                   const std::vector<std::pair<size_t, size_t>>& l,
                   const std::vector<std::array<size_t, 4>>& f,
                   const std::vector<bool>& b)
        : N_pt(c.size()),
          N_el(l.size()),
          coordinates(c),
          lines(l),
          faces(f),
          boundary(b) {
        num_boundary_points =
            std::count(boundary.begin(), boundary.end(), true);
#ifdef __DEBUG__
        std::cout << "Mesh2D Created" << std::endl;
#endif
    }

    Mesh2D::~Mesh2D() {
#ifdef __DEBUG__
        std::cout << "Mesh2D Deleted" << std::endl;
#endif
    }

    const std::vector<std::pair<PrecisionType, PrecisionType>>& Mesh2D::get_coordinates()
        const {
        return coordinates;
    }

    const std::vector<std::pair<size_t, size_t>>& Mesh2D::get_lines() const {
        return lines;
    }

    const std::vector<std::array<size_t, 4>>& Mesh2D::get_faces() const {
        return faces;
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

        // Check if faces are present
        if (faces.size() > 0) {
            size_t num_cells = faces.size();
            file << "CELLS " << num_cells << " " << 5 * num_cells << "\n";
            for (const auto& conn : faces) {
                file << "4 " << conn[0] << " " << conn[1] << " " << conn[2]
                     << " " << conn[3] << "\n";
            }

            file << "CELL_TYPES " << num_cells << "\n";
            for (size_t i = 0; i < num_cells; ++i) {
                file << "9\n";  // VTK_QUAD
            }
        } else {
            size_t num_cells = lines.size();
            file << "CELLS " << num_cells << " " << 3 * num_cells << "\n";
            for (const auto& conn : lines) {
                file << "2 " << conn.first << " " << conn.second << "\n";
            }

            file << "CELL_TYPES " << num_cells << "\n";
            for (size_t i = 0; i < num_cells; ++i) {
                file << "3\n";  // VTK_LINE
            }
        }

        file.close();
#ifdef __DEBUG__
        std::cout << "Data written to " << filename << std::endl;
#endif
    }

    void Mesh2D::write_field_to_vtk(const std::string& filename,
                                    const std::vector<PrecisionType>& field,
                                    const std::string field_name) const {
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

        if (faces.size() > 0) {
            size_t num_cells = faces.size();
            file << "CELLS " << num_cells << " " << 5 * num_cells << "\n";
            for (const auto& conn : faces) {
                file << "4 " << conn[0] << " " << conn[1] << " " << conn[2]
                     << " " << conn[3] << "\n";
            }

            file << "CELL_TYPES " << num_cells << "\n";
            for (size_t i = 0; i < num_cells; ++i) {
                file << "9\n";  // VTK_QUAD
            }
        } else {
            size_t num_cells = lines.size();
            file << "CELLS " << num_cells << " " << 3 * num_cells << "\n";
            for (const auto& conn : lines) {
                file << "2 " << conn.first << " " << conn.second << "\n";
            }

            file << "CELL_TYPES " << num_cells << "\n";
            for (size_t i = 0; i < num_cells; ++i) {
                file << "3\n";  // VTK_LINE
            }
        }

        file << "POINT_DATA " << N_pt << "\n";
        file << "SCALARS " << field_name << " float 1\n";
        file << "LOOKUP_TABLE default\n";
        for (const auto& val : field) {
            file << val << "\n";
        }

        file.close();
    }

    std::vector<PrecisionType> Mesh2D::get_field_from_expr(
        const std::string& expr) const {
#ifdef __DEBUG__
        std::cout << "Generating field from expression: " << expr << std::endl;
#endif
        PrecisionType x, y;

        exprtk::symbol_table<PrecisionType> symbol_table;
        symbol_table.add_variable("x", x);
        symbol_table.add_variable("y", y);
        symbol_table.add_constants();

        exprtk::expression<PrecisionType> expression;
        expression.register_symbol_table(symbol_table);

        exprtk::parser<PrecisionType> parser;
        parser.compile(expr, expression);

        std::vector<PrecisionType> field(N_pt, 0.0);
        for (size_t i = 0; i < N_pt; i++) {
            x = coordinates[i].first;
            y = coordinates[i].second;
            field[i] = expression.value();
        }

        return field;
    };

    RectMesh2D::RectMesh2D(size_t Nx, size_t Ny)
        : Mesh2D(generate_coordinates(Nx, Ny), generate_lines(Nx, Ny),
                 generate_faces(Nx, Ny)),
          Nx(Nx),
          Ny(Ny),
          dx(1.0 / (Nx - 1)),
          dy(1.0 / (Ny - 1)) {
        this->boundary.reserve(Nx * Ny);
        for (size_t i = 0; i < Nx * Ny; i++) {
            bool is_boundary_pt = (i / Nx == 0 || i % Nx == 0 ||
                                   i % Nx == Nx - 1 || i / Nx == Ny - 1);
            this->boundary.push_back(is_boundary_pt);

            if (is_boundary_pt) {
                num_boundary_points++;
            }
        }

        this->N_pt = Nx * Ny;
        this->N_el = this->lines.size();
#ifdef __DEBUG__
        std::cout << "RectMesh2D Created" << std::endl;
#endif
    }

    RectMesh2D::~RectMesh2D() {
#ifdef __DEBUG__
        std::cout << "RectMesh2D Deleted" << std::endl;
#endif
    }

    const std::vector<std::pair<PrecisionType, PrecisionType>>& RectMesh2D::get_coordinates()
        const {
        return coordinates;
    }

    std::vector<std::pair<PrecisionType, PrecisionType>> RectMesh2D::generate_coordinates(
        size_t Nx, size_t Ny) {
        // Throw error if Nx or Ny is less than 2
        if (Nx < 2 || Ny < 2) {
            throw std::invalid_argument("Nx and Ny must be greater than 1");
        }

        std::vector<std::pair<PrecisionType, PrecisionType>> coords;
        coords.reserve(Nx * Ny);

        PrecisionType dx = 1.0 / (Nx - 1);
        PrecisionType dy = 1.0 / (Ny - 1);

        for (size_t j = 0; j < Ny; j++) {
            for (size_t i = 0; i < Nx; i++) {
                coords.emplace_back(i * dx, j * dy);
            }
        }
        return coords;
    }

    std::vector<std::pair<size_t, size_t>> RectMesh2D::generate_lines(
        size_t Nx, size_t Ny) {
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

    std::vector<std::array<size_t, 4>> RectMesh2D::generate_faces(size_t Nx,
                                                                  size_t Ny) {
        std::vector<std::array<size_t, 4>> faces;
        faces.reserve((Nx - 1) * (Ny - 1));

        for (size_t j = 0; j < Ny - 1; j++) {
            for (size_t i = 0; i < Nx - 1; i++) {
                size_t p0 = j * Nx + i;
                size_t p1 = p0 + 1;
                size_t p2 = p0 + Nx;
                size_t p3 = p2 + 1;

                faces.push_back({p0, p1, p3, p2});
            }
        }

        return faces;
    }

}  // namespace fin_diff