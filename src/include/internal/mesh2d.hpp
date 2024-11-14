#ifndef MESH2D_HPP
#define MESH2D_HPP

#include <cstddef>
#include <string>
#include <utility>
#include <vector>

#include "exprtk.hpp"

namespace fin_diff {

    class Mesh2D {
        // Assume square

       public:
        Mesh2D(const std::vector<std::pair<double, double>> &coordinates,
               const std::vector<std::pair<size_t, size_t>> &lines) : Mesh2D(coordinates, lines, std::vector<std::array<size_t, 4>>()) {};
        Mesh2D(const std::vector<std::pair<double, double>> &coordinates,
               const std::vector<std::pair<size_t, size_t>> &lines,
               const std::vector<bool> &boundary) : Mesh2D(coordinates, lines, std::vector<std::array<size_t, 4>>(), boundary) {};
        Mesh2D(const std::vector<std::pair<double, double>> &coordinates,
               const std::vector<std::pair<size_t, size_t>> &lines,
               const std::vector<std::array<size_t, 4>>& faces);
        Mesh2D(const std::vector<std::pair<double, double>> &coordinates,
               const std::vector<std::pair<size_t, size_t>> &lines,
               const std::vector<std::array<size_t, 4>>& faces,
               const std::vector<bool> &boundary);
        virtual ~Mesh2D();

        size_t get_num_points() const { return N_pt; }
        size_t get_num_elements() const { return N_el; }

        std::vector<bool> get_boundary() const { return boundary; }
        bool is_boundary(size_t i) const { return boundary[i]; }

        virtual const std::vector<std::pair<double, double>> &get_coordinates() const;
        virtual const std::vector<std::pair<size_t, size_t>> &get_lines() const;
        virtual const std::vector<std::array<size_t, 4>> &get_faces() const;

        void write_to_vtk(const std::string &filename) const;

        std::vector<double> get_empty_field() const {
            return std::vector<double>(N_pt, 0.0);
        }

        std::vector<double> get_field_from_expr(const std::string &expr) const;

        void write_field_to_vtk(const std::string &filename,
                                const std::vector<double> &field,
                                const std::string field_name) const;

       protected:
        size_t N_pt;
        size_t N_el;

        std::vector<std::pair<double, double>> coordinates;
        std::vector<std::pair<size_t, size_t>> lines;
        std::vector<std::array<size_t, 4>> faces;

        std::vector<bool> boundary;
        size_t num_boundary_points = 0;
    };

    class RectMesh2D : public Mesh2D {
       public:
        RectMesh2D(size_t Nx, size_t Ny);
        ~RectMesh2D();

        size_t get_Nx() const { return Nx; }
        size_t get_Ny() const { return Ny; }

        double get_dx() const { return dx; }
        double get_dy() const { return dy; }

        const std::vector<std::pair<double, double>> &get_coordinates()
            const override;

       private:
        size_t Nx, Ny;
        double dx, dy;

        std::vector<std::pair<double, double>> generate_coordinates(size_t Nx,
                                                                    size_t Ny);
        std::vector<std::pair<size_t, size_t>> generate_lines(size_t Nx, size_t Ny);
        std::vector<std::array<size_t, 4>> generate_faces(size_t Nx, size_t Ny);
    };

}  // namespace fin_diff

#endif