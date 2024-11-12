#ifndef MESH2D_HPP
#define MESH2D_HPP

#include <cstddef>
#include <string>
#include <utility>
#include <vector>

namespace fin_diff {

class Mesh2D {
    // Assume square

   public:
    Mesh2D(const std::vector<std::pair<double, double>> &coordinates,
           const std::vector<std::pair<size_t, size_t>> &lines);
    Mesh2D(const std::vector<std::pair<double, double>> &coordinates,
           const std::vector<std::pair<size_t, size_t>> &lines,
           const std::vector<bool> &boundary);
    virtual ~Mesh2D();

    size_t get_num_points() const;
    size_t get_num_elements() const;

    std::vector<bool> get_boundary() const { return boundary; }
    bool is_boundary(size_t i) const { return boundary[i]; }

    virtual const std::vector<std::pair<double, double>> &get_coordinates() const;
    virtual const std::vector<std::pair<size_t, size_t>> &get_lines() const;
    // virtual const std::vector<std::pair<size_t, size_t>>& get_faces() const; //
    // TODO Face Support for visualization

    // const std::vector<double>& get_line_dx() const { return line_dx; }
    // const std::vector<double>& get_line_dy() const { return line_dy; }

    void write_to_vtk(const std::string &filename) const;


   protected:
    size_t N_pt;
    size_t N_el;

    std::vector<std::pair<double, double>> coordinates;
    std::vector<std::pair<size_t, size_t>> lines;
    // std::vector<std::pair<size_t, size_t>> faces;

    // std::vector<double> line_dx;
    // std::vector<double> line_dy;

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
    // std::vector<std::pair<size_t, size_t>> generate_faces(size_t Nx, size_t
    // Ny);
};

}  // namespace fin_diff

#endif