#ifndef MESH2D_HPP
#define MESH2D_HPP

#include <cstddef>
#include <string>
#include <utility>
#include <vector>

namespace fin_diff {

class Mesh2D {
   public:
    Mesh2D(const std::vector<std::pair<double, double>>& coordinates,
           const std::vector<std::pair<size_t, size_t>>& lines);
    virtual ~Mesh2D();

    size_t get_num_points() const;
    size_t get_num_elements() const;

    virtual const std::vector<std::pair<double, double>>& get_coordinates()
        const;
    virtual const std::vector<std::pair<size_t, size_t>>& get_lines() const;
    // virtual const std::vector<std::pair<size_t, size_t>>& get_faces() const;
    virtual const std::vector<double>& get_line_distance() const;

    void write_to_vtk(const std::string& filename) const;

   protected:
    size_t N_pt;
    size_t N_el;

    std::vector<std::pair<double, double>> coordinates;
    std::vector<std::pair<size_t, size_t>> lines;
    // std::vector<std::pair<size_t, size_t>> faces;

    std::vector<double> line_distance;
};

class RectMesh2D : public Mesh2D {
   public:
    RectMesh2D(size_t Nx, size_t Ny);
    ~RectMesh2D();

    size_t get_Nx() const;
    size_t get_Ny() const;

    const std::vector<std::pair<double, double>>& get_coordinates()
        const override;

   private:
    size_t Nx, Ny;

    std::vector<std::pair<double, double>> generate_coordinates(size_t Nx,
                                                                size_t Ny);
    std::vector<std::pair<size_t, size_t>> generate_lines(size_t Nx, size_t Ny);
    // std::vector<std::pair<size_t, size_t>> generate_faces(size_t Nx, size_t
    // Ny);
};

}  // namespace fin_diff

#endif