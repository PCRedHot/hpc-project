#define CATCH_CONFIG_MAIN
#include "catch2/catch_all.hpp"

#include "matrix.hpp" 
#include "mesh2d.hpp"
#include "discrete_2d.hpp"

using namespace fin_diff;

TEST_CASE("Discretisation2D Creation 1", "[discretisation2d_1]") {
    int Nx = 10;
    int Ny = 10;

    RectMesh2D mesh(Nx, Ny);
    Discretisation2D disc(&mesh);

    auto lhs = disc.get_lhs();
    
    double dx = mesh.get_dx();
    double dy = mesh.get_dy();
    double mul = dx * dx * dy * dy;

    CHECK(lhs.get_n_rows() == Nx * Ny);
    CHECK(lhs.get_n_cols() == Nx * Ny);
    CHECK(disc.get_rhs().get_n_rows() == Nx * Ny);
    CHECK(disc.get_rhs().get_n_cols() == 1);

    CHECK(disc.get_lhs()(0, 0) == 1.0);                         // Boundary Point
    CHECK(disc.get_lhs()(14, 14) == 2 * (dx * dx + dy * dy));   // Point Connected to one boundary point 
    
    CHECK(disc.get_lhs()(14, 15) == -1.0 * dy * dy);            // vertical line
    CHECK(disc.get_lhs()(14, 14 + Nx) == -1.0 * dy * dy);       // horizontal line

    // Forcing Term: val * dx^2 * dy^2
    disc.add_forcing_term(1.0);

    auto rhs = disc.get_rhs();
    CHECK(rhs(0, 0) == 0.0);
    CHECK(rhs(14, 0) == mul);

    
    // Dirichlet BC contribution to line: - val * dx^2 or - val * dy^2
    disc.clear_forcing_term();

    double g = 1.0;
    disc.add_dirichlet_bc(g);

    rhs = disc.get_rhs();

    CHECK(rhs(0, 0) == g);
    CHECK(rhs(14, 0) == -g * dx * dx);

}