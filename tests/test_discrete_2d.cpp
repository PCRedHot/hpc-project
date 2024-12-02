#define CATCH_CONFIG_MAIN
#include "catch2/catch_all.hpp"
#include <catch2/catch_approx.hpp>

#include "discrete_2d.hpp"
#include "matrix.hpp"
#include "mesh2d.hpp"

#ifdef PRECISION_FLOAT
using PrecisionType = float;
#else
using PrecisionType = double;
#endif

using namespace fin_diff;

TEST_CASE("Discretisation2D Creation 1", "[discretisation2d_1]") {
    int Nx = 10;
    int Ny = 10;

    auto mesh = std::make_shared<RectMesh2D>(Nx, Ny);
    auto disc = std::make_shared<Discretisation2D>(mesh);

    auto lhs = disc->get_lhs();

    PrecisionType dx = mesh->get_dx();
    PrecisionType dy = mesh->get_dy();
    PrecisionType mul = dx * dx * dy * dy;

    CHECK(lhs.get_n_rows() == Nx * Ny);
    CHECK(lhs.get_n_cols() == Nx * Ny);
    CHECK(disc->get_rhs().get_n_rows() == Nx * Ny);
    CHECK(disc->get_rhs().get_n_cols() == 1);

    CHECK(disc->get_lhs()(0, 0) == Catch::Approx(1.0));  // Boundary Point
    CHECK(
        disc->get_lhs()(14, 14) ==
        Catch::Approx(
            2 * (dx * dx + dy * dy)));  // Point Connected to one boundary point

    CHECK(disc->get_lhs()(14, 15) ==
          Catch::Approx(-1.0 * dy * dy));  // vertical line
    CHECK(disc->get_lhs()(14, 14 + Nx) ==
          Catch::Approx(-1.0 * dy * dy));  // horizontal line

    // Forcing Term: val * dx^2 * dy^2
    disc->add_forcing_term(1.0);

    auto rhs = disc->get_rhs();

    CHECK(rhs(0, 0) == Catch::Approx(0.0));
    CHECK(rhs(14, 0) == Catch::Approx(mul));

    // Dirichlet BC contribution to line: - val * dx^2 or - val * dy^2
    disc->clear_forcing_term();

    PrecisionType g = 1.0;
    disc->add_dirichlet_bc(g);

    rhs = disc->get_rhs();

    CHECK(rhs(0, 0) == Catch::Approx(g));
    CHECK(rhs(14, 0) == Catch::Approx(-g * dx * dx));
}