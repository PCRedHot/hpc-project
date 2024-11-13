#define CATCH_CONFIG_MAIN
#include "catch2/catch_all.hpp"

#include <vector>
#include <random>

#include "solver.hpp"
#include "solver_config.hpp"
#include "discrete_2d.hpp"
#include "mesh2d.hpp"



using namespace fin_diff;

TEST_CASE("Test Problem 1", "[test_problem_1]") {
    // Random
    // Initialize random number generator
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 1.0);


    // Define the problem
    auto f = [](double x, double y) { return 0.0; };
    auto g = [](double x, double y) { return 0.0; };


    // Define the mesh
    int Nx = 10;
    int Ny = 10;
    RectMesh2D mesh(Nx, Ny);

    std::vector<double> u_init = std::vector<double>(Nx * Ny);
    for (size_t i = 0; i < u_init.size(); i++) {
        u_init[i] = dis(gen);
    }

    // Define the discretisation
    Discretisation2D disc(&mesh);

    disc.add_forcing_term(f);
    disc.add_dirichlet_bc(g);

    // Define the solver
    JacobiSolver solver(&disc);
    SolverConfig config;
    config.set_max_iter(1000);
    config.set_convergence_tol(1e-6);
    solver.set_config(config);

    // Solve the problem
    auto u = solver.solve();

    std::cout << "u: " << std::endl;
    for (size_t i = 0; i < u.size(); i++) {
        std::cout << u[i] << " ";
    }
    std::cout << std::endl;
}

TEST_CASE("Test Problem 2", "[test_problem_2]") {
    // Define the problem
    // auto f = "8 * pi * pi * sin(2 * pi * x) * sin(2 * pi * y)";
    // auto g = "sin(2 * pi * x) * sin(2 * pi * y)";

    auto f = [](double x, double y) { return 8 * M_PI * M_PI * std::sin(2 * M_PI * x) * std::sin(2 * M_PI * y); };
    auto g = [](double x, double y) { return std::sin(2 * M_PI * x) * std::sin(2 * M_PI * y); };


    // Define the mesh
    int Nx = 10;
    int Ny = 10;
    RectMesh2D mesh(Nx, Ny);

    std::vector<double> u_init = std::vector<double>(Nx * Ny, 0.0);

    // Define the discretisation
    Discretisation2D disc(&mesh);

    disc.add_forcing_term(f);
    disc.add_dirichlet_bc(g);

    // Define the solver
    JacobiSolver solver(&disc);
    SolverConfig config;
    config.set_max_iter(1000);
    config.set_convergence_tol(1e-6);
    solver.set_config(config);

    // Solve the problem
    auto u = solver.solve();

    std::cout << "u: " << std::endl;
    for (size_t i = 0; i < u.size(); i++) {
        std::cout << u[i] << " ";
    }
    std::cout << std::endl;
}

TEST_CASE("Test Problem 2 - expr", "[test_problem_2_expr]") {
    // Define the problem
    auto f = "8 * pi * pi * sin(2 * pi * x) * sin(2 * pi * y)";
    auto g = "sin(2 * pi * x) * sin(2 * pi * y)";

    // auto f = [](double x, double y) { return 8 * M_PI * M_PI * std::sin(2 * M_PI * x) * std::sin(2 * M_PI * y); };
    // auto g = [](double x, double y) { return std::sin(2 * M_PI * x) * std::sin(2 * M_PI * y); };


    // Define the mesh
    int Nx = 10;
    int Ny = 10;
    RectMesh2D mesh(Nx, Ny);

    std::vector<double> u_init = std::vector<double>(Nx * Ny, 0.0);

    // Define the discretisation
    Discretisation2D disc(&mesh);

    disc.add_forcing_term(f);
    disc.add_dirichlet_bc(g);

    // Define the solver
    JacobiSolver solver(&disc);
    SolverConfig config;
    config.set_max_iter(1000);
    config.set_convergence_tol(1e-6);
    solver.set_config(config);

    // Solve the problem
    auto u = solver.solve();

    std::cout << "u: " << std::endl;
    for (size_t i = 0; i < u.size(); i++) {
        std::cout << u[i] << " ";
    }
    std::cout << std::endl;
}

TEST_CASE("Test Problem 3", "[test_problem_3]") {
    // Define the problem
    // auto f = "3 / 2 * sin(2 * x) * sinh(y)";
    // auto g = "1 / 2 * sin(2 * x) * sinh(y)";

    auto f = [](double x, double y) { return 3.0 / 2.0 * std::sin(2 * x) * std::sinh(y); };
    auto g = [](double x, double y) { return 1.0 / 2.0 * std::sin(2 * x) * std::sinh(y); };


    // Define the mesh
    int Nx = 10;
    int Ny = 10;
    RectMesh2D mesh(Nx, Ny);

    std::vector<double> u_init = std::vector<double>(Nx * Ny, 0.0);

    // Define the discretisation
    Discretisation2D disc(&mesh);

    disc.add_forcing_term(f);
    disc.add_dirichlet_bc(g);

    // Define the solver
    JacobiSolver solver(&disc);
    SolverConfig config;
    config.set_max_iter(1000);
    config.set_convergence_tol(1e-6);
    solver.set_config(config);

    // Solve the problem
    auto u = solver.solve();

    std::cout << "u: " << std::endl;
    for (size_t i = 0; i < u.size(); i++) {
        std::cout << u[i] << " ";
    }
    std::cout << std::endl;
}

TEST_CASE("Test Problem 3 - expr", "[test_problem_3_expr]") {
    // Define the problem
    auto f = "3 / 2 * sin(2 * x) * sinh(y)";
    auto g = "1 / 2 * sin(2 * x) * sinh(y)";

    // auto f = [](double x, double y) { return 3.0 / 2.0 * std::sin(2 * x) * std::sinh(y); };
    // auto g = [](double x, double y) { return 1.0 / 2.0 * std::sin(2 * x) * std::sinh(y); };


    // Define the mesh
    int Nx = 10;
    int Ny = 10;
    RectMesh2D mesh(Nx, Ny);

    std::vector<double> u_init = std::vector<double>(Nx * Ny, 0.0);

    // Define the discretisation
    Discretisation2D disc(&mesh);

    disc.add_forcing_term(f);
    disc.add_dirichlet_bc(g);

    // Define the solver
    JacobiSolver solver(&disc);
    SolverConfig config;
    config.set_max_iter(1000);
    config.set_convergence_tol(1e-6);
    solver.set_config(config);

    // Solve the problem
    auto u = solver.solve();

    std::cout << "u: " << std::endl;
    for (size_t i = 0; i < u.size(); i++) {
        std::cout << u[i] << " ";
    }
    std::cout << std::endl;
}

TEST_CASE("Test Problem 3 - heavy", "[test_problem_3_heavy]") {
    // Define the problem
    // auto f = "3 / 2 * sin(2 * x) * sinh(y)";
    // auto g = "1 / 2 * sin(2 * x) * sinh(y)";

    auto f = [](double x, double y) { return 3.0 / 2.0 * std::sin(2 * x) * std::sinh(y); };
    auto g = [](double x, double y) { return 1.0 / 2.0 * std::sin(2 * x) * std::sinh(y); };


    // Define the mesh
    int Nx = 20;
    int Ny = 20;
    RectMesh2D mesh(Nx, Ny);

    std::vector<double> u_init = std::vector<double>(Nx * Ny, 0.0);

    // Define the discretisation
    Discretisation2D disc(&mesh);

    disc.add_forcing_term(f);
    disc.add_dirichlet_bc(g);

    // Define the solver
    JacobiSolver solver(&disc);
    SolverConfig config;
    config.set_max_iter(1000);
    config.set_convergence_tol(1e-6);
    solver.set_config(config);

    // Solve the problem
    auto u = solver.solve();

    std::cout << "u: " << std::endl;
    for (size_t i = 0; i < u.size(); i++) {
        std::cout << u[i] << " ";
    }
    std::cout << std::endl;
}

