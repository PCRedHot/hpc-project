#include <cmath>
#include <iostream>
#include <memory>
#include <string>

#include "arg.hpp"
#include "discrete_2d.hpp"
#include "mesh2d.hpp"
#include "param.hpp"
#include "solver.hpp"

using namespace fin_diff;

int main(int argc, char* argv[]) {
#ifdef __DEBUG__
    std::cout << "Debug mode is enabled" << std::endl;
#endif

    // Augment Processing
    AugumentParser parser(argc, argv);

    Parameter param(parser);
    param.print();

    if (param.get_dimension() == 2) {
        // 2D problem
        // Create the mesh
        auto mesh =
            std::make_shared<RectMesh2D>(param.get_num_x(), param.get_num_y());

        // Create the discretisation
        auto disc = std::make_shared<Discretisation2D>(mesh);

        // Add forcing term
        disc->add_forcing_term(param.get_forcing_term_expr());

        // Add Dirichlet BC
        disc->add_dirichlet_bc(param.get_dirichlet_bc_expr());

        std::unique_ptr<Solver> solver = nullptr;
        std::string solver_type = param.get_solver();
        if (solver_type == "jacobi")
            solver = std::make_unique<JacobiSolver>(disc);
        else
            throw std::invalid_argument("Invalid solver: " + solver_type);

        // Create the solver
        solver->set_config(SolverConfig(param));

        // Solve the problem
        auto result = solver->solve();

        // Check if exact solution given
        if (!param.get_exact_sol_expr().empty()) {
            // Calculate the error
            auto exact_sol_expr = param.get_exact_sol_expr();
            auto exact_sol = mesh->get_field_from_expr(exact_sol_expr);
            double inf_norm = 0.0;
            for (size_t i = 0; i < mesh->get_num_points(); i++) {
                double diff = abs(result[i] - exact_sol[i]);
#ifdef __DEBUG__
                std::cout << "Exact: " << exact_sol[i]
                          << " Approx: " << result[i] << " Diff: " << diff
                          << std::endl;
#endif
                if (diff > inf_norm) inf_norm = diff;
            }

            std::cout << "Exact Solution: " << exact_sol_expr << std::endl;
            std::cout << "Exact Error: " << inf_norm << std::endl;
        }

        mesh->write_field_to_vtk(param.get_output_file(), result, "u");

    } else if (param.get_dimension() == 3) {
        throw std::invalid_argument("3D problem not supported yet");
    } else {
        throw std::invalid_argument("Invalid dimension: " +
                                    std::to_string(param.get_dimension()));
    }

    return 0;
}