#include <iostream>
#include <memory>
#include <string>

#include "arg.hpp"
#include "param.hpp"
#include "discrete_2d.hpp"	
#include "solver.hpp"
#include "mesh2d.hpp"

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
        auto mesh = std::make_shared<RectMesh2D>(param.get_num_x(), param.get_num_y());
        
        // Create the discretisation
        auto disc = std::make_shared<Discretisation2D>(mesh);

        // Add forcing term
        disc->add_forcing_term(param.get_forcing_term_expr());

        // Add Dirichlet BC
        disc->add_dirichlet_bc(param.get_dirichlet_bc_expr());
        
        std::unique_ptr<Solver> solver = nullptr;
        std::string solver_type = param.get_solver();
        if (solver_type == "jacobi") solver = std::make_unique<JacobiSolver>(disc);
        else throw std::invalid_argument("Invalid solver: " + solver_type);
        
        // Create the solver
        solver->set_config(SolverConfig(param));

        // Solve the problem
        auto result = solver->solve();

        mesh->write_field_to_vtk(param.get_output_file(), result, "u");

    } else if (param.get_dimension() == 3) {
        throw std::invalid_argument("3D problem not supported yet");
    } else {
        throw std::invalid_argument("Invalid dimension: " + std::to_string(param.get_dimension()));
    }

    return 0;
}