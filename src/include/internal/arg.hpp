#ifndef ARG_HPP
#define ARG_HPP

/*
Arguments for the command line
*/

#include "argparse.hpp"

namespace fin_diff {
    class AugumentParser {
       public:
        AugumentParser(int argc, char* argv[]) : parser("Finite Difference Solver") {
            parser.add_argument("-p", "--param")
                .help("Path to parameter file")
                .default_value(std::string(""));

            parser.add_argument("-d", "--dimension")
                .help("Dimension of the problem")
                .choices("2")
                .default_value(2)
                .action([](const std::string& value) { return std::stoi(value); });
            parser.add_argument("-x", "--num_x")
                .help("Number of grid points in x direction")
                .default_value(5)
                .action([](const std::string& value) { return std::stoi(value); });
            parser.add_argument("-y", "--num_y")
                .help("Number of grid points in y direction")
                .default_value(5)
                .action([](const std::string& value) { return std::stoi(value); });
            parser.add_argument("-z", "--num_z")
                .help("Number of grid points in z direction (Only for 3D)")
                .default_value(5)
                .action([](const std::string& value) { return std::stoi(value); });

            parser.add_argument("-s", "--solver")
                .help("Solver to use")
                .choices("jacobi")
                .default_value(std::string("jacobi"));
            parser.add_argument("-t", "--tolerance")
                .help("Convergence tolerance")
                .default_value(1e-6)
                .action([](const std::string& value) { return std::stod(value); });
            parser.add_argument("-m", "--max_iter")
                .help("Maximum number of iterations")
                .default_value(100)
                .action([](const std::string& value) { return std::stoi(value); });

                
            parser.add_argument("-i", "--init-guess")
                .help("Initial guess value function")
                .default_value(std::string("0.0"));
            parser.add_argument("-f", "--forcing-term")
                .help("Forcing term function f(x)")
                .default_value(std::string("0.0"));
            parser.add_argument("-g", "--dirichlet-bc")
                .help("Dirichlet boundary condition function g(x)")
                .default_value(std::string("0.0"));

            parser.add_argument("-o", "--output")
                .help("Output folder")
                .default_value(std::string("./output"));
                
            try {
                parser.parse_args(argc, argv);
            } catch (const std::exception& err) {
                std::cerr << err.what() << std::endl;
                std::cerr << parser;
                std::exit(1);
            }
            
        }

        std::string get_param_file() { return parser.get<std::string>("--param"); }
        int get_dimension() { return parser.get<int>("--dimension"); }
        int get_num_x() { return parser.get<int>("--num_x"); }
        int get_num_y() { return parser.get<int>("--num_y"); }
        int get_num_z() { return parser.get<int>("--num_z"); }
        std::string get_solver() { return parser.get<std::string>("--solver"); }
        double get_tolerance() { return parser.get<double>("--tolerance"); }
        int get_max_iter() { return parser.get<int>("--max_iter"); }
        std::string get_init_guess_func_expr() { return parser.get<std::string>("--init-guess"); }
        std::string get_forcing_term_func_expr() { return parser.get<std::string>("--forcing-term"); }
        std::string get_dirichlet_bc_func_expr() { return parser.get<std::string>("--dirichlet-bc"); }
        std::string get_output_folder() { return parser.get<std::string>("--output"); }

    private:
        argparse::ArgumentParser parser;
    };
}   // namespace fin_diff

#endif