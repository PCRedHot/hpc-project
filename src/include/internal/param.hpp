#ifndef PARAM_HPP
#define PARAM_HPP

#include <fstream>
#include <iostream>
#include <string>

#include "arg.hpp"
#include "helper.hpp"

#define _PARAM_PRINT_KEY_WIDTH 21
#define _PARAM_PRINT_VAL_WIDTH 31

#define _PARAM_DIMENSION_KEY "DIMENSION"
#define _PARAM_NUM_X_KEY "NUM_X"
#define _PARAM_NUM_Y_KEY "NUM_Y"
#define _PARAM_NUM_Z_KEY "NUM_Z"
#define _PARAM_SOLVER_KEY "SOLVER"
#define _PARAM_TOLERANCE_KEY "TOLERANCE"
#define _PARAM_MAX_ITER_KEY "MAX_ITER"
#define _PARAM_FORCING_TERM_KEY "FORCING_TERM"
#define _PARAM_DIRICHLET_BC_KEY "DIRICHLET_BC"
#define _PARAM_INIT_GUESS_KEY "INIT_GUESS"
#define _PARAM_OUTPUT_FILE_KEY "OUTPUT"

using namespace fin_diff;

namespace fin_diff {

    class Parameter {
       public:
        Parameter(AugumentParser& parser) {
            // Check if the parameter file is provided
            std::string param_file = parser.get_param_file();
            if (!param_file.empty()) {
                *this = Parameter(param_file);
                return;
            }

            // Parse the arguments
            dimension = parser.get_dimension();
            num_x = parser.get_num_x();
            num_y = parser.get_num_y();
            num_z = parser.get_num_z();
            solver = parser.get_solver();
            tolerance = parser.get_tolerance();
            max_iter = parser.get_max_iter();
            init_guess_expr = parser.get_init_guess_func_expr();

            output_file = parser.get_output_file();

            check_valid();
        }

        Parameter(const std::string& param_file) {
#ifdef __DEBUG__
            show_cwd();
#endif
            // Parse the parameter file

            std::ifstream file(param_file);

            if (!file.is_open()) {
                std::cerr << "Error: Parameter file not found: " << param_file << std::endl;
                throw std::invalid_argument("Parameter file not found");
            }

            std::string line;
            while (std::getline(file, line)) {
                // Erase leading and trailing whitespaces
                line.erase(line.begin(), std::find_if(line.begin(), line.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
                line.erase(std::find_if(line.rbegin(), line.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), line.end());

                // Erase comments with #
                size_t pos = line.find("#");
                if (pos != std::string::npos) {
                    line.erase(pos);
                }

                // Erase comments with //
                pos = line.find("//");
                if (pos != std::string::npos) {
                    line.erase(pos);
                }

                // Erase trailing whitespaces again
                line.erase(std::find_if(line.rbegin(), line.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), line.end());

                // Skip empty lines
                if (line.empty()) {
                    continue;
                }

#ifdef __DEBUG__
                std::cout << "Line: " << line << std::endl;
#endif

                // Split the line into key and value
                pos = line.find(" ");

                if (pos == std::string::npos) {
                    std::cerr << "Error: Invalid parameter file format" << std::endl;
                    std::cerr << "Line: " << line << std::endl;
                }

                std::string key = line.substr(0, pos);
                std::string value = line.substr(pos + 1);

                try {
                    if (key == _PARAM_DIMENSION_KEY)
                        dimension = std::stoi(value);
                    else if (key == _PARAM_NUM_X_KEY)
                        num_x = std::stoi(value);
                    else if (key == _PARAM_NUM_Y_KEY)
                        num_y = std::stoi(value);
                    else if (key == _PARAM_NUM_Z_KEY)
                        num_z = std::stoi(value);
                    else if (key == _PARAM_SOLVER_KEY)
                        solver = value;
                    else if (key == _PARAM_TOLERANCE_KEY)
                        tolerance = std::stod(value);
                    else if (key == _PARAM_MAX_ITER_KEY)
                        max_iter = std::stoi(value);
                    else if (key == _PARAM_INIT_GUESS_KEY)
                        init_guess_expr = value;
                    else if (key == _PARAM_FORCING_TERM_KEY)
                        forcing_term_expr = value;
                    else if (key == _PARAM_DIRICHLET_BC_KEY)
                        dirichlet_bc_expr = value;
                    else if (key == _PARAM_OUTPUT_FILE_KEY)
                        output_file = value;
                    else {
                        std::cerr << "Error: Invalid key in parameter file" << std::endl;
                        throw std::invalid_argument("Invalid key in parameter file");
                    }
                } catch (const std::invalid_argument& err) {
                    std::cerr << "Error: Invalid parsing value in parameter file" << std::endl;
                    std::cerr << "Key: " << key << ", Value: " << value << std::endl;
                    file.close();

                    throw std::invalid_argument("Invalid parsing value in parameter file");
                }
            }

            file.close();

            try {
                check_valid();
            } catch (const std::invalid_argument& err) {
                std::cerr << "Error: Invalid parameter file" << std::endl;
                throw std::invalid_argument("Invalid parameter file");
            }
        }

        int get_dimension() const { return dimension; }
        int get_num_x() const { return num_x; }
        int get_num_y() const { return num_y; }
        int get_num_z() const { return num_z; }
        std::string get_solver() const { return solver; }
        double get_tolerance() const { return tolerance; }
        int get_max_iter() const { return max_iter; }
        std::string get_init_guess_expr() const { return init_guess_expr; }
        std::string get_forcing_term_expr() const { return forcing_term_expr; }
        std::string get_dirichlet_bc_expr() const { return dirichlet_bc_expr; }
        std::string get_output_file() const { return output_file; }

        void print() const {
            std::cout << '|' << std::setw(_PARAM_PRINT_KEY_WIDTH) << centered("Parameter") << "| " << std::setw(_PARAM_PRINT_VAL_WIDTH) << centered("Value") << std::endl;
            std::cout << '|' << std::setw(_PARAM_PRINT_KEY_WIDTH) << std::setfill('-') << "" << "| " << std::setw(_PARAM_PRINT_VAL_WIDTH) << "" << std::endl
                      << std::setfill(' ');
            std::cout << '|' << std::setw(_PARAM_PRINT_KEY_WIDTH) << centered("Dimension") << "| " << dimension << std::endl;
            std::cout << '|' << std::setw(_PARAM_PRINT_KEY_WIDTH) << centered("Nx") << "| " << num_x << std::endl;
            std::cout << '|' << std::setw(_PARAM_PRINT_KEY_WIDTH) << centered("Ny") << "| " << num_y << std::endl;
            std::cout << '|' << std::setw(_PARAM_PRINT_KEY_WIDTH) << centered("Nz") << "| " << num_z << std::endl;
            std::cout << '|' << std::setw(_PARAM_PRINT_KEY_WIDTH) << centered("Solver") << "| " << solver << std::endl;
            std::cout << '|' << std::setw(_PARAM_PRINT_KEY_WIDTH) << centered("Tolerance") << "| " << tolerance << std::endl;
            std::cout << '|' << std::setw(_PARAM_PRINT_KEY_WIDTH) << centered("Maximum iterations") << "| " << max_iter << std::endl;
            std::cout << '|' << std::setw(_PARAM_PRINT_KEY_WIDTH) << centered("Initial guess") << "| " << init_guess_expr << std::endl;
            std::cout << '|' << std::setw(_PARAM_PRINT_KEY_WIDTH) << centered("Forcing term") << "| " << forcing_term_expr << std::endl;
            std::cout << '|' << std::setw(_PARAM_PRINT_KEY_WIDTH) << centered("Dirichlet BC") << "| " << dirichlet_bc_expr << std::endl;
            std::cout << '|' << std::setw(_PARAM_PRINT_KEY_WIDTH) << centered("Output file") << "| " << output_file << std::endl;
        }

       private:
        int dimension;
        int num_x;
        int num_y;
        int num_z = 0;
        std::string solver;
        double tolerance;
        int max_iter;

        std::string init_guess_expr;
        std::string forcing_term_expr;
        std::string dirichlet_bc_expr;

        std::string output_file;

        void check_valid() {
            // Check if the dimension is valid
            if (dimension != 2 && dimension != 3) {
                std::cerr << "Error: Invalid dimension: " << dimension << std::endl;
                throw std::invalid_argument("Invalid dimension");
            }

            if ((dimension == 3 && num_z == 0) || (dimension == 2 && num_z != 0)) {
                std::cerr << "Error: Invalid number of z: " << num_z << " in " << dimension << "D" << std::endl;
                throw std::invalid_argument("Invalid number of z");
            }

            // Check if the solver is valid
            if (solver != "jacobi") {
                std::cerr << "Error: Invalid solver: " << solver << std::endl;
                throw std::invalid_argument("Invalid solver");
            }

            // Check if the tolerance is valid
            if (tolerance <= 0) {
                std::cerr << "Error: Invalid tolerance: " << tolerance << std::endl;
                throw std::invalid_argument("Invalid tolerance");
            }

            // Check if the maximum number of iterations is valid
            if (max_iter <= 0) {
                std::cerr << "Error: Invalid maximum number of iterations: " << max_iter << std::endl;
                throw std::invalid_argument("Invalid maximum number of iterations");
            }

            // Check if the output folder is valid
            if (output_file.empty()) {
                std::cerr << "Error: Empty output folder" << std::endl;
                throw std::invalid_argument("Empty output folder");
            }
        }
    };

}  // namespace fin_diff

#endif