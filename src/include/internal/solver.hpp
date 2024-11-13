#ifndef SOLVER_HPP
#define SOLVER_HPP

#include <vector>

#include "discrete_2d.hpp"
#include "matrix.hpp"
#include "solver_config.hpp"

namespace fin_diff {

    class Solver {
       public:
        Solver() {};
        ~Solver() {};

        virtual std::vector<double> solve() = 0;
        virtual std::vector<double> solve(std::vector<double> init_u) = 0;

        void set_config(SolverConfig config) { this->config = config; }
        SolverConfig* get_config() { return &config; }

       protected:
        SolverConfig config = SolverConfig();
    };

    class JacobiSolver : public Solver {
       public:
        JacobiSolver(Discretisation* disc) : disc(disc) {}

        std::vector<double> solve() override {
            Mesh2D mesh = disc->get_mesh();

            std::vector<double> init_u(mesh.get_num_points(), 0.0);
            return this->solve(init_u);
        };

        std::vector<double> solve(std::vector<double> init_u) override {
            auto A = disc->get_lhs();
            auto b = disc->get_rhs();

#ifdef __DEBUG__
            std::cout << "A: " << std::endl;
            A.print();

            std::cout << "b: " << std::endl;
            b.print();
#endif

            MatrixDiagonal<double> D = MatrixDiagonal<double>(A);

#ifdef __DEBUG__
            std::cout << "D: " << std::endl;
            D.print();
#endif
            MatrixDiagonal<double> D_inv = MatrixDiagonal<double>(D);
            D_inv.inv();

            MatrixCRS<double> LU = A - D;

            Matrix<double> u = Matrix<double>(b.get_n_rows(), 1, init_u);
            Matrix<double> u_new = Matrix<double>(u.get_n_rows(), 1);

            // Solve the system
            const double tol = config.get_convergence_tol();
            const size_t max_iter = config.get_max_iter();

            size_t n_iter = 0;
            double curr_tol = 1e20;
#ifdef __DEBUG__
            std::cout << "D_inv: " << std::endl;
            D_inv.print();

            std::cout << "LU: " << std::endl;
            LU.print();
#endif
            while (n_iter < max_iter && curr_tol > tol) {
                u_new = D_inv * (b - LU * u);

                u = u_new;
                n_iter++;

                // Calculate the residual
                Matrix<double> residual = A * u - b;

                // Calculate the norm of the residual
                curr_tol = residual.norm() / residual.get_n_rows();

#ifdef __DEBUG__
                std::cout << "Iteration: " << n_iter
                          << ", Residual: " << curr_tol << std::endl;
#endif
            }

            return u.to_vector();
        }

       private:
        Discretisation* disc;
    };
}  // namespace fin_diff

#endif