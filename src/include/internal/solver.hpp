#ifndef SOLVER_HPP
#define SOLVER_HPP

#include <memory>
#include <vector>

#include "discrete_2d.hpp"
#include "matrix.hpp"
#include "solver_config.hpp"

#ifdef PRECISION_FLOAT
using PrecisionType = float;
#else
using PrecisionType = double;
#endif

namespace fin_diff {

    class Solver {
       public:
        // explicit Solver(Discretisation* d) : disc(std::shared_ptr<Discretisation>(d)) {};
        explicit Solver(std::shared_ptr<Discretisation> d) : disc(std::move(d)) {};
        virtual ~Solver() = default;

        virtual std::vector<PrecisionType> solve() = 0;
        virtual std::vector<PrecisionType> solve(std::vector<PrecisionType> init_u) = 0;

        void set_config(SolverConfig config) { this->config = config; }
        SolverConfig* get_config() { return &config; }

       protected:
        SolverConfig config = SolverConfig();

        std::shared_ptr<Discretisation> disc;
    };

    class JacobiSolver : public Solver {
       public:
        // explicit JacobiSolver(Discretisation* d) : Solver(d) {}
        explicit JacobiSolver(std::shared_ptr<Discretisation> d) : Solver(std::move(d)) {}

        ~JacobiSolver() {};

        std::vector<PrecisionType> solve() override {
            Mesh2D mesh = disc->get_mesh();

            std::vector<PrecisionType> init_u(mesh.get_num_points(), 0.0);
            return this->solve(init_u);
        };

        std::vector<PrecisionType> solve(std::vector<PrecisionType> init_u) override {
            auto A = disc->get_lhs();
            auto b = disc->get_rhs();

#ifdef __DEBUG__
            std::cout << "A: " << std::endl;
            A.print();

            std::cout << "b: " << std::endl;
            b.print();
#endif

            MatrixDiagonal<PrecisionType> D = MatrixDiagonal<PrecisionType>(A);

#ifdef __DEBUG__
            std::cout << "D: " << std::endl;
            D.print();
#endif
            MatrixDiagonal<PrecisionType> D_inv = MatrixDiagonal<PrecisionType>(D);
            D_inv.inv();

            MatrixCRS<PrecisionType> LU = A - D;

            Matrix<PrecisionType> u = Matrix<PrecisionType>(b.get_n_rows(), 1, init_u);
            Matrix<PrecisionType> u_new = Matrix<PrecisionType>(u.get_n_rows(), 1);

            // Solve the system
            const auto tol = config.get_convergence_tol();
            const size_t max_iter = config.get_max_iter();

            size_t n_iter = 0;
            PrecisionType curr_tol = 1e20;
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
                Matrix<PrecisionType> residual = A * u - b;

                // Calculate the norm of the residual
                // curr_tol = residual.norm() / residual.get_n_rows();
                curr_tol = residual.weighted_norm();

                std::cout << "Iteration: " << n_iter
                          << ", Residual: " << curr_tol << std::endl;
            }

            return u.to_vector();
        }
    };
}  // namespace fin_diff

#endif