#ifndef SOLVER_CONFIG_HPP
#define SOLVER_CONFIG_HPP

#include <vector>

#include "discrete_2d.hpp"
#include "matrix.hpp"
#include "param.hpp"

namespace fin_diff {

    class SolverConfig {
       public:
        SolverConfig() = default;
        SolverConfig(const SolverConfig& o) = default;
        SolverConfig(const Parameter& param) {
            set_convergence_tol(param.get_tolerance());
            set_max_iter(param.get_max_iter());
        }

        ~SolverConfig() = default;

        const double get_convergence_tol() const { return convergence_tol; }
        const size_t get_max_iter() const { return max_iter; }

        void set_convergence_tol(double tol) { convergence_tol = tol; }
        void set_max_iter(size_t iter) { max_iter = iter; }

       private:
        double convergence_tol = 1e-6;
        size_t max_iter = 1000;
    };

}  // namespace fin_diff

#endif