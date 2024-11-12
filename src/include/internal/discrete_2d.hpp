#ifndef LHS2D_HPP
#define LHS2D_HPP

#include <functional>

#include "matrix.hpp"
#include "mesh2d.hpp"

namespace fin_diff {

    class Discretisation {
    public:
        virtual ~Discretisation() = default;

        virtual MatrixCRS<double> get_lhs() = 0;
        virtual Matrix<double> get_rhs() = 0;
    };

    class Discretisation2D : public Discretisation {
       public:
        Discretisation2D(const RectMesh2D* mesh);

        ~Discretisation2D();

        MatrixCRS<double> get_lhs() override{ 
            _calculate_lhs();
            return lhs;
        }
        Matrix<double> get_rhs() override{ 
            _calculate_rhs();
            return rhs;
        }

        void add_forcing_term(double val);
        void add_forcing_term(const std::function<double(double, double)>& func);
        void add_forcing_term(const std::string& expr);

        void add_dirichlet_bc(double val);
        void add_dirichlet_bc(const std::function<double(double, double)>& func);
        void add_dirichlet_bc(const std::string& expr);

        void clear_forcing_term();
        void clear_dirichlet_bc();

       private:
        const RectMesh2D* mesh;

        MatrixCRS<double> lhs;
        Matrix<double> rhs;

        std::vector<double> forcing_terms_val;
        std::vector<std::function<double(double, double)>> forcing_terms_func;
        std::vector<std::string> forcing_terms_expr;

        std::vector<double> dirichlet_bcs_val;
        std::vector<std::function<double(double, double)>> dirichlet_bcs_func;
        std::vector<std::string> dirichlet_bcs_expr;

        void _calculate() {
            // Check if mesh is set
            if (mesh->get_Nx() == 0 || mesh->get_Ny() == 0) {
                throw std::runtime_error("Mesh is not set");
            }
            
            _calculate_rhs();
            _calculate_lhs();
        };

        void _calculate_rhs(const bool force = false);
        void _calculate_lhs(const bool force = false);

        bool rhs_ready = false;
        bool lhs_ready = false;
    };

}  // namespace fin_diff

#endif