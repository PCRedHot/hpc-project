#ifndef LHS2D_HPP
#define LHS2D_HPP

#ifdef PRECISION_FLOAT
using PrecisionType = float;
#else
using PrecisionType = double;
#endif

#include <functional>
#include <memory>

#include "matrix.hpp"
#include "mesh2d.hpp"

namespace fin_diff {

    class Discretisation {
       public:
        // explicit Discretisation(Mesh2D* m) :
        // mesh(std::shared_ptr<Mesh2D>(m)), lhs(mesh->get_num_points(),
        // mesh->get_num_points()), rhs(mesh->get_num_points(), 1) {}
        explicit Discretisation(std::shared_ptr<Mesh2D> m)
            : mesh(std::move(m)),
              lhs(mesh->get_num_points(), mesh->get_num_points()),
              rhs(mesh->get_num_points(), 1) {}

        virtual ~Discretisation() = default;

        virtual MatrixCRS<PrecisionType> get_lhs() { return lhs; }
        virtual Matrix<PrecisionType> get_rhs() { return rhs; }

        virtual const Mesh2D& get_mesh() const { return *mesh; };

       protected:
        std::shared_ptr<Mesh2D> mesh;

        MatrixCRS<PrecisionType> lhs;
        Matrix<PrecisionType> rhs;
    };

    class Discretisation2D : public Discretisation {
       public:
        // explicit Discretisation2D(RectMesh2D* m);
        explicit Discretisation2D(std::shared_ptr<RectMesh2D> m);

        ~Discretisation2D();

        MatrixCRS<PrecisionType> get_lhs() override {
            _calculate_lhs();
            return lhs;
        }
        Matrix<PrecisionType> get_rhs() override {
            _calculate_rhs();
            return rhs;
        }

        void add_forcing_term(PrecisionType val);
        void add_forcing_term(
            const std::function<PrecisionType(PrecisionType, PrecisionType)>&
                func);
        void add_forcing_term(const std::string& expr);

        void add_dirichlet_bc(PrecisionType val);
        void add_dirichlet_bc(
            const std::function<PrecisionType(PrecisionType, PrecisionType)>&
                func);
        void add_dirichlet_bc(const std::string& expr);

        void clear_forcing_term();
        void clear_dirichlet_bc();

        const RectMesh2D& get_mesh() const {
            return static_cast<const RectMesh2D&>(*mesh);
        }

       private:
        std::vector<PrecisionType> forcing_terms_val;
        std::vector<std::function<PrecisionType(PrecisionType, PrecisionType)>>
            forcing_terms_func;
        std::vector<std::string> forcing_terms_expr;

        std::vector<PrecisionType> dirichlet_bcs_val;
        std::vector<std::function<PrecisionType(PrecisionType, PrecisionType)>>
            dirichlet_bcs_func;
        std::vector<std::string> dirichlet_bcs_expr;

        void _calculate() {
            // Check if mesh is set
            if (get_mesh().get_Nx() == 0 || get_mesh().get_Ny() == 0) {
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