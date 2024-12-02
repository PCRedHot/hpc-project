
#include "discrete_2d.hpp"

#include <functional>
#include <vector>

#include "exprtk.hpp"
#include "matrix.hpp"
#include "mesh2d.hpp"

namespace fin_diff {

//     Discretisation2D::Discretisation2D(RectMesh2D *m) : Discretisation(m) {
//         _calculate_lhs();
// #ifdef __DEBUG__
//         std::cout << "Discretisation2D Created" << std::endl;
// #endif
//     }

    Discretisation2D::Discretisation2D(std::shared_ptr<RectMesh2D> m) : Discretisation(std::move(m)) {
        _calculate_lhs();
#ifdef __DEBUG__
        std::cout << "Discretisation2D Created" << std::endl;
#endif
    }

    Discretisation2D::~Discretisation2D() {
#ifdef __DEBUG__
        std::cout << "Discretisation2D Deleted" << std::endl;
#endif
    }

    void Discretisation2D::add_forcing_term(double val) {
        this->forcing_terms_val.push_back(val);
        this->rhs_ready = false;
    }

    void Discretisation2D::add_forcing_term(const std::function<double(double, double)> &func) {
        this->forcing_terms_func.push_back(func);
        this->rhs_ready = false;
    }

    void Discretisation2D::add_forcing_term(const std::string &expr) {
        this->forcing_terms_expr.push_back(expr);
        this->rhs_ready = false;
    }

    void Discretisation2D::add_dirichlet_bc(double val) {
        this->dirichlet_bcs_val.push_back(val);
        this->rhs_ready = false;
    }

    void Discretisation2D::add_dirichlet_bc(const std::function<double(double, double)> &func) {
        this->dirichlet_bcs_func.push_back(func);
        this->rhs_ready = false;
    }

    void Discretisation2D::add_dirichlet_bc(const std::string &expr) {
        this->dirichlet_bcs_expr.push_back(expr);
        this->rhs_ready = false;
    }

    void Discretisation2D::clear_forcing_term() {
        this->forcing_terms_val.clear();
        this->forcing_terms_func.clear();
        this->forcing_terms_expr.clear();

        this->rhs_ready = false;
    }

    void Discretisation2D::clear_dirichlet_bc() {
        this->dirichlet_bcs_val.clear();
        this->dirichlet_bcs_func.clear();
        this->dirichlet_bcs_expr.clear();

        this->rhs_ready = false;
    }

    void Discretisation2D::_calculate_lhs(const bool force) {
        if (!force && lhs_ready) return;

        lhs.clear();

        auto lines = get_mesh().get_lines();
        auto coords = get_mesh().get_coordinates();

        double dx_sq = get_mesh().get_dx() * get_mesh().get_dx();
        double dy_sq = get_mesh().get_dy() * get_mesh().get_dy();

        for (int i = 0; i < lines.size(); i++) {
            auto line = lines[i];

            size_t p1 = line.first;
            size_t p2 = line.second;

            bool is_boundary_p1 = get_mesh().is_boundary(p1);
            bool is_boundary_p2 = get_mesh().is_boundary(p2);

            if (is_boundary_p1 && is_boundary_p2) {
                // Both points are on the boundary
                continue;
            }

            double d_sq = 0.0;

            if (p2 - p1 == 1) {  // p2 is always greater than p1
                // Horizontal line
                d_sq = dy_sq;
            } else {
                // Vertical line
                d_sq = dx_sq;
            }

            if (!is_boundary_p1) {
                // p2 contribution on p1
                lhs(p1, p1) += d_sq;

                if (!is_boundary_p2) {
                    lhs(p1, p2) -= d_sq;
                }  // else -> Move to RHS
            }

            if (!is_boundary_p2) {
                // p1 contribution on p2
                lhs(p2, p2) += d_sq;

                if (!is_boundary_p1) {
                    lhs(p2, p1) -= d_sq;
                }  // else -> Move to RHS
            }
        }

        auto boundary = get_mesh().get_boundary();
        for (int i = 0; i < get_mesh().get_num_points(); i++) {
            if (boundary[i]) {
                lhs(i, i) = 1.0;
            }
        }

        lhs_ready = true;
    }

    void Discretisation2D::_calculate_rhs(const bool force) {
        if (!force && rhs_ready) return;

        rhs.clear();

        auto lines = get_mesh().get_lines();
        auto coords = get_mesh().get_coordinates();

        // preprocess string expressions
        typedef exprtk::symbol_table<double> symbol_table_t;
        typedef exprtk::expression<double> expression_t;
        typedef exprtk::parser<double> parser_t;

        double x, y;

        symbol_table_t symbol_table;
        symbol_table.add_variable("x", x);
        symbol_table.add_variable("y", y);
        symbol_table.add_constants();

        std::vector<expression_t> forcing_terms_expr_parsed;
        std::vector<expression_t> dirichlet_bcs_expr_parsed;

        for (auto expr : forcing_terms_expr) {
            expression_t expression;
            expression.register_symbol_table(symbol_table);

            parser_t parser;
            parser.compile(expr, expression);

            forcing_terms_expr_parsed.push_back(expression);
        }

        for (auto expr : dirichlet_bcs_expr) {
            expression_t expression;
            expression.register_symbol_table(symbol_table);

            parser_t parser;
            parser.compile(expr, expression);

            dirichlet_bcs_expr_parsed.push_back(expression);
        }

        // Apply forcing terms and set Dirichlet BCs
        double dx_sq = get_mesh().get_dx() * get_mesh().get_dx();
        double dy_sq = get_mesh().get_dy() * get_mesh().get_dy();

        double f_mul = dx_sq * dy_sq;

        for (int i = 0; i < get_mesh().get_num_points(); i++) {
            bool is_boundary = get_mesh().is_boundary(i);

            x = coords[i].first;
            y = coords[i].second;

            if (is_boundary) {
                // Apply Dirichlet BCs
                for (auto val : dirichlet_bcs_val) {
                    rhs(i, 0) += val;
#ifdef __DEBUG__
                    std::cout << "Point " << i << " (" << x << ", " << y << "): Add val " << val << std::endl;
#endif
                }

                for (auto func : dirichlet_bcs_func) {
                    rhs(i, 0) += func(x, y);
#ifdef __DEBUG__
                    std::cout << "Point " << i << " (" << x << ", " << y << "): Add func val " << func(x, y) << std::endl;
#endif
                }

                for (auto expr : dirichlet_bcs_expr_parsed) {
                    rhs(i, 0) += expr.value();
#ifdef __DEBUG__
                    std::cout << "Point " << i << " (" << x << ", " << y << "): Add expr val " << expr.value() << std::endl;
#endif
                }
            } else {
                // Apply forcing terms
                for (auto val : forcing_terms_val) {
                    rhs(i, 0) += val * f_mul;
                }

                for (auto func : forcing_terms_func) {
                    rhs(i, 0) += func(x, y) * f_mul;
                }

                for (auto expr : forcing_terms_expr_parsed) {
                    rhs(i, 0) += expr.value() * f_mul;
                }
            }
        }

        // Add Elimated Boundary Point Influences
        for (auto line : lines) {
            size_t p1 = line.first;
            size_t p2 = line.second;

            bool is_boundary_p1 = get_mesh().is_boundary(p1);
            bool is_boundary_p2 = get_mesh().is_boundary(p2);

            if (!(is_boundary_p1 ^ is_boundary_p2)) {
                // Both points are on the boundary or both are not on the boundary
                continue;
            }

            double d_sq = 0.0;
            if (coords[p1].first == coords[p2].first) {
                // Horizontal line
                d_sq = dy_sq;
            } else {
                // Vertical line
                d_sq = dx_sq;
            }

            if (!is_boundary_p1) {
                // p2 bc contribution on p1
                rhs(p1, 0) -= d_sq * rhs(p2, 0);
            } else {
                // p1 bc contribution on p2
                rhs(p2, 0) -= d_sq * rhs(p1, 0);
            }
        }

        rhs_ready = true;
    }
}  // namespace fin_diff
