#define CATCH_CONFIG_MAIN
#include "catch2/catch_all.hpp"

#include "param.hpp"

using namespace fin_diff;

TEST_CASE("Read Parameter File 1", "[read_param_file_1]") {
    Parameter param("./data/param1.param");

    CHECK(param.get_dimension() == 2);
    CHECK(param.get_num_x() == 5);
    CHECK(param.get_num_y() == 5);
    CHECK(param.get_num_z() == 0);
    CHECK(param.get_solver() == "jacobi");
    CHECK(param.get_init_guess_expr() == "0.0");
    CHECK(param.get_tolerance() == 1e-6);
    CHECK(param.get_max_iter() == 100);
    CHECK(param.get_output_file() == "./output");
}

TEST_CASE("Read Parameter File 2", "[read_param_file_2]") {
    CHECK_THROWS_AS(Parameter("./data/param2.param"), std::invalid_argument);
}

TEST_CASE("Read Parameter File 3", "[read_param_file_3]") {
    CHECK_THROWS_AS(Parameter("./data/param_not_exist.param"), std::invalid_argument);
}
