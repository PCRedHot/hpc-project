#define CATCH_CONFIG_MAIN
#include "catch2/catch_all.hpp"

#include "matrix.hpp" // Include your matrix header

using namespace fin_diff;

TEST_CASE("Dense Matrix get and set 1", "[dense_matrix_get_set_1]") {
    Matrix<int> mat(3, 3);

    mat.set(1, 1, 42);
    REQUIRE(mat.get(1, 1) == 42);
}

TEST_CASE("Dense Matrix get and set 2", "[dense_matrix_get_set_2]") {
    Matrix<double> mat(2, 6);

    mat.set(1, 4, 42.5);
    REQUIRE(mat.get(1, 4) == 42.5);
}

TEST_CASE("Dense Matrix Out of Bound 1", "[dense_matrix_out_of_bound_1]") {
    Matrix<int> mat(3, 3);

    REQUIRE_THROWS_AS(mat.get(3, 3), std::out_of_range);
    REQUIRE_THROWS_AS(mat.get(-1, 3), std::out_of_range);
    REQUIRE_THROWS_AS(mat.get(-1, 2), std::out_of_range);
    REQUIRE_THROWS_AS(mat.get(0, 7), std::out_of_range);
}