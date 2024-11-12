#define CATCH_CONFIG_MAIN
#include "catch2/catch_all.hpp"

#include "matrix.hpp" // Include your matrix header

using namespace fin_diff;

TEST_CASE("Dense Matrix get and set 1", "[dense_matrix_get_set_1_legacy]") {
    Matrix<int> mat(3, 3);

    mat.set(1, 1, 42);
    REQUIRE(mat.get(1, 1) == 42);
    REQUIRE(mat.get(1, 0) == 0);
}

TEST_CASE("Dense Matrix get and set 2", "[dense_matrix_get_set_2_legacy]") {
    Matrix<double> mat(2, 6);

    mat.set(1, 4, 42.5);
    REQUIRE(mat.get(1, 4) == 42.5);
}

TEST_CASE("Dense Matrix Out of Bound 1", "[dense_matrix_out_of_bound_1_legacy]") {
    Matrix<int> mat(3, 3);

    REQUIRE_THROWS_AS(mat.get(3, 3), std::out_of_range);
    REQUIRE_THROWS_AS(mat.get(-1, 3), std::out_of_range);
    REQUIRE_THROWS_AS(mat.get(-1, 2), std::out_of_range);
    REQUIRE_THROWS_AS(mat.get(0, 7), std::out_of_range);
}

TEST_CASE("Sparse Matrix get and set 1", "[sparse_matrix_get_set_1_legacy]") {
    MatrixCRS<int> mat(3);

    mat.set(1, 1, 42);
    mat.set(1, 2, 43);
    REQUIRE(mat.get(1, 1) == 42);
    REQUIRE(mat.get(1, 2) == 43);
    REQUIRE(mat.get(1, 0) == 0);
}

TEST_CASE("Sparse Matrix get and set 2", "[sparse_matrix_get_set_2_legacy]") {
    Matrix<int> mat(3, 3);
    mat.set(0, 0, 1);
    mat.set(0, 2, 2);
    mat.set(2, 0, 3);
    mat.set(2, 1, 4);

    MatrixCRS<int> crs = MatrixCRS<int>(mat);

    crs.set(1, 1, 42);
    crs.set(2, 1, 0);
    crs.set(0, 0, 0);
    crs.set(0, 2, 24);

    REQUIRE(crs.get(0, 0) == 0);
    REQUIRE(crs.get(0, 2) == 24);
    REQUIRE(crs.get(2, 0) == 3);
    REQUIRE(crs.get(2, 1) == 0);
    REQUIRE(crs.get(1, 1) == 42);
}

TEST_CASE("Sparsen Matrix from Dense Matrix 1", "[spase_from_dense_matrix_1_legacy]") {
    Matrix<int> mat(3, 3);
    mat.set(0, 0, 1);
    mat.set(0, 2, 2);
    mat.set(2, 0, 3);
    mat.set(2, 1, 4);

    MatrixCRS<int> crs = MatrixCRS<int>(mat);

    REQUIRE(crs.get(0, 0) == 1);
    REQUIRE(crs.get(0, 2) == 2);
    REQUIRE(crs.get(2, 0) == 3);
    REQUIRE(crs.get(2, 1) == 4);
    REQUIRE(crs.get(1, 1) == 0);
}

