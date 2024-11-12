#define CATCH_CONFIG_MAIN
#include "catch2/catch_all.hpp"

#include "matrix.hpp"

using namespace fin_diff;

TEST_CASE("Dense Matrix get and set 1", "[dense_matrix_get_set_1]") {
    Matrix<int> mat(3, 3);

    mat(1, 1) = 42;
    CHECK(mat(1, 1) == 42);
    CHECK(mat(1, 0) == 0);
}

TEST_CASE("Dense Matrix get and set 2", "[dense_matrix_get_set_2]") {
    Matrix<double> mat(2, 6);

    mat(1, 4) = 42.5;
    CHECK(mat(1, 4) == 42.5);
}

TEST_CASE("Dense Matrix Out of Bound 1", "[dense_matrix_out_of_bound_1]") {
    Matrix<int> mat(3, 3);

    CHECK_THROWS_AS(mat(3, 3), std::out_of_range);
    CHECK_THROWS_AS(mat(-1, 3), std::out_of_range);
    CHECK_THROWS_AS(mat(-1, 2), std::out_of_range);
    CHECK_THROWS_AS(mat(0, 7), std::out_of_range);
}

TEST_CASE("Sparse Matrix get and set 1", "[sparse_matrix_get_set_1]") {
    MatrixCRS<int> mat(3);

    mat(1, 1) = 42;
    CHECK(mat(1, 1) == 42);
}

TEST_CASE("Sparse Matrix Out of Bound 1", "[sparse_matrix_out_of_bound_1]") {
    MatrixCRS<int> mat(3);

    CHECK_THROWS_AS(mat(3, 3), std::out_of_range);
    CHECK_THROWS_AS(mat(-1, 3), std::out_of_range);
    CHECK_THROWS_AS(mat(-1, 2), std::out_of_range);
    CHECK_THROWS_AS(mat(0, 7), std::out_of_range);
}

TEST_CASE("Matrix Manipulation 1", "[matrix_man_1]") {
    Matrix<int> mat(3, 3);

    mat(1, 1) = 42;
    mat(1, 2) = 43;
    mat(2, 1) = 44;
    mat(2, 2) = 45;

    mat += 1;

    CHECK(mat(1, 1) == 43);
    CHECK(mat(1, 2) == 44);
    CHECK(mat(2, 1) == 45);
    CHECK(mat(2, 2) == 46);
}

TEST_CASE("Matrix Manipulation 2", "[matrix_man_2]") {
    Matrix<int> mat(3, 3);

    mat(1, 1) = 42;
    mat(1, 2) = 43;
    mat(2, 1) = 44;
    mat(2, 2) = 45;

    mat(1, 1) += 1;
    mat(1, 2) += 20;
    mat(2, 1) += 3;
    mat(2, 2) -= 4;

    CHECK(mat(1, 1) == 43);
    CHECK(mat(1, 2) == 63);
    CHECK(mat(2, 1) == 47);
    CHECK(mat(2, 2) == 41);
}

TEST_CASE("Matrix Manipulation 3", "[matrix_man_3]") {
    MatrixCRS<int> mat(3, 3);

    mat(1, 1) = 42;
    mat(1, 2) = 43;
    mat(2, 1) = 44;
    mat(2, 2) = 45;

    mat(1, 1) += 1;
    mat(1, 2) += 20;
    mat(2, 1) += 3;
    mat(2, 2) -= 4;

    CHECK(mat(1, 1) == 43);
    CHECK(mat(1, 2) == 63);
    CHECK(mat(2, 1) == 47);
    CHECK(mat(2, 2) == 41);
}

TEST_CASE("Matrix Manipulation 4", "[matrix_man_4]") {
    Matrix<int> mat1(3, 3, 10);
    Matrix<int> mat2(3, 3, 3);
    Matrix<int> mat3(3, 4, -2);

    mat1 += mat2;

    CHECK(mat1(1, 1) == 13);
    CHECK(mat1(1, 2) == 13);
    CHECK(mat1(2, 1) == 13);

    CHECK_THROWS_AS(mat1 += mat3, std::invalid_argument);

}


TEST_CASE("Matrix Manipulation 5", "[matrix_man_5]") {
    MatrixCRS<int> mat1 = MatrixCRS<int>::from_dense(Matrix<int>(3, 3, 10));
    MatrixCRS<int> mat2 = MatrixCRS<int>::from_dense(Matrix<int>(3, 3, 3));
    MatrixCRS<int> mat3 = MatrixCRS<int>::from_dense(Matrix<int>(3, 4, -2));

    mat1 += mat2;

    CHECK(mat1(1, 1) == 13);
    CHECK(mat1(1, 2) == 13);
    CHECK(mat1(2, 1) == 13);

    CHECK_THROWS_AS(mat1 += mat3, std::invalid_argument);
}

TEST_CASE("Matrix Manipulation 6", "[matrix_man_6]") {
    MatrixCRS<int> mat1 = MatrixCRS<int>::from_dense(Matrix<int>(3, 3, std::vector<int>{10, 20, -30, 40, 50}));
    MatrixCRS<int> mat2 = MatrixCRS<int>::from_dense(Matrix<int>(3, 3, std::vector<int>{0, 0, 20, 10, 0, 0, 0}));
    MatrixCRS<int> mat3 = MatrixCRS<int>::from_dense(Matrix<int>(3, 4, -2));

    mat1 += mat2;

    CHECK(mat1(0, 0) == 10);
    CHECK(mat1(0, 1) == 20);
    CHECK(mat1(0, 2) == -10);
    CHECK(mat1(1, 0) == 50);
    CHECK(mat1(1, 1) == 50);
    CHECK(mat1(1, 2) == 0);
    CHECK(mat1(2, 0) == 0);
    CHECK(mat1(2, 1) == 0);
    CHECK(mat1(2, 2) == 0);

    CHECK_THROWS_AS(mat1 += mat3, std::invalid_argument);
}