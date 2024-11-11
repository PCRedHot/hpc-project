#include "matrix.hpp"

#include <cstddef>
#include <stdlib.h>
#include <vector>

namespace fin_diff {

template <typename _T>
void MatrixBase<_T>::print() const {
    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            std::cout << get(i, j) << " ";
        }
        std::cout << std::endl;
    }
}

template <typename _T>
Matrix<_T>::Matrix(size_t rows, size_t cols) : rows(rows), cols(cols) {
    data = new _T[rows * cols];
}

template <typename _T>
Matrix<_T>::~Matrix() {
    delete[] data;
}

template <typename _T>
size_t MatrixBase<_T>::get_rows() const {
    return rows;
}

template <typename _T>
size_t MatrixBase<_T>::get_cols() const {
    return cols;
}

template <typename _T>
_T Matrix<_T>::get(size_t i, size_t j) const {
    this->validate(i, j);

    return data[i * cols + j];
}

template <typename _T>
void Matrix<_T>::set(size_t i, size_t j, _T val) {
    this->validate(i, j);

    data[i * cols + j] = val;
}

// MatrixCRS class implementation
template<typename _T>
MatrixCRS<_T>::MatrixCRS(size_t n) : rows(n), cols(n) {
    this->construct(n, n);
}

template <typename _T>
MatrixCRS<_T>::MatrixCRS(size_t rows, size_t cols) : rows(rows), cols(cols) {
    values = std::vector<_T>();
    col_indices = std::vector<size_t>();

    row_ptrs = malloc((rows + 1) * sizeof(size_t));

    for (size_t i = 0; i < rows; i++) {
        row_ptrs[i] = 0;
    }

    row_ptrs[rows] = 1;
}

template <typename _T>
MatrixCRS<_T>::~MatrixCRS() {
    free(row_ptrs);
}

template <typename _T>
_T MatrixCRS<_T>::get(size_t i, size_t j) const{
    this->validate(i, j);
    
    return _T();  // Placeholder
}

template <typename _T>
void MatrixCRS<_T>::set(size_t i, size_t j, _T val) {
    this->validate(i, j);

    // Find the correct position to insert the value
    // Update the row_ptrs array
    size_t cur_row_index = row_ptrs[i];

    if (cur_row_index == row_ptrs[i + 1]) {
        // No elements in the row

    }
}

template <typename _T>
MatrixCRS<_T> MatrixCRS<_T>::from_dense(const Matrix<_T>& mat) {
    // Implement conversion from dense matrix to CRS format
    return MatrixCRS<_T>(mat.get_rows(), mat.get_cols());
}

}  // namespace fin_diff