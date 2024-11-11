#include "matrix.hpp"

#include <cstddef>
#include <stdlib.h>
#include <vector>
#include <iostream>

namespace fin_diff {

template <typename _T>
void Matrix<_T>::print() const {
    for (size_t i = 0; i < this->rows; ++i) {
        for (size_t j = 0; j < this->cols; ++j) {
            std::cout << get(i, j) << " ";
        }
        std::cout << std::endl;
    }
}

template <typename _T>
Matrix<_T>::Matrix(size_t rows, size_t cols) : MatrixBase<_T>(rows, cols) {
    std::cout << rows << "x" << cols << " Matrix Created" << std::endl;
    
    data = static_cast<_T*>(std::calloc(rows * cols, sizeof(_T)));
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

    return data[i * this->cols + j];
}

template <typename _T>
void Matrix<_T>::set(size_t i, size_t j, _T val) {
    this->validate(i, j);

    data[i * this->cols + j] = val;
}

// MatrixCRS class implementation
template <typename _T>
MatrixCRS<_T>::MatrixCRS(size_t n) : MatrixBase<_T>(n, n) {
    this->construct(n, n);
}

template <typename _T>
MatrixCRS<_T>::MatrixCRS(size_t rows, size_t cols) : MatrixBase<_T>(rows, cols) {
    this->construct(rows, cols);
}

template <typename _T>
void MatrixCRS<_T>::construct(size_t rows, size_t cols) {
    this->row_ptrs = new size_t[rows + 1];
    for (size_t i = 0; i <= rows; ++i) {
        row_ptrs[i] = 0;
    }
}


template <typename _T>
MatrixCRS<_T>::~MatrixCRS() {
    free(row_ptrs);
}

template <typename _T>
_T MatrixCRS<_T>::get(size_t i, size_t j) const{
    // this->validate(i, j); // This will be handled by the _find method

    size_t index = this->_find(i, j);

    if (index != -1) {
        return values[index];
    }

    return _T{};
}

template <typename _T>
void MatrixCRS<_T>::set(size_t i, size_t j, _T val) {
    if (val == _T{}) {
        this->_remove(i, j);
        return;
    }

    this->validate(i, j);

    // Find the correct position to insert the value
    size_t cur_ptr = row_ptrs[i];
    size_t nxt_row_start_index = row_ptrs[i + 1];

    while (cur_ptr < nxt_row_start_index) {
        if (col_indices[cur_ptr] == j) {
            values[cur_ptr] = val;
            return;
        } else if (col_indices[cur_ptr] > j) {
            break;
        }
        cur_ptr++;
    }

    // Insert the value at the correct position
    values.insert(values.begin() + cur_ptr, val);
    col_indices.insert(col_indices.begin() + cur_ptr, j);

    // Update the row pointers
    for (size_t k = i + 1; k <= this->rows; k++) {
        row_ptrs[k]++;
    }
}

template <typename _T>
MatrixCRS<_T> MatrixCRS<_T>::from_dense(const Matrix<_T>& mat) {
    MatrixCRS<_T> crs(mat.get_rows(), mat.get_cols());

    for (size_t i = 0; i < mat.get_rows(); i++) {
        bool has_non_zero = false;
        for (size_t j = 0; j < mat.get_cols(); j++) {
            _T val = mat.get(i, j);
            if (val != _T{}) {
                crs.values.push_back(val);
                crs.col_indices.push_back(j);
                
                if (!has_non_zero) {
                    has_non_zero = true;
                    crs.row_ptrs[i] = crs.values.size() - 1;
                }
            }
        }

        if (!has_non_zero) {
            crs.row_ptrs[i] = crs.values.size();
        }
    }

    crs.row_ptrs[crs.get_rows()] = crs.values.size();

    return crs;
}

template <typename _T>
void MatrixCRS<_T>::print() const {
    // TODO: Better way to print the matrix
    for (size_t i = 0; i < this->rows; ++i) {
        for (size_t j = 0; j < this->cols; ++j) {
            std::cout << get(i, j) << " ";
        }
        std::cout << std::endl;
    }
}

template <typename _T>
void MatrixCRS<_T>::_remove(size_t i, size_t j) {
    size_t index = this->_find(i, j);

    if (index != -1) {
        values.erase(values.begin() + index);
        col_indices.erase(col_indices.begin() + index);

        for (size_t k = i + 1; k <= this->rows; k++) {
            row_ptrs[k]--;
        }
    }
}

template <typename _T>
size_t MatrixCRS<_T>::_find(size_t i, size_t j) const {
    this->validate(i, j);

    size_t cur_ptr = row_ptrs[i];
    size_t nxt_row_start_index = row_ptrs[i + 1];

    while (cur_ptr < nxt_row_start_index) {
        if (col_indices[cur_ptr] == j) {
            return cur_ptr;
        } else if (col_indices[cur_ptr] > j) {
            break;
        }
        cur_ptr++;
    }

    return -1;
}

}  // namespace fin_diff