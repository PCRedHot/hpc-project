#include <stdlib.h>

#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <vector>

#include "matrix.hpp"

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
    Matrix<_T>::Matrix() : MatrixBase<_T>(0, 0) {
#ifdef __DEBUG__
        std::cout << "Null Matrix Created" << std::endl;
#endif
    }

    template <typename _T>
    Matrix<_T>::Matrix(size_t rows, size_t cols) : MatrixBase<_T>(rows, cols) {
        data = std::shared_ptr<_T[]>(new _T[rows * cols]);
        for (size_t i = 0; i < rows * cols; i++) {
            data[i] = _T{};
        }
#ifdef __DEBUG__
        std::cout << rows << "x" << cols << " Matrix Created" << std::endl;
#endif
    }

    template <typename _T>
    Matrix<_T>::Matrix(size_t rows, size_t cols, _T val)
        : MatrixBase<_T>(rows, cols) {
        data = std::shared_ptr<_T[]>(new _T[rows * cols]);
        for (size_t i = 0; i < rows * cols; i++) {
            data[i] = val;
        }
#ifdef __DEBUG__
        std::cout << rows << "x" << cols << " Matrix Created" << std::endl;
#endif
    }

    template <typename _T>
    Matrix<_T>::Matrix(size_t rows, size_t cols, std::vector<_T> vals)
        : MatrixBase<_T>(rows, cols) {
        size_t n_vals = vals.size();

        data = std::shared_ptr<_T[]>(new _T[rows * cols]);
        for (size_t i = 0; i < rows * cols; i++) {
            if (i < n_vals) {
                data[i] = vals[i];
            } else {
                data[i] = _T{};
            }
        }
#ifdef __DEBUG__
        std::cout << rows << "x" << cols << " Matrix Created" << std::endl;
#endif
    }

    template <typename _T>
    Matrix<_T>::Matrix(const Matrix<_T>& o)
        : MatrixBase<_T>(o.get_n_rows(), o.get_n_cols()) {
        data = std::shared_ptr<_T[]>(new _T[o.get_n_rows() * o.get_n_cols()]);
        for (size_t i = 0; i < o.get_n_rows() * o.get_n_cols(); i++) {
            data[i] = o.data[i];
        }
    }

    // MatrixCRS class implementation
    template <typename _T>
    MatrixCRS<_T>::MatrixCRS() : MatrixBase<_T>(0, 0) {
        this->construct(0, 0);
    }

    template <typename _T>
    MatrixCRS<_T>::MatrixCRS(size_t n) : MatrixBase<_T>(n, n) {
        this->construct(n, n);
    }

    template <typename _T>
    MatrixCRS<_T>::MatrixCRS(size_t rows, size_t cols)
        : MatrixBase<_T>(rows, cols) {
        this->construct(rows, cols);
    }

    template <typename _T>
    MatrixCRS<_T>::MatrixCRS(const Matrix<_T>& o)
        : MatrixBase<_T>(o.get_n_rows(), o.get_n_cols()) {
        this->construct(o.get_n_rows(), o.get_n_cols());

        for (size_t i = 0; i < o.get_n_rows(); i++) {
            bool has_non_zero = false;
            for (size_t j = 0; j < o.get_n_cols(); j++) {
                _T val = o.get(i, j);
                if (val != _T{}) {
                    this->values.push_back(val);
                    this->col_indices.push_back(j);

                    if (!has_non_zero) {
                        has_non_zero = true;
                        this->row_ptrs[i] = this->values.size() - 1;
                    }
                }
            }

            if (!has_non_zero) {
                this->row_ptrs[i] = this->values.size();
            }
        }

        row_ptrs[o.get_n_rows()] = values.size();
    }

    template <typename _T>
    MatrixCRS<_T>::MatrixCRS(const MatrixCRS<_T>& o)
        : MatrixBase<_T>(o.get_n_rows(), o.get_n_cols()) {
        this->construct(o.get_n_rows(), o.get_n_cols());

        // Copy the values
        this->values.assign(o.values.begin(), o.values.end());
        this->col_indices.assign(o.col_indices.begin(), o.col_indices.end());

        for (size_t i = 0; i <= this->rows; i++) {
            this->row_ptrs[i] = o.row_ptrs[i];
        }
    }

    template <typename _T>
    void MatrixCRS<_T>::construct(size_t rows, size_t cols) {
        this->row_ptrs = std::shared_ptr<size_t[]>(new size_t[rows + 1]);
        for (size_t i = 0; i <= rows; i++) {
            this->row_ptrs[i] = 0;
        }
    }

    // new get and set methods using Proxy class
    template <typename _T>
    _T MatrixCRS<_T>::_unsave_get(size_t i, size_t j) const {
        size_t index = this->_find(i, j);

        if (index != -1) {
            return values[index];
        }

        return _T{};
    }

    template <typename _T>
    void MatrixCRS<_T>::_unsave_set(size_t i, size_t j, _T val) {
        if (val == _T{}) {
            this->_remove(i, j);
            return;
        }

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
    void MatrixCRS<_T>::_unsave_set(size_t i, size_t j, size_t index, _T val) {
        if (val == _T{}) {
            return this->_remove(i, j, index);
        }

        if (index != -1) {
            values[index] = val;
        } else {
            _unsave_set(i, j, val);
        }
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
    void MatrixCRS<_T>::_remove(size_t i, size_t j, size_t index) {
        if (index == -1) {
            return this->_remove(i, j);
        }

        values.erase(values.begin() + index);
        col_indices.erase(col_indices.begin() + index);

        for (size_t k = i + 1; k <= this->rows; k++) {
            row_ptrs[k]--;
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

        return static_cast<size_t>(-1);
    }

    // Diagonal Matrix class implementation
    template <typename _T>
    MatrixDiagonal<_T>::MatrixDiagonal(size_t n) : MatrixBase<_T>(n, n) {
        data = std::shared_ptr<_T[]>(new _T[n]);
        for (size_t i = 0; i < n; i++) {
            data[i] = _T{};
        }
#ifdef __DEBUG__
        std::cout << n << "x" << n << " Diagonal Matrix Created" << std::endl;
#endif
    }

    template <typename _T>
    MatrixDiagonal<_T>::MatrixDiagonal(size_t n, _T val)
        : MatrixBase<_T>(n, n) {
        data = std::shared_ptr<_T[]>(new _T[n]);
        for (size_t i = 0; i < n; i++) {
            data[i] = val;
        }
#ifdef __DEBUG__
        std::cout << n << "x" << n << " Diagonal Matrix Created" << std::endl;
#endif
    }

    template <typename _T>
    MatrixDiagonal<_T>::MatrixDiagonal(size_t n, std::vector<_T> vals)
        : MatrixBase<_T>(n, n) {
        size_t n_vals = vals.size();

        data = std::shared_ptr<_T[]>(new _T[n]);
        for (size_t i = 0; i < n; i++) {
            if (i < n_vals) {
                data[i] = vals[i];
            } else {
                data[i] = _T{};
            }
        }
#ifdef __DEBUG__
        std::cout << n << "x" << n << " Diagonal Matrix Created" << std::endl;
#endif
    }

    template <typename _T>
    MatrixDiagonal<_T>::MatrixDiagonal(const MatrixBase<_T>& o)
        : MatrixBase<_T>(o.get_n_rows(), o.get_n_cols()) {
        data = std::shared_ptr<_T[]>(new _T[o.get_n_rows()]);
        for (size_t i = 0; i < o.get_n_rows(); i++) {
            data[i] = o.get(i, i);
        }
    }

    template <typename _T>
    MatrixDiagonal<_T>::MatrixDiagonal(const MatrixDiagonal<_T>& o)
        : MatrixBase<_T>(o.get_n_rows(), o.get_n_cols()),
          data(new _T[o.get_n_rows()]) {
        for (size_t i = 0; i < o.get_n_rows(); i++) {
            data[i] = o.data[i];
        }
    }

    template <typename _T>
    MatrixDiagonal<_T>& MatrixDiagonal<_T>::operator=(
        const MatrixDiagonal<_T>& o) {
        if (this != &o) {
            if (this->get_n_rows() != o.get_n_rows() ||
                this->get_n_cols() != o.get_n_cols()) {
                throw std::invalid_argument("Matrix dimensions do not match");
            }
            for (size_t i = 0; i < o.get_n_rows(); i++) {
                data[i] = o.data[i];
            }
        }
        return *this;
    }
}  // namespace fin_diff