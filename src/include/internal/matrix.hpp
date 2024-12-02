#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace fin_diff {

    // Forward declaration
    template <typename _T>
    class Matrix;

    template <typename _T>
    class MatrixCRS;

    template <typename _T>
    class MatrixDiagonal;

    template <typename _T>
    class MatrixBase {
       public:
        MatrixBase(size_t rows, size_t cols) : rows(rows), cols(cols) {}
        // virtual ~MatrixBase() = default;

        size_t get_n_rows() const { return rows; }
        size_t get_n_cols() const { return cols; }

        virtual _T get(size_t i, size_t j) const = 0;
        virtual void set(size_t i, size_t j, _T val) = 0;

        template <typename _T2>
        MatrixBase<_T>& operator+=(MatrixBase<_T2> o) {
            if (this->get_n_rows() != o.get_n_rows() ||
                this->get_n_cols() != o.get_n_cols()) {
                throw std::invalid_argument("Matrix dimensions do not match");
            }

            for (int i = 0; i < this->get_n_cols() * this->get_n_rows(); i++) {
                this->set(
                    i / this->get_n_cols(), i % this->get_n_cols(),
                    this->get(i / this->get_n_cols(), i % this->get_n_cols()) +
                        o.get(i / this->get_n_cols(), i % this->get_n_cols()));
            }

            return *this;
        }

        template <typename _T2>
        Matrix<_T> operator+(const MatrixBase<_T2>& o) {
            if (this->get_n_rows() != o.get_n_rows() ||
                this->get_n_cols() != o.get_n_cols()) {
                throw std::invalid_argument("Matrix dimensions do not match");
            }

            Matrix<_T> res(this->get_n_rows(), this->get_n_cols());

            for (size_t i = 0; i < this->get_n_rows(); i++) {
                for (size_t j = 0; j < this->get_n_cols(); j++) {
                    res(i, j) = this->get(i, j) + o.get(i, j);
                }
            }

            return res;
        }

        template <typename _T2>
        Matrix<_T> operator-(const MatrixBase<_T2>& o) {
            if (this->get_n_rows() != o.get_n_rows() ||
                this->get_n_cols() != o.get_n_cols()) {
                throw std::invalid_argument("Matrix dimensions do not match");
            }

            Matrix<_T> res(this->get_n_rows(), this->get_n_cols());

            for (size_t i = 0; i < this->get_n_rows(); i++) {
                for (size_t j = 0; j < this->get_n_cols(); j++) {
                    res(i, j) = this->get(i, j) - o.get(i, j);
                }
            }

            return res;
        }

        template <typename _T2>
        Matrix<_T> operator*(const MatrixBase<_T2>& o) {
            this->validate_matrix_mul(o);

            size_t n_rows = this->get_n_rows();
            size_t n_cols = o.get_n_cols();
            size_t n_inner = this->get_n_cols();

            Matrix<_T> res(n_rows, n_cols);

            for (size_t i = 0; i < n_rows; i++) {
                for (size_t j = 0; j < n_cols; j++) {
                    _T sum = 0;
                    for (size_t k = 0; k < n_inner; k++) {
                        sum += this->get(i, k) * o.get(k, j);
                    }
                    res(i, j) = sum;
                }
            }

            return res;
        }

        void print() const {
            for (size_t i = 0; i < this->rows; ++i) {
                for (size_t j = 0; j < this->cols; ++j) {
                    std::cout << get(i, j) << " ";
                }
                std::cout << std::endl;
            }
        };

       protected:
        size_t rows, cols;

        void validate(size_t i, size_t j) const {
            if (i < 0 || j < 0 || i >= get_n_rows() || j >= get_n_cols()) {
                std::ostringstream oss;
                oss << "Index out of range: i=" << i << ", j=" << j
                    << ", rows=" << get_n_rows() << ", cols=" << get_n_cols();
                throw std::out_of_range(oss.str());
            }
        }

        void validate_matrix_mul(const MatrixBase<_T>& o) const {
            if (this->get_n_cols() != o.get_n_rows()) {
                std::cerr << "Matrix dimensions do not match: "
                          << this->get_n_rows() << "x" << this->get_n_cols()
                          << " times " << o.get_n_rows() << "x"
                          << o.get_n_cols() << std::endl;
                throw std::invalid_argument("Matrix dimensions do not match");
            }
        }
    };

    template <typename _T>
    class Matrix : public MatrixBase<_T> {
       public:
        Matrix();
        Matrix(size_t rows, size_t cols);
        Matrix(size_t rows, size_t cols, _T val);
        Matrix(size_t rows, size_t cols, std::vector<_T> vals);

        Matrix(const Matrix<_T>& o);

        _T get(size_t i, size_t j) const override { return (*this)(i, j); }
        void set(size_t i, size_t j, _T val) override { (*this)(i, j) = val; }

        // Overload the operator() for two-dimensional access
        _T& operator()(size_t row, size_t col) {
            this->validate(row, col);
            return data[row * this->cols + col];
        }
        const _T& operator()(size_t row, size_t col) const {
            this->validate(row, col);
            return data[row * this->cols + col];
        }

        Matrix<_T>& operator+=(_T val) {
            for (int i = 0; i < this->get_n_cols() * this->get_n_rows(); i++) {
                data[i] += val;
            }

            return *this;
        }

        Matrix<_T>& operator-=(_T val) {
            for (int i = 0; i < this->get_n_cols() * this->get_n_rows(); i++) {
                data[i] -= val;
            }

            return *this;
        }

        Matrix<_T>& operator*=(_T val) {
            for (int i = 0; i < this->get_n_cols() * this->get_n_rows(); i++) {
                data[i] *= val;
            }

            return *this;
        }

        Matrix<_T>& operator/=(_T val) {
            for (int i = 0; i < this->get_n_cols() * this->get_n_rows(); i++) {
                data[i] /= val;
            }

            return *this;
        }

        template <typename _T2>
        Matrix<_T>& operator+=(const Matrix<_T2>& o) {
            if (this->get_n_rows() != o.get_n_rows() ||
                this->get_n_cols() != o.get_n_cols()) {
                throw std::invalid_argument("Matrix dimensions do not match");
            }

            for (int i = 0; i < this->get_n_cols() * this->get_n_rows(); i++) {
                data[i] += o.data[i];
            }

            return *this;
        }

        template <typename _T2>
        Matrix<_T> operator*(const MatrixBase<_T2>& o) {
            this->validate_matrix_mul(o);

            size_t n_rows = this->get_n_rows();
            size_t n_cols = o.get_n_cols();
            size_t n_inner = this->get_n_cols();

            Matrix<_T> res(n_rows, n_cols);

            for (size_t i = 0; i < n_rows; i++) {
                for (size_t j = 0; j < n_cols; j++) {
                    _T sum = 0;
                    for (size_t k = 0; k < n_inner; k++) {
                        sum += this->get(i, k) * o.get(k, j);
                    }
                    res(i, j) = sum;
                }
            }

            return res;
        }

        double norm() const {
            if (this->get_n_rows() != 1 && this->get_n_cols() != 1) {
                throw std::invalid_argument(
                    "Norm can only be calculated for vectors (1D matrices)");
            }

            double sum = 0;
            for (int i = 0; i < this->get_n_cols() * this->get_n_rows(); i++) {
                sum += data[i] * data[i];
            }

            return std::sqrt(sum);
        }

        double weighted_norm() const {
            if (this->get_n_rows() != 1 && this->get_n_cols() != 1) {
                throw std::invalid_argument(
                    "Norm can only be calculated for vectors (1D matrices)");
            }

            double sum = 0;
            for (int i = 0; i < this->get_n_cols() * this->get_n_rows(); i++) {
                sum += data[i] * data[i];
            }

            sum /= this->get_n_cols();

            return std::sqrt(sum);
        }

        std::vector<_T> to_vector() const {
            if (this->get_n_cols() != 1 && this->get_n_rows() != 1) {
                throw std::invalid_argument(
                    "Matrix cannot be converted to a vector");
            }

            std::vector<_T> vec(this->get_n_cols() * this->get_n_rows());
            for (int i = 0; i < this->get_n_cols() * this->get_n_rows(); i++) {
                vec[i] = data[i];
            }

            return vec;
        }

        void clear() {
            for (int i = 0; i < this->get_n_cols() * this->get_n_rows(); i++) {
                data[i] = _T{};
            }
        };

        void print() const;

       private:
        std::shared_ptr<_T[]> data;
    };

    template <typename _T>
    class MatrixCRS : public MatrixBase<_T> {
       public:
        MatrixCRS();
        MatrixCRS(size_t n);
        MatrixCRS(size_t rows, size_t cols);
        MatrixCRS(const Matrix<_T>& o);
        MatrixCRS(const MatrixCRS<_T>& o);

        void construct(size_t rows, size_t cols);

        // static MatrixCRS<_T> from_dense(const Matrix<_T>& mat);

        _T get(size_t i, size_t j) const override { return (*this)(i, j); }
        void set(size_t i, size_t j, _T val) override { (*this)(i, j) = val; }

        // Proxy class for handling element access and modification
        class Proxy {
           public:
            Proxy(MatrixCRS<_T>& matrix, size_t row, size_t col)
                : matrix(matrix), row(row), col(col) {}

            operator _T() const { return matrix._unsave_get(row, col); }

            Proxy& operator=(_T val) {
                matrix._unsave_set(row, col, val);
                return *this;
            }

            Proxy& operator+=(_T val) {
                size_t index = matrix._find(row, col);
                _T entry_val = index != -1 ? matrix.values[index] : 0.0;

                matrix._unsave_set(row, col, index, entry_val + val);
                return *this;
            }

            Proxy& operator-=(_T val) {
                size_t index = matrix._find(row, col);
                _T entry_val = index != -1 ? matrix.values[index] : 0.0;

                matrix._unsave_set(row, col, index, entry_val - val);
                return *this;
            }

           private:
            MatrixCRS<_T>& matrix;
            size_t row, col;
        };

        Proxy operator()(size_t row, size_t col) {
            this->validate(row, col);

            return Proxy(*this, row, col);
        }

        const _T operator()(size_t row, size_t col) const {
            this->validate(row, col);

            return _unsave_get(row, col);
        }

        template <typename _T2>
        MatrixCRS<_T>& operator+=(MatrixCRS<_T2> o) {
            if (this->get_n_rows() != o.get_n_rows() ||
                this->get_n_cols() != o.get_n_cols()) {
                throw std::invalid_argument("Matrix dimensions do not match");
            }

            size_t r = 0;
            for (int i = 0; i < o.row_ptrs[this->get_n_rows()]; i++) {
                if (i == o.row_ptrs[r + 1]) {
                    r++;
                }

                size_t c = o.col_indices[i];

                size_t index = this->_find(r, c);

                _T o_val = static_cast<_T>(o.values[i]);
                _T val = index != -1 ? this->values[index] : 0.0;
                this->_unsave_set(r, c, index, val + o_val);
            }
            return *this;
        }

        template <typename _T2>
        MatrixCRS<_T>& operator+=(Matrix<_T2> o) {
            if (this->get_n_rows() != o.get_n_rows() ||
                this->get_n_cols() != o.get_n_cols()) {
                throw std::invalid_argument("Matrix dimensions do not match");
            }

            for (size_t i = 0; i < o.get_n_rows(); i++) {
                for (size_t j = 0; j < o.get_n_cols(); j++) {
                    _T val = o.get(i, j);
                    if (val != _T{}) {
                        size_t index = this->_find(i, j);
                        _T this_val = index != -1 ? this->values[index] : 0.0;
                        this->_unsave_set(i, j, index, this_val + val);
                    }
                }
            }

            return *this;
        }

        template <typename _T2>
        MatrixCRS<_T>& operator+=(const MatrixDiagonal<_T2> o) {
            if (this->get_n_rows() != o.get_n_rows() ||
                this->get_n_cols() != o.get_n_cols()) {
                throw std::invalid_argument("Matrix dimensions do not match");
            }

            for (size_t i = 0; i < o.get_n_rows(); i++) {
                size_t index = this->_find(i, i);
                _T val = index != -1 ? this->values[index] : 0.0;
                this->_unsave_set(i, i, index, val + o(i, i));
            }

            return *this;
        }

        template <typename _T2>
        MatrixCRS<_T> operator-(const MatrixDiagonal<_T2> o) {
            if (this->get_n_rows() != o.get_n_rows() ||
                this->get_n_cols() != o.get_n_cols()) {
                throw std::invalid_argument("Matrix dimensions do not match");
            }

            MatrixCRS<_T> res(*this);

            for (size_t i = 0; i < this->get_n_rows(); i++) {
                size_t index = res._find(i, i);
                _T val = index != -1 ? res.values[index] : 0.0;
                res._unsave_set(i, i, index, val - o(i, i));
            }

            return res;
        }

        template <typename _T2>
        Matrix<_T> operator-(const MatrixBase<_T2>& o) {
            if (this->get_n_rows() != o.get_n_rows() ||
                this->get_n_cols() != o.get_n_cols()) {
                throw std::invalid_argument("Matrix dimensions do not match");
            }

            Matrix<_T> res(this->get_n_rows(), this->get_n_cols());

            for (size_t i = 0; i < this->get_n_rows(); i++) {
                for (size_t j = 0; j < this->get_n_cols(); j++) {
                    res(i, j) = this->get(i, j) - o.get(i, j);
                }
            }

            return res;
        }

        template <typename _T2>
        MatrixCRS<_T> operator*(const MatrixBase<_T2>& o) {
            this->validate_matrix_mul(o);

            size_t n_rows = this->get_n_rows();
            size_t n_cols = o.get_n_cols();
            size_t n_inner = this->get_n_cols();

            MatrixCRS<_T> res(n_rows, n_cols);

            for (size_t i = 0; i < n_rows; i++) {
                for (size_t j = 0; j < n_cols; j++) {
                    _T sum = 0;
                    for (size_t k = 0; k < n_inner; k++) {
                        sum += this->get(i, k) * o.get(k, j);
                    }
                    res(i, j) = sum;
                }
            }

            return res;
        }

        void clear() {
            values.clear();
            col_indices.clear();
            for (size_t i = 0; i <= this->rows; i++) {
                row_ptrs[i] = 0;
            }
        }

        void print() const;

       private:
        std::vector<_T> values;
        std::vector<size_t> col_indices;
        std::shared_ptr<size_t[]> row_ptrs;

        void _remove(size_t i, size_t j);
        void _remove(size_t i, size_t j, size_t index);
        size_t _find(size_t i, size_t j) const;

        _T _unsave_get(size_t i, size_t j) const;
        void _unsave_set(size_t i, size_t j, _T val);
        void _unsave_set(size_t i, size_t j, size_t index, _T val);
    };

    template <typename _T>
    class MatrixDiagonal : public MatrixBase<_T> {
       public:
        MatrixDiagonal(size_t n);
        MatrixDiagonal(size_t n, _T val);
        MatrixDiagonal(size_t n, std::vector<_T> vals);
        MatrixDiagonal(const MatrixBase<_T>& o);
        MatrixDiagonal(const MatrixDiagonal<_T>& o);
        MatrixDiagonal<_T>& operator=(const MatrixDiagonal<_T>& o);

        _T get(size_t i, size_t j) const override {
            this->validate(i, j);

            if (i != j) {
                return _T{};
            }

            return data[i];
        }

        void set(size_t i, size_t j, _T val) override {
            this->diagonal_validate(i, j);

            data[i] = val;
        }

        // Overload the operator() for two-dimensional access
        _T& operator()(size_t i, size_t j) {
            this->diagonal_validate(i, j);

            return data[i];
        }

        const _T operator()(size_t i, size_t j) const {
            this->validate(i, j);

            if (i != j) {
                return _T{};
            }

            return data[i];
        }

        MatrixDiagonal<_T>& operator+=(_T val) {
            for (int i = 0; i < this->get_n_cols(); i++) {
                data[i] += val;
            }

            return *this;
        }

        MatrixDiagonal<_T>& operator-=(_T val) {
            for (int i = 0; i < this->get_n_cols(); i++) {
                data[i] -= val;
            }

            return *this;
        }

        MatrixDiagonal<_T>& operator*=(_T val) {
            for (int i = 0; i < this->get_n_cols(); i++) {
                data[i] *= val;
            }

            return *this;
        }

        MatrixDiagonal<_T>& operator/=(_T val) {
            for (int i = 0; i < this->get_n_cols(); i++) {
                data[i] /= val;
            }

            return *this;
        }

        template <typename _T2>
        MatrixDiagonal<_T>& operator+=(MatrixDiagonal<_T2> o) {
            if (this->get_n_rows() != o.get_n_rows()) {
                throw std::invalid_argument("Matrix dimensions do not match");
            }

            for (int i = 0; i < this->get_n_cols(); i++) {
                data[i] += o.data[i];
            }

            return *this;
        }

        template <typename _T2>
        Matrix<_T> operator*(const MatrixBase<_T2>& o) {
            this->validate_matrix_mul(o);

            size_t n_rows = this->get_n_rows();
            size_t n_cols = o.get_n_cols();

            Matrix<_T> res(n_rows, n_cols);

            for (size_t i = 0; i < n_rows; i++) {
                for (size_t j = 0; j < n_cols; j++) {
                    res(i, j) = this->get(i, i) * o.get(i, j);
                }
            }

            return res;
        }

        void inv() {
            for (int i = 0; i < this->get_n_cols(); i++) {
                data[i] = 1 / data[i];
            }
        }

        void clear() {
            for (int i = 0; i < this->get_n_cols(); i++) {
                data[i] = _T{};
            }
        }

       private:
        std::shared_ptr<_T[]> data;

        void diagonal_validate(size_t i, size_t j) const {
            this->validate(i, j);
            if (i != j) {
                throw std::invalid_argument("Matrix is a diagonal matrix");
            }
        }
    };
}  // namespace fin_diff

#include "matrix.tpp"

#endif