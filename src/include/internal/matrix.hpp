#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <vector>

namespace fin_diff {

    template <typename _T>
    class MatrixBase {
       public:
        MatrixBase(size_t rows, size_t cols) : rows(rows), cols(cols) {}
        // virtual ~MatrixBase() = default;

        size_t get_n_rows() const;
        size_t get_n_cols() const;

        
        virtual _T get(size_t i, size_t j) const = 0;
        virtual void set(size_t i, size_t j, _T val) = 0;



        template <typename _T2>
        MatrixBase<_T>& operator+=(MatrixBase<_T2> o) {
            if (this->get_n_rows() != o.get_n_rows() || this->get_n_cols() != o.get_n_cols()) {
                throw std::invalid_argument("Matrix dimensions do not match");
            }

            for (int i = 0; i < this->get_n_cols() * this->get_n_rows(); i++) {
                this->set(i / this->get_n_cols(), i % this->get_n_cols(), this->get(i / this->get_n_cols(), i % this->get_n_cols()) + o.get(i / this->get_n_cols(), i % this->get_n_cols()));
            }

            return *this;
        }

        virtual void print() const = 0;

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
    };

    template <typename _T>
    class Matrix : public MatrixBase<_T> {
       public:
        Matrix();
        Matrix(size_t rows, size_t cols);
        Matrix(size_t rows, size_t cols, _T val);
        Matrix(size_t rows, size_t cols, std::vector<_T> vals);
        // ~Matrix();

        
        _T get(size_t i, size_t j) const override;
        void set(size_t i, size_t j, _T val) override;

        // Overload the operator() for two-dimensional access
        _T& operator()(size_t row, size_t col);
        const _T& operator()(size_t row, size_t col) const;

        Matrix<_T>& operator+=(_T val) {
            for (int i = 0; i < this->get_n_cols() * this->get_n_rows(); i++) {
                data[i] += val;
            }

            return *this;
        }

        template <typename _T2>
        Matrix<_T>& operator+=(Matrix<_T2> o) {
            if (this->get_n_rows() != o.get_n_rows() || this->get_n_cols() != o.get_n_cols()) {
                throw std::invalid_argument("Matrix dimensions do not match");
            }

            for (int i = 0; i < this->get_n_cols() * this->get_n_rows(); i++) {
                data[i] += o.data[i];
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

        void clear() {
            for (int i = 0; i < this->get_n_cols() * this->get_n_rows(); i++) {
                data[i] = _T{};
            }
        };

        void print() const override;

       private:
        std::shared_ptr<_T[]> data;
    };

    template <typename _T>
    class MatrixCRS : public MatrixBase<_T> {
       public:
        MatrixCRS();
        MatrixCRS(size_t n);
        MatrixCRS(size_t rows, size_t cols);
        // ~MatrixCRS();

        void construct(size_t rows, size_t cols);

        _T get(size_t i, size_t j) const override;
        void set(size_t i, size_t j, _T val) override;

        // Proxy class for handling element access and modification
        class Proxy {
           public:
            Proxy(MatrixCRS<_T>& matrix, size_t row, size_t col)
                : matrix(matrix), row(row), col(col) {}

            operator _T() const {
                return matrix._unsave_get(row, col);
            }

            Proxy& operator=(_T val) {
                matrix._unsave_set(row, col, val);
                return *this;
            }

            Proxy& operator+=(_T val) {
                size_t index = matrix._find(row, col);
                double entry_val = index != -1 ? matrix.values[index] : 0.0;

                matrix._unsave_set(row, col, index, entry_val + val);
                return *this;
            }

            Proxy& operator-=(_T val) {
                size_t index = matrix._find(row, col);
                double entry_val = index != -1 ? matrix.values[index] : 0.0;

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
            if (this->get_n_rows() != o.get_n_rows() || this->get_n_cols() != o.get_n_cols()) {
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
                this->_unsave_set(r, c, index, val+o_val);
            }
            return *this;
        }

        void clear() {
            values.clear();
            col_indices.clear();
            for (size_t i = 0; i <= this->rows; i++) {
                row_ptrs[i] = 0;
            }
        }

        static MatrixCRS<_T> from_dense(const Matrix<_T>& mat);

        void print() const override;

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


        _T get(size_t i, size_t j) const override;
        void set(size_t i, size_t j, _T val) override;

        void print() const override;

       private:
        std::shared_ptr<_T[]>  data;

        void diagonal_validate(size_t i, size_t j) const;
    };
}  // namespace fin_diff

#include "matrix.tpp"

#endif