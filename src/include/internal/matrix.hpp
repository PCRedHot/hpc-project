#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <cstddef>
#include <iostream>
#include <vector>

namespace fin_diff {

template <typename _T>
class MatrixBase {
   public:
    virtual ~MatrixBase() = default;

    size_t get_rows() const;
    size_t get_cols() const;

    virtual _T get(size_t i, size_t j) const = 0;
    virtual void set(size_t i, size_t j, _T val) = 0;
    void print() const = 0;

   protected:
    size_t rows, cols;

    void validate(size_t i, size_t j) const {
        if (i >= get_rows() || j >= get_cols()) {
            throw std::out_of_range("Index out of range");
        }
    }
};

template <typename _T>
class Matrix : public MatrixBase<_T> {
   public:
    Matrix(size_t rows, size_t cols);
    ~Matrix();

    _T get(size_t i, size_t j) const override;
    void set(size_t i, size_t j, _T val) override;

   private:
    size_t rows, cols;
    _T *data;
};

template <typename _T>
class MatrixCRS : public MatrixBase<_T> {
   public:
    MatrixCRS(size_t n);
    MatrixCRS(size_t rows, size_t cols);
    ~MatrixCRS();

    _T get(size_t i, size_t j) const override;
    void set(size_t i, size_t j, _T val) override;

    static MatrixCRS<_T> from_dense(const Matrix<_T> &mat);

   private:
    size_t rows, cols;
    std::vector<_T> values;
    std::vector<size_t> col_indices;
    size_t *row_ptrs;
};

}  // namespace fin_diff

#endif  // MATRIX_HPP