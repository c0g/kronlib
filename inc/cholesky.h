//
// Created by Thomas Nickson on 9/07/15.
//

#ifndef KRONMAT_CHOLESKY_H
#define KRONMAT_CHOLESKY_H

#include <vector>
#include "blas/blas.h"
#include "kronlib.h"
#include "matrix.h"
#include <thrust/copy.h>
#include <thrust/iterator/counting_iterator.h>
namespace kronlib {

template <typename Matrix, typename = void>
class Cholesky {
    using Storage = typename Matrix::Storage;
    using T = typename Storage::value_type;
public:
    Storage L;
    size_t n;
    const Matrix matrix;
public:
    Cholesky(const Matrix & mat) : n{mat.nR()}, matrix{mat} {
        L = matrix.getConstData();
        potrf(matrix.getContext(), CblasColMajor, 'L', n, L, n);
        thrust::counting_iterator<size_t> idx{0};
        zero_if_not_lower_triangular<T> zinlt{n};
        thrust::transform(idx, idx + n*n, L.begin(), L.begin(), zinlt);
    }
    Matrix inv() const{
        Matrix ans = cholmat();
        potri(matrix.getContext(), CblasColMajor, 'L', n, ans.getMutableData(), n);
        return ans;
    }
    Matrix solve(const Matrix & other) const{
        Matrix ans = other;
        int n = other.nR();
        int m = other.nC();
        potrs(matrix.getContext(), CblasColMajor, 'L', n, m, L, n, ans.getMutableData(), n);
        return ans;
    }
    T logdet() const {
        thrust::counting_iterator<size_t> idx(0);
        auto diag_idxs = thrust::make_transform_iterator(idx, mult_by<size_t>(n + 1));
        auto log_vals = thrust::make_transform_iterator(L.begin(), logarithm<T>());
        auto log_diags = thrust::make_permutation_iterator(log_vals, diag_idxs);
        T half_ld = thrust::reduce(log_diags, log_diags + n);
        return 2 * half_ld;
    }
    Matrix cholmat() const {
        // Purpose is: return a matrix that is the lower triangular
        return Matrix(matrix.getContext(), L, n, n);
    }
    Matrix operator*(const Matrix & other) {
        return matrix * other;
    }
    long N() const {
        return n;
    }

};


//////////////////////////////////////// 
// IF CUDA < 7 CHOLESKY NOT SUPPORTED //
// COPY MATRIX TO HOST                //
////////////////////////////////////////
template <typename T>
class Cholesky<Matrix<backend::CUDAContext<T>>, typename std::enable_if< (backend::CUDAContext<T>::version < 7) >::type> {
    private:
        Cholesky();
};

};
#endif //KRONMAT_CHOLESKY_H
