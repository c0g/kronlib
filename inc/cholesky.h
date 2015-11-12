//
// Created by Thomas Nickson on 9/07/15.
//

#ifndef KRONMAT_CHOLESKY_H
#define KRONMAT_CHOLESKY_H

#include <vector>
#include "matrix.h"
#include <thrust/copy.h>

template <typename Storage>
class Cholesky {
    using T = typename Storage::value_type;
public:
    Storage L;
    size_t n;
public:
    Cholesky(const Matrix<Storage> & mat) {
        assert(mat.nR() == mat.nC());
        n = mat.nR();
        L = mat.getConstData();
        potrf(CblasColMajor, 'L', n, L, n);
    }
    /*
    Matrix<Storage> inv() const { 
        Matrix<Storage> ans{L, n, n};
        potri(CblasColMajor, 'L', n, ans.data.data(), n);
        //Fill upper triangle
        for (int r = 1; r < n; ++r) {
            for (int c = 0; c < r; ++c) {
                ans(c, r) = ans(r, c);
            }
        }
        return ans;
    };
    */
    Matrix<Storage> solve(const Matrix<Storage> & other) const{
        Matrix<Storage> ans = other;
        int n = other.nR();
        int m = other.nC();
        potrs(CblasColMajor, 'L', n, m, L, n, ans.getMutableData(), n);
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
    Matrix<Storage> cholmat() const {
        // Initalize ans as a n,n matrix filled with L
        Matrix<Storage> ans(L, n, n);
        // zero if not lower triangular either returns the value if inside the lower triangle or 0
        zero_if_not_lower_triangular<T> zinlt{n};
        // provide counting iterator to give index
        thrust::counting_iterator<size_t> idx(0);
        thrust::transform(idx, idx + n*n, ans.getConstData().begin(), ans.getMutableData().begin(), zinlt);
        return ans;
    }
    long N() const {
        return n;
    }
};

template <typename Storage> 
bool operator== (const Matrix<Storage> & mat, const Cholesky<Storage> & chol) {
    assert(mat.nR() == mat.nC());
    assert(mat.nR() == chol.N());
    return chol.cholmat() == mat;
}

template <typename Storage> 
bool operator== (const Cholesky<Storage> & chol, const Matrix<Storage> & mat) {
    return (mat == chol);
}

#endif //KRONMAT_CHOLESKY_H
