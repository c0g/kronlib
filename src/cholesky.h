//
// Created by Thomas Nickson on 9/07/15.
//

#ifndef KRONMAT_CHOLESKY_H
#define KRONMAT_CHOLESKY_H

#include <vector>
#include "matrix.h"

template <typename T>
class Cholesky {
public:
    CBLAS_ORDER order = CblasRowMajor;
    std::vector<T> L;
    long n;
public:
    Cholesky(const Matrix<T> & mat) {
        assert(mat.nR() == mat.nC());
        n = mat.nR();
        L = mat.getConstData();
        order = mat.getOrder();
        potrf(mat.getOrder(), 'L', n, L.data(), n);
    }
    Matrix<T> inv() const { 
        std::vector<T> data = L;
        potri(CblasRowMajor, 'L', n, data.data(), n);
        Matrix<T> mat(std::move(data), n, n, n, false);

        //Fill upper triangle
        for (int r = 1; r < n; ++r) {
            for (int c = 0; c < r; ++c) {
                mat(c, r) = mat(r, c);
            }
        }
        return mat;
    };
    Matrix<T> solve(const Matrix<T> & other) const{
        std::vector<T> B(other.getConstData());
        int n = other.nR();
        int m = other.nC();
        potrs(order, 'L', n, m, L.data(), n, B.data(), m);
        return Matrix<T>(std::move(B), n, m, m, false);
    }
    T logdet() const {
        T ld = 0;
        for (int r = 0; r < n; ++r) {
            ld += std::log(L[r * n + r]);
        }
        return 2 * ld;
    }
    Matrix<T> cholmat() const {
        return Matrix<T>(L, n, n, n, false);
    }
    long N() const {
        return n;
    }
    T operator()(long r, long c) const {
        if (c > r) {
            return 0;
        } else {
            return L[r * n + c];
        }
    }
};

template <typename T> 
bool operator== (const Matrix<T> & mat, const Cholesky<T> & chol) {
    assert(mat.nR() == mat.nC());
    assert(mat.nR() == chol.N());
    bool match = true;
    for (int r = 0; r < chol.N(); ++r){
        for (int c = 0; c <= r; ++c) {
            match &= (mat(r, c) == chol(r, c));
        }
    }
    return match;
}

template <typename T> 
bool operator== (const Cholesky<T> & chol, const Matrix<T> & mat) {
    return (mat == chol);
}

#endif //KRONMAT_CHOLESKY_H