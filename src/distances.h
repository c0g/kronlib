//
// Created by Thomas Nickson on 23/06/15.
//

#ifndef KRONMAT_DISTANCES_H
#define KRONMAT_DISTANCES_H

#include <iostream>
#include "matrix.h"

template <typename T>
inline Matrix<T> pdist2(const Matrix<T> & A, const Matrix<T> & B)
{
    Matrix<T> ans(A.nR(), B.nR());
    for (int r = 0; r < A.nR(); ++(++r)) {
        for (int c = 0; c < B.nR(); ++(++c)) {
            ans(r, c) = A(r, 0) * A(r, 0) + B(c, 0) * B(c, 0);
            ans(r, c+1) = A(r, 0) * A(r, 0) + B(c+1, 0) * B(c+1, 0);
            ans(r+1, c) = A(r+1, 0) * A(r+1, 0) + B(c, 0) * B(c, 0);
            ans(r+1, c+1) = A(r+1, 0) * A(r+1, 0) + B(c+1, 0) * B(c+1, 0);
        }
    }
    blas_gemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
              A.nR(), B.nR(), 1,
              -2.0, A.getConstData().data(), 1, B.getConstData().data(), B.nR(),
              1.0, ans.getMutableData().data(), B.nR());
    
    return ans;
}

#endif //KRONMAT_DISTANCES_H
