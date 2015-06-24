//
// Created by Thomas Nickson on 23/06/15.
//

#ifndef KRONMAT_DISTANCES_H
#define KRONMAT_DISTANCES_H

#include <iostream>
#include "matrix.h"

template <typename T>
inline matrix<T> pdist2(const matrix<T> & A, const matrix<T> & B)
{
    matrix<T> ans(A.nR(), B.nR());
    for (int r = 0; r < A.nR(); r+=2) {
        for (int c = 0; c < B.nR(); c+=2) {
            ans(r, c) = A(r, 0) * A(r, 0) + B(c, 0) * B(c, 0);
            ans(r, c+1) = A(r, 0) * A(r, 0) + B(c+1, 0) * B(c+1, 0);
            ans(r+1, c) = A(r+1, 0) * A(r+1, 0) + B(c, 0) * B(c, 0);
            ans(r+1, c+1) = A(r+1, 0) * A(r+1, 0) + B(c+1, 0) * B(c+1, 0);
        }
    }
    blas_gemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
              A.nR(), B.nR(), 1,
              -2.0, A.dataptr(), 1, B.dataptr(), B.nR(),
              1.0, ans.mutable_dataptr(), ans.nC());
    return ans;
}

#endif //KRONMAT_DISTANCES_H
