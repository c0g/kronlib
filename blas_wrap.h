// Created by Thomas Nickson on 12/06/15.
//

#ifndef KRONMAT_BLAS_WRAP_H
#define KRONMAT_BLAS_WRAP_H

namespace tommat {
    extern "C" {
#include <cblas.h>
    };


    void blas_gemm(const enum CBLAS_ORDER __Order,
                   const enum CBLAS_TRANSPOSE __TransA,
                   const enum CBLAS_TRANSPOSE __TransB, const int __M, const int __N,
                   const int __K, const double __alpha, const double *__A,
                   const int __lda, const double *__B, const int __ldb,
                   const double __beta, double *__C, const int __ldc) {
        cblas_dgemm(__Order, __TransA, __TransB, __M, __N, __K,
                    __alpha, __A, __lda,
                    __B, __ldb,
                    __beta, __C, __ldc);
    }

    void blas_gemm(const enum CBLAS_ORDER __Order,
                   const enum CBLAS_TRANSPOSE __TransA,
                   const enum CBLAS_TRANSPOSE __TransB, const int __M, const int __N,
                   const int __K, const float __alpha, const float *__A,
                   const int __lda, const float *__B, const int __ldb,
                   const float __beta, float *__C, const int __ldc) {
        cblas_sgemm(__Order, __TransA, __TransB, __M, __N, __K,
                    __alpha, __A, __lda,
                    __B, __ldb,
                    __beta, __C, __ldc);
    }
}
#endif //KRONMAT_BLAS_WRAP_H
