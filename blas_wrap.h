// Created by Thomas Nickson on 12/06/15.
//

#ifndef KRONMAT_BLAS_WRAP_H
#define KRONMAT_BLAS_WRAP_H


    extern "C" {
#include <cblas.h>
    };


    inline void blas_gemm(const enum CBLAS_ORDER __Order,
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

    inline void blas_gemm(const enum CBLAS_ORDER __Order,
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

    inline float blas_dot(int N, const float * A, int strideA, const float * B, int strideB) {
        return cblas_sdot(N, A, strideA, B, strideB);
    }

    inline double blas_dot(int N, const double * A, int strideA, const double * B, int strideB) {
        return cblas_ddot(N, A, strideA, B, strideB);
    }

    inline void blas_ger(int M, int N, float alpha,
            const float * X, int incX,
            const float * Y, int incY,
            float * A, int LDA)
    {
        cblas_sger(CblasRowMajor, M, N, alpha, X, incX, Y, incY, A, LDA);
    }

    inline void blas_ger(int M, int N, double alpha,
                         const double * X, int incX,
                         const double * Y, int incY,
                         double * A, int LDA)
    {
        cblas_dger(CblasRowMajor, M, N, alpha, X, incX, Y, incY, A, LDA);
    }
    inline void blas_copy(int N, const float * X, int incX, float * Y, int incY)
    {
        cblas_scopy(N, X, incX, Y, incY);
    }
    inline void blas_axpy(int N, float mult, const float * X, int incX, float * Y, int incY)
    {
        cblas_saxpy(N, mult, X, incX, Y, incY);
    }
#endif //KRONMAT_BLAS_WRAP_H
