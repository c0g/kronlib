// Created by Thomas Nickson on 12/06/15.
//

#ifndef KRONMAT_CPU_BLAS_H
#define KRONMAT_CPU_BLAS_H

#include "storage.h"
#include <cblas.h>
#include "blas_ops.h"

template <typename Storage, typename State> // 'state' is used in CUDA
inline void blas_gemm(State, const enum BlasOrder Order,
                      const enum BlasTranspose TransA,
                      const enum BlasTranspose TransB, const int M, const int N,
                      const int K, const typename Storage::value_type * alpha, const Storage & vecA,
                      const int lda, const Storage & vecB, const int ldb, 
                      const typename Storage::value_type * beta, Storage & vecC, const int ldc)
{
    using T = typename Storage::value_type;
    const T * A = thrust::raw_pointer_cast(vecA.data());
    const T * B = thrust::raw_pointer_cast(vecB.data());


    T * C = thrust::raw_pointer_cast(vecC.data());
    cpu_blas_gemm(cblasOrder(Order), cblasTranspose(TransA), cblasTranspose(TransB), 
                M, N, K,
                *alpha, A, lda,
                B, ldb,
                *beta, C, ldc);
}
inline void cpu_blas_gemm(const enum CBLAS_ORDER Order,
                      const enum CBLAS_TRANSPOSE TransA,
                      const enum CBLAS_TRANSPOSE TransB, const int M, const int N,
                      const int K, const double alpha, const double *A,
                      const int lda, const double *B, const int ldb,
                      const double beta, double *C, const int ldc)
{
    cblas_dgemm(Order, TransA, TransB, M, N, K,
                alpha, A, lda,
                B, ldb,
                beta, C, ldc);
}

inline void cpu_blas_gemm(const enum CBLAS_ORDER Order,
                      const enum CBLAS_TRANSPOSE TransA,
                      const enum CBLAS_TRANSPOSE TransB, const int M, const int N,
                      const int K, const float alpha, const float *A,
                      const int lda, const float *B, const int ldb,
                      const float beta, float *C, const int ldc)
{
    cblas_sgemm(Order, TransA, TransB, M, N, K,
                alpha, A, lda,
                B, ldb,
                beta, C, ldc);
}

template <typename Storage>
inline typename Storage::value_type blas_dot(int N, const Storage & vecA, int strideA, const Storage &  vecB, int strideB)
{
    using T = typename Storage::value_type;
    const T * A = thrust::raw_pointer_cast(vecA.data());
    const T * B = thrust::raw_pointer_cast(vecB.data());
    return nt_blas_dot(N, A, strideA, B, strideB);
}

inline float nt_blas_dot(int N, const float * A, int strideA, const float * B, int strideB)
{
    return cblas_sdot(N, A, strideA, B, strideB);
}

inline double nt_blas_dot(int N, const double * A, int strideA, const double * B, int strideB)
{
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
#endif
