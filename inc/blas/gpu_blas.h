// Created by Thomas Nickson on 12/06/15.
//

#ifndef KRONMAT_CUDA_BLAS_H
#define KRONMAT_CUDA_BLAS_H


#include "blas_ops.h"
#include <cublas_v2.h>
#include "storage.h"
namespace kronlib {

inline void cuda_blas_gemm(const cublasHandle_t handle,
                      const cublasOperation_t TransA,
                      const cublasOperation_t TransB, const int M, const int N,
                      const int K, const double * alpha, const double *A,
                      const int lda, const double *B, const int ldb,
                      const double * beta, double *C, const int ldc)
{
    cublasDgemm(handle, TransA, TransB, M, N, K,
                alpha, A, lda,
                B, ldb,
                beta, C, ldc);
}

inline void cuda_blas_gemm(const cublasHandle_t handle,
                      const cublasOperation_t TransA,
                      const cublasOperation_t TransB, const int M, const int N,
                      const int K, const float * alpha, const float *A,
                      const int lda, const float *B, const int ldb,
                      const float * beta, float *C, const int ldc)
{
    cublasSgemm(handle, TransA, TransB, M, N, K,
                alpha, A, lda,
                B, ldb,
                beta, C, ldc);
}

template <typename T>
inline void blas_gemm(const std::shared_ptr<kronlib::backend::CUDAContext<T>> context,
                      const enum BlasOrder Order,
                      const enum BlasTranspose TransA,
                      const enum BlasTranspose TransB, const int M, const int N,
                      const int K, const T * alpha, const typename kronlib::backend::CUDAContext<T>::Storage & vecA,
                      const int lda, const typename kronlib::backend::CUDAContext<T>::Storage & vecB, const int ldb,
                      const T * beta, typename kronlib::backend::CUDAContext<T>::Storage & vecC, const int ldc)
{
    const T * A = thrust::raw_pointer_cast(vecA.data());
    const T * B = thrust::raw_pointer_cast(vecB.data());
    T * C = thrust::raw_pointer_cast(vecC.data());
    cuda_blas_gemm(context->blas, cublasTranspose(TransA), cublasTranspose(TransB), M, N, K,
                alpha, A, lda,
                B, ldb,
                beta, C, ldc);
}
};
#endif // CUDA_BLAS_H

