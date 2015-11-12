// Created by Thomas Nickson on 12/06/15.
//

#ifndef KRONMAT_BLAS_WRAP_H
#define KRONMAT_BLAS_WRAP_H


#include <cblas.h>
#include <lapacke.h>


inline void blas_gemm(const enum CBLAS_ORDER Order,
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

inline void blas_gemm(const enum CBLAS_ORDER Order,
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

inline float blas_dot(int N, const float * A, int strideA, const float * B, int strideB)
{
    return cblas_sdot(N, A, strideA, B, strideB);
}

inline double blas_dot(int N, const double * A, int strideA, const double * B, int strideB)
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

inline void potrf(int order, char UPLO, int N, float * a, int lda) {
    LAPACKE_spotrf(order, UPLO, N, a, lda);
}

inline void potrf(int order, char UPLO, int N, double * a, int lda) {
    LAPACKE_dpotrf(order, UPLO, N, a, lda);
}
inline void potri(int order, char UPLO, lapack_int n, float * a, lapack_int lda) {
    LAPACKE_spotri(order, UPLO, n, a, lda);
}
inline void potri(int order, char UPLO, lapack_int n, double * a, lapack_int lda) {
    LAPACKE_dpotri(order, UPLO, n, a, lda);
}
inline void potrs(int order, char UPLO, lapack_int n, lapack_int nrhs, const float * a, lapack_int lda, float * b, lapack_int ldb) {
    LAPACKE_spotrs(order, UPLO, n, nrhs, a, lda, b, ldb);
}
inline void potrs(int order, char UPLO, lapack_int n, lapack_int nrhs, const double * a, lapack_int lda, double * b, lapack_int ldb) {
    LAPACKE_dpotrs(order, UPLO, n, nrhs, a, lda, b, ldb);
}
#endif //KRONMAT_BLAS_WRAP_H
