// Created by Thomas Nickson on 12/06/15.
//

#ifndef KRONMAT_CPU_LAPACKE_H
#define KRONMAT_CPU_LAPACKE_H


#include <cblas.h>
#include <lapacke.h>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>


inline void LAPACKE_potrf(int order, char UPLO, int N, float * a, int lda) {
    LAPACKE_spotrf(order, UPLO, N, a, lda);
}

inline void LAPACKE_potrf(int order, char UPLO, int N, double * a, int lda) {
    LAPACKE_dpotrf(order, UPLO, N, a, lda);
}
template <typename Storage, typename Handle> //handle used for CUDA versions
void potrf(Handle, int order, char UPLO, int N, Storage & a, int lda) {
    std::cout << "In CPU potrf" << std::endl;
    using T = typename Storage::value_type;
    T * aptr = thrust::raw_pointer_cast(a.data());
    LAPACKE_potrf(order, UPLO, N, aptr, lda);
}


inline void LAPACKE_potri(int order, char UPLO, lapack_int n, float * a, lapack_int lda) {
    LAPACKE_spotri(order, UPLO, n, a, lda);
}
inline void LAPACKE_potri(int order, char UPLO, lapack_int n, double * a, lapack_int lda) {
    LAPACKE_dpotri(order, UPLO, n, a, lda);
}
template <typename Storage, typename Handle>
inline void potri(Handle, int order, char UPLO, lapack_int n, Storage & a, lapack_int lda) {
    using T = typename Storage::value_type;
    T * aptr = thrust::raw_pointer_cast(a.data());
    LAPACKE_potri(order, UPLO, n, aptr, lda);
}

inline void LAPACKE_potrs(int order, char UPLO, lapack_int n, lapack_int nrhs, const float * a, lapack_int lda, float * b, lapack_int ldb) {
    LAPACKE_spotrs(order, UPLO, n, nrhs, a, lda, b, ldb);
}
inline void LAPACKE_potrs(int order, char UPLO, lapack_int n, lapack_int nrhs, const double * a, lapack_int lda, double * b, lapack_int ldb) {
    LAPACKE_dpotrs(order, UPLO, n, nrhs, a, lda, b, ldb);
}
template <typename Storage, typename Handle>
inline void potrs(Handle, int order, char UPLO, lapack_int n, lapack_int nrhs, const Storage & a, lapack_int lda, Storage & b, lapack_int ldb) {
    using T = typename Storage::value_type;
    const T * aptr = thrust::raw_pointer_cast(a.data());
    T * bptr = thrust::raw_pointer_cast(b.data());
    LAPACKE_potrs(order, UPLO, n, nrhs, aptr, lda, bptr, ldb);
}
#endif //KRONMAT_CPU_LAPACKE
