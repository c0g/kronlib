// Created by Thomas Nickson on 12/06/15.
//

#ifndef KRONMAT_GPU_LAPACKE
#define KRONMAT_GPU_LAPACKE

#include <thrust/device_vector.h>
#include <thrust/host_vector.h>

#if KRONMAT_CHOL_VERSION == CUSTOM
// TODO: get rid of copy back to host for Cholesky

#include <lapacke.h>


#elif KRONMAT_CHOL_VERSION == CUDA

inline void LAPACKE_potrf(int order, char UPLO, int N, float * a, int lda) {
    LAPACKE_spotrf(order, UPLO, N, a, lda);
}

inline void CUSOLVER_potrf(cusolverDnHandle_t handle, int order, char UPLO, int N, double * a, int lda) {
    LAPACKE_dpotrf(order, UPLO, N, a, lda);
}
template <typename Storage> //handle used for CUDA versions
inline void potrf(cusolverDnHandle_t handle, int order, char UPLO, int N, Storage & a, int lda) {
    using T = typename Storage::value_type;
    T * aptr = thrust::raw_pointer_cast(a.data());
    CUSOLVER_potrf(order, UPLO, N, aptr, lda);
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

#endif  // KRONMAT_CHOL_VERSION

#endif // KRONMAT_GPU_LAPACKE
