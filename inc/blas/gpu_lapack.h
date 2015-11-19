// Created by Thomas Nickson on 12/06/15.
//

#ifndef KRONMAT_GPU_LAPACKE
#define KRONMAT_GPU_LAPACKE

#include <thrust/device_vector.h>
#include <thrust/host_vector.h>

#if CUSOLVER

#include <cusolverDn.h>


    void CUSOLVER_potrf(cusolverDnHandle_t handle, cublasFillMode_t UPLO, 
            int n, double * Aptr, size_t lda, int * minor) {
        int Lwork;
        cusolverStatus_t status;
        status =  cusolverDnDpotrf_bufferSize(handle, UPLO, n, Aptr, lda, &Lwork); // TODO: Check these errors (As if).
        if (status != CUSOLVER_STATUS_SUCCESS) 
        {
            std::cout << "Buffer size failed at " << __LINE__ << " in " << __FILE__ << " with code " << std::endl;
        }
        double * workPtr;
        cudaMalloc(&workPtr, Lwork * sizeof(double));
        status = cusolverDnDpotrf(handle, UPLO, n, Aptr, lda, workPtr, Lwork, minor); 
        if (status != CUSOLVER_STATUS_SUCCESS) 
        {
            std::cout << "POTRF failed at " << __LINE__ << " in " << __FILE__ << " with code " << status << std::endl;
        }
        cudaFree(workPtr);
    }
    void CUSOLVER_potrf(cusolverDnHandle_t handle, cublasFillMode_t UPLO, 
            int n, float * Aptr, size_t lda, int * minor) {
        int Lwork;
        cusolverStatus_t status;
        status =  cusolverDnSpotrf_bufferSize(handle, UPLO, n, Aptr, lda, &Lwork); // TODO: Check these errors (As if).
        if (status != CUSOLVER_STATUS_SUCCESS) 
        {
            std::cout << "Buffer size failed at " << __LINE__ << " in " << __FILE__ << " with code " << std::endl;
        }
        float * workPtr;
        cudaMalloc(&workPtr, Lwork * sizeof(float));
        status = cusolverDnSpotrf(handle, UPLO, n, Aptr, lda, workPtr, Lwork, minor); 
        if (status != CUSOLVER_STATUS_SUCCESS) 
        {
            std::cout << "POTRF failed at " << __LINE__ << " in " << __FILE__ << " with code " << status << std::endl;
        }
        cudaFree(workPtr);
    }

    template <typename T>
    inline void potrf(const std::shared_ptr<kronlib::backend::CUDAContext<T>> context,
            int order, char UPLO, int N, typename kronlib::backend::CUDAContext<T>::Storage & a, int lda)
    {
        T * Aptr = thrust::raw_pointer_cast(a.data());
        CUSOLVER_potrf(context->solver, cublasFillMode(UPLO), N, Aptr, lda, context->minor);
        if (context->getMinor() != 0) 
        {
            std::cout << "Matrix not PD, failing minor " << context->getMinor() << std::endl;
        }
    }

        

    void CUSOLVER_potrs(cusolverDnHandle_t handle, cublasFillMode_t UPLO, 
            int n, int nrhs, const double * Aptr, size_t lda, double * Bptr, size_t ldb, int * minor) {
        cusolverStatus_t status;
        status = cusolverDnDpotrs(handle, UPLO, n, nrhs, Aptr, lda, Bptr, ldb, minor); 
        if (status != CUSOLVER_STATUS_SUCCESS) 
        {
            std::cout << "POTRS failed at " << __LINE__ << " in " << __FILE__ << " with code " << status << std::endl;
        }
    }
    void CUSOLVER_potrs(cusolverDnHandle_t handle, cublasFillMode_t UPLO, 
            int n, int nrhs, const float * Aptr, size_t lda, float * Bptr, size_t ldb, int * minor) {
        cusolverStatus_t status;
        status = cusolverDnSpotrs(handle, UPLO, n, nrhs, Aptr, lda, Bptr, ldb, minor); 
        if (status != CUSOLVER_STATUS_SUCCESS) 
        {
            std::cout << "POTRS failed at " << __LINE__ << " in " << __FILE__ << " with code " << status << std::endl;
        }
    }

    template <typename T>
    inline void potrs(const std::shared_ptr<kronlib::backend::CUDAContext<T>> context,
            int order, char UPLO, int N, int nrhs, const typename kronlib::backend::CUDAContext<T>::Storage & a, int lda,
            typename kronlib::backend::CUDAContext<T>::Storage & b, int ldb)
    {
        const T * Aptr = thrust::raw_pointer_cast(a.data());
        T * Bptr = thrust::raw_pointer_cast(b.data());
        CUSOLVER_potrs(context->solver, cublasFillMode(UPLO), N, nrhs, Aptr, lda, Bptr, ldb, context->minor);
        if (context->getMinor() != 0) 
        {
            std::cout << "Matrix not PD, failing minor " << context->getMinor() << std::endl;
        }
    }

        

            


#else  // Not CUSOLVER
#endif  // CUSOLVER

#endif // KRONMAT_GPU_LAPACKE
