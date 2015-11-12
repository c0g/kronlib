#ifndef KRONMAT_BACKEND_H
#define KRONMAT_BACKEND_H

#include <memory>

#include <thrust/system/cpp/execution_policy.h>
#include <thrust/system/tbb/execution_policy.h>
#include <thrust/system/cuda/execution_policy.h>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/system/cuda/vector.h>
#include <thrust/system/tbb/vector.h>

#include <cuda_runtime.h>
#include <cublas_v2.h>
#include <cusolverDn.h>


namespace kronlib {
    namespace backend {
        template <typename T>
        class HostContext {
            private:
                T zero_val = 0;
                T one_val = 1;
            public:
                using Storage = thrust::host_vector<T>;
                thrust::system::cpp::detail::par_t exec() const { return thrust::cpp::par; }
                T * one() { return &one_val; };
                T * zero() { return &zero_val; };
        };

        template <typename T>
        class TBBContext {
            private:
                T zero_val = 0;
                T one_val = 1;
            public:
                using Storage = thrust::system::tbb::vector<T>;
                thrust::system::tbb::detail::par_t exec() const { return thrust::tbb::par; }
                T * one() { return &one_val; };
                T * zero() { return &zero_val; };
        };

        template <typename T>
        class CUDAContext {
            private: 
                T * devMemory;
            public:
                cudaStream_t stream;
                cublasHandle_t blas;
#ifdef CUSOLVER
                cusolverDnHandle_t solver;
#endif

                using Storage = thrust::system::cuda::vector<T>;
                CUDAContext(cudaStream_t _stream = 0) {
                    stream = _stream;
                    T vals[2];
                    vals[0] = 0; vals[1] = 1;
                    cudaMalloc(&devMemory, sizeof(T) * 2);
                    cudaMemcpy(devMemory, vals, sizeof(T) * 2, cudaMemcpyHostToDevice);
                    cublasCreate(&blas);
                    cublasSetStream(blas, stream);
                    cublasSetPointerMode(blas, CUBLAS_POINTER_MODE_DEVICE);
                    cublasSetAtomicsMode(blas, CUBLAS_ATOMICS_ALLOWED);
#ifdef CUSOLVER
                    cusolverDnCreate(&solver);
                    cusolverDnSetStream(solver, stream);
#endif
                }
                ~CUDAContext() {
#ifdef CUSOLVER
                    cusolverDnDestroy(solver);
#endif
                    cublasDestroy(blas);
                    cudaFree(devMemory);
                }
                thrust::system::cuda::detail::execute_on_stream exec() const { return thrust::cuda::par.on(stream); }
                T * one() { return &devMemory[1]; };
                T * zero() { return &devMemory[0]; };
        };
        template <typename T>
        using Host = std::shared_ptr<HostContext<T>>;
        template <typename T>
        using TBB = std::shared_ptr<TBBContext<T>>;
        template <typename T>
        using CUDA = std::shared_ptr<CUDAContext<T>>;
        template <typename T>
        Host<T> getHostContext() { return std::make_shared<HostContext<T>>(); };
        template <typename T>
        TBB<T> getTBBContext() { return std::make_shared<TBBContext<T>>(); };
        template <typename T>
        CUDA<T> getCUDAContext() { return std::make_shared<CUDAContext<T>>(); };
        template <typename T>
        CUDA<T> getCUDAContext(cudaStream_t stream) { return std::make_shared<CUDAContext<T>>(stream); };
    };
};


#endif // KRONMAT_BACKEND_H
