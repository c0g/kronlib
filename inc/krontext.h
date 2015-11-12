// Copyright: Tom Nickson 2015
#ifndef KRONMAT_KRONTEXT_H
#define KRONMAT_KRONTEXT_H

#include "storage.h"
#if THRUST_DEVICE_SYSTEM == THRUST_DEVICE_SYSTEM_CUDA
#include <cublas_v2.h>
  #if !(KRONMAT_CHOL_VERSION == CUSTOM)
    #include <cusolverDn.h>
  #endif
#endif


namespace kronmat {

    template<typename Backend>
    class Context<Backend
};
/*
template <typename T>
class Context { // Context for the KronMat library
    T * hzero;
    T * hone;
    T * dzero;
    T * done;

    T * hvals;
    T * dvals;
    
#if THRUST_DEVICE_SYSTEM == THRUST_DEVICE_SYSTEM_CUDA
    cublasHandle_t blasHandle;
  #if !(KRONMAT_CHOL_VERSION == CUSTOM)
    cusolverDnHandle_t solverHandle;
  #else
    int solverHandle = 0;
   #endif
#else 
    int solverHandle = 0;
    int blasHandle = 0;
#endif

    Context() {
        hvals = (T * ) malloc(sizeof(T) * 2);
        hvals[0] = 0; hvals[1] = 1;
        hzero = &hvals[0];
        hone = &hvals[1];
#if THRUST_DEVICE_SYSTEM == THRUST_DEVICE_SYSTEM_CUDA
        cudaMalloc(sizeof(T) * 2, &dvals);
        cudaMemCpy(dvals, hvals, sizeof(T) * 2, cudaMemCpyHostToDevice);
        dzero = &dvals[0];
        done = &dvals[1];
        cublasCreate(&blasHandle);
        cublasSetPointerMode(&blasHandle, CUBLAS_POINTER_MODE_DEVICE);
        cublasSetAtomicsMode(&blasHandle, CUBLAS_ATOMICS_ALLOWED);
  #if !(KRONMAT_CHOL_VERSION == CUSTOM)
        cusolverDnCreate(&solverHandle);      
  #endif 
#else
        dzero = hzero;
        done = hone;
#endif
                
    }
    ~Context() {
        free(vals);
#if THRUST_DEVICE_SYSTEM == THRUST_DEVICE_SYSTEM_CUDA
        cudaFree(dvals);
        cublasDestroy(&blasHandle);
  #if !(KRONMAT_CHOL_VERSION == CUSTOM)
        cusolverDnDestroy(&solverHandle);      
  #endif 
#endif 
    }
}
*/
#endif 
}
