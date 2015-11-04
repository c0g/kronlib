#ifndef KRONMAT_BLAS_H
#define KRONMAT_BLAS_H

#include "blas/blas_ops.h"
#include "storage.h"
#include "blas/cpu_blas.h"
#include "blas/cpu_lapack.h"

#if THRUST_DEVICE_SYSTEM == THRUST_DEVICE_SYSTEM_CUDA // we're on the GPU, add template specialisations for device_vector
#include "blas/gpu_blas.h"
#include "blas/gpu_lapack.h"
#endif // THRUST_DEVICE_SYSTEM_CUDA

#endif // KRONMAT_BLAS_H
