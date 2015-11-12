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

#endif  // KRONMAT_CHOL_VERSION

#endif // KRONMAT_GPU_LAPACKE
