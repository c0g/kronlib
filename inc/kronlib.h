#ifndef KRONLIB_H
#define KRONLIB_H

#include "backend.h"
#include "matrix.h"
#include "util.h"
#include "cholesky.h"
#include "kronecker_matrix.h"

namespace kronlib {

template <typename T>
using HostMatrix = Matrix<backend::HostContext<T>>;

template <typename T>
using CUDAMatrix = Matrix<backend::CUDAContext<T>>;

template <typename T>
using TBBMatrix = Matrix<backend::TBBContext<T>>;


};


#endif 
