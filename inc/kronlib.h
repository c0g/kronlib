#ifndef KRONLIB_H
#define KRONLIB_H

#include "backend.h"
#include "matrix.h"
namespace kronlib {

template <typename T>
using HostMatrix = Matrix<backend::HostContext<T>>;

template <typename T>
using CUDAMatrix = Matrix<backend::CUDAContext<T>>;

template <typename T>
using TBBMatrix = Matrix<backend::TBBContext<T>>;

};


#endif 
