#ifndef STORAGE_H
#define STORAGE_H

#include <thrust/device_vector.h>
#include <thrust/host_vector.h>

template <typename T>
using device = thrust::device_vector<T>;

template <typename T>
using host = thrust::host_vector<T>;

#endif
