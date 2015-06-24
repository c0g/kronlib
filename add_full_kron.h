//
// Created by Thomas Nickson on 17/06/15.
//

#ifndef KRONMAT_ADD_FULL_KRON_H
#define KRONMAT_ADD_FULL_KRON_H

#include <vector>
#include "blas_wrap.h"
#include "basic_kronmat.h"

using namespace tommat;

template <typename T>
std::vector<T> vec_dot_kvs(const basic_kvs<T> & kvs, const std::vector<T> & obs)
{
    std::vector<long> shapes = kvs.shapes;
    long full_width = 1;
    for (auto width : shapes) full_width *= width;


    std::vector<T> kroned_vecs;
    kroned_vecs.resize(full_width, 0);


    for (int n = 0; n < kvs.height; ++n) {
        T alpha = -2 * obs[n];
        long current_width = shapes[0];
        std::vector<T> kronme(
            kvs.data[0].begin() + n * shapes[0],
            kvs.data[0].begin() + (n+1) * shapes[0]
        );
        for (int d = 1; d < kvs.dim - 1; ++d) {
            std::vector<T> data(kvs.data[d].begin() + n * shapes[d],kvs.data[d].begin() + (n  + 1) * shapes[d] );
            std::vector<T> destination;
            destination.resize(current_width * shapes[d], 0);

            blas_ger(current_width, shapes[d], alpha,
                     kronme.data(), 1,
                     data.data(), 1,
                     destination.data(), shapes[d]);
            alpha = 1.0;
            kronme = destination;
            current_width *= shapes[d];
        }
        std::vector<T> data (
            kvs.data[kvs.dim-1].begin() + n * shapes[kvs.dim-1],
            kvs.data[kvs.dim-1].begin() + (n  + 1) * shapes[kvs.dim-1]
        );
        blas_ger(current_width, shapes[kvs.dim - 1], 1.0,
                kronme.data(), 1,
                data.data(), 1,
                kroned_vecs.data(), shapes[kvs.dim-1]);
    }
    return kroned_vecs;
}

#endif //KRONMAT_ADD_FULL_KRON_H
