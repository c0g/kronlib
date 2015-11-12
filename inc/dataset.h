//
// Created by Thomas Nickson on 23/06/15.
//

#ifndef KRONMAT_DATASET_H
#define KRONMAT_DATASET_H

#include "matrix.h"

template <typename T>

struct dataset {
    dataset(const matrix<T> &X, const matrix<T> &Y) : X(X), Y(Y), n(X.nR()), d(X.nC()), out_d(Y.nC())
    {
        assert(X.nR() == Y.nR());
    }
    long n;
    long d;
    long out_d;
    const matrix<T> X;
    const matrix<T> Y;
};

#endif //KRONMAT_DATASET_H
