//
// Created by Thomas Nickson on 23/06/15.
//

#ifndef KRONMAT_KERNEL_H
#define KRONMAT_KERNEL_H

#include "matrix.h"
#include "distances.h"

template <typename T>
struct sqexp_hyp {
    T log_lengthscale;
    T log_outputscale;
    T lengthscale2() {
        return exp(2 * log_lengthscale);
    }
    T outputscale2() {
        return exp(2 * log_outputscale);
    }
    T dlengthscale2_dloglengthscale() {
        return 2 * lengthscale2();
    }
    T dlengthscale2_dloglengthscale() {
        return 2 * lengthscale2();
    }
};

template <typename T>
class sqexp1d {
private:

public:
    matrix operator()(const & sqexp_hyp, const matrix & X1, const matrix & X2)
    {
        auto sq_dist = pdist2(X1, X2);
    }
};

#endif //KRONMAT_KERNEL_H
