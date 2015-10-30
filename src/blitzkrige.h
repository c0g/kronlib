//
// Created by Thomas Nickson on 13/07/15.
//

#ifndef BLITZKRIGE_LIKELIHOOD_H
#define BLITZKRIGE_LIKELIHOOD_H

#include "kernel.h"
#include "matrix.h"


template <typename T>
class Blitzkrige {
public:
    T operator(const sqexp_hyp<T> & hyp) const;
    Matrix<T> d(const sqexp_hyp<T> & hyp) const;
    Matrix<T> predict(const sqexp_hyp<T> & hyp, const Matrix<T> & Xp) const;
private:

    Matrix<T> dlik_dkern(const sqexp_hyp<T> & hyp);
    
}