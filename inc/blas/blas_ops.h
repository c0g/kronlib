#ifndef BLAS_OPS_H
#define BLAS_OPS_H
#include "cblas.h"
#include "cublas_v2.h"

enum BlasTranspose {N, T, C};
enum BlasOrder {COL, ROW};

CBLAS_ORDER cblasOrder(BlasOrder order) {
    switch (order) {
        case COL:
            return CblasColMajor;
        case ROW:
            return CblasRowMajor;
        default:
            return CblasColMajor;
    }
}

CBLAS_TRANSPOSE cblasTranspose(BlasTranspose trans) {
    switch (trans) {
        case N:
            return CblasNoTrans;
        case T:
            return CblasTrans;
        case C:
            return CblasConjTrans;
        default:
            return CblasNoTrans;
    }
}

cublasOperation_t cublasTranspose(BlasTranspose trans) {
    switch(trans) {
        case N:
            return CUBLAS_OP_N;
        case T:
            return CUBLAS_OP_T;
        case C:
            return CUBLAS_OP_C;
        default:
            return CUBLAS_OP_N;
    }
}

#endif
