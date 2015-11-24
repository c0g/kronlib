#ifndef BLAS_OPS_H
#define BLAS_OPS_H
#include "cblas.h"
#include "cublas_v2.h"

enum BlasTranspose {None, Trans, Conj};
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
        case None:
            return CblasNoTrans;
        case Trans:
            return CblasTrans;
        case Conj:
            return CblasConjTrans;
        default:
            return CblasNoTrans;
    }
}

cublasOperation_t cublasTranspose(BlasTranspose trans) {
    switch(trans) {
        case None:
            return CUBLAS_OP_N;
        case Trans:
            return CUBLAS_OP_T;
        case Conj:
            return CUBLAS_OP_C;
        default:
            return CUBLAS_OP_N;
    }
}

cublasFillMode_t cublasFillMode(char UPLO) {
    switch (UPLO) {
        case 'L': 
            return CUBLAS_FILL_MODE_LOWER;
        case 'U':
            return CUBLAS_FILL_MODE_UPPER;
        default:
            return CUBLAS_FILL_MODE_LOWER;
    }
}

#endif
