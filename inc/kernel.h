//
// Created by Thomas Nickson on 23/06/15.
//

#ifndef KRONMAT_KERNEL_H
#define KRONMAT_KERNEL_H

#include "matrix.h"
#include "distances.h"
#include <stdexcept>

template <typename T>
struct sqexp_hyp {
    T log_lengthscale;
    T log_outputscale;
    T lengthscale2() const {
        return exp(2 * log_lengthscale);
    }
    T outputscale2() const {
        return exp(2 * log_outputscale);
    }
    T dlengthscale2_dloglengthscale() const {
        return 2 * lengthscale2();
    }
    T dlengthscale2_dlogoutputscale() const {
        return 2 * outputscale2();
    }
    void set_log(int idx, T val) {
        switch (idx) {
            case 0: 
                log_lengthscale = val;
                break;
            case 1: 
                log_outputscale = val;
                break;
            default:
                throw std::out_of_range("Only 0, 1 valid idx");
        } 
    }
    T get_log(int idx) const {
        switch (idx) {
            case 0: 
                return log_lengthscale;
            case 1: 
                return log_outputscale;
            default:
                throw std::out_of_range("Only 0, 1 valid idx");
        } 
    }
    void set_real(int idx, T val) {
        switch (idx) {
            case 0: 
                log_lengthscale = std::log(val);
                break;
            case 1: 
                log_outputscale = std::log(val);
                break;
            default:
                throw std::out_of_range("Only 0, 1 valid idx");
        } 
    }
    T get_real(int idx) const {
        switch (idx) {
            case 0: 
                return std::exp(log_lengthscale);
            case 1: 
                return std::exp(log_outputscale);
            default:
                throw std::out_of_range("Only 0, 1 valid idx");
        } 
    }
};

// template <typename T>
// struct SqExpKandGrad {
//     Matrix<T> K;

// }

template <typename T>
class sqexp1d {
private:
public:
    Matrix<T> operator()(const sqexp_hyp<T> & sqexp_hyp, const Matrix<T> & X1, const Matrix<T> & X2)
    {
        auto K = pdist2(X1, X2);
        K.minus_inplace();
        K /= sqexp_hyp.lengthscale2();
        K *= 0.5;
        K.exp_inplace();
        K *= sqexp_hyp.outputscale2();
        return K;
    }

    std::vector<Matrix<T>> dhyp(const sqexp_hyp<T> & sqexp_hyp, const Matrix<T> & X1, const Matrix<T> & X2)
    {
        auto dist = pdist2(X1, X2);
        auto K = dist;
        K.minus_inplace();
        K /= sqexp_hyp.lengthscale2();
        K *= 0.5;
        K.exp_inplace();
        K *= sqexp_hyp.outputscale2();

        std::vector<Matrix<T>> ans;
        auto KhadDist = K.hadamard(dist);
        KhadDist /= sqexp_hyp.lengthscale2();
        ans.push_back(KhadDist);

        K *= 2;
        ans.push_back(K); 
        return ans;
    }
};

#endif //KRONMAT_KERNEL_H
