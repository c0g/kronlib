//
// Created by Thomas Nickson on 12/06/15.
//

#ifndef KRONMAT_MATRIX_H
#define KRONMAT_MATRIX_H

#include <vector>
#include <iostream>
#include <assert.h>
#include "blas_wrap.h"
#include <cmath>


template <typename T>
class Matrix;

template <typename T>
class KroneckerMatrix;

template <typename T>
class KroneckerVectorStack;

template <typename T>
class Cholesky;

template <typename T>
Matrix<T> exp(const Matrix<T> &);

template<typename T>
class Matrix  {
    friend KroneckerMatrix<T>;
    friend KroneckerVectorStack<T>;
    friend Cholesky<T>;
public:
    Matrix() : data{}, nr{0}, nc{0}, r_stride{0}, r0{0}, c0{0}, trans{false} {};
    Matrix(std::vector<T> data_, long r_, long c_,
           long r_stride_, bool trans_) :
        data{data_}, nr{r_}, nc{c_}, r_stride{r_stride_}, trans{trans_} {};
    Matrix(long r_, long c_) : Matrix()
    {
        set_size(r_, c_);
    };
    void set_size(long r, long c)
    {
        if (trans) {
            nr = c;
            nc = r;
        } else {
            nr = r;
            nc = c;
        }
        r_stride = c;
        data.resize(nr * nc, 0);

    }
private:
    std::vector<T> data;

    long nr;
    long nc;
    long r_stride;
    long r0 = 0;
    long c0 = 0;
    CBLAS_ORDER order = CblasRowMajor;
    bool trans;
    const long vecidx(long r, long c) const
    {
        if (trans) {
            return (c) * r_stride + (r);
        } else {
            return (r) * r_stride + (c);
        }
    }
    const long offset() const
    {
        return 0;
    }
    const T * dataptr() const
    {
        return data.data();
    }

    T * mutable_dataptr()
    {
        return data.data();
    }

public:

    const T * begindata() const
    {
        return data.data();
    }
    CBLAS_ORDER getOrder() const
    {
        return order;
    }
    const std::vector<T> & getConstData() const
    {
        return data;
    }
    std::vector<T> & getMutableData()
    {
        return data;
    }
    long nR() const
    {
        long ans = trans ? nc : nr;
        return ans;
    }
    long nC() const
    {
        long ans = trans ? nr : nc;
        return ans;
    }

    bool isTrans() const
    {
        return trans;
    }


    T trace()
    {
        assert(nR() == nC());
        T trace = 0;
        for (int idx = 0; idx < nR(); ++idx) {
            trace += (*this)(idx, idx);
        }
        return trace;
    }

    Matrix<T> transpose() const
    {
        return Matrix(data, nr, nc, r_stride, !trans);
    }

    void mutable_transpose()
    {
        trans = !trans;
    }

    Matrix<T> reshape(long newr, long newc) const
    {
        //TODO: return a 'view' if possible;
        assert(nR() * nC() == newr * newc);
        std::vector<T> newdata;
        for (int r = 0; r < nR(); ++r) {
            for (int c = 0; c < nC(); ++c) {
                newdata.push_back((*this)(r, c));
            }
        }
        return Matrix(newdata, newr, newc, newc, false);
    }

    void mutable_reshape(long newr, long newc) {
        assert(nR() * nC() == newr * newc);
        assert(!trans);
        nr = newr;
        nc = newc;
        r_stride = newc;
    }



    T& operator()(long ridx, long cidx)
    {
        assert(ridx < nR());
        assert(cidx < nC());
        return data[vecidx(ridx, cidx)];
    }
    T operator()(long ridx, long cidx) const
    {
        assert(ridx < nR());
        assert(cidx < nC());
        return data[vecidx(ridx, cidx)];
    }



    Matrix<T> operator*(const Matrix<T> &other) const
    {
        CBLAS_TRANSPOSE meTrans;
        CBLAS_TRANSPOSE otherTrans;
        int M, N, K, lda, ldb, ldc;

        assert(nC() == other.nR());

        M = nR();  //results matrices rows
        N = other.nC(); //results matrices cols
        K = nC(); // this matrices columns == other matrices rows

        ldc = other.nC();

        // check if this Matrix is transposed:
        if (isTrans()) {
            meTrans = CblasTrans;
            lda = nR();
        } else {
            meTrans = CblasNoTrans;
            lda = nC();
        }

        if (other.isTrans()) {
            otherTrans = CblasTrans;
            ldb = other.nR();
        } else {
            otherTrans = CblasNoTrans;
            ldb = other.nC();
        }

        std::vector<T> new_data;
        new_data.resize(nR() * other.nC());

        blas_gemm(CblasRowMajor, meTrans, otherTrans, M, N, K, 1.0,
                  dataptr(), lda, other.dataptr(), ldb,
                  0.0, new_data.data(), ldc);

        return Matrix(new_data, nR(), other.nC(), other.nC(), false);
    }

    Matrix<T> Tdot(const Matrix<T> &other)
    {
        mutable_transpose();
        auto ans = (*this) * other;
        mutable_transpose();
        return ans;
    }

    Matrix<T> operator*(const KroneckerVectorStack<T> kvs) {
        assert(nR() == 1 && "Only works when matrix is a row vector!");
        return vec_dot_kvs((*this), kvs);
    }

    Matrix<T> operator+(const Matrix<T> &other) const
    {
        assert(nR() == other.nR());
        assert(nC() == other.nC());

        std::vector<T> new_data;
        for (long r = 0; r < nR(); ++r) {
            for (long c = 0; c < nC(); ++c) {
                new_data.push_back((*this)(r, c) + other(r, c));
            }
        }
        return Matrix(new_data, nR(), other.nC(), other.nC(), false);
    }

    Matrix<T> operator-(const Matrix<T> &other) const
    {
        auto ans = other;
        ans.minus_inplace();
        ans += (*this);
        return ans;
    }

    Matrix<T> hadamard(const Matrix<T> &other) const
    {
        assert(nR() == other.nR());
        assert(nC() == other.nC());

        std::vector<T> new_data;
        for (long r = 0; r < nR(); ++r) {
            for (long c = 0; c < nC(); ++c) {
                new_data.push_back((*this)(r, c) * other(r, c));
            }
        }
        return Matrix(new_data, nR(), other.nC(), other.nC(), false);
    }
    Matrix<T> elemwise_div(const Matrix<T> &other) const
    {
        assert(nR() == other.nR());
        assert(nC() == other.nC());

        std::vector<T> new_data;
        for (long r = 0; r < nR(); ++r) {
            for (long c = 0; c < nC(); ++c) {
                new_data.push_back((*this)(r, c) / other(r, c));
            }
        }
        return Matrix(new_data, nR(), other.nC(), other.nC(), false);
    }

    // defined in kronecker Matrix;
    // Matrix<T> Matrix<T>::operator*(const KroneckerMatrix<T> &other);


    void operator+=(const Matrix<T> &other)
    {
        assert(nR() == other.nR());
        assert(nC() == other.nC());
        for (long r = 0; r < nR(); ++r) {
            for (long c = 0; c < nC(); ++c) {
                (*this)(r, c) += other(r, c);
            }
        }
    }

    void operator+=(const T val)
    {
        for (auto & el : data) el += val;
    }
    void operator-=(const T val)
    {
        (*this) += -val;
    }
    void operator*=(const T val)
    {
        for (auto & el : data) el *= val;
    }
    void operator/=(const T val)
    {
        (*this) *= 1.0 / val;
    }

    Matrix operator+(const T val) const
    {
        Matrix<T> ans = (*this);
        ans += val;
        return ans;
    }
    Matrix operator-(const T val) const
    {
        Matrix<T> ans = (*this);
        ans -= val;
        return ans;
    }
    Matrix operator*(const T val) const
    {
        Matrix<T> ans = (*this);
        ans *= val;
        return ans;
    }
    Matrix operator/(const T val) const
    {
        Matrix<T> ans = (*this);
        ans /= val;
        return ans;
    }

    Matrix<T> operator-() const
    {
        Matrix<T> ans = (*this);
        ans.minus_inplace();
        return ans;
    }

    bool operator==(const Matrix & other) const
    {
        assert(nR() == other.nR());
        assert(nC() == other.nC());
        bool match = true;
        for (int r = 0; r < nR(); ++r) {
            for (int c = 0; c < nC(); ++c) {
                match &= (*this)(r, c) == other(r, c);
            }
        }
        return match;
    }

    bool operator!=(const Matrix & other) const
    {
        return !((*this) == other);
    }

    void minus_inplace()
    {
        for (T & el : data) el = -el;
    }

    void exp_inplace()
    {
        for (T & el : data) el = exp(el);
    }

    Matrix<T> solve(const Matrix<T> & other)
    {

    }

private:
    struct literal_assign_helper {

        explicit literal_assign_helper(Matrix * m_): m(m_)
        {
            next();
        }
        ~literal_assign_helper()
        {
        }

        const literal_assign_helper& operator, (
            const T& val
        ) const
        {
            (*m)(r, c) = val;
            next();
            return *this;
        }

    private:
        void next() const
        {
            ++c;
            if (c == m->nC()) {
                c = 0;
                ++r;
            }
            if (r == m->nR() && c == m->nC()) {
                done = true;
            }
        }
        friend class Matrix;
        mutable bool done = false;
        mutable long r = 0;
        mutable long c = 0;
        Matrix * m;

    };

public:

    Matrix& operator = (
        const literal_assign_helper& val
    )
    {
        *this = *val.m;
        return *this;
    }

    const literal_assign_helper operator = (
        const T& val
    )
    {
        for (T& vec_val : data) {
            vec_val = val;
        }
        return literal_assign_helper(this);
    }

    void apply_lambda();

    friend Matrix exp<T>(const Matrix & m);
};
template <typename T>
std::ostream& operator<<(
    std::ostream& out, Matrix<T> M
)
{
    for (long row = 0; row < M.nR(); ++row) {
        for (long col = 0; col < M.nC(); ++col) {
            out << M(row, col) << " ";
        }
        out << std::endl;
    }
    return out;
}

template <typename T>
Matrix<T> exp(const Matrix<T> & m)
{
    Matrix<T> ans = m;
    ans.exp_inplace();
    return ans;
}

template <typename T, typename N>
Matrix<T> operator*(N val, const Matrix<T> & m)
{
    return m * val;
}



template <typename T>
class LU {
private:
    std::vector<T> LU;
    std::vector<int> pivots;
public:
    Matrix<T> inv() const;
    Matrix<T> solve(const Matrix<T> & other) const;
};

#endif //KRONMAT_Matrix_H
