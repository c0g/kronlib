
// Created by Thomas Nickson on 12/06/15.
//

#ifndef KRONMAT_MATRIX_H
#define KRONMAT_MATRIX_H

#include <vector>
#include <iostream>
#include <assert.h>
#include "blas_wrap.h"
#include <cmath>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/fill.h>

/*
template<typename Storage, typename T>
class Matrix;

template<typename Storage, typename T>
class KroneckerMatrix;

template<typename Storage, typename T>
class KroneckerVectorStack;

template<typename Storage, typename T>
class Cholesky;

template<typename Storage, typename T>
Matrix<Storage, T> exp(const Matrix<T, Storage> &);
*/
template <typename T>
using device = thrust::device_vector<T>;

template <typename T>
using host = thrust::host_vector<T>;

template<typename Storage>
class Matrix  {
	using T = typename Storage::value_type;
    //friend KroneckerMatrix<Storage, T>;
    //friend KroneckerVectorStack<Storage, T>;
    //friend Cholesky<Storage, T>;
public:
    Matrix() : data{}, nr{0}, nc{0}, r0{0}, c0{0} {};
    Matrix(Storage data_, long r_, long c_) :
        data{data_}, nr{r_}, nc{c_} {
	};
    Matrix(long r_, long c_) : Matrix()
    {
        set_size(r_, c_);
    };
    void set_size(long r, long c)
    {
        nr = r;
    	nc = c;
        data.resize(nr * nc, 0);

    }
private:

    Storage data;
    long nr;
    long nc;
    long r0 = 0;
    long c0 = 0;
    long vecidx(long r, long c) const
    {
        return r + nr * c;
    }
    long offset() const
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
    const Storage & getConstData() const
    {
        return data;
    }
    Storage & getMutableData()
    {
        return data;
    }
    long nR() const
    {
        return nr;
    }
    long nC() const
    {
        return nc;
    }

    void operator=(const T& val) {
        thrust::fill(data.begin(), data.end(), val);
    }

    T trace() const 
    {
        assert(nR() == nC());
        T trace = 0;
        for (int idx = 0; idx < nR(); ++idx) {
            trace += (*this)(idx, idx);
        }
        return trace;
    }

    Matrix<Storage> transpose() const
    {
        Storage newdata;
        for (int r = 0; r < nR(); ++r) {
            for (int c = 0; c < nC(); ++c) {
                newdata.push_back((*this)(r, c));
            }
        }
        return Matrix(newdata, nC(), nR());
    }

    Matrix<Storage> reshape(long newr, long newc) const
    {
        assert(nR() * nC() == newr * newc);
        Storage newdata;
        for (int c = 0; c < nC(); ++c) {
            for (int r = 0; r < nR(); ++r) {
                newdata.push_back((*this)(r, c));
            }
        }
        return Matrix(newdata, newr, newc);
    }

    T operator()(long ridx, long cidx) const
    {
        assert(ridx < nR());
        assert(cidx < nC());
        return data[vecidx(ridx, cidx)];
    }



    Matrix<Storage> operator*(const Matrix<Storage> &other) const
    {
        CBLAS_TRANSPOSE meTrans = CblasNoTrans;
        CBLAS_TRANSPOSE otherTrans = CblasNoTrans;
        int M, N, K, lda, ldb, ldc;

        assert(nC() == other.nR());

        M = nR();  //results matrices rows
        N = other.nC(); //results matrices cols
        K = nC(); // this matrices columns == other matrices rows

        lda = nR(); // size of leading dim of this matrix
        ldb = other.nR(); // size of leading dim of other matrix
        ldc = nR(); // size of leading dim of answer matrix

        Storage new_data;
        new_data.resize(nR() * other.nC());


        blas_gemm(CblasColMajor, meTrans, otherTrans, M, N, K, 1.0,
                  dataptr(), lda, other.dataptr(), ldb,
                  0.0, new_data.data(), ldc);

        return Matrix(new_data, nR(), other.nC());
    }

    /* Removed until I can get clearer idea of how to make mutable transpose work
    Matrix<Storage> Tdot(const Matrix<Storage> &other)
    {
        mutable_transpose();
        auto ans = (*this) * other;
        mutable_transpose();
        return ans;
    }
    */
/*
    Matrix<Storage> operator*(const KroneckerVectorStack<T> kvs) {
        assert(nR() == 1 && "Only works when matrix is a row vector!");
        return vec_dot_kvs((*this), kvs);
    }
*/
    Matrix operator+(const Matrix<Storage> &other) const
    {
        assert(nR() == other.nR());
        assert(nC() == other.nC());

        Storage new_data;
        for (long c = 0; c < nC(); ++c) {
            for (long r = 0; r < nR(); ++r) {
                new_data.push_back((*this)(r, c) + other(r, c));
            }
        }
        return Matrix(new_data, nR(), other.nC());
    }

    Matrix operator-(const Matrix<Storage> &other) const
    {
        Matrix<Storage> ans = other;
        ans.minus_inplace();
        ans += (*this);
        return ans;
    }

    Matrix hadamard(const Matrix<Storage> &other) const
    {
        assert(nR() == other.nR());
        assert(nC() == other.nC());

        Storage new_data;
        for (long c = 0; c < nC(); ++c) {
            for (long r = 0; r < nR(); ++r) {
                new_data.push_back((*this)(r, c) * other(r, c));
            }
        }
        return Matrix(new_data, nR(), other.nC());
    }
    Matrix elemwise_div(const Matrix<Storage> &other) const
    {
        assert(nR() == other.nR());
        assert(nC() == other.nC());

        Storage new_data;
        for (long c = 0; c < nC(); ++c) {
            for (long r = 0; r < nR(); ++r) {
                new_data.push_back((*this)(r, c) / other(r, c));
            }
        }
        return Matrix(new_data, nR(), other.nC());
    }

    // defined in kronecker Matrix;
    // Matrix<Storage> Matrix<Storage>::operator*(const KroneckerMatrix<Storage> &other);


    void operator+=(const Matrix<Storage> &other)
    {
        assert(nR() == other.nR());
        assert(nC() == other.nC());
        for (long c = 0; c < nC(); ++c) {
            for (long r = 0; r < nR(); ++r) {
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
        Matrix<Storage> ans = (*this);
        ans += val;
        return ans;
    }
    Matrix operator-(const T val) const
    {
        Matrix<Storage> ans = (*this);
        ans -= val;
        return ans;
    }
    Matrix operator*(const T val) const
    {
        Matrix<Storage> ans = (*this);
        ans *= val;
        return ans;
    }
    Matrix operator/(const T val) const
    {
        Matrix<Storage> ans = (*this);
        ans /= val;
        return ans;
    }

    Matrix operator-() const
    {
        Matrix<Storage> ans = (*this);
        ans.minus_inplace();
        return ans;
    }

    bool operator==(const Matrix & other) const
    {
        assert(nR() == other.nR());
        assert(nC() == other.nC());
        bool match = true;
        for (long c = 0; c < nC(); ++c) {
            for (long r = 0; r < nR(); ++r) {
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

    Matrix<Storage> solve(const Matrix<Storage> & other)
    {
        return *this;
    }
    void apply_lambda();

    //friend Matrix exp<Storage>(const Matrix & m);
};

template <typename Storage>
std::ostream& operator<<(
    std::ostream& out, Matrix<Storage> M
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

template <typename Storage>
Matrix<Storage> exp(const Matrix<Storage> & m)
{
    Matrix<Storage> ans = m;
    ans.exp_inplace();
    return ans;
}

template <typename Storage, typename N>
Matrix<Storage> operator*(N val, const Matrix<Storage> & m)
{
    return m * val;
}

#endif //KRONMAT_Matrix_H
