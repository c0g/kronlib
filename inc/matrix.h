
// Created by Thomas Nickson on 12/06/15.
//

#ifndef KRONMAT_MATRIX_H
#define KRONMAT_MATRIX_H

#include <vector>
#include <iostream>
#include <assert.h>
#include "backend.h"
#include "blas/blas.h"
#include "functors.h"
#include <cmath>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/fill.h>
#include <thrust/inner_product.h>
#include <thrust/functional.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/iterator/transform_iterator.h>
#include <thrust/gather.h>
#include <cublas_v2.h>

namespace kronlib {
template<typename Backend>
class Matrix;

template<typename Backend>
class KroneckerMatrix;

/*
template<typename Backend, typename T>
class KroneckerVectorStack;


template<typename Backend, typename T>
Matrix<Backend, T> exp(const Matrix<T, Storage> &);
*/

template<typename Backend>
class Matrix  {
	using T = typename Backend::Storage::value_type;
    friend KroneckerMatrix<Backend>;
    //friend KroneckerVectorStack<Backend, T>;
public:
    using Storage = typename Backend::Storage;
    Matrix(std::shared_ptr<Backend> context_, Storage data_, size_t r_,size_t c_) : context{context_}, data{data_}, nr{r_}, nc{c_} {}
    Matrix(size_t r_, size_t c_) : Matrix(std::make_shared<Backend>(), Storage(r_ * c_), r_, c_) {}
    template <typename OtherBackend>
        Matrix(const Matrix<OtherBackend> & other) : Matrix(std::make_shared<Backend>(), other.getConstData(), other.nR(), other.nC()) {}

private:
    std::shared_ptr<Backend> context;
    Storage data;
    size_t nr;
    size_t nc;
    size_t r0 = 0;
    size_t c0 = 0;
    size_t vecidx(size_t r, size_t c) const
    {
        return r + nr * c;
    }
    size_t offset() const
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

    std::shared_ptr<Backend> getContext() const { return context; };
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
    size_t nR() const
    {
        return nr;
    }
    size_t nC() const
    {
        return nc;
    }

    T trace() const 
    {
        // Diagonal is at 0, nr + 1, 2 * nr + 2, etc. Sum those items.
        // Only care about leading square matrix:
        size_t n = std::min(nr, nc);
        thrust::counting_iterator<size_t> idxs(0);
        auto map = thrust::make_transform_iterator(idxs, mult_by<size_t>(nr + 1));
        return thrust::reduce(thrust::make_permutation_iterator(data.begin(), map),
                                         thrust::make_permutation_iterator(data.begin(), map + n));
    }

    Matrix<Backend> transpose() const
    {
        Matrix ans(nc, nr);
        thrust::counting_iterator<size_t> indices(0);
        auto transposed = thrust::make_transform_iterator(indices, transpose_index(nr, nc));
        Storage & ansData = ans.getMutableData();
        thrust::gather(context->exec(), transposed, transposed + ansData.size(), data.begin(), ansData.begin());
        return ans;
    }
    void reshape_inplace(size_t newr, size_t newc) {
        assert(nR() * nC() == newr * newc);
        nr = newr;
        nc = newc;
    }
    Matrix<Backend> reshape(size_t newr, size_t newc) const
    {
        assert(nR() * nC() == newr * newc);
        Storage newdata = data;
        return Matrix(context, newdata, newr, newc);
    }
    void setat(size_t ridx, size_t cidx, T val) {
        assert(ridx < nR());
        assert(cidx < nC());
        data[vecidx(ridx, cidx)] = val;     
    }
    T operator()(size_t ridx, size_t cidx) const
    {
        assert(ridx < nR());
        assert(cidx < nC());
        return data[vecidx(ridx, cidx)];
    }

    Matrix<Backend> operator*(const Matrix<Backend> &other) const
    {
        BlasTranspose meTrans = N;
        BlasTranspose otherTrans = N;
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

        blas_gemm(context, COL, meTrans, otherTrans, M, N, K, context->one(),
                  data, lda, other.getConstData(), ldb,
                  context->zero(), new_data, ldc);

        return Matrix(context, new_data, nR(), other.nC());
    }

    /* Removed until I can get clearer idea of how to make mutable transpose work
    Matrix<Backend> Tdot(const Matrix<Storage> &other)
    {
        mutable_transpose();
        auto ans = (*this) * other;
        mutable_transpose();
        return ans;
    }
    */
/*
    Matrix<Backend> operator*(const KroneckerVectorStack<T> kvs) {
        assert(nR() == 1 && "Only works when matrix is a row vector!");
        return vec_dot_kvs((*this), kvs);
    }
*/

    //Element wise operations with other matrices
    bool operator==(const Matrix & other) const
    {
        assert(nR() == other.nR());
        assert(nC() == other.nC());
        bool init = true;
        auto ans = thrust::inner_product(context->exec(), data.begin(), data.end(), other.getConstData().begin(), init, and_reduce(), thrust::equal_to<T>());
        return ans;
    }
    bool tol_eq(const Matrix & other, T tol) const
    {
        assert(nR() == other.nR());
        assert(nC() == other.nC());
        bool init = true;
        auto ans = thrust::inner_product(context->exec(), data.begin(), data.end(), other.getConstData().begin(), init, and_reduce(), tol_equal<T>(tol));
        return ans;
    }


    bool operator!=(const Matrix & other) const
    {
        return !((*this) == other);
    }
    Matrix operator+(const Matrix<Backend> &other) const
    {
        assert(nR() == other.nR());
        assert(nC() == other.nC());
        Matrix ans = *this;
        ans += other;
        return ans;
    }

    Matrix operator-(const Matrix<Backend> &other) const
    {
        assert(nR() == other.nR());
        assert(nC() == other.nC());
        Matrix ans = *this;
        ans -= other;
        return ans;
    }

    Matrix elemwise_mult(const Matrix<Backend> &other) const
    {
        assert(nR() == other.nR());
        assert(nC() == other.nC());
        Matrix ans = *this;
        ans.elemwise_mult_inplace(other);
        return ans;

    }
    Matrix elemwise_div(const Matrix<Backend> &other) const
    {
        assert(nR() == other.nR());
        assert(nC() == other.nC());
        Matrix ans = *this;
        ans.elemwise_div_inplace(other);
        return ans;
    }
    void operator+=(const Matrix<Backend> &other)
    {
        assert(nR() == other.nR());
        assert(nC() == other.nC());
        thrust::transform(context->exec(), data.begin(), data.end(), other.getConstData().begin(), data.begin(), thrust::plus<T>());
    }
    void operator-=(const Matrix<Backend> &other)
    {
        assert(nR() == other.nR());
        assert(nC() == other.nC());
        thrust::transform(context->exec(), data.begin(), data.end(), other.getConstData().begin(), data.begin(), thrust::minus<T>());
    }
    void elemwise_div_inplace(const Matrix<Backend> &other)
    {
        assert(nR() == other.nR());
        assert(nC() == other.nC());
        thrust::transform(context->exec(), data.begin(), data.end(), other.getConstData().begin(), data.begin(), thrust::divides<T>());
    }
    void elemwise_mult_inplace(const Matrix<Backend> &other)
    {
        assert(nR() == other.nR());
        assert(nC() == other.nC());
        thrust::transform(context->exec(), data.begin(), data.end(), other.getConstData().begin(), data.begin(), thrust::multiplies<T>());
    }
    
    // Operations with scalars
    Matrix operator+(T val) const
    {
        Matrix<Backend> ans = (*this);
        ans += val;
        return ans;
    }
    Matrix operator-(T val) const
    {
        Matrix<Backend> ans = (*this);
        ans -= val;
        return ans;
    }
    Matrix operator*(T val) const
    {
        Matrix<Backend> ans = (*this);
        ans *= val;
        return ans;
    }
    Matrix operator/(T val) const
    {
        Matrix<Backend> ans = (*this);
        ans /= val;
        return ans;
    }
    void operator+=(T val)
    {
        thrust::transform(context->exec(), data.begin(), data.end(), data.begin(), add_by<T>(val));
    }
    void operator-=(T val)
    {
        thrust::transform(context->exec(), data.begin(), data.end(), data.begin(), minus_by<T>(val));
    }
    void operator*=(T val)
    {
        thrust::transform(context->exec(), data.begin(), data.end(), data.begin(), mult_by<T>(val));
    }
    void operator/=(T val)
    {
        thrust::transform(context->exec(), data.begin(), data.end(), data.begin(), div_by<T>(val));
    }

    Matrix operator-() const
    {
        Matrix<Backend> ans = (*this);
        ans.negate_inplace();
        return ans;
    }
    void negate_inplace()
    {
        auto exec = context->exec();
        thrust::transform(exec, data.begin(), data.end(), data.begin(), thrust::negate<T>()); 
    }
    void exp_inplace()
    {
        thrust::transform(context->exec(), data.begin(), data.end(), data.begin(), exponentiate<T>()); 
    }

    void apply_lambda();
    // defined in kronecker Matrix;
    // Matrix<Backend> Matrix<Storage>::operator*(const KroneckerMatrix<Storage> &other);

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
        )  const 
        {
            m->setat(r, c, val);
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
        mutable size_t r = 0;
        mutable size_t c = 0;
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
        thrust::fill(data.begin(), data.end(), val);
        return literal_assign_helper(this);
    }
};

template <typename Backend>
std::ostream& operator<<(
    std::ostream& out, Matrix<Backend> M
)
{
    for (size_t row = 0; row < M.nR(); ++row) {
        for (size_t col = 0; col < M.nC(); ++col) {
            out << M(row, col) << " ";
        }
        out << std::endl;
    }
    return out;
}

template <typename Backend>
Matrix<Backend> exp(const Matrix<Backend> & m)
{
    Matrix<Backend> ans = m;
    ans.exp_inplace();
    return ans;
}

template <typename Backend, typename N>
Matrix<Backend> operator*(N val, const Matrix<typename Backend::Storage> & m)
{
    return m * val;
}

}
#endif //KRONMAT_Matrix_H
