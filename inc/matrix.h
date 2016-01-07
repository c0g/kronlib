// Created by Thomas Nickson on 12/06/15.
//

#ifndef KRONMAT_MATRIX_H
#define KRONMAT_MATRIX_H

#include <vector>
#include <functional>
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
    Matrix() : Matrix(std::make_shared<Backend>(), Storage{0}, 0, 0) {};
    Matrix(std::shared_ptr<Backend> context_, Storage data_, size_t r_,size_t c_) : context{context_}, data{data_}, nr{r_}, nc{c_} {}
    Matrix(size_t r_, size_t c_) : Matrix(std::make_shared<Backend>(), Storage(r_ * c_), r_, c_) {}
    Matrix(std::shared_ptr<Backend> context, size_t r_, size_t c_) : Matrix(context, Storage(r_ * c_), r_, c_) {}
    //Matrix(const Matrix & other) : Matrix(other.getContext(), other.getConstData(), other.nR(), other.nC()) {}
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
    void swapStorage(Storage & newData) { std::swap(data, newData); } 
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

    Matrix<Backend> operator*(const Matrix<Backend> & other) const {return dot(other); }
    Matrix<Backend> dot(const Matrix<Backend> & other) const { return dot(None, other, None); }
    Matrix<Backend> dot(BlasTranspose meTrans, const Matrix<Backend> &other, BlasTranspose otherTrans) const
    {
        int Mme, Nme, Mother, Nother, Mc, Nc, K;
        switch (meTrans) 
        {
            case None:
                Mme = nR();
                Nme = nC();
                break;
            case Trans:
            case Conj:
                Mme = nC();
                Nme = nR();
                break;
        }
        switch (otherTrans)
        {
            case None:
                Mother = other.nR();
                Nother = other.nC();
                break;
            case Trans:
            case Conj:
                Mother = other.nC();
                Nother = other.nR();
                break;
        }
        if (Mother != Nme)
        { 
            std::cout << "Matrices must align \n";
            std::cout << "You gave a matrix of " << Mme << " x " << Nme << " and " << Mother << " x " << Nother << std::endl;
            std::cout << "This is an error" << std::endl;
            exit(-1);
        }

        Mc = Mme;
        Nc = Nother;
        K = Nme;
        return dot(Mc, K, Nc, meTrans, other, otherTrans);
    }

    Matrix<Backend> dot(size_t M, size_t K, size_t N, BlasTranspose meTrans, const Matrix<Backend> &other, BlasTranspose otherTrans) const
    {
        int lda, ldb, ldc;

        switch (meTrans)
        {
            case None:
                lda = M; // Matrix is col-major and not trans .: leading dim is M
                break;
            case Trans:
            case Conj:
                lda = K; // Matrix is col-major and trans .: leading dim is K
                break;
        }
        switch (otherTrans) 
        {
            case None:
                ldb = K; // Matrix is col-major and not trans .: leading dim is K
                break;
            case Trans:
            case Conj:
                ldb = N; // Matrix is col-major and trans .: leading dim is N
                break;
        }
        ldc = M; // size of leading dim of answer matrix - as col-major in memory
        Matrix<Backend> ans(getContext(), M, N);
        dot_into(M, K, N, meTrans, lda, other, ldb, otherTrans, ans, ldc);
        return ans;
    }
    void dot_into(size_t M, size_t K, size_t N, BlasTranspose meTrans, int lda, const Matrix<Backend> &other, int ldb, BlasTranspose otherTrans, Matrix<Backend> &dest, int ldc) const {
        blas_gemm(context, COL, meTrans, otherTrans, M, N, K, context->one(),
                  data, lda, other.getConstData(), ldb,
                  context->zero(), dest.getMutableData(), ldc);
    }

    struct nop_indexer : std::unary_function<size_t, size_t> {
        __host__ __device__
        size_t operator()(size_t idx) { return idx; }
    };
    struct row_indexer : std::unary_function<size_t, size_t> {
        size_t leading;
        __host__ __device__
            row_indexer( size_t leading ) : leading{leading} {}
        __host__ __device__
            size_t operator()( size_t idx ) { return idx / leading; }
    };
    void row_wise_tiled_add_inplace( const Matrix<Backend> & row )
    {
        row_indexer odx{nR()};
        nop_indexer mdx;
        thrust::plus<T> op;
        permuted_op_inplace(mdx, odx, op, row);
    }

    struct col_indexer : std::unary_function<size_t, size_t> {
        size_t leading;
        __host__ __device__
            col_indexer( size_t leading ) : leading{leading} {}
        __host__ __device__
            size_t operator()( size_t idx ) { return idx % leading; }
    };
    void col_wise_tiled_add_inplace( const Matrix<Backend> & col )
    {
        col_indexer odx{nR()};
        nop_indexer mdx;
        thrust::plus<T> op;
        permuted_op_inplace(mdx, odx, op, col);
    }

    void row_wise_tiled_mult_inplace( const Matrix<Backend> & row )
    {
        row_indexer odx{nR()};
        nop_indexer mdx;
        thrust::multiplies<T> op;
        permuted_op_inplace(mdx, odx, op, row);
    }

    void col_wise_tiled_mult_inplace( const Matrix<Backend> & col )
    {
        col_indexer odx{nR()};
        nop_indexer mdx;
        thrust::multiplies<T> op;
        permuted_op_inplace(mdx, odx, op, col);
    }

    template <typename MeIndex, typename OtherIndex, typename Op>
    void permuted_op_inplace(MeIndex mdx, OtherIndex odx, Op op, const Matrix<Backend> & other)
    {
        permuted_op_into(mdx, odx, mdx, op, other, (*this));
    }
    template <typename MeIndex, typename OtherIndex, typename TargetIndex, typename Op>
    void permuted_op_into(MeIndex mdx, OtherIndex odx, TargetIndex tdx, Op op, const Matrix<Backend> & other, Matrix<Backend> & target)
    {
        auto counter = thrust::make_counting_iterator(0);
        auto mdxer = thrust::make_transform_iterator(counter, mdx);
        auto odxer = thrust::make_transform_iterator(counter, odx);
        auto tdxer = thrust::make_transform_iterator(counter, tdx);
        auto mbegin = thrust::make_permutation_iterator(data.begin(), mdxer);
        auto mend = mbegin + data.size();
        auto otherbegin = thrust::make_permutation_iterator(other.getConstData().begin(), odxer);
        auto tbegin = thrust::make_permutation_iterator(target.getMutableData().begin(), tdxer);
        thrust::transform(mbegin, mend, otherbegin, tbegin, op);
    }
    template <typename T>
    struct square : std::unary_function<T,T> {
        __host__ __device__
            T operator()(T val) { return val * val; }
    };
    Matrix<Backend> sumsq_cols() const 
    {
        Matrix<Backend> squared = elemwise_mult(*this);
        Matrix<Backend> ones{nC(), 1};
        ones = 1;
        return squared.dot(ones);
    }

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

} // kronlib namespace
#endif //KRONMAT_Matrix_H
