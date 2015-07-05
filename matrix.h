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
#include "google_float.h"


    template <typename T>
    class matrix;

    template <typename T>
    class kronecker_matrix;

    template <typename T>
    matrix<T> exp(const matrix<T> &);

    template<typename T>
    class matrix {
    public:
        matrix() : data{}, nr{}, nc{}, r_stride{}, r0{}, c0{}, trans{} {};
        matrix(std::vector<T> data_, long r_, long c_,
               long r_stride_, long r0_, long c0_, bool trans_) :
                data{data_}, nr{r_}, nc{c_}, r_stride{r_stride_}, r0{r0_}, c0{c0_}, trans{trans_} {};
        matrix(long r_, long c_) : matrix()
        {
            set_size(r_, c_);
        };
        void set_size(long r, long c)
        {
            nr = r;
            nc = c;
            r_stride = c;
            data.resize(nr * nc, 0);
        }
    private:
        std::vector<T> data;
        long nr;
        long nc;
        long r_stride;
        long r0;
        long c0;
        CBLAS_ORDER order = CblasRowMajor;
        bool trans;
        const long vecidx(long r, long c) const
        {
            if (trans)
            {
                return (c0 + c) * r_stride + (r0 + r);
            } else {
                return (r0 + r) * r_stride + (c0 + c);
            }
        }
        const long offset() const {
            return vecidx(0,0);
        }

    public:
        long stride() const {
            return r_stride;
        }
        long use_count() const {
            return data.use_count();
        }
        long nR() const {
            long ans = trans ? nc : nr;
            return ans;
        }
        long nC() const {
            long ans = trans ? nr : nc;
            return ans;
        }

        bool isTrans() const {
            return trans;
        }


        const T * dataptr() const {
            return &(data.data()[offset()]);
        }

        T * mutable_dataptr() {
            return &(data.data()[offset()]);
        }

        matrix<T> transpose() const {
            return matrix(data, nr, nc, r_stride, r0, c0, !trans);
        }

        T& operator()(long ridx, long cidx) {
            assert(ridx<nR());
            assert(cidx<nC());
            return data[vecidx(ridx, cidx)];
        }
        T operator()(long ridx, long cidx) const {
            assert(ridx<nR());
            assert(cidx<nC());
            return data[vecidx(ridx, cidx)];
        }

        long getOffset() {
            return r_stride * r0 + c0;
        }
        matrix<T> operator*(const kronecker_matrix<T> &other);


        matrix<T> operator*(const matrix<T> &other) const {
            CBLAS_TRANSPOSE meTrans;
            CBLAS_TRANSPOSE otherTrans;
            int M, N, K, lda, ldb, ldc;

            assert(nC() == other.nR());

            M = nR();  //results matrices rows
            N = other.nC(); //results matrices cols
            K = nC(); // this matrices columns == other matrices rows

            ldc = other.nC();

            // check if this matrix is transposed:
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

            return matrix<T>(new_data, nR(), other.nC(), other.nC(), 0, 0, false);
        }

        matrix<T> Tdot(const matrix<T> &other)  {
            //TODO: This is ugly!
            // Avoid copying for switching the transpose switch on for this version. Ideally want to be able to make
            // read only views into a matrix?
            trans = !trans;
            auto ans = (*this) * other;
            trans = !trans;
            return ans;
        }

        matrix<T> operator+(const matrix<T> &other) const {
            assert(nR() == other.nR());
            assert(nC() == other.nC());

            std::vector<T> new_data;
            for (long r = 0; r < nR(); ++r)
            {
                for (long c = 0; c < nC(); ++c)
                {
                    new_data.push_back((*this)(r, c) + other(r, c));
                }
            }
            return matrix<T>(new_data, nR(), other.nC(), other.nC(), 0, 0, false);
        }

        matrix<T> hadamard(const matrix<T> &other) const {
            assert(nR() == other.nR());
            assert(nC() == other.nC());

            std::vector<T> new_data;
            for (long r = 0; r < nR(); ++r)
            {
                for (long c = 0; c < nC(); ++c)
                {
                    new_data.push_back((*this)(r, c) * other(r, c));
                }
            }
            return matrix<T>(new_data, nR(), other.nC(), other.nC(), 0, 0, false);
        }
        matrix<T> elemwise_div(const matrix<T> &other) const {
            assert(nR() == other.nR());
            assert(nC() == other.nC());

            std::vector<T> new_data;
            for (long r = 0; r < nR(); ++r) {
                for (long c = 0; c < nC(); ++c) {
                    new_data.push_back((*this)(r, c) / other(r, c));
                }
            }
            return matrix<T>(new_data, nR(), other.nC(), other.nC(), 0, 0, false);
        }


        void operator+=(const matrix<T> &other) {
            assert(nR() == other.nR());
            assert(nC() == other.nC());
            for (long r = 0; r < nR(); ++r) {
                for (long c = 0; c < nC(); ++c) {
                    (*this)(r, c) += other(r, c);
                }
            }
        }

        void operator+=(const T val) {
            for (auto & el : data) el += val;
        }
        void operator-=(const T val) {
            for (auto & el : data) el -= val;
        }
        void operator*=(const T val) {
            for (auto & el : data) el *= val;
        }
        void operator/=(const T val) {
            for (auto & el : data) el /= val;
        }

        matrix operator+(const T val) const
        {
            matrix<T> ans = (*this);
            ans += val;
            return ans;
        }
        matrix operator-(const T val) const
        {
            matrix<T> ans = (*this);
            ans -= val;
            return ans;
        }
        matrix operator*(const T val) const
        {
            matrix<T> ans = (*this);
            ans *= val;
            return ans;
        }
        matrix operator/(const T val) const
        {
            matrix<T> ans = (*this);
            ans /= val;
            return ans;
        }

        matrix<T> operator-() const
        {
            matrix<T> ans = (*this);
            for (T & el : (*ans.data)) el = -el;
            return ans;
        }

        bool operator==(const matrix & other) const {
            assert(nR() == other.nR());
            assert(nC() == other.nC());
            bool match = true;
            for (int r = 0; r < nR(); ++r) {
                for (int c = 0; c < nC(); ++c) {
                    match &= FloatingPoint<T>((*this)(r,c)).AlmostEquals(FloatingPoint<T>(other(r,c)));
                }
            }
            return match;
        }

        void minus_inplace()
        {
            for (T & el : (*data)) el = -el;
        }

        void exp_inplace()
        {
            for (T & el : (*data)) el = exp(el);
        }

        private:
        struct literal_assign_helper
        {

            explicit literal_assign_helper(matrix * m_): m(m_)
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
                (*m)(r,c) = val;
                next();
                return *this;
            }

        private:
            void next() const
            {
                ++c;
                if (c == m->nC())
                {
                    c = 0;
                    ++r;
                }
                if (r == m->nR() && c == m->nC()) {
                    done = true;
                }
            }
            friend class matrix;
            mutable bool done = false;
            mutable long r = 0;
            mutable long c = 0;
            matrix * m;

        };

    public:

        matrix& operator = (
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

        friend matrix exp<T>(const matrix & m);
    };
    template <typename T>
    std::ostream& operator<<(
            std::ostream& out, matrix<T> M
    )
    {
        for (long row = 0; row < M.nR(); ++row) {
            for (long col = 0; col < M.nC(); ++col) {
                out << M(row,col) << " ";
            }
            out << std::endl;
        }
        return out;
    }

template <typename T>
matrix<T> exp(const matrix<T> & m)
{
    matrix<T> ans = m;
    for (T & el : (*ans.data)) el = exp(el);
    return ans;
}

template <typename T>
class kronecker_matrix {
private:
    std::shared_ptr<std::vector<matrix<T>>> data;
};
#endif //KRONMAT_MATRIX_H
