//
// Created by Thomas Nickson on 12/06/15.
//

#ifndef KRONMAT_MATRIX_H
#define KRONMAT_MATRIX_H

#include <vector>
#include <assert.h>
#include "blas_wrap.h"

namespace tommat {
    template<typename T>
    class matrix {
    public:
        matrix(std::shared_ptr<std::vector<T>> data_, long r_, long c_,
               long r_stride_, long r0_, long c0_, bool trans_) :
                data{data_}, nr{r_}, nc{c_}, r_stride{r_stride_}, r0{r0_}, c0{c0_}, trans{trans_} { };
        matrix(long r_, long c_) : nr{r_}, nc{c_}, r_stride{c_}, r0{0}, c0{0}, trans{false}
        {
            data = std::make_shared<std::vector<T>>();
            data->resize(nR() * nC());
        };
    private:
        std::shared_ptr<std::vector<T>> data; // Stored as shared pointer so we can share data for transposes
        long nr;
        long nc;
        long r_stride;
        long r0;
        long c0;
        CBLAS_ORDER order = CblasRowMajor;
        bool trans;
        const long offset() const {
            return r0 * r_stride + c0;
        }

    public:
        long nR() const {
            return nr;
        }

        bool isTrans() const {
            return trans;
        }

        long nC() const {
            return nc;
        }
        const T * dataptr() const {
            return &(data->data()[offset()]);
        }

        matrix<T> transpose() {
            return matrix(data, nr, nc, r_stride, r0, c0, !trans);
        }

        T& operator()(long ridx, long cidx) {
            return (*data)[(r0 + ridx) * r_stride + (c0 + cidx)];
        }
        T operator()(long ridx, long cidx) const {
            return (*data)[(r0 + ridx) * r_stride + (c0 + cidx)];
        }

        long getOffset() {
            return r_stride * r0 + c0;
        }

        matrix<T> operator*(const matrix &other) {
            CBLAS_TRANSPOSE meTrans;
            CBLAS_TRANSPOSE otherTrans;



            = isTrans() ? CblasTrans : CblasNoTrans;
            = other.isTrans() ? CblasTrans : CblasNoTrans;
            assert(nC() == other.nR());

            std::shared_ptr<std::vector<T>> newdata;
            newdata->reserve(nR() * other.nC());

            blas_gemm(order, meTrans, otherTrans, nR(), nC(), other.nC(), 1.0,
                            dataptr(), 1, other.dataptr(), 1,
                            1.0, newdata->data(), 1);
            return matrix<T>(newdata, nR(), other.nC(), other.nC(), 0, 0, false);
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
            for (T& vec_val : (*data)) {
                vec_val = val;
            }
            return literal_assign_helper(this);
        }
    };
    template <typename T>
    std::ostream& operator<<(
            std::ostream& out, matrix<T> M
    )
    {
        for (int row = 0; row < M.nR(); ++row) {
            for (int col = 0; col < M.nC(); ++col) {
                out << M(row,col) << " ";
            }
            out << std::endl;
        }
        return out;
    }

}

#endif //KRONMAT_MATRIX_H
