//
// Created by Thomas Nickson on 11/06/15.
//

#ifndef KRONMAT_STATIC_FUNCTIONS_H
#define KRONMAT_STATIC_FUNCTIONS_H

#include <vector>
#include "matrix.h"
#include <deque>
#include "kronvecstack.h"
#include "cholesky.h"
#include <thrust/functional.h>


template <typename T>
std::vector<Matrix<T>> kronmat_dot_kronmat (
                        const std::vector<Matrix<T>> & M1, const std::vector<Matrix<T>> & M2
                    )
{
    std::vector<Matrix<T>> ans;
    assert(M1.size() == M2.size());
    for (int i = 0; i < M1.size(); ++i) {
        ans.push_back(M1[i] * M2[i]);
    }
    return ans;
}

template <typename Storage>
Matrix<Storage> kronmat_dot_fullvec ( const std::vector<Matrix<Storage>> & K, const Matrix<Storage> & V )
{
    assert(V.nC() == 1);
    size_t nmats = K.size();
    Matrix<Storage> x = V;

    for (int n = nmats - 1; n > -1; --n) {
        size_t thisC = K[n].nC();
        size_t xSize = x.nR() * x.nC();
        x = x.reshape(thisC, xSize / thisC);
        x = (K[n] * x).transpose().reshape(xSize,1);
    }
    return x;
}
template <typename Storage>
Matrix<Storage> kronmat_solve_fullvec ( const std::vector<Matrix<Storage>> & K, const Matrix<Storage> & V )
{
    assert(V.nC() == 1);
    size_t nmats = K.size();
    Matrix<Storage> x = V;

    for (int n = nmats - 1; n > -1; --n) {
        size_t thisC = K[n].nC();
        size_t xSize = x.nR() * x.nC();
        x = x.reshape(thisC, xSize / thisC);
        Cholesky<Storage> Ln(K[n]);
        x = (Ln.solve(x)).transpose().reshape(xSize, 1);
    }
    return x;
}
template <typename Storage>
Matrix<Storage> kron_full(
    const std::vector<Matrix<Storage>> & K
)
{
    size_t nr = 1;
    size_t nc = 1;
    size_t nmats = K.size();
    for (int i = 0; i < nmats; ++i) {
        nr *= K[i].nR();
        nc *= K[i].nC();
    }
    Matrix<Storage> full_mat(nr, nc);
    full_mat = 1;

    size_t row_acc = 1;
    size_t col_acc = 1;
    std::deque<size_t> rowstrides;
    std::deque<size_t> colstrides;

    rowstrides.push_front(1);
    colstrides.push_front(1);
    //for (long n = nmats - 1; n >= 0; --n) {
    for (size_t offset = 1; offset <= nmats; ++offset) {
        size_t n = nmats - offset;
        row_acc *= K[n].nR();
        col_acc *= K[n].nC();
        rowstrides.push_front(row_acc);
        colstrides.push_front(col_acc);
    }
    Storage & target = full_mat.getMutableData();
    thrust::fill(target.begin(), target.end(), 1);
    thrust::counting_iterator<size_t> idx_start(0);
    for (size_t d = 0; d < nmats; ++d) 
    {
        kronecker_index idxer{row_acc, col_acc, K[d].nR(), K[d].nC(), rowstrides[d], rowstrides[d+1], colstrides[d], colstrides[d+1]};
        auto submat_idx = thrust::make_transform_iterator(idx_start, idxer);
        auto submat_shuffle = thrust::make_permutation_iterator(K[d].getConstData().begin(), submat_idx);
        thrust::transform(target.begin(), target.end(), submat_shuffle, target.begin(), thrust::multiplies<typename Storage::value_type>());
    }
    return full_mat;
}
/*
template <typename T>
Matrix<T> kvs_full( const std::vector<Matrix<T>> & KVS,
                    size_t start, size_t end)
{
    size_t nr = end - start;
    size_t nc = 1;
    size_t nmats = KVS.size();
    for (int i = 0; i < nmats; ++i) {
        nc *= KVS[i].nC();
    }
    Matrix<T> full_mat(nr, nc);
    full_mat = 1;

    size_t col_acc = 1;
    std::deque<size_t> colstrides;

    colstrides.push_front(1);
    for (long n = nmats - 1; n >= 0; --n) {
        col_acc *= KVS[n].nC();
        colstrides.push_front(col_acc);
    }

    for (size_t rout = start; rout < end; ++rout) {
        for (size_t cout = 0; cout < nc; ++cout) {
            for (size_t d = 0; d < nmats; ++d) {
                size_t colm = (cout % colstrides[d]) / colstrides[d + 1];
                full_mat(rout - start, cout) *= KVS[d](rout, colm);
            }
        }
    }
    return full_mat;
}
*/
template <typename T>
Matrix<T> kvs_full(const std::vector<Matrix<T>> & KVS)
{
    return kvs_full(KVS, 0, KVS[0].nR());
}

template <typename T>
Matrix<T> vec_dot_kvs(const Matrix<T> & m, const KroneckerVectorStack<T> & kvs)
{
    assert(m.nR() == 1 && "Must be row vector");
    std::vector<size_t> shapes_mut;
    size_t full_width = 1;

    for (const auto & m : kvs.sub_matrices) shapes_mut.push_back(m.nC());
    for (const auto & m : kvs.sub_matrices) full_width *= m.nC();
    const std::vector<size_t> shapes = shapes_mut;
    size_t kvs_height = m.nC();
    size_t kvs_dim = kvs.sub_matrices.size();

    size_t final_r = 1;
    size_t final_c = full_width;


    std::vector<T> kroned_vecs;
    kroned_vecs.resize(full_width, 0);


    for (int n = 0; n < kvs_height; ++n) {
        T alpha = m(0, n);
        size_t current_width = shapes[0];
        std::vector<T> kronme(
            kvs.sub_matrices[0].begindata() + n * shapes[0],
            kvs.sub_matrices[0].begindata() + (n + 1) * shapes[0]
        );
        for (auto & val : kronme) val *= alpha;
        for (int d = 1; d < kvs_dim - 1; ++d) {
            std::vector<T> data(
                kvs.sub_matrices[d].begindata() + n * shapes[d],
                kvs.sub_matrices[d].begindata() + (n + 1) * shapes[d]
            );
            std::vector<T> destination;
            destination.resize(current_width * shapes[d], 0);

            blas_ger(current_width, shapes[d], 1.0,
                     kronme.data(), 1,
                     data.data(), 1,
                     destination.data(), shapes[d]);
            alpha = 1.0;
            kronme = destination;
            current_width *= shapes[d];
        }
        std::vector<T> data (
            kvs.sub_matrices[kvs_dim - 1].begindata() + n * shapes[kvs_dim - 1],
            kvs.sub_matrices[kvs_dim - 1].begindata() + (n  + 1) * shapes[kvs_dim - 1]
        );
        blas_ger(current_width, shapes[kvs_dim - 1], 1.0,
                 kronme.data(), 1,
                 data.data(), 1,
                 kroned_vecs.data(), shapes[kvs_dim - 1]);
    }
    return Matrix<T>(kroned_vecs, final_r, final_c);
}

/*
template <typename T>
std::vector<T> kvs_dot_vec(const basic_kvs<T> & kvs, const std::vector<T> & vec) {
    // Size of initial reshape of vec
    size_t init_width = 1;
    size_t last_idx = kvs.dim - 1;
    for (int i = 0; i < last_idx; ++i) {
        init_width *= kvs.shapes[i];
    }

    std::vector<T> alpha;
    alpha.resize(kvs.height * init_width);

    //    We can use BLAS for the first part. Gives us a N x \prod_{i=1}^{D-1} dim_i matrix.
    blas_gemm(CblasRowMajor, CblasNoTrans, CblasTrans,
              kvs.height, init_width, kvs.shapes[last_idx],
              1.0,
              kvs.data[last_idx].data(), kvs.shapes[last_idx],
              vec.data(), kvs.shapes[last_idx],
              1.0,
              alpha.data(), init_width);

    std::vector<float> beta;
    size_t next_width = init_width / kvs.shapes[last_idx];
    size_t next_idx = last_idx - 1;
    beta.resize(kvs.height * next_width);

    size_t N = kvs.height;
    std::vector<size_t> shapes = kvs.shapes;
    size_t full_width = 1;
    for (auto width : shapes) full_width *= width;
    size_t current_width = full_width / shapes.back();


    for (size_t d = kvs.dim - 2; d >= 0; --d)
    {
        current_width = current_width / shapes[d];
        std::vector<T> tmp;
        int this_width = shapes[d];

        for (size_t n = 0; n < kvs.height; ++n)
        {
            for (size_t idx = 0; idx < current_width; ++idx)
            {
                T val = blas_dot(this_width,
                                 &kvs.data[d][n * this_width], 1,
                                 &alpha[n * current_width * this_width + idx * this_width ], 1);
                tmp.push_back(val);
            }
        }
        alpha = tmp;
    }
    return alpha;
}
*/

#endif //KRONMAT_STATIC_FUNCTIONS_H
