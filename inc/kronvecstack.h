//
// Created by Thomas Nickson on 06/07/2015.
//

#ifndef KRONMAT_KRONVECSTACK_H
#define KRONMAT_KRONVECSTACK_H

#include "matrix.h"
#include <functional>
#include <numeric>
namespace kronlib {
template <typename Matrix>
class KroneckerVectorStack {

private:
    std::vector<Matrix> sub_matrices;
public:
    KroneckerVectorStack(const std::vector<Matrix> & sub_matrices) : sub_matrices{sub_matrices} {}
    KroneckerVectorStack() {}

public:
    const std::vector<Matrix> & getSubMatrices() const
    {
        return sub_matrices;
    }

    void add_dimension(const Matrix & mat) { push_back(mat); }
    void push_back(const Matrix & mat)
    {
        sub_matrices.push_back(mat);
    }
    Matrix full() const {
        return kvs_full(sub_matrices);
    }
    void print_submatrices(std::ostream& out) const
    {
        for (const auto & m : sub_matrices) out << "Matrix: " << std::endl << m << std::endl;;
    }
    bool operator==(const KroneckerVectorStack & other) const
    {
        return std::inner_product( getSubMatrices().begin(),
                getSubMatrices().end(),
                other.getSubMatrices().begin(),
                true,
                std::logical_or<bool>(), // reduction operator
                std::equal_to<Matrix>() // comparison operator
                );
    }
};
template <typename T>
std::ostream& operator<<(
    std::ostream& out, KroneckerVectorStack<T> K
)
{
    K.print_submatrices(out);
    return out;
}
struct key_iterator : std::unary_function<size_t, size_t> {
    size_t h_kvs_el;
    __host__ __device__
        key_iterator(size_t h_kvs_el) : h_kvs_el{h_kvs_el} {}

    __host__ __device__ 
        size_t operator()(size_t tmp_idx) 
        { 
            size_t idx = tmp_idx / h_kvs_el;
            //std::cout << "Key: " << idx << "   ";
            return idx;
        }
};

template<typename T>
using zipped = thrust::tuple<T, T>;

template <typename T>
struct mult_zip : std::unary_function<zipped<T>,  T> {
    __host__ __device__
        T operator()(const zipped<T> & vals) { return thrust::get<0>(vals) * thrust::get<1>(vals); }
};


struct kvs_iterator : std::unary_function<size_t, size_t> {
    size_t h_kvs_el, h_tmp;
    __host__ __device__ 
        kvs_iterator(size_t h_kvs_el, size_t h_tmp) : h_kvs_el{h_kvs_el}, h_tmp{h_tmp} {}

    __host__ __device__
        size_t operator()(size_t tmp_idx) 
        {
            size_t idx_pt1 = h_kvs_el * (tmp_idx / h_tmp);
            size_t idx_pt2 = tmp_idx % h_kvs_el;
            //std::cout << "Idx: " << idx_pt1 + idx_pt2 << " (pt 1: "<< idx_pt1 << " pt 2: " << idx_pt2 << ") \n";
            return idx_pt1 + idx_pt2;
        }
};


template<typename MatrixType>
MatrixType dot(const MatrixType & vec, const KroneckerVectorStack<MatrixType> & kvs) 
{
    assert(vec.nR() == 1); // must be a row vector!
    // We do one stage like a normal vec * kronmat -
    // vec.reshape(m, k).transpose() * kvs_0, where kvs_0 has shape k * k. m = \prod_{d != 0} k_d
    // This uses the inplace-reshape and traspose facility of the matrix class to perform this without copying
    // that function sig is:
    // M, K, N, meTrans, otherMatrix, otherTrans
    size_t kvsM = 1; 
    for (const auto & mat : kvs.getSubMatrices()) kvsM *= mat.nR();
    size_t kvsN = kvs.getSubMatrices().back().nC();
    size_t kvs0M = kvs.getSubMatrices().back().nR();
    size_t kvsNot0M = kvsM / kvs0M;
    MatrixType tmpmat = vec.dot(kvsNot0M, kvs0M, kvsN, Trans, kvs.getSubMatrices().back(), None); 
    typename MatrixType::Storage tmp = tmpmat.getConstData();
    size_t h_tmp = tmpmat.nR();
    size_t tmp_count = tmpmat.nR() * tmpmat.nC();
    size_t dim = kvs.getSubMatrices().size();
    size_t out_count = tmp_count / kvs.getSubMatrices()[dim - 2].nR();

    // Next, working backwards through dims, we create a new matrix where each 
    // element is the dot products along columns of step1 and that kvs member
    // we then swap this matrix and step1 and continue in the loop
    //for (const auto & val : tmp) std::cout << val << ", ";
    for (size_t d = dim - 2; d != 0; --d) 
    {
        // src_it indexes tmp, which is a dense matrix. 
        // kron_it indexes kvs_el, which is a kvs matrix and needs to tile over tmp.
        // Visually:
        // +-----------+  +-----------+
        // |           |  |   kvs_el  |
        // |    tmp    |  +-----------+
        // |           |  |   kvs_el  |
        // +-----------+  +-----------+
        // In this example, we reduce tmp to two rows, by dotting each column 
        // segment with the adjacent one from kvs_el.
        //
        // This is implemented by: 
        // * tmp_it, an iterator to count the elements in tmp
        // * kvs_idx, an iterator to select the correct element of kvs_el for the given element of tmp.
        //      This is: tmp_it % height(kvs_el) + height(kvs_el) * tmp_it % height(tmp)
        // * kvs_it : an iterator looking up indexes in kvs using kvs_it
        // * key, an iterator defined as tmp_it / height(kve_el)
        // * out, a pre-allocated vector of size (height(tmp) / height(kvs_el)) x width(kvs_el)
        // 
        // We used thrust::reduce_by_key to sum the product of the values pointed to by kvs_it and tmp_it into 
        // an output array. The index into the output array is given by key - hence, the first element in 
        // the output is the dot of part of the first column of tmp with the first column of kvs_el, the next is
        // the dot of the second part of the first column of tmp with the (same) first column of kvs_el, etc.

        const MatrixType & kvs_el = kvs.getSubMatrices()[d];
        size_t h_kvs_el = kvs_el.nR();
        size_t w_kvs_el = kvs_el.nC();

        typename MatrixType::Storage out{out_count};

        thrust::counting_iterator<size_t> tmp_idx(0);
        auto kvs_idx = thrust::make_transform_iterator(tmp_idx, kvs_iterator(h_kvs_el, h_tmp));
        auto kvs_it = thrust::make_permutation_iterator(kvs_el.getConstData().begin(), kvs_idx);
        auto key_first = thrust::make_transform_iterator(tmp_idx, key_iterator(h_kvs_el));
        auto key_last = key_first + tmp_count;
        auto tmp_zip_kvs = thrust::make_zip_iterator(thrust::make_tuple(tmp.begin(), kvs_it));
        auto tmp_times_kvs = thrust::make_transform_iterator(tmp_zip_kvs, mult_zip<typename MatrixType::Storage::value_type>());
        thrust::reduce_by_key(key_first, key_last, tmp_times_kvs, thrust::make_discard_iterator(), out.begin());
        /*
        for (size_t idx = 0; idx < tmp_count; idx++) 
        {
            std::cout << kvs_idx[idx] << ", ";
            std::cout << kvs_it[idx] << ", ";
            std::cout << tmp[idx] << ", ";
            std::cout << key_first[idx] << ", ";
            std::cout << tmp_times_kvs[idx] << std::endl;
        }
        */

        tmp_count = out_count;
        out_count /= kvs.getSubMatrices()[d - 1].nR();
        h_tmp /= h_kvs_el;
        thrust::swap( out, tmp );
        //for (const auto & val : tmp) std::cout << val << ", ";
    }
    // Finally, data is in tmp. Allocate a matrix and do the same algo into that:
    const MatrixType & kvs_el = kvs.getSubMatrices()[0];
    size_t h_kvs_el = kvs_el.nR();
    size_t w_kvs_el = kvs_el.nC();
    MatrixType ans(1, kvs_el.nC());
    thrust::counting_iterator<size_t> tmp_idx(0);
    auto kvs_idx = thrust::make_transform_iterator(tmp_idx, kvs_iterator(h_kvs_el, h_tmp));
    auto kvs_it = thrust::make_permutation_iterator(kvs_el.getConstData().begin(), kvs_idx);
    auto key_first = thrust::make_transform_iterator(tmp_idx, key_iterator(h_kvs_el));
    auto key_last = key_first + tmp_count;
    auto tmp_zip_kvs = thrust::make_zip_iterator(thrust::make_tuple(tmp.begin(), kvs_it));
    auto tmp_times_kvs = thrust::make_transform_iterator(tmp_zip_kvs, mult_zip<typename MatrixType::Storage::value_type>());
    thrust::reduce_by_key(key_first, key_last, tmp_times_kvs, thrust::make_discard_iterator(), ans.getMutableData().begin());
    return ans;
}
} // kronlib
#endif //KRONMAT_KRONVECSTACK_H
