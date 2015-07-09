//
// Created by Thomas Nickson on 11/06/15.
//

#ifndef KRONMAT_KRONMAT_H
#define KRONMAT_KRONMAT_H

#include <vector>
#include <deque>
#include "matrix.h"
#include "static_functions.h"
#include "util.h"


template <typename T>
class kronecker_matrix  {
    friend matrix<T>;

public:


private:
public:
    kronecker_matrix(const std::vector<matrix<T>> &sub_matrices) : sub_matrices(sub_matrices) { }
    kronecker_matrix() {}

private:
    std::vector<matrix<T>> sub_matrices;
    T multiplier = 1;
    T sign = 1;

public:
    void push_matrix(const matrix<T> & mat)
    {
        sub_matrices.push_back(mat);
    }

    kronecker_matrix<T> operator*(const kronecker_matrix<T> & other)
    {
        auto newsub_matrices = kronmat_dot_kronmat(sub_matrices, other.sub_matrices);
        return kronecker_matrix(newsub_matrices);
    }

    matrix<T> operator*(const matrix<T> & other) {
        // assert(other.nC() == 1 && "The full matrix can only be a column vector");
        return kronmat_dot_fullvec(sub_matrices, other);
    }
    matrix<T> Tdot(const matrix<T> & other) {
        // assert(other.nC() == 1 && "The full matrix can only be a column vector");
        for (matrix<T> & m : sub_matrices) m.mutable_transpose();
        auto ans = kronmat_dot_fullvec((*this), other);
        for (matrix<T> & m : sub_matrices) m.mutable_transpose();
    }


    long nR()
    {
        long rows = 1;
        for (auto m : sub_matrices) {
            rows *= m.nR();
        }
        return rows;
    }
    long nC()
    {
        long cols = 1;
        for (auto m : sub_matrices) {
            cols *= m.nC();
        }
        return cols;
    }

    matrix<T> full()
    {
        return kron_full(sub_matrices);
    }



    bool operator==(const kronecker_matrix<T> & other) const
    {
        bool match = true;
        assert(sub_matrices.size() == other.sub_matrices.size());
        for (size_t idx = 0; idx < sub_matrices.size(); ++idx) {
            match &= (sub_matrices[idx] == other.sub_matrices[idx]);
        }
        return match;
    }
    void print_submatrices(std::ostream & out)
    {
        int matnum = 1;
        for (auto m : sub_matrices) {
            out << "Submatrix: " << matnum++ << std::endl;
            out << m << std::endl;
        }
    }
};
template <typename T>
std::ostream& operator<<(
    std::ostream& out, kronecker_matrix<T> K
)
{
    K.print_submatrices(out);
    return out;
}

template <typename T, typename N>
kronecker_matrix<T> operator*(N val, kronecker_matrix<T> m)
{
    return m * val;
}



// template <typename T>
// matrix<T> matrix<T>::operator*(const kronecker_matrix<T> &other) {

// //     kronecker_matrix has the operator kron_mat * matrix
// //     we can use that by noting:
// //         a * b = (b' * a')'

//     mutable_transpose();
//     other.mutable_transpose();
//     auto ans = other * (*this);
//     ans.mutable_transpose();
//     other.mutable_transpose();
//     mutable_transpose();
//     return ans;
// }

#endif //KRONMAT_KRONMAT_H
