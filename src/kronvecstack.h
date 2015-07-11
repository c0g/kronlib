//
// Created by Thomas Nickson on 06/07/2015.
//

#ifndef KRONMAT_KRONVECSTACK_H
#define KRONMAT_KRONVECSTACK_H

#include "matrix.h"
#include "kronecker_matrix.h"

template <typename T>
class KroneckerVectorStack {

    friend Matrix<T>;
    friend KroneckerMatrix<T>;
public:
    KroneckerVectorStack(std::vector<Matrix<T>> sub_matrices) : sub_matrices{sub_matrices}, trans{false} {}
    KroneckerVectorStack() {}
    std::vector<Matrix<T>> sub_matrices;

public:
    bool isTrans() const {
        return sub_matrices[0].isTrans();
    }
    void push_matrix(const Matrix<T> & mat)
    {
        sub_matrices.push_back(mat);
    }
    KroneckerVectorStack operator*(const KroneckerMatrix<T> & other) const
    {
        assert(!trans); // can only do this with the kronecker dimension
        auto newsub_matrices = kronmat_dot_kronmat(sub_matrices, other.sub_matrices);
        return KroneckerVectorStack(std::move(newsub_matrices));
    }

    Matrix<T> full() const {
        return kvs_full(sub_matrices);
    }

    void mutable_transpose()
    {
        for (Matrix<T> & m : sub_matrices) m.mutable_transpose();
    }

    KroneckerVectorStack transpose() const
    {
        KroneckerVectorStack ans = (*this);
        ans.mutable_transpose();
        return ans;
    }

};


#endif //KRONMAT_KRONVECSTACK_H
