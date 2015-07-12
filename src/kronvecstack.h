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
    KroneckerVectorStack(std::vector<Matrix<T>> sub_matrices) : sub_matrices{sub_matrices} {}
    KroneckerVectorStack() {}
    std::vector<Matrix<T>> sub_matrices;

public:
    bool isTrans() const {
        //TODO: better way to keep transpose in sync!
        return sub_matrices[0].isTrans();
    }
    void push_matrix(const Matrix<T> & mat)
    {
        sub_matrices.push_back(mat);
    }
    KroneckerVectorStack operator*(const KroneckerMatrix<T> & other) const
    {
        assert(!isTrans()); // can only do this with the kronecker dimension
        auto newsub_matrices = kronmat_dot_kronmat(sub_matrices, other.sub_matrices);
        return KroneckerVectorStack(newsub_matrices);
    }

    Matrix<T> full() const {
        auto ans = kvs_full(sub_matrices);
        if (isTrans()) {
            ans.mutable_transpose();
        }
        return ans;
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

    void print_submatrices(std::ostream& out)
    {
        for (const auto & m : sub_matrices) out << "Matrix: " << std::endl << m << std::endl;;
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
#endif //KRONMAT_KRONVECSTACK_H
