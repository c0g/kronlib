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
class KroneckerVectorStack;

template <typename T>
class KroneckerMatrix  {
    friend Matrix<T>;
    friend KroneckerVectorStack<T>;

public:


private:
public:
    KroneckerMatrix(const std::vector<Matrix<T>> &sub_matrices) : sub_matrices(sub_matrices) { }
    KroneckerMatrix() {}

private:
    std::vector<Matrix<T>> sub_matrices;
public:
    void push_matrix(const Matrix<T> & mat)
    {
        sub_matrices.push_back(mat);
    }

    KroneckerMatrix<T> operator*(const KroneckerMatrix<T> & other)
    {
        auto newsub_matrices = kronmat_dot_kronmat(sub_matrices, other.sub_matrices);
        return KroneckerMatrix(std::move(newsub_matrices));
    }

    KroneckerVectorStack<T> operator*(const KroneckerVectorStack<T> & other)
    {
        auto newsub_matrices = kronmat_dot_kronmat(sub_matrices, other.sub_matrices);
        return KroneckerVectorStack<T>(std::move(newsub_matrices));
    }

    Matrix<T> operator*(const Matrix<T> & other)
    {
        assert(other.nC() == 1 && "The full Matrix can only be a column vector");
        return kronmat_dot_fullvec(sub_matrices, other);
    }
    Matrix<T> Tdot(const Matrix<T> & other)
    {
        assert(other.nC() == 1 && "The full Matrix can only be a column vector");
        for (Matrix<T> & m : sub_matrices) m.mutable_transpose();
        auto ans = kronmat_dot_fullvec((*this), other);
        for (Matrix<T> & m : sub_matrices) m.mutable_transpose();
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

    Matrix<T> full()
    {
        return kron_full(sub_matrices);
    }

    void mutable_transpose()
    {
        for (Matrix<T> & m : sub_matrices) m.mutable_transpose();
    }

    KroneckerMatrix transpose() const
    {
        KroneckerMatrix ans = (*this);
        ans.mutable_transpose();
        return ans;
    }

    bool operator==(const KroneckerMatrix<T> & other) const
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
            out << "SubMatrix: " << matnum++ << std::endl;
            out << m << std::endl;
        }
    }
};
template <typename T>
std::ostream& operator<<(
    std::ostream& out, KroneckerMatrix<T> K
)
{
    K.print_submatrices(out);
    return out;
}



// template <typename T>
// Matrix<T> Matrix<T>::operator*(const KroneckerMatrix<T> &other) {

// //     KroneckerMatrix has the operator kron_mat * Matrix
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
