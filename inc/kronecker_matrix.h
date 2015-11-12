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



//template <typename T>
//class KroneckerVectorStack;

template <typename Storage>
class KroneckerMatrix  {
    friend Matrix<Storage>;
//    friend KroneckerVectorStack<T>;

public:


private:
public:
    KroneckerMatrix(const std::vector<Matrix<Storage>> &sub_matrices) : sub_matrices(sub_matrices) { }
    KroneckerMatrix() {}

private:
    std::vector<Matrix<Storage>> sub_matrices;
public:

    void push_matrix(const Matrix<Storage> & mat)
    {
        sub_matrices.push_back(mat);
    }

    KroneckerMatrix<Storage> operator*(const KroneckerMatrix<Storage> & other) const
    {
        auto newsub_matrices = kronmat_dot_kronmat(sub_matrices, other.sub_matrices);
        return KroneckerMatrix(newsub_matrices);
    }

    /*
    KroneckerVectorStack<T> operator*(const KroneckerVectorStack<T> & other) const
    {
        assert(other.isTrans()); // kronecker dimensions need to line up
        auto newsub_matrices = kronmat_dot_kronmat(sub_matrices, other.sub_matrices);
        KroneckerVectorStack<T> ans(newsub_matrices);
        return ans.transpose();
    }
    */
    Matrix<Storage> solve(const Matrix<Storage> & other) const
    {
        assert(other.nC() == 1 && "The full Matrix can only be a column vector");
        return kronmat_solve_fullvec(sub_matrices, other);
    }
    Matrix<Storage> operator*(const Matrix<Storage> & other) const
    {
        assert(other.nC() == 1 && "The full Matrix can only be a column vector");
        return kronmat_dot_fullvec(sub_matrices, other);
    }

    long nR() const
    {
        long rows = 1;
        for (auto m : sub_matrices) {
            rows *= m.nR();
        }
        return rows;
    }
    long nC() const
    {
        long cols = 1;
        for (auto m : sub_matrices) {
            cols *= m.nC();
        }
        return cols;
    }

    Matrix<Storage> full() const 
    {
        auto ans = kron_full(sub_matrices);
        return ans;

    }

    KroneckerMatrix transpose() const
    {
        KroneckerMatrix<Storage> ans;
        for (const auto & m : sub_matrices)
        {
            ans.push_matrix(m.transpose());
        }
        return ans;
    }

    bool operator==(const KroneckerMatrix<Storage> & other) const
    {
        bool match = true;
        assert(sub_matrices.size() == other.sub_matrices.size());
        for (size_t idx = 0; idx < sub_matrices.size(); ++idx) {
            match &= (sub_matrices[idx] == other.sub_matrices[idx]);
        }
        return match;
    }
    void print_submatrices(std::ostream & out) const
    {
        int matnum = 1;
        for (auto m : sub_matrices) {
            out << "SubMatrix: " << matnum++ << std::endl;
            out << m << std::endl;
        }
    }
};
template <typename Storage>
std::ostream& operator<< (
    std::ostream& out, KroneckerMatrix<Storage> K
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
