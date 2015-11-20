//
// Created by Thomas Nickson on 11/06/15.
//

#ifndef KRONMAT_KRONMAT_H
#define KRONMAT_KRONMAT_H

#include <vector>
#include <deque>
#include "kronlib.h"

namespace kronlib {

template <typename MatrixType>
class Kronecker  {
private:
    std::vector<MatrixType> subMatrices;
public:
    Kronecker(const std::vector<MatrixType> &subMatrices) : subMatrices{subMatrices} { }
    Kronecker() {}

    typename std::vector<MatrixType>::size_type dim() const { return subMatrices.size(); }

    const std::vector<MatrixType> & getSubMatrices() const 
    {
        return subMatrices;
    }

    std::vector<MatrixType> & getMutableSubMatrices()
    {
        return subMatrices;
    }

    void push(const MatrixType & mat)
    {
        subMatrices.push_back(mat);
    }

    Kronecker<MatrixType> operator*(const Kronecker<MatrixType> & other) const
    {
        auto newsubMatrices = kronmat_dot_kronmat(subMatrices, other.subMatrices);
        return Kronecker(newsubMatrices);
    }

    KroneckerVectorStack<MatrixType> operator*(const KroneckerVectorStack<MatrixType> & other) const
    {
        auto newsubMatrices = kronmat_dot_kronmat(subMatrices, other.getSubMatrices());
        return KroneckerVectorStack{newsubMatrices};
    }
    MatrixType operator*(const MatrixType & other) const
    {
        assert(other.nC() == 1 && "The full Matrix can only be a column vector");
        return kronmat_dot_fullvec(subMatrices, other);
    }

    long nR() const
    {
        long rows = 1;
        for (auto m : subMatrices) {
            rows *= m.nR();
        }
        return rows;
    }
    long nC() const
    {
        long cols = 1;
        for (auto m : subMatrices) {
            cols *= m.nC();
        }
        return cols;
    }

    MatrixType full() const 
    {
        auto ans = kron_full(subMatrices);
        return ans;

    }

    Kronecker transpose() const
    {
        Kronecker<MatrixType> ans;
        for (const auto & m : subMatrices)
        {
            ans.push(m.transpose());
        }
        return ans;
    }

    bool operator==(const Kronecker<MatrixType> & other) const
    {
        bool match = true;
        assert(subMatrices.size() == other.subMatrices.size());
        for (size_t idx = 0; idx < subMatrices.size(); ++idx) {
            match &= (subMatrices[idx] == other.subMatrices[idx]);
        }
        return match;
    }
    void print_submatrices(std::ostream & out) const
    {
        int matnum = 1;
        for (auto m : subMatrices) {
            out << "SubMatrix: " << matnum++ << std::endl;
            out << m << std::endl;
        }
    }
};


template <typename MatrixType>
class Kronecker<Cholesky<MatrixType>>  {
private:
    std::vector<Cholesky<MatrixType>> subChols;
public:
    Kronecker(const std::vector<Cholesky<MatrixType>> &subMatrices) 
    {
        for (const auto & mat : subMatrices) { push(mat); }
    }
    explicit Kronecker(const Kronecker<MatrixType> & kronMat) : Kronecker{ kronMat.getSubMatrices } { }
    Kronecker() {}
    typename std::vector<Cholesky<MatrixType>>::size_type dim() const { return subChols.size(); }
    void push(const Cholesky<MatrixType> & chol) { subChols.push_back(chol); } 

    Kronecker<MatrixType> solve(const Kronecker<MatrixType> & other) 
    {
        Kronecker<MatrixType> ans;
        assert(dim() == other.dim());
        auto chol_it = std::begin(subChols);
        auto mat_it = std::begin(other.getSubMatrices());
        auto mat_end = std::end(other.getSubMatrices());
        for (/* chol_it and mat_it */ ; 
                (chol_it < subChols.end()) && (mat_it < mat_end); 
                ++chol_it, ++mat_it) 
        {
            ans.push( chol_it->solve( *mat_it ));
        }
        return ans;
    }
    MatrixType solve(const MatrixType & other)
    {
        assert(other.nC() == 1); // Matrix must be a column vector
        return kronchol_solve_fullvec(subChols, other);
    }


};


template <typename MatrixType>
std::ostream& operator<< (
    std::ostream& out, Kronecker<MatrixType> K
) 
{
    K.print_submatrices(out);
    return out;
}



// template <typename T>
// Matrix<T> Matrix<T>::operator*(const Kronecker<T> &other) {

// //     Kronecker has the operator kron_mat * Matrix
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
};
#endif //KRONMAT_KRONMAT_H


