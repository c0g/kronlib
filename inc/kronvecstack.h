//
// Created by Thomas Nickson on 06/07/2015.
//

#ifndef KRONMAT_KRONVECSTACK_H
#define KRONMAT_KRONVECSTACK_H

#include "matrix.h"
#include <functional>
namespace kronlib {
template <typename Matrix>
class KroneckerVectorStack {

private:
    std::vector<Matrix> sub_matrices;
public:
    KroneckerVectorStack(const std::vector<Matrix> & sub_matrices) : sub_matrices{sub_matrices} {}
    KroneckerVectorStack() {}

public:
    const std::vector<Matrix> & getSubMatrices() 
    {
        return sub_matrices;
    }

    void add_dimension(const Matrix & mat)
    {
        sub_matrices.push_back(mat);
    }
    Matrix full() const {
        return kvs_full(sub_matrices);
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

template<typename MatrixType>
MatrixType dot(const KroneckerVectorStack<MatrixType> & kvs, const MatrixType & vec) 
{
    assert(vec.nC() == 1); // Must be a column vector!

}

template<typename MatrixType>
MatrixType dot(const MatrixType & vec, const KroneckerVectorStack<MatrixType> & kvs) 
{
    assert(vec.nR() == 1); // must be a row vector!
    // We do one stage like a normal vec * kronmat -
    // vec.reshape(m, k).transpose() * kvs_0

}
} // kronlib
#endif //KRONMAT_KRONVECSTACK_H
