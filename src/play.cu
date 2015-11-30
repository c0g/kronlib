#include <iostream>
#include "kronlib.h"

using namespace kronlib;

int main() {

    HostMatrix<float> mat1(2,2);
    mat1 = 1, 2,
           3, 4;

    auto ans1 = mat1 * mat1;
    auto ans2 = mat1.elemwise_mult(ans1);
    auto ans3 = mat1.dot(Trans, mat1, None);
    
    HostMatrix<double> dmat = ans3;
    Cholesky<HostMatrix<float>> chol{ dmat };
    std::cout << chol.solve(dmat);

    Kronecker<HostMatrix<float>> kron{ {ans1, ans2, ans3} };

    CUDAMatrix<float> cuda{ mat1 };

    CUDAMatrix<float> one(100, 100);
    CUDAMatrix<float> two(500, 100);
    CUDAMatrix<float> three(100, 100);
    KroneckerVectorStack<CUDAMatrix<float>> kvs{ std::vector<CUDAMatrix<float>>{one, two, three} };
    CUDAMatrix<float> leftvec(1, 100 * 500 * 100);
    std::cout << dot(leftvec, kvs).nR() << std::endl;
    std::cout << dot(leftvec, kvs).nC() << std::endl;

    return 0;
}
