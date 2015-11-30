#include <iostream>
#include "kronlib.h"

using namespace kronlib;

int main() {
    CUDAMatrix<float> one(100, 100);
    CUDAMatrix<float> two(100, 100);
    CUDAMatrix<float> three(100, 100);
    KroneckerVectorStack<CUDAMatrix<float>> kvs{ std::vector<CUDAMatrix<float>>{one, two, three} };
    CUDAMatrix<float> leftvec(1, 100 * 100 * 100);
    std::cout << dot(leftvec, kvs).nR() << std::endl;
    std::cout << dot(leftvec, kvs).nC() << std::endl;

    return 0;
}
