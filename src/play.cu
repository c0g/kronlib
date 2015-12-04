#include <iostream>
#include <boost/timer/timer.hpp>
#include "kronlib.h"

using namespace kronlib;
using timer = boost::timer::auto_cpu_timer;

int main() {

    CUDAMatrix<float> _cmat1(3, 3);
    auto cans_ = _cmat1 * _cmat1;

    timer ct1;
    for (int i = 0; i < 50; ++i) {
        CUDAMatrix<float> cmat1(2000, 2000);
        auto cans1 = cmat1 * cmat1;
    }
    std::cout << ct1.format() << std::endl;
    ct1.stop();

    timer ht1;
    for (int i = 0; i < 50; ++i) {
        HostMatrix<float> mat1(2000, 2000);
        auto ans1 = mat1 * mat1;
    }
    std::cout << ht1.format() << std::endl;
    ht1.stop();

    timer ct2;
    for (int i = 0; i < 50; ++i) {
        CUDAMatrix<float> cone(200, 500);
        CUDAMatrix<float> ctwo(200, 500);
        CUDAMatrix<float> cthree(200, 500);
        KroneckerVectorStack<CUDAMatrix<float>> ckvs;
        ckvs.push_back( cone );
        ckvs.push_back( ctwo );
        ckvs.push_back( cthree );
        CUDAMatrix<float> cleftvec(1, 200 * 200 * 200);
        dot(cleftvec, ckvs);
    }
    std::cout << ct2.format() << std::endl;
    ct2.stop();

    timer ht2;
    for (int i = 0; i < 50; ++i) {
        HostMatrix<float> one(200, 500);
        HostMatrix<float> two(200, 500);
        HostMatrix<float> three(200, 500);
        KroneckerVectorStack<HostMatrix<float>> kvs;
        kvs.push_back( one );
        kvs.push_back( two );
        kvs.push_back( three );
        HostMatrix<float> leftvec(1, 200 * 200 * 200);
        dot(leftvec, kvs);
    }
    std::cout << ht2.format() << std::endl;
    ht2.stop();

    return 0;
}
