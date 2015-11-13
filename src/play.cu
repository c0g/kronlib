#include <iostream>
#include <dlib/timing.h>
#include "kronlib.h"

using namespace kronlib;

int main() {
    TBBMatrix<float> tbb(3000, 3000);
    CUDAMatrix<float> cuda(3000, 3000);
    HostMatrix<float> host(3000, 3000);

    dlib::timing::start(1, "tbb");
    auto tbbp = tbb * tbb.transpose();
    tbbp.negate_inplace();

    dlib::timing::stop(1);
    dlib::timing::start(2, "cuda");
    auto cudap = cuda * cuda.transpose();
    cudap.negate_inplace();
    dlib::timing::stop(2);
    dlib::timing::start(3, "host");
    auto hostp = host * host.transpose();
    hostp.negate_inplace();
    dlib::timing::stop(3);
    dlib::timing::print();
    std::cout << tbbp(0,0) << cudap(0,0) << hostp(0,0) <<std::endl;


    return 0;
}
