#include <iostream>
#include "kronlib.h"

using namespace kronlib;

int main() {
    TBBMatrix<float> tbb(3000, 3000);
    CUDAMatrix<float> cuda(3000, 3000);
    HostMatrix<float> host(3000, 3000);

    auto tbbp = tbb * tbb.transpose();
    tbbp.negate_inplace();

    auto cudap = cuda * cuda.transpose();
    cudap.negate_inplace();
    auto hostp = host * host.transpose();
    hostp.negate_inplace();

    return 0;
}
