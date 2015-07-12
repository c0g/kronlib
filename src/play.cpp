#include <iostream>
#include <vector>
#include <random>
#include <boost/timer/timer.hpp>
#include <functional>
#include <stdlib.h>
#include <stdio.h>
//#include "add_full_kron.h"
#include "matrix.h"
#include "kernel.h"


#include "kronvecstack_dot_vec.h"
#include "distances.h"
#include "parameter.h"
#include "kronecker_matrix.h"
#include "cholesky.h"

using ntype = double;


//Matrix<ntype> make_mat()
//{
//    return Matrix<ntype>(3,3);
//}

int main() {

    Matrix<float> x{10, 1};
    x = 1,2,3,4,5,6,7,8,9,10;

    const float delta = 1e-5;
//    const float invdelta = 1e4;
    const float onepdelt = 1+delta;
    sqexp_hyp<float> hyp{1, 1};
    sqexp_hyp<float> hyp_ll{onepdelt, 1};
    sqexp_hyp<float> hyp_ls{1, onepdelt};

    sqexp1d<float> k;

    auto normal = k(hyp, x, x);
    auto dll = k(hyp_ll, x, x);
    auto dls = k(hyp_ls, x, x);

    auto derivs = k.dhyp(hyp, x, x);
    std::cout << (dll - normal) / delta << std::endl;
    std::cout << derivs[0] << std::endl;

    std::cout << (dls - normal) / delta << std::endl;
    std::cout << derivs[1] << std::endl;





};


