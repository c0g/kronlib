#include <iostream>
#include <vector>
#include <random>
#include <boost/timer/timer.hpp>
#include <functional>
#include <stdlib.h>
#include <stdio.h>
//#include "add_full_kron.h"
#include "matrix.h"


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


    Matrix<float> left(1, 3);
    left = 1, 2, 3;

    Matrix<float> kvs1(3, 2);
    kvs1 = 1, 2, 3, 4, 5, 6;

    Matrix<float> kvs2(3, 2);
    kvs2 = 2, 3, 4, 5, 6, 7;

    Matrix<float> kvs3(3, 3);
    kvs3 = 1, 2, 3, 4, 5, 6, 7, 8, 9;

    Matrix<float> kvs4(3, 3);
    kvs4 = 1, 2, 3, 4, 5, 6, 7, 8, 9;


    KroneckerVectorStack<float> kvs;
    kvs.push_matrix(kvs1);
    kvs.push_matrix(kvs2);
    kvs.push_matrix(kvs3);
    kvs.push_matrix(kvs4);

    auto ans = left * kvs;

    Matrix<float> real_ans(1, 60);
    real_ans =
            0.,   2.,   4.,   6.,   8.,   0.,   4.,   8.,  12.,  16.,   0.,
            6.,  12.,  18.,  24.,   0.,   3.,   6.,   9.,  12.,   0.,   6.,
            12.,  18.,  24.,   0.,   9.,  18.,  27.,  36.,   0.,   4.,   8.,
            12.,  16.,   0.,   8.,  16.,  24.,  32.,   0.,  12.,  24.,  36.,
            48.,   0.,   6.,  12.,  18.,  24.,   0.,  12.,  24.,  36.,  48.,
            0.,  18.,  36.,  54.,  72.,
            240.,   288.,   336.,   384.,   432.,   300.,   360.,   420.,
            480.,   540.,   360.,   432.,   504.,   576.,   648.,   300.,
            360.,   420.,   480.,   540.,   375.,   450.,   525.,   600.,
            675.,   450.,   540.,   630.,   720.,   810.,   320.,   384.,
            448.,   512.,   576.,   400.,   480.,   560.,   640.,   720.,
            480.,   576.,   672.,   768.,   864.,   400.,   480.,   560.,
            640.,   720.,   500.,   600.,   700.,   800.,   900.,   600.,
            720.,   840.,   960.,   1080.,
            2100.,  2310.,  2520.,  2730.,  2940.,  2400.,  2640.,  2880.,
            3120.,  3360.,  2700.,  2970.,  3240.,  3510.,  3780.,  2450.,
            2695.,  2940.,  3185.,  3430.,  2800.,  3080.,  3360.,  3640.,
            3920.,  3150.,  3465.,  3780.,  4095.,  4410.,  2520.,  2772.,
            3024.,  3276.,  3528.,  2880.,  3168.,  3456.,  3744.,  4032.,
            3240.,  3564.,  3888.,  4212.,  4536.,  2940.,  3234.,  3528.,
            3822.,  4116.,  3360.,  3696.,  4032.,  4368.,  4704.,  3780.,
            4158.,  4536.,  4914.,  5292.;



};


