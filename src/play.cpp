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

using ntype = double;


//matrix<ntype> make_mat()
//{
//    return matrix<ntype>(3,3);
//}

int main() {


    matrix<float> mat1(2, 2);
    mat1 = 1,2,3,4;
    matrix<float> mat2(2, 2);
    mat2 = 3,4,5,6;
    matrix<float> mat3(2, 3);
    mat3 = 3,4,5,6,7,8;



    kronecker_matrix<float> kron_mat;
    kron_mat.push_matrix(mat1);
    kron_mat.push_matrix(mat2);
    kron_mat.push_matrix(mat3);

    auto mat = kron_mat.full();



};


