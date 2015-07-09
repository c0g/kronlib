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


//matrix<ntype> make_mat()
//{
//    return matrix<ntype>(3,3);
//}

int main() {


    matrix<float> mat(3, 3);
    mat = 1, 3, 4, 5, 7, 9, 4 ,5 ,6;

    Cholesky<float> chol(mat);
    matrix<float> cholmat = chol.cholmat();



};


