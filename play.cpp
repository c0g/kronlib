#include <iostream>
#include <vector>
#include <random>
#include <ctime>
#include <math.h>
#include <boost/timer.hpp>
//#include <dlib/matrix.h>

#include "matrix.h"
#include "kronecker_matrix.h"
#include "static_functions.h"

using namespace tommat;
int MAT_SIZE = 10;
int MAT_R=5000;
int MAT_C=5000;
int VEC_LEN = MAT_SIZE * MAT_SIZE;
int main() {

    matrix<double> m(3, 2);
    m = 1,2,
        3,4,
        5,6;
    std::cout << m << std::endl;
    std::cout << m.transpose() << std::endl;

    auto mm = m.transpose() * m;
    std::cout << mm << std::endl;


}

