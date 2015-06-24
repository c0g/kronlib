//#include <iostream>
//#include <vector>
//#include <random>
//#include <ctime>
//#include <math.h>
//#include <boost/timer.hpp>
//#include <dlib/matrix.h>
#include <stdlib.h>

//#include "matrix.h"
//#include "kronecker_matrix.h"
//#include "static_functions.h"

#include <Accelerate/Accelerate.h>
#include <stdio.h>
//#include <vDSP.h>

//using namespace tommat;
//int MAT_SIZE = 10;
//int MAT_R=5000;
//int MAT_C=5000;
//int VEC_LEN = MAT_SIZE * MAT_SIZE;

inline void mul(float *__restrict__ a, float *__restrict__ b, float *__restrict__ c) {
    (*c) = (*a) * (*b);
}

int main() {

    float * v1;
    float * v2;
    float * v3;

    v1 = (float *)(malloc(1000000000 * sizeof(float)));
    v2 = (float *)(malloc(1000000000 * sizeof(float)));
    v3 = (float *)(malloc(1000000000 * sizeof(float)));

    for (int i = 0; i < 1000000000; ++i) {
        v1[i] = 2.0;
        v2[i] = 1.0;
    }


    vDSP_vmul(v1, 1, v2, 1, v3, 1, 1000000000);
    

    free(v1);
    free(v2);
    free(v3);


}

