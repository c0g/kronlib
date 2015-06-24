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

using ntype = double;


//matrix<ntype> make_mat()
//{
//    return matrix<ntype>(3,3);
//}

int main() {


    matrix<ntype> m(3000,2000);
    matrix<ntype> m4(3000,2000);
    for (int i = 0; i < 3000; ++i) {
        for (int j = 0; j < 2000; ++j) {
            m(i,j) = i-j;
            m4(i,j) = j-i;
        }
    }

    std::cout << std::endl;
    boost::timer::auto_cpu_timer t1;
    auto mres = m4.transpose() * m ;
    std::cout << t1.format() << " " << mres(0,0) << std::endl;

    boost::timer::auto_cpu_timer t2;
    mres = m4.Tdot(m);
    std::cout << t2.format() << " " << mres(0,0) << std::endl;



};


