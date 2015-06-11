#include <iostream>
#include <vector>
#include <random>
#include <ctime>
#include <math.h>
#include <boost/timer.hpp>

extern "C" {
    #include <cblas.h>
}

using namespace std;
int MAT_SIZE = 1000;
int VEC_LEN = MAT_SIZE * MAT_SIZE;
int main() {

    auto normal = bind( normal_distribution<double>(0, 1), default_random_engine());
    vector<double> v1, v2, mm1, mm2;
    cout << "Hello, Wordl!" << endl;

    for (int i = 0; i < VEC_LEN; ++i) {
        v1.push_back(normal());
        v2.push_back(normal());
        mm2.push_back(0);
    }
    boost::timer t1;
    double ans = cblas_ddot(VEC_LEN, v1.data(), 1, v2.data(), 1);
    cout << "CBlas dot product is " << ans <<  " in " << t1.elapsed() << endl;
    double ans2 = 0;

    boost::timer t4;
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, MAT_SIZE, MAT_SIZE, MAT_SIZE, 1.0, v1.data(), MAT_SIZE, v2.data(), MAT_SIZE, 1.0, mm2.data(), MAT_SIZE);
    cout << "CLBAS DGEMM product is " << mm2[0] <<  " in " << t4.elapsed() << endl;


    boost::timer t2;
    for (int i = 0; i < VEC_LEN; ++i) {
        ans2 += v1[i] * v2[i];
    }
    cout << "My dot product is " << ans <<  " in " << t2.elapsed() << endl;

    boost::timer t3;
    for (int r = 0; r < MAT_SIZE; ++r) {
        for (int mid = 0; mid < MAT_SIZE; ++mid) {
            double tmp = 0;
            for (int c = 0; c < MAT_SIZE; ++c) {
                tmp += v2[MAT_SIZE * r + mid] * v1[MAT_SIZE * mid + c];
            }
            mm1.push_back(tmp);
        }
    }
    cout << "My mm product is " << mm1[0] <<  " in " << t3.elapsed() << endl;




    return 0;
}