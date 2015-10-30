#include <iostream>

#include <cstdlib>
#include <ctime>

extern "C" {
#include <cblas.h>
}

using ntype = double;


int main() {

    int M = 2;
    int N = 3;
    int D = 4;
    double * kvs1 = new double[M * N]{1};
    double * kvs2 = new double[M * N]{1};
    double * kvs3 = new double[M * N]{1};
    double * kvs4 = new double[M * N]{1};

    double * vec = new double[N*N*N*N]{2};

    long * shapes = new long[D]{M, M, M, M};

    double ** kvs = new double*[D]{kvs1, kvs2, kvs3, kvs4};



    double * test_ans = new double[N*N*N*N]{0};
    for (int nidx = 0; nidx < N; ++ nidx)
    {
        for (int midx = 0; midx < M*M*M*M; ++ midx)
        {
            long idx0 = (midx / M*M*M) % 2;
            long idx1 = (midx / M*M) % 2;
            long idx2 = (midx / M) % 2;
            long idx3 = (midx) % 2;
            test_ans[midx] += vec[nidx] * kvs1[idx0 + M * nidx] * kvs2[idx1 + M * nidx] * kvs2[idx2 + M * nidx] * kvs3[idx3 + M * nidx];
        }
    }


    std::cout << std::endl;
    for (int el = 0; el < 16; ++el)
        std::cout << test_ans[el] << " ";
    std::cout << std::endl;



    //Now for the real deal
    double * src = new double[N*N*N*N]{0};
    double * dst = new double[N*N*N*N]{0};

    double * ans = new double[N*N*N*N]{0};

    std::cout << std::endl;
    for (int el = 0; el < 16; ++el)
        std::cout << ans[el] << " ";
    std::cout << std::endl;

    delete[] src;
    delete[] dst;
    delete[] kvs;
    delete[] kvs1;
    delete[] kvs2;
    delete[] kvs3;
    delete[] kvs4;
}


