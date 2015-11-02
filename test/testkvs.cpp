//
// Created by Thomas Nickson on 05/07/2015.
//

#include "matrix.h"
#include "kronecker_matrix.h"
#include "kronvecstack.h"
#include "iostream"
#include "gtest/gtest.h"

TEST(KVS, Full) 
{
    Matrix<float> kvs1(3, 2);
    kvs1 = 1, 2, 3, 4, 5, 6;

    Matrix<float> kvs2(3, 2);
    kvs2 = 2, 3, 4, 5, 6, 7;

    KroneckerVectorStack<float> kvs;
    kvs.push_matrix(kvs1);
    kvs.push_matrix(kvs2);

    Matrix<float> real_ans(3, 4);
    real_ans = 2.,   3.,   4.,   6.,  12.,  15.,  16.,  20.,  30.,  35.,  36.,
        42.;

    EXPECT_EQ(kvs.full(), real_ans);


}

TEST(KVS, LDotVecKVS_2)
{
    Matrix<float> left(1, 3);
    left = 1, 2, 3;

    Matrix<float> kvs1(3, 2);
    kvs1 = 1, 2, 3, 4, 5, 6;

    Matrix<float> kvs2(3, 2);
    kvs2 = 2, 3, 4, 5, 6, 7;

    KroneckerVectorStack<float> kvs;
    kvs.push_matrix(kvs1);
    kvs.push_matrix(kvs2);

    auto ans = left * kvs;

    Matrix<float> real_ans(1, 4);
    real_ans = 116.,  138.,  144.,  172.;

    EXPECT_EQ(ans, real_ans);
}

TEST(KVS, LDotVecKVS_3)
{
    Matrix<float> left(1, 3);
    left = 1, 2, 3;

    Matrix<float> kvs1(3, 2);
    kvs1 = 1, 2, 3, 4, 5, 6;

    Matrix<float> kvs2(3, 2);
    kvs2 = 2, 3, 4, 5, 6, 7;

    Matrix<float> kvs3(3, 3);
    kvs3 = 1, 2, 3, 4, 5, 6, 7, 8, 9;

    KroneckerVectorStack<float> kvs;
    kvs.push_matrix(kvs1);
    kvs.push_matrix(kvs2);
    kvs.push_matrix(kvs3);

    auto ans = left * kvs;

    Matrix<float> real_ans(1, 12);
    real_ans =  728.,   844.,   960.,   858.,   996.,  1134.,   888.,  1032.,
    1176.,  1048.,  1220.,  1392.;

    EXPECT_EQ(ans, real_ans);
}

TEST(KVS, LDotVecKVS_4)
{
    Matrix<float> left(1, 3);
    left = 1, 2, 3;

    Matrix<float> kvs1(3, 3);
    kvs1 = 1, 2, 3, 4, 5, 6, 7, 8, 9;

    Matrix<float> kvs2(3, 2);
    kvs2 = 2, 3, 4, 5, 6, 7;

    Matrix<float> kvs3(3, 3);
    kvs3 = 1, 2, 3, 4, 5, 6, 7, 8, 9;

    Matrix<float> kvs4(3, 4);
    kvs4 = 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12;


    KroneckerVectorStack<float> kvs;
    kvs.push_matrix(kvs1);
    kvs.push_matrix(kvs2);
    kvs.push_matrix(kvs3);
    kvs.push_matrix(kvs4);

    auto ans = left * kvs;

    Matrix<float> real_ans(1, 3 * 2 * 3 * 4);
    real_ans =  8580.,   9592.,  10604.,  11616.,   9876.,  11048.,  12220.,
         13392.,  11172.,  12504.,  13836.,  15168.,  10064.,  11256.,
         12448.,  13640.,  11590.,  12972.,  14354.,  15736.,  13116.,
         14688.,  16260.,  17832.,   9876.,  11048.,  12220.,  13392.,
         11376.,  12736.,  14096.,  15456.,  12876.,  14424.,  15972.,
         17520.,  11590.,  12972.,  14354.,  15736.,  13358.,  14964.,
         16570.,  18176.,  15126.,  16956.,  18786.,  20616.,  11172.,
         12504.,  13836.,  15168.,  12876.,  14424.,  15972.,  17520.,
         14580.,  16344.,  18108.,  19872.,  13116.,  14688.,  16260.,
         17832.,  15126.,  16956.,  18786.,  20616.,  17136.,  19224.,
         21312.,  23400.;

    EXPECT_EQ(ans, real_ans);
}

TEST(KVS, LDotWithKron)
{
    Matrix<float> kron1(2,2);
    kron1 = 1,2,3,4;
    Matrix<float> kron2(2,2);
    kron2 = 2,3,4,5;

    KroneckerMatrix<float> kron;
    kron.push_matrix(kron1);
    kron.push_matrix(kron2);

    Matrix<float> kvs1(2,2);
    kvs1 = 1,2,3,4;
    Matrix<float> kvs2(2,2);
    kvs2 = 1,2,3,4;

    KroneckerVectorStack<float> kvs;
    kvs.push_matrix(kvs1);
    kvs.push_matrix(kvs2);


    EXPECT_EQ(kvs.full() * kron.full(), (kvs * kron).full());
}


int main(int argc, char **argv)
{
    std::cout << "Running test" << std::endl;
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
