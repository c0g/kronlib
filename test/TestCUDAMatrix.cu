//
// Created by Thomas Nickson on 05/07/2015.
//

#include "kronlib.h"
#include "iostream"
#include "gtest/gtest.h"
using namespace kronlib;
TEST(Matrix, SetScalarCUDA)
{
    CUDAMatrix<float> mat(2, 3);
    mat = 1;
    for (int r = 0; r < 2; ++r) {
        for (int c = 0; c < 3; ++c) {
            EXPECT_EQ(mat(r, c), 1);
        }
    }
}
TEST(Matrix, SetSeqCUDA)
{
    CUDAMatrix<float> mat(2, 3);
    mat = 1, 2, 3, 4, 5, 6;
    int ans = 1;
    for (int r = 0; r < 2; ++r) {
        for (int c = 0; c < 3; ++c) {
            EXPECT_EQ(mat(r, c), ans);
            ++ans;
        }
    }
}

TEST(Matrix, MinusInplaceCUDA)
{
    CUDAMatrix<float> mat(2, 2);
    mat = 1, 2, 3, 4;
    
    CUDAMatrix<float> mat_minus(2, 2);
    mat_minus = -1, -2, -3, -4;

    mat.negate_inplace();

    EXPECT_EQ(mat, mat_minus);
}

TEST(Matrix, MinusNotInplace)
{
    CUDAMatrix<float> mat(2, 2);
    mat = 1, 2, 3, 4;
    
    CUDAMatrix<float> mat_minus(2, 2);
    mat_minus = -1, -2, -3, -4;

    auto mat_minus2 = -mat;

    EXPECT_EQ(mat_minus2, mat_minus);
}

TEST(Matrix, ExpInplace)
{
    CUDAMatrix<float> mat(2, 2);
    mat = 1, 2, 3, 4;
    
    CUDAMatrix<float> mat_exp(2, 2);
    mat_exp = std::exp(1), std::exp(2), std::exp(3), std::exp(4);

    mat.exp_inplace();

    EXPECT_EQ(mat, mat_exp);
}

TEST(Matrix, ExpNotInplace)
{
    CUDAMatrix<float> mat(2, 2);
    mat = 1, 2, 3, 4;
    
    CUDAMatrix<float> mat_exp(2, 2);
    mat_exp = std::exp(1), std::exp(2), std::exp(3), std::exp(4);

    auto mat_exp2 = exp(mat);

    EXPECT_EQ(mat_exp2, mat_exp);
}

TEST(Matrix, EqualsAssignment)
{
    CUDAMatrix<float> mat(2, 3);
    mat = 1, 2, 3, 4, 5, 6;
    CUDAMatrix<float> mat2 = mat + 1;
    for (int r = 0; r < 2; ++r) {
        for (int c = 0; c < 3; ++c) {
            if ((r == 0) && (c == 0)) {
                EXPECT_EQ(mat(r, c), 1);
            } else if ((r == 0) && (c == 1)) {
                EXPECT_EQ(mat(r, c), 2);
            } else if ((r == 0) && (c == 2)) {
                EXPECT_EQ(mat(r, c), 3);
            } else if ((r == 1) && (c == 0)) {
                EXPECT_EQ(mat(r, c), 4);
            } else if ((r == 1) && (c == 1)) {
                EXPECT_EQ(mat(r, c), 5);
            } else if ((r == 1) && (c == 2)) {
                EXPECT_EQ(mat(r, c), 6);
            }
        }
    }
}

TEST(Matrix, Equality)
{
    CUDAMatrix<float> mat1(2, 3);
    mat1 = 1, 2, 3, 4, 5, 6;
    CUDAMatrix<float> mat2(2, 3);
    mat2 = 1, 2, 3, 4, 5, 6;
    ASSERT_TRUE(mat1 == mat2);
}
TEST(Matrix, Inequality)
{
    CUDAMatrix<float> mat1(2, 3);
    mat1 = 1, 2, 3, 4, 5, 6;
    CUDAMatrix<float> mat2(2, 3);
    mat2 = 1, 2, 3, 4, 5, 7;
    ASSERT_FALSE(mat1 == mat2);
}

TEST(Matrix, TransposeIndexing)
{
    CUDAMatrix<float> mat(20, 30);
    for (int r = 0; r < 20; ++r) {
        for (int c = 0; c < 30; ++c) {
            mat.setat(r, c, r + c);
        }
    }
    auto mat2 = mat.transpose();
    for (int r = 0; r < 20; ++r) {
        for (int c = 0; c < 30; ++c) {
            EXPECT_EQ(mat(r, c), mat2(c, r));
        }
    }
}

TEST(Matrix, Negate)
{
    CUDAMatrix<float> mat(2, 3);
    mat = 1, 2, 3, 4, 5, 6;
    auto mat2 = -mat;
    for (int r = 0; r < 2; ++r) {
        for (int c = 0; c < 3; ++c) {
            EXPECT_EQ(mat(r, c), -mat2(r, c));
        }
    }
}

TEST(Matrix, ScalarAdd)
{
    CUDAMatrix<float> mat(2, 3);
    mat = 1, 2, 3, 4, 5, 6;
    CUDAMatrix<float> mat2 = mat + 1;
    for (int r = 0; r < 2; ++r) {
        for (int c = 0; c < 3; ++c) {
            EXPECT_EQ(mat(r, c) + 1, mat2(r, c));
        }
    }
}

TEST(Matrix, ScalarMinus)
{
    CUDAMatrix<float> mat(2, 3);
    mat = 1, 2, 3, 4, 5, 6;
    CUDAMatrix<float> mat2 = mat - 1;
    for (int r = 0; r < 2; ++r) {
        for (int c = 0; c < 3; ++c) {
            EXPECT_EQ(mat(r, c) - 1, mat2(r, c));
        }
    }
}

TEST(Matrix, MatrixAdd)
{
    CUDAMatrix<float> mat1(2, 3);
    mat1 = 1, 2, 3, 4, 5, 6;
    CUDAMatrix<float> mat2(2, 3);
    mat2 = 0, 1, 2, 3, 4, 5;

    auto mat_sum = mat1 + mat2;

    for (int r = 0; r < 2; ++r) {
        for (int c = 0; c < 3; ++c) {
            EXPECT_EQ(mat1(r, c) + mat2(r, c), mat_sum(r, c));
        }
    }
}

TEST(Matrix, MatrixMinus)
{
    CUDAMatrix<float> mat1(2, 3);
    mat1 = 1, 2, 3, 4, 5, 6;
    CUDAMatrix<float> mat2(2, 3);
    mat2 = 0, 1, 2, 3, 4, 5;

    auto mat_sum = mat1 - mat2;

    for (int r = 0; r < 2; ++r) {
        for (int c = 0; c < 3; ++c) {
            EXPECT_EQ(mat1(r, c) - mat2(r, c), mat_sum(r, c));
        }
    }
}

TEST(Matrix, ScalarMul)
{
    CUDAMatrix<float> mat(2, 3);
    mat = 1, 2, 3, 4, 5, 6;
    CUDAMatrix<float> mat2 = mat * 2;
    for (int r = 0; r < 2; ++r) {
        for (int c = 0; c < 3; ++c) {
            EXPECT_EQ(mat(r, c) * 2, mat2(r, c));
        }
    }
}

TEST(Matrix, ScalarDiv)
{
    CUDAMatrix<float> mat(2, 3);
    mat = 1, 2, 3, 4, 5, 6;
    CUDAMatrix<float> mat2 = mat / 2;
    for (int r = 0; r < 2; ++r) {
        for (int c = 0; c < 3; ++c) {
            EXPECT_EQ(mat(r, c) / 2, mat2(r, c));
        }
    }
}



TEST(Matrix, Hadamard)
{
    CUDAMatrix<float> mat1(2, 3);
    mat1 = 1, 2, 3, 4, 5, 6;
    CUDAMatrix<float> mat2(2, 3);
    mat2 = 0, 1, 2, 3, 4, 5;
    auto mat_hadam = mat1.elemwise_mult(mat2);
    for (int r = 0; r < 2; ++r) {
        for (int c = 0; c < 3; ++c) {
            EXPECT_EQ(mat1(r, c) * mat2(r, c), mat_hadam(r, c));
        }
    }
}
TEST(Matrix, SimpleProduct)
{
    CUDAMatrix<float> mat1(2, 3);
    mat1 = 1, 2, 3, 4, 5, 6;
    CUDAMatrix<float> mat2(3, 2);
    mat2 = 0, 1, 2, 3, 4, 5;

    CUDAMatrix<float> mat_ans12(2, 2);
    mat_ans12 = 16, 22, 34, 49;

    CUDAMatrix<float> mat_ans21(3, 3);
    mat_ans21 = 4, 5, 6, 14, 19, 24, 24, 33, 42;

    auto mat_prod12 = mat1 * mat2;
    auto mat_prod21 = mat2 * mat1;

    EXPECT_EQ(mat_prod12, mat_ans12);
    EXPECT_EQ(mat_prod21, mat_ans21);
}

TEST(Matrix, ProductWithTrans)
{
    CUDAMatrix<float> mat1(2, 3);
    mat1 = 1, 2, 3, 4, 5, 6;
    auto mat2 = mat1.transpose();

    CUDAMatrix<float> mat_ans(2, 2);
    mat_ans = 14, 32, 32, 77;

    auto mat_prod = mat1 * mat2;

    EXPECT_EQ(mat_prod, mat_ans);
}

TEST(Matrix, ProductShape)
{
    CUDAMatrix<float> mat1(2, 3);
    mat1 = 1, 2, 3, 4, 5, 6;

    CUDAMatrix<float> mat2(3, 1);
    mat2 = 1, 2, 3;


    CUDAMatrix<float> mat_ans(2, 1);
    mat_ans = 14, 32;

    auto mat_prod = mat1 * mat2;

    EXPECT_EQ(mat_prod, mat_ans);
    EXPECT_EQ(mat_prod.nR(), mat_ans.nR());
    EXPECT_EQ(mat_prod.nC(), mat_ans.nC());

}

TEST(Matrix, BigProductSingleCUDA)
{
    int M = 200;
    int N = 200;
    int K = 200;
    CUDAMatrix<float> hmat1(M, N);
    CUDAMatrix<float> hmat2(N, K);
    for (int r = 0; r < M; ++r) {
        for (int c = 0; c < N; ++c) {
            hmat1.setat(r, c, (r + c)/1000);
        }
    }
    for (int r = 0; r < N; ++r) {
        for (int c = 0; c < K; ++c) {
            hmat2.setat(r, c, (r + c)/1000);
        }
    }
    // Manual Matrix mult
    CUDAMatrix<float> hmat_ans(M, K);
    hmat_ans = 0;
    for (int r = 0; r < M; ++r) {
        for (int c = 0; c < K; ++c) {
            for (int el = 0; el < N; ++el) {
                hmat_ans.setat(r, c, hmat_ans(r,c) + hmat1(r, el) * hmat2(el, c));
            }
        }
    }
    auto hmat_prod = hmat1 * hmat2;

    EXPECT_EQ(hmat_prod, hmat_ans);
    EXPECT_EQ(hmat_ans.nR(), hmat_prod.nR());
    EXPECT_EQ(hmat_ans.nC(), hmat_prod.nC());

}

TEST(Matrix, Reshape)
{
    CUDAMatrix<float> mat1(2, 3);
    mat1 = 1, 2, 3, 
           4, 5, 6;

    CUDAMatrix<float> mat2(3, 2);
    mat2 = 1, 5,
           4, 3,
           2, 6;

    EXPECT_EQ(mat2.reshape(2, 3), mat1);
}

TEST(Matrix, TransposeReshape)
{
    CUDAMatrix<float> mat(2, 3);
    mat = 1, 3, 4, 5, 7, 9;

    CUDAMatrix<float> ans(2, 3);
    ans = 1, 4, 7, 3, 5, 9;

    auto mat_trans = mat.transpose();
    auto mat_res = mat_trans.reshape(2, 3);

    EXPECT_EQ(ans, mat_res);
}

TEST(Matrix, ReshapeTranspose)
{
    CUDAMatrix<float> mat(2, 3);
    mat = 1, 3, 4, 5, 7, 9;

    CUDAMatrix<float> ans(2, 3);
    ans = 1, 5, 3, 7, 4, 9;

    auto mat_trans = mat.reshape(3, 2).transpose();

    EXPECT_EQ(ans, mat_trans);
}

TEST(Matrix, TransposeReshapeTranspose)
{
    CUDAMatrix<float> mat(2, 3);
    mat = 1, 3, 4, 5, 7, 9;

    CUDAMatrix<float> ans(3, 2);
    ans = 1, 3, 4, 5, 7, 9;

    auto mat_trans = mat.transpose().reshape(2, 3).transpose();

    EXPECT_EQ(ans, mat_trans);
}

TEST(Matrix, ReshapeTransposeReshape)
{
    CUDAMatrix<float> mat(2, 3);
    mat = 1, 3, 4, 5, 7, 9;

    CUDAMatrix<float> ans(2, 3);
    ans = 1, 4, 7, 3, 5 ,9;

    auto mat_trans = mat.reshape(2, 3).transpose().reshape(2, 3);

    EXPECT_EQ(ans, mat_trans);
}

TEST(Matrix, TraceEqual)
{
    CUDAMatrix<float> mat(2, 2);
    mat = 1, 2, 3, 4;
    EXPECT_EQ(mat.trace(), 5);
}
TEST(Matrix, TraceTall)
{
    CUDAMatrix<float> mat(3, 2);
    mat = 1, 2, 3, 4, 5, 6;
    EXPECT_EQ(mat.trace(), 5);
}
TEST(Matrix, TraceWide)
{
    CUDAMatrix<float> mat(2, 3);
    mat = 1, 2, 3, 4, 5, 6;
    EXPECT_EQ(mat.trace(), 6);
}
int main(int argc, char **argv)
{
    std::cout << "Running test" << std::endl;
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
