//
// Created by Thomas Nickson on 05/07/2015.
//

#include "testmatrix.h"
#include "../matrix.h"
#include "iostream"
#include "gtest/gtest.h"

TEST(Matrix, EqualsAssignment) {
    matrix<float> mat(2, 3);
    mat = 1, 2, 3, 4, 5, 6;
    matrix<float> mat2 = mat + 1;
    for (int r = 0; r < 2; ++r) {
        for (int c = 0; c < 3; ++c) {
            if ((r == 0) && (c == 0)) {
                EXPECT_EQ(mat(r,c), 1);
            }
            else if ((r == 0) && (c == 1)) {
                EXPECT_EQ(mat(r,c), 2);
            }
            else if ((r == 0) && (c == 2)) {
                EXPECT_EQ(mat(r,c), 3);
            }
            else if ((r == 1) && (c == 0)) {
                EXPECT_EQ(mat(r,c), 4);
            }
            else if ((r == 1) && (c == 1)) {
                EXPECT_EQ(mat(r,c), 5);
            }
            else if ((r == 1) && (c == 2)) {
                EXPECT_EQ(mat(r,c), 6);
            }
        }
    }
}

TEST(Matrix, Equality) {
    matrix<float> mat1(2, 3);
    mat1 = 1, 2, 3, 4, 5, 6;
    matrix<float> mat2(2, 3);
    mat2 = 1, 2, 3, 4, 5, 6;
    ASSERT_TRUE(mat1 == mat2); 
}
TEST(Matrix, Inequality) {
    matrix<float> mat1(2, 3);
    mat1 = 1, 2, 3, 4, 5, 6;
    matrix<float> mat2(2, 3);
    mat2 = 1, 2, 3, 4, 5, 7;
    ASSERT_FALSE(mat1 == mat2); 
}

TEST(Matrix, TransposeIndexing) {
    matrix<float> mat(20, 30);
    for (int r = 0; r < 20; ++r) {
        for (int c = 0; c < 30; ++c) {
            mat(r, c) = r + c;
        }
    }
    auto mat2 = mat.transpose();
    for (int r = 0; r < 20; ++r) {
        for (int c = 0; c < 30; ++c) {
            EXPECT_EQ(mat(r,c), mat2(c, r));
        }
    }
}

TEST(Matrix, ScalarAdd) {
    matrix<float> mat(2, 3);
    mat = 1, 2, 3, 4, 5, 6;
    matrix<float> mat2 = mat + 1;
    for (int r = 0; r < 2; ++r) {
        for (int c = 0; c < 3; ++c) {
            EXPECT_EQ(mat(r,c) + 1, mat2(r, c));
        }
    }
}

TEST(Matrix, MatrixAdd) {
    matrix<float> mat1(2, 3);
    mat1 = 1, 2, 3, 4, 5, 6;
    matrix<float> mat2(2, 3);
    mat2 = 0, 1, 2, 3, 4, 5;

    auto mat_sum = mat1 + mat2;

    for (int r = 0; r < 2; ++r) {
        for (int c = 0; c < 3; ++c) {
            EXPECT_EQ(mat1(r,c) + mat2(r,c), mat_sum(r, c));
        }
    }
}

TEST(Matrix, ScalarMul) {
    matrix<float> mat(2, 3);
    mat = 1, 2, 3, 4, 5, 6;
    matrix<float> mat2 = mat * 2;
    for (int r = 0; r < 2; ++r) {
        for (int c = 0; c < 3; ++c) {
            EXPECT_EQ(mat(r,c) * 2, mat2(r, c));
        }
    }
}

TEST(Matrix, Hadamard) {
    matrix<float> mat1(2, 3);
    mat1 = 1, 2, 3, 4, 5, 6;
    matrix<float> mat2(2, 3);
    mat2 = 0, 1, 2, 3, 4, 5;
    auto mat_hadam = mat1.hadamard(mat2);
    for (int r = 0; r < 2; ++r) {
        for (int c = 0; c < 3; ++c) {
            EXPECT_EQ(mat1(r,c) * mat2(r,c), mat_hadam(r,c));
        }
    }
}

TEST(Matrix, SimpleProduct) {
    matrix<float> mat1(2, 3);
    mat1 = 1, 2, 3, 4, 5, 6;
    matrix<float> mat2(3, 2);
    mat2 = 0, 1, 2, 3, 4, 5;
    
    matrix<float> mat_ans12(2,2);
    mat_ans12 = 16, 22, 34, 49;

    matrix<float> mat_ans21(3,3);
    mat_ans21 = 4, 5, 6, 14, 19, 24, 24, 33, 42;

    auto mat_prod12 = mat1 * mat2;
    auto mat_prod21 = mat2 * mat1;

    EXPECT_EQ(mat_prod12, mat_ans12);
    EXPECT_EQ(mat_prod21, mat_ans21);
}

TEST(Matrix, ProductWithTrans) {
    matrix<float> mat1(2, 3);
    mat1 = 1, 2, 3, 4, 5, 6;
    auto mat2 = mat1.transpose();
    
    matrix<float> mat_ans(2,2);
    mat_ans = 14, 32, 32, 77;

    auto mat_prod = mat1 * mat2;
    
    EXPECT_EQ(mat_prod, mat_ans);
}

TEST(Matrix, ProductShape) {
    matrix<float> mat1(2, 3);
    mat1 = 1, 2, 3, 4, 5, 6;
    
    matrix<float> mat2(3, 1);
    mat2 = 1, 2, 3;
    
    
    matrix<float> mat_ans(2,1);
    mat_ans = 14, 32;

    auto mat_prod = mat1 * mat2;

    EXPECT_EQ(mat_prod, mat_ans);
    EXPECT_EQ(mat_prod.nR(), mat_ans.nR());
    EXPECT_EQ(mat_prod.nC(), mat_ans.nC());

}

TEST(Matrix, BigProductSingle) {
    int M = 100;
    int N = 100;
    int K = 100;
    matrix<float> mat1(M, N);
    matrix<float> mat2(N, K);
    for (int r = 0; r < M; ++r) {
        for (int c = 0; c < N; ++c) {
            mat1(r, c) = (r + c);
        }
    }
    for (int r = 0; r < N; ++r) {
        for (int c = 0; c < K; ++c) {
            mat2(r, c) = (r + c);
        }
    }
    // Manual matrix mult
    matrix<float> mat_ans(M, K);
    mat_ans = 0;
    for (int r = 0; r < M; ++r) {
        for (int c = 0; c < K; ++c) {
            for (int el = 0; el < N; ++el) {
                mat_ans(r, c) += mat1(r, el) * mat2(el, c);
            }
        }
    }
    
    auto mat_prod = mat1 * mat2;

    EXPECT_EQ(mat_prod, mat_ans);
    EXPECT_EQ(mat_ans.nR(), mat_prod.nR());
    EXPECT_EQ(mat_ans.nC(), mat_prod.nC());

}

TEST(Matrix, BigProductDouble) {
    int M = 500;
    int N = 500;
    int K = 100;
    matrix<double> mat1(M, N);
    matrix<double> mat2(N, K);
    for (int r = 0; r < M; ++r) {
        for (int c = 0; c < N; ++c) {
            mat1(r, c) = (r + c);
        }
    }
    for (int r = 0; r < N; ++r) {
        for (int c = 0; c < K; ++c) {
            mat2(r, c) = (r + c);
        }
    }
    // Manual matrix mult
    matrix<double> mat_ans(M, K);
    mat_ans = 0;
    for (int r = 0; r < M; ++r) {
        for (int c = 0; c < K; ++c) {
            for (int el = 0; el < N; ++el) {
                mat_ans(r, c) += mat1(r, el) * mat2(el, c);
            }
        }
    }
    
    auto mat_prod = mat1 * mat2;

    EXPECT_EQ(mat_prod, mat_ans);
    EXPECT_EQ(mat_ans.nR(), mat_prod.nR());
    EXPECT_EQ(mat_ans.nC(), mat_prod.nC());

}




int main(int argc, char **argv) {
    std::cout << "Running test" << std::endl;
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}