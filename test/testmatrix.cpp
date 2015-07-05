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




int main(int argc, char **argv) {
    std::cout << "Running test" << std::endl;
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}