//
// Created by Thomas Nickson on 05/07/2015.
//

#include "matrix.h"
#include "cholesky.h"
#include "iostream"
#include "gtest/gtest.h"

TEST(Cholesky, CalcAndEquality)
{
    Matrix<host<float>> mat(3, 3);
    mat = 1, 1, 0.5, 1, 10, 3.5, 0.5, 3.5, 2.25;
    Cholesky<host<float>> chol(mat);
    Matrix<host<float>> mat_ans(3, 3);
    mat_ans = 1, 0, 0, 
	      1, 3, 0, 
              0.5, 1, 1;
    EXPECT_TRUE(mat_ans == chol);
}
/*
TEST(Cholesky, Inverse)
{
    Matrix<host<float>> mat(3, 3);
    mat = 1, 1, 0.5, 1, 10, 3.5, 0.5, 3.5, 2.25;
    Cholesky<host<float>> chol(mat);
    auto mat_inv = chol.inv();
    Matrix<host<float>> mat_ans(3, 3);
    mat_ans = 1.13888889, -0.05555556, -0.16666667, -0.05555556,  0.22222222,
    -0.33333333, -0.16666667, -0.33333333,  1.;
    for (int r = 0; r < mat.nR(); ++r) {
        for (int c = 0; c < mat.nC(); ++c) {
            EXPECT_TRUE(std::abs(mat_ans(r, c) - mat_inv(r, c)) < 1e-5);
        }
    }
}
*/

TEST(Cholesky, Solve)
{
    Matrix<host<float>> mat(3, 3);
    mat = 1, 1, 0.5, 1, 10, 3.5, 0.5, 3.5, 2.25;
    Cholesky<host<float>> chol(mat);

    Matrix<host<float>> solve_with(3, 4);
    solve_with = 1, 1, 1, 	
		1, 2, 3, 
		1, 3, 5,
		0, 3, 6;

    Matrix<host<float>> ans(3, 4);
    ans = 0.19444444,  0.97222222,  0.58333333,
 	-0.02777778, -1.27777778,    0.61111111,
 	-0.83333333, -1.38888889,  4.16666667,
 	-1.16666667,    2.5       ,  4.83333333;
    Matrix<host<float>> solved = chol.solve(solve_with);

    for (int r = 0; r < solved.nR(); ++r) {
        for (int c = 0; c < solved.nC(); ++c) {
            EXPECT_TRUE(std::abs(ans(r, c) - solved(r, c)) < 1e-5);
        }
    }
}

TEST(Cholesky, Solve2)
{
    Matrix<host<double>> mat(4, 4);
    mat = 125,  195,  275,  429,  195,  305,  429,  671,  275,  429,  625,
        975,  429,  671,  975, 1525;
    Cholesky<host<double>> chol(mat);

    Matrix<host<double>> solve_with(4, 1);
    solve_with = 0, 1, 2, 3;

    Matrix<host<double>> ans(4, 1);
    ans = -64.375,  41.125,  28.375, -18.125;

    auto solved = chol.solve(solve_with);

    for (int r = 0; r < solved.nR(); ++r) {
        for (int c = 0; c < solved.nC(); ++c) {
            EXPECT_TRUE(std::abs(ans(r, c) - solved(r, c)) < 1e-5);
        }
    }

}
TEST(Cholesky, LogDet)
{
    Matrix<host<float>> mat(3, 3);
    mat = 1, 1, 0.5, 1, 10, 3.5, 0.5, 3.5, 2.25;
    Cholesky<host<float>> chol(mat);

    float ld = 2.1972245773362196;

    EXPECT_FLOAT_EQ(ld, chol.logdet());

}

int main(int argc, char **argv)
{
    std::cout << "Running test" << std::endl;
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
