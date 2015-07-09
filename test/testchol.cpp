//
// Created by Thomas Nickson on 05/07/2015.
//

#include "matrix.h"
#include "cholesky.h"
#include "iostream"
#include "gtest/gtest.h"

TEST(Cholesky, CalcAndEquality)
{
    matrix<float> mat(3, 3);
    mat = 1, 1, 0.5, 1, 10, 3.5, 0.5, 3.5, 2.25;
    Cholesky<float> chol(mat);
    matrix<float> mat_ans(3, 3);
    mat_ans = 1, 0, 0, 1, 3, 0, 0.5, 1, 1;
    EXPECT_EQ(mat_ans, chol);
}
TEST(Cholesky, Inverse)
{
    matrix<float> mat(3, 3);
    mat = 1, 1, 0.5, 1, 10, 3.5, 0.5, 3.5, 2.25;
    Cholesky<float> chol(mat);
    auto mat_inv = chol.inv();
    matrix<float> mat_ans(3, 3);
    mat_ans = 1.13888889, -0.05555556, -0.16666667, -0.05555556,  0.22222222,
    -0.33333333, -0.16666667, -0.33333333,  1.;

    for (int r = 0; r < mat.nR(); ++r) {
        for (int c = 0; c < mat.nC(); ++c) {
            EXPECT_TRUE(std::abs(mat_ans(r, c) - mat_inv(r, c)) < 1e-5);
        }
    }
}


TEST(Cholesky, Solve)
{
    matrix<float> mat(3, 3);
    mat = 1, 1, 0.5, 1, 10, 3.5, 0.5, 3.5, 2.25;
    Cholesky<float> chol(mat);

    matrix<float> solve_with(3, 4);
    solve_with = 1, 1, 1, 1, 2, 3, 1, 3, 5, 0, 3, 6;

    matrix<float> ans(3, 4);
    ans = 0.19444444,  0.97222222,  0.58333333, -0.02777778, -1.27777778,
    0.61111111, -0.83333333, -1.38888889,  4.16666667, -1.16666667,
    2.5       ,  4.83333333;

    auto solved = chol.solve(solve_with);

    for (int r = 0; r < solved.nR(); ++r) {
        for (int c = 0; c < solved.nC(); ++c) {
            EXPECT_TRUE(std::abs(ans(r, c) - solved(r, c)) < 1e-5);
        }
    }

}

TEST(Cholesky, LogDet)
{
    matrix<float> mat(3, 3);
    mat = 1, 1, 0.5, 1, 10, 3.5, 0.5, 3.5, 2.25;
    Cholesky<float> chol(mat);

    float ld = 2.1972245773362196;

    EXPECT_FLOAT_EQ(ld, chol.logdet());

}

int main(int argc, char **argv)
{
    std::cout << "Running test" << std::endl;
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}