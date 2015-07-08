//
// Created by Thomas Nickson on 05/07/2015.
//

#include "matrix.h"
#include "kronecker_matrix.h"
#include "iostream"
#include "gtest/gtest.h"

TEST(KronMatrix, RowsAndColumn)
{
    matrix<float> mat1(2, 3);
    matrix<float> mat2(2, 3);
    kronecker_matrix<float> kron_mat;
    kron_mat.push_matrix(mat1);
    kron_mat.push_matrix(mat2);
    EXPECT_EQ(mat1.nR() * mat2.nR(), kron_mat.nR());
    EXPECT_EQ(mat1.nC() * mat2.nC(), kron_mat.nC());
}

TEST(KronMatrix, KronFull2)
{
    matrix<float> mat1(2, 2);
    mat1 = 1, 2, 3, 4;
    matrix<float> mat2(2, 2);
    mat2 = 3, 4, 5, 6;

    matrix<float> mat_ans(4, 4);
    mat_ans = 3,  4,  6,  8,  5,  6, 10, 12,  9, 12, 12, 16, 15, 18, 20, 24;

    kronecker_matrix<float> kron_mat;
    kron_mat.push_matrix(mat1);
    kron_mat.push_matrix(mat2);

    EXPECT_EQ(kron_mat.full(), mat_ans);
}

TEST(KronMatrix, KronFull3)
{
    matrix<float> mat1(2, 2);
    mat1 = 1, 2, 3, 4;
    matrix<float> mat2(2, 2);
    mat2 = 3, 4, 5, 6;
    matrix<float> mat3(2, 3);
    mat3 = 3, 4, 5, 6, 7, 8;

    matrix<float> mat_ans(8, 12);
    mat_ans = 9,  12,  15,  12,  16,  20,  18,  24,  30,  24,  32,  40,  18,
    21,  24,  24,  28,  32,  36,  42,  48,  48,  56,  64,  15,  20,
    25,  18,  24,  30,  30,  40,  50,  36,  48,  60,  30,  35,  40,
    36,  42,  48,  60,  70,  80,  72,  84,  96,  27,  36,  45,  36,
    48,  60,  36,  48,  60,  48,  64,  80,  54,  63,  72,  72,  84,
    96,  72,  84,  96,  96, 112, 128,  45,  60,  75,  54,  72,  90,
    60,  80, 100,  72,  96, 120,  90, 105, 120, 108, 126, 144, 120,
    140, 160, 144, 168, 192;

    kronecker_matrix<float> kron_mat;
    kron_mat.push_matrix(mat1);
    kron_mat.push_matrix(mat2);
    kron_mat.push_matrix(mat3);

    EXPECT_EQ(kron_mat.full(), mat_ans);
}

TEST(KronMatrix, KronDotKron)
{
    matrix<float> mat1a(2, 2);
    mat1a = 1, 2, 3, 4;
    matrix<float> mat2a(2, 2);
    mat2a = 3, 4, 5, 6;

    kronecker_matrix<float> kron_mat1;
    kron_mat1.push_matrix(mat1a);
    kron_mat1.push_matrix(mat2a);

    matrix<float> mat1b(2, 2);
    mat1b = 3, 4, 5, 6;
    matrix<float> mat2b(2, 2);
    mat2b = 4, 5, 6, 7;

    kronecker_matrix<float> kron_mat2;
    kron_mat2.push_matrix(mat1b);
    kron_mat2.push_matrix(mat2b);

    auto kron_mat = kron_mat1 * kron_mat2;

    matrix<float> mat_ans(4, 4);
    mat_ans = 468,  559,  576,  688,  728,  871,  896, 1072, 1044, 1247, 1296,
    1548, 1624, 1943, 2016, 2412;

    EXPECT_EQ(kron_mat.full(), mat_ans);
}


TEST(KronMatrix, TestScalarMult)
{
    matrix<float> mat1a(2, 2);
    mat1a = 1, 2, 3, 4;
    matrix<float> mat1b(2, 2);
    mat1b = 2, 4, 6, 8;

    matrix<float> mat2a(2, 2);
    mat2a = 2, 4, 6, 8;
    matrix<float> mat2b(2, 2);
    mat2b = 4, 8, 12, 16;

    kronecker_matrix<float> mat1;
    mat1.push_matrix(mat1a);
    mat1.push_matrix(mat1b);

    kronecker_matrix<float> mat2;
    mat2.push_matrix(mat2a);
    mat2.push_matrix(mat2b);
    mat1 *= 4;

    EXPECT_EQ(mat1, mat2);
}

TEST(KronMatrix, TestNegativeScalarMult)
{
    matrix<float> mat1a(2, 2);
    mat1a = 1, 2, 3, 4;
    matrix<float> mat1b(2, 2);
    mat1b = 2, 4, 6, 8;

    matrix<float> mat2a(2, 2);
    mat2a = -2, -4, -6, -8;
    matrix<float> mat2b(2, 2);
    mat2b = -4, -8, -12, -16;

    kronecker_matrix<float> mat1;
    mat1.push_matrix(mat1a);
    mat1.push_matrix(mat1b);

    kronecker_matrix<float> mat2;
    mat2.push_matrix(mat2a);
    mat2.push_matrix(mat2b);
    mat1 *= -4;

    EXPECT_EQ(mat1, mat2);
}

int main(int argc, char **argv)
{
    std::cout << "Running test" << std::endl;
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}