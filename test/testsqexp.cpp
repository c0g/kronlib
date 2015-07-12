//
// Created by Thomas Nickson on 05/07/2015.
//

#include "matrix.h"
#include "kernel.h"
#include "iostream"
#include "gtest/gtest.h"

TEST(SqExp, GetHyp)
{
    sqexp_hyp<float> hyp{0.1, 0.2};
    EXPECT_FLOAT_EQ(hyp.get_real(0), std::exp(0.1));
    EXPECT_FLOAT_EQ(hyp.get_real(1), std::exp(0.2));
    EXPECT_FLOAT_EQ(hyp.get_log(0), 0.1);
    EXPECT_FLOAT_EQ(hyp.get_log(1), 0.2);
    EXPECT_FLOAT_EQ(hyp.get_real(0)*hyp.get_real(0), hyp.lengthscale2());
    EXPECT_FLOAT_EQ(hyp.get_real(1)*hyp.get_real(1), hyp.outputscale2());
}

TEST(SqExp, SetHyp)
{
    sqexp_hyp<float> hyp{0.1, 0.1};
    hyp.set_real(0, 10);
    EXPECT_FLOAT_EQ(hyp.get_real(0), 10);
    EXPECT_FLOAT_EQ(hyp.get_log(0), std::log(10));

    hyp.set_real(1, 11);
    EXPECT_FLOAT_EQ(hyp.get_real(1), 11);
    EXPECT_FLOAT_EQ(hyp.get_log(1), std::log(11));
}

TEST(SqExp, K)
{
    double log_lengthscale = 0.3;
    double log_outputscale = 0.2;
    sqexp_hyp<double> hyp{log_lengthscale, log_outputscale};

    sqexp1d<double> kernel;

    Matrix<double> x(10, 1);
    x = 1, 2, 3, 4, 5, 6, 7, 8, 9, 10;

    auto k = kernel(hyp, x, x);

    double l2 = std::exp( 2 * log_lengthscale);
    double s2 = std::exp( 2 * log_outputscale);
    for (int r = 0; r < 10; ++r) {
        for (int c = 0; c < 10; ++c) {
            float val = s2 * std::exp( - std::pow(x(r, 0) - x(c, 0), 2) / (2 * l2));
            EXPECT_TRUE(std::abs(val - k(r, c)) < 1e-5);
        }
    }
}

TEST(SqExp, Grads)
{
    Matrix<double> x{10, 1};
    x = 1,2,3,4,5,6,7,8,9,10;

    const double delta = 1e-8;
    const double onepdelt = 1+delta;
    sqexp_hyp<double> hyp{1, 1};
    sqexp_hyp<double> hyp_ll{onepdelt, 1};
    sqexp_hyp<double> hyp_ls{1, onepdelt};

    sqexp1d<double> k;

    auto normal = k(hyp, x, x);
    auto dll = k(hyp_ll, x, x);
    auto dls = k(hyp_ls, x, x);

    auto grad1 = (dll - normal) / delta;
    auto grad2 = (dls - normal) / delta;
    auto derivs = k.dhyp(hyp, x, x);

    for (int r = 0; r < 10; ++r) {
        for (int c = 0; c < 10; ++c) {
            EXPECT_TRUE(std::abs(derivs[0](r,c) - grad1(r, c)) < 1e-5);
            EXPECT_TRUE(std::abs(derivs[1](r,c) - grad2(r, c)) < 1e-5);
        }
    }
}

int main(int argc, char **argv)
{
    std::cout << "Running test" << std::endl;
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}