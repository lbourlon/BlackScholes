#include "gtest/gtest.h"
#include "1.baseline.hpp"

TEST(FactorialTest, HandlesZeroInput) {
    ftype S = 300.0;
    ftype K = 250.0;
    ftype tau = 1;
    ftype sig = 0.15;
    ftype r = 0.03;

    BlackScholesResult res = black_scholes(S, K, sig, r, tau);
    BlackScholesResult expected = {
        .call = 58.81977,
        .put = 1.43116,
    };

    EXPECT_NEAR(expected.call, res.call, 0.00001);
    EXPECT_NEAR(expected.put, res.put, 0.00001);
}

int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
