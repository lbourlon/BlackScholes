#include "1.baseline.hpp"
#include "gtest/gtest.h"

TEST(FactorialTest, HandlesZeroInput) {
    const BsInput input = {
        .S = 300.0, .K = 250.0, .tau = 1, .sig = 0.15, .r = 0.03,
    };

    BsOutput res = black_scholes(input);
    BsOutput expected = {
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
