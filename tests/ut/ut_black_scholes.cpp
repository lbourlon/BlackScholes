#include "0.common.hpp"
#include "black_scholes.hpp"
#include "gtest/gtest.h"
#include <print>
#include <random>

using namespace testing;

using bs_implem = std::function<BsOutput(BsInput)>;
class TestRealValues : public TestWithParam<bs_implem> {};
class RandomizedControl : public TestWithParam<bs_implem> {};

BsInput get_random_bs_input() {
    static std::mt19937 rng(298739);
    std::uniform_real_distribution<double> price(1, 10000);
    std::uniform_real_distribution<double> percentage(0.001, 1); // 0.1% -> 50%
    std::uniform_real_distribution<double> variation(-5, 5);
    std::uniform_real_distribution<double> tau(0, 20);

    double spot_price = price(rng);
    double strike_price = spot_price + variation(rng) * spot_price / 2.0;
    return {
        strike_price, spot_price, tau(rng), percentage(rng), percentage(rng),
    };
}

void bs_expect_absolute_err(BsOutput expected, BsOutput computed) {
    auto absolute_err = 0.01;
    EXPECT_NEAR(expected.call, computed.call, absolute_err);
    EXPECT_NEAR(expected.put, computed.put, absolute_err);
}

TEST_P(TestRealValues, TestWithRealValues) {
    auto func = GetParam();
    // Spot price, strike price, maturity, volatility, risk free rate
    bs_expect_absolute_err(func({300.0, 250.0, 1, 0.15, 0.03}), {58.8198, 1.4312});
    bs_expect_absolute_err(func({9877, 2871, 1.9, 0.21, 0.0411}), {7221.6670, 0.0008});
    bs_expect_absolute_err(func({821, 781, 3, 0.33, 0.049}), {251.1488, 104.3814});
    bs_expect_absolute_err(func({287, 190, 299, 0.13, .029}), {286.9674, 0.0});
    bs_expect_absolute_err(func({912, 821, 1.2, .1263, .109}), {193.6817, 2.0215});
    bs_expect_absolute_err(func({2, 1.2, 4.98, .10, 0.0005}), {0.8042, 0.0013});
}

TEST_P(RandomizedControl, WithControlFunc) {
    auto func = GetParam();
    int fail_count = 0;
    int total_count = 10000;
    for (size_t i = 0; i < total_count; ++i) {
        BsInput input = get_random_bs_input();
        BsOutput expected = black_scholes_numerical_approx_N(input);
        BsOutput computed = func(input);

        // Using the spread as a testing metric makes intuitive sense in my head
        // TODO : look at what real papers do to test these
        double expected_spread = std::abs(expected.call - expected.put);
        double computed_spread = std::abs(computed.call - computed.put);

        double error = std::abs(expected_spread - computed_spread);
        double absolute_err = 1.1e-3;

        if (error > absolute_err) {
            std::println("============= BAD SPREAD ==============");
            std::println("Iteration {} errored:", i);
            std::println("expected({}) != computed({})", computed_spread,
                         expected_spread);
            std::println("error({}) > max_error({})", error, absolute_err);
            std::println("{}\n-----------------------------", input);
            std::println("Expected output {}", expected);
            std::println("Computed output {}\n", computed);
            std::println("=============-------------==============");
            fail_count += 1;
        }
    }
    if (fail_count != 0) {
        std::println("FAILED {}/{} times", fail_count, total_count);
        EXPECT_TRUE(fail_count == 0);
    }
}

const auto all_implementations = Values( //
    black_scholes,                       //
    black_scholes_alt,                   //
    black_scholes_numerical_approx_N,    //
    black_scholes_fast_math,             //
    black_scholes_pipelined,             //
    black_scholes_simd_oriented          //
);

const auto fast_implementations = Values( //
    black_scholes_numerical_approx_N,     // spread is accurate up to 1e-12
    black_scholes_fast_math,              //
    black_scholes_pipelined,              //
    black_scholes_simd_oriented           //
);

INSTANTIATE_TEST_SUITE_P(Bs1, TestRealValues, all_implementations);
INSTANTIATE_TEST_SUITE_P(Bs2, RandomizedControl, fast_implementations);

int main(int argc, char **argv) {
    InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
