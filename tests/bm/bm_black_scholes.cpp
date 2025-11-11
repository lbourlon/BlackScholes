#include <random>
#include "black_scholes.hpp"
#include "benchmark/benchmark.h"

BsInput get_random_bs_input() {

    std::uniform_real_distribution<double> price_dist(0, 10'000);
    std::uniform_real_distribution<double> percent(0, 1);
    std::uniform_real_distribution<double> time(1, 25);
    std::default_random_engine re;

    return {price_dist(re),
        price_dist(re),
        time(re),
        percent(re),
        percent(re),
    };
}

static void BM_black_scholes(benchmark::State &state, std::function<BsOutput(BsInput)> bs_implementation) {
    const BsInput input = get_random_bs_input();
    for (auto _ : state) {
        bs_implementation(input);
    }
}


constexpr size_t iter = 2000000;
BENCHMARK_CAPTURE(BM_black_scholes, baseline, black_scholes);
BENCHMARK_CAPTURE(BM_black_scholes, nuperical_approx_N, black_scholes_numerical_approx_N)->Iterations(iter);
BENCHMARK_CAPTURE(BM_black_scholes, fast_math, black_scholes_fast_math)->Iterations(iter);
BENCHMARK_CAPTURE(BM_black_scholes, pipelined, black_scholes_pipelined)->Iterations(iter);
BENCHMARK_CAPTURE(BM_black_scholes, simd_oriented, black_scholes_simd_oriented)->Iterations(iter);
BENCHMARK_MAIN();
