#include <random>
#include "black_scholes.hpp"
#include "benchmark/benchmark.h"

// std::array<BsInput, 512>
BsInput get_random_bs_input() {

    std::uniform_real_distribution<double> price_dist(0, 10'000);
    std::uniform_real_distribution<double> percent(0, 1);
    std::uniform_real_distribution<double> time(1, 25);
    std::default_random_engine re;

    return {
        .S = price_dist(re),
        .K = price_dist(re),
        .tau = time(re),
        .sig = percent(re),
        .r = percent(re),
    };
}

static void BM_black_scholes(benchmark::State &state, std::function<BsOutput(BsInput)> bs_implementation) {
    for (auto _ : state) {
        state.PauseTiming();
        const BsInput input = get_random_bs_input();
        state.ResumeTiming();
        bs_implementation(input);
    }
}


BENCHMARK_CAPTURE(BM_black_scholes, baseline, black_scholes);
BENCHMARK_CAPTURE(BM_black_scholes, bs_black76, black_scholes_alt);
BENCHMARK_CAPTURE(BM_black_scholes, no_math, black_scholes_no_math);
BENCHMARK_MAIN();
