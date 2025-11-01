#include "1.baseline.hpp"
#include "2.bs_black76.hpp"
#include "benchmark/benchmark.h"
#include <random>

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

static void BM_black_scholes(benchmark::State &state) {
    for (auto _ : state) {
        state.PauseTiming();
        const BsInput input = get_random_bs_input();
        state.ResumeTiming();
        black_scholes(input);
    }
}

static void BM_black_scholes_alt(benchmark::State &state) {
    for (auto _ : state) {
        state.PauseTiming();
        const BsInput input = get_random_bs_input();
        state.ResumeTiming();
        black_scholes_alt(input);
    }
}

BENCHMARK(BM_black_scholes);
BENCHMARK(BM_black_scholes_alt);
BENCHMARK_MAIN();
