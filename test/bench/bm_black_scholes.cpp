#include "1.baseline.hpp"
#include "benchmark/benchmark.h"

static void BM_black_scholes(benchmark::State &state) {
    const BsInput input = {
        .S = 300.0, .K = 250.0, .tau = 1, .sig = 0.15, .r = 0.03,
    };

    for (auto _ : state) {
        black_scholes(input);
    }
}

BENCHMARK(BM_black_scholes);
BENCHMARK_MAIN();
