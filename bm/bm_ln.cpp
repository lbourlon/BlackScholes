#include <random>
#include <cmath>
#include "benchmark/benchmark.h"

#include "0.fast_math.hpp"

static void BM_ln(benchmark::State &state, std::function<double(double)> func_impl) {
    std::uniform_real_distribution<double> rand_gen(-10'000, 10'000);
    std::default_random_engine re;

    double x = rand_gen(re);
    for (auto _ : state) {
        benchmark::DoNotOptimize(func_impl(x));
    }
}

BENCHMARK_CAPTURE(BM_ln, std_log, fm::ln_std);
BENCHMARK_CAPTURE(BM_ln, ln_sqrt2, fm::ln_sqrt2);
BENCHMARK_CAPTURE(BM_ln, ln_sqrt2_twice, fm::ln_sqrt2_twice);
BENCHMARK_CAPTURE(BM_ln, ln_mercator, fm::ln_mercator);
BENCHMARK_CAPTURE(BM_ln, ln_power_of_2, fm::ln_power_of_2);
BENCHMARK_CAPTURE(BM_ln, ln_newton_halley, fm::ln_newton_halley);
BENCHMARK_CAPTURE(BM_ln, ln_newton_halley_fast, fm::ln_newton_halley_fast);
BENCHMARK_MAIN();
