#include <random>
#include <cmath>
#include "benchmark/benchmark.h"

#include "0.fast_math.hpp"

static void BM_sqrt(benchmark::State &state, std::function<double(double)> func_impl) {
    std::uniform_real_distribution<double> rand_gen(-10'000, 10'000);
    std::default_random_engine re;

    double x = rand_gen(re);
    for (auto _ : state) {
        benchmark::DoNotOptimize(func_impl(x));
    }
}

BENCHMARK_CAPTURE(BM_sqrt, sqrt_std, fm::sqrt_std);
BENCHMARK_CAPTURE(BM_sqrt, sqrt_exp_ln, fm::sqrt_exp_ln);
BENCHMARK_CAPTURE(BM_sqrt, sqrt_pow_ln, fm::sqrt_pow_ln);
BENCHMARK_CAPTURE(BM_sqrt, sqrt_newton_raphson, fm::sqrt_newton_raphson);
BENCHMARK_CAPTURE(BM_sqrt, sqrt_kahan_ng, fm::sqrt_kahan_ng);
BENCHMARK_MAIN();
