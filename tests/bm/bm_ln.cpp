#include <random>
#include <cmath>
#include "benchmark/benchmark.h"

#include "0.fast_math.hpp"

static void BM_log(benchmark::State &state, std::function<double(double)> func_impl) {
    std::uniform_real_distribution<double> rand_gen(-10'000, 10'000);
    std::default_random_engine re;

    double x = rand_gen(re);
    for (auto _ : state) {
        benchmark::DoNotOptimize(func_impl(x));
    }
}

BENCHMARK_CAPTURE(BM_log, std_log, fm::log_std);
BENCHMARK_CAPTURE(BM_log, log_sqrt2, fm::log_sqrt2);
BENCHMARK_CAPTURE(BM_log, log_sqrt2_twice, fm::log_sqrt2_twice);
BENCHMARK_CAPTURE(BM_log, log_mercator, fm::log_mercator);
BENCHMARK_CAPTURE(BM_log, log_power_of_2, fm::log_power_of_2);
BENCHMARK_CAPTURE(BM_log, log_newton_halley, fm::log_newton_halley);
BENCHMARK_CAPTURE(BM_log, log_newton_halley_fast, fm::log_newton_halley_fast);
BENCHMARK_MAIN();
