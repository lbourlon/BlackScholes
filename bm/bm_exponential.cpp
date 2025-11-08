#include <random>
#include "benchmark/benchmark.h"

#include "0.fast_math.hpp"

static void BM_exponential(benchmark::State &state, std::function<double(double)> func_impl) {
    std::uniform_real_distribution<double> rand_gen(-10'000, 10'000);
    std::default_random_engine re;

    for (auto _ : state) {

        double x = rand_gen(re); // <- bias upwards but improves reproducibility
        benchmark::DoNotOptimize(func_impl(x));
    }
}

// constexpr size_t iter = 2'000;


BENCHMARK_CAPTURE(BM_exponential, std_exp, fm::exp_std);
BENCHMARK_CAPTURE(BM_exponential, exp_reduced, fm::exp_reduced);
BENCHMARK_CAPTURE(BM_exponential, exp_reduced_fast, fm::exp_reduced_fast);
BENCHMARK_CAPTURE(BM_exponential, exp_reduced_const_array, fm::exp_reduced_const_array);
BENCHMARK_MAIN();
