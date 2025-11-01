#include "benchmark/benchmark.h"
#include "1.baseline.hpp"

static void BM_black_scholes(benchmark::State& state) {
  for (auto _ : state) {
    black_scholes(1.1,2.2,3.3,4.4,5.5);
  }
}

BENCHMARK(BM_black_scholes);
BENCHMARK_MAIN();
