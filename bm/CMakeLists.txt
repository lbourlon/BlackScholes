include(CTest)
add_executable(bm_black_scholes bm_black_scholes.cpp)

target_include_directories(bm_black_scholes
  PRIVATE
  ${CMAKE_SOURCE_DIR}/black_scholes/src
)
target_link_libraries(bm_black_scholes
  PRIVATE
  benchmark::benchmark
  black_scholes
)

add_compile_options(bm_black_scholes -O3 -march=native)
gtest_discover_tests(bm_black_scholes)
add_test( NAME bench COMMAND bm_black_scholes --benchmark_min_time=0.1)
