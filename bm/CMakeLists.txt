include(CTest)

function(create_benchmark target)
  add_executable(${target} ${target}.cpp)

  target_include_directories(${target} PRIVATE
    ${CMAKE_SOURCE_DIR}/black_scholes/src
  )
  target_link_libraries(${target}
    PRIVATE
    benchmark::benchmark
    black_scholes
  )

  add_compile_options(-O3 -march=native -mavx2 -mavx)
  gtest_discover_tests(${target})
  add_test( NAME bench-${target} COMMAND ${target} --benchmark_min_time=1
  --benchmark_enable_random_interleaving=true)
endfunction()

create_benchmark(bm_black_scholes)
create_benchmark(bm_exponential)
