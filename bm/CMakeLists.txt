add_executable(bm_black_scholes bm_black_scholes.cpp)

target_include_directories(bm_black_scholes
  PRIVATE
  ${CMAKE_SOURCE_DIR}/src
)
target_link_libraries(bm_black_scholes
  PRIVATE
  benchmark::benchmark
)

add_compile_options(bm_black_scholes -O3 -march=native)
