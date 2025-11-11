
add_compile_options(-std=c++23 -O3 -march=native -mavx2 -mavx -msse2 -msse -g)
add_library(black_scholes
   1.baseline.cpp
   2.bs_numerical_N.cpp
   3.bs_fast_math.cpp
   4.bs_pipelined.cpp
   5.bs_less_operations.cpp
)

add_compile_options(-O3 -march=native -mavx2 -mavx)
add_executable(bs main.cpp)
target_link_libraries(bs PRIVATE black_scholes)
