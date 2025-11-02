
add_compile_options(-std=c++23 -O3 -march=native -mavx2 -mavx -msse2 -msse)
add_library(black_scholes
   src/1.baseline.cpp
   src/2.bs_black76.cpp
   src/3.bs_fast_math.cpp
)

add_compile_options(-O3 -march=native -mavx2 -mavx)
add_executable(bs_bin main.cpp)

target_link_libraries(bs_bin PRIVATE black_scholes)
