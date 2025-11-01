add_library(black_scholes
   src/1.baseline.cpp
   src/2.bs_black76.cpp
   src/3.no_math.cpp
)

add_compile_options(-O3 -march=native)
add_executable(bs_bin main.cpp)

target_link_libraries(bs_bin PRIVATE black_scholes)
