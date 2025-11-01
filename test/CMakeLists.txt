include(FetchContent)

FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG        52eb8108c5bdec04579160ae17225d66034bd723 # v1.17.0
  GIT_SHALLOW    TRUE
  TLS_VERIFY ${CMAKE_TLS_VERIFY}
)
set(BUILD_GTEST OFF)
set(BUILD_GMOCK OFF)
set(INSTALL_GTEST OFF)

FetchContent_Declare(
  benchmark
  GIT_REPOSITORY https://github.com/google/benchmark.git
  GIT_TAG        eddb0241389718a23a42db6af5f0164b6e0139af # tag v1.9.4
  GIT_SHALLOW    TRUE
  TLS_VERIFY ${CMAKE_TLS_VERIFY}
)
set(BENCHMARK_ENABLE_TESTING OFF)
set(BENCHMARK_ENABLE_GTEST_TESTS OFF)
set(BENCHMARK_ENABLE_EXAMPLES OFF)
set(BENCHMARK_DOWNLOAD_DEPENDENCIES OFF)

FetchContent_MakeAvailable(googletest)
FetchContent_MakeAvailable(benchmark)

option(BUILD_TESTING "Build the tests" ON)
add_subdirectory(unit)
add_subdirectory(bench)

