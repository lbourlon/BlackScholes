cmake_minimum_required(VERSION 3.22)

if (CMAKE_VERSION VERSION_GREATER_EQUAL "3.24.0")
   cmake_policy(SET CMP0135 NEW)
endif()

project(BlackScholes)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

find_program(CCACHE_FOUND ccache)
if(CCACHE_FOUND)
   set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE ccache)
endif()

set(CMAKE_TLS_VERIFY TRUE)
set(CMAKE_BUILD_TYPE Release)

add_subdirectory(black_scholes)

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
if (${BUILD_TESTING})
   add_subdirectory(ut)
   add_subdirectory(bm)
endif()
