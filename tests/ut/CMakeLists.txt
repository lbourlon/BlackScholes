include(GoogleTest)

add_executable(ut_black_scholes ut_black_scholes.cpp)

target_include_directories(ut_black_scholes
  PRIVATE
  ${CMAKE_SOURCE_DIR}/src
)
target_link_libraries(ut_black_scholes
  PRIVATE
  gtest
  black_scholes
)

# gtest_discover_tests(ut_black_scholes)
