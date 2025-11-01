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

add_compile_options(-O3 -march=native)
add_executable(blackscholes src/main.cpp)

add_subdirectory(test)
