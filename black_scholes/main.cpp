#include <bitset>
#include <cmath>
#include <print>

#include "src/black_scholes.hpp"
#include "src/0.fast_math.hpp"

int main() {
    const BsInput input = { 300.0, 250.0, 1, 0.15, 0.03 };
    BsOutput res = black_scholes(input);
    std::println("Call: {}/58.81977, Put: {}/1.43116", res.call, res.put);
}
