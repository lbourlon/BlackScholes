#include <cmath>
#include <print>

#include "0.common.hpp"
#include "1.baseline.hpp"

int main() {
    const BsInput input = {
        .S = 300.0, .K = 250.0, .tau = 1, .sig = 0.15, .r = 0.03,
    };

    BsOutput res = black_scholes(input);
    std::println("Call: {}/58.81977, Put: {}/1.43116", res.call, res.put);

    const BsInput input2 = {.S = 821.0, .K = 781.0, .tau = 3, .sig = 0.33, .r = 0.049};
    res = black_scholes(input2);
    std::println("Call: {}/58.81977, Put: {}/1.43116", res.call, res.put);
}
