#include <cmath>
#include <print>
#include "0.math.hpp"
#include "0.common.hpp"
#include "1.baseline.hpp"

int main() {
    ftype S = 300.0;
    ftype K = 250.0;
    ftype tau = 1;
    ftype sig = 0.15;
    ftype r = 0.03;

    BlackScholesResult res = black_scholes(S, K, sig, r, tau);
    std::println("Call: {}/58.81977, Put: {}/1.43116", res.call, res.put);

    S = 821.0; K = 781.0; tau = 3; sig = 0.33; r = 0.049;
    res = black_scholes(S, K, sig, r, tau);
    std::println("Call: {}/58.81977, Put: {}/1.43116", res.call, res.put);
}
