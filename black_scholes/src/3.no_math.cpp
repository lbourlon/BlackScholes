#include "0.common.hpp"
#include <cmath>
#include <numbers>
#include <algorithm>

using std::numbers::pi;
// using std::numbers::e;
inline auto ln(auto x) { return std::log(x); }
inline auto sqrt(auto x) { return std::sqrt(x); }
inline auto exp(auto x) {
    return std::exp(x); }

inline auto p2(auto x) { return x * x; }

inline double N(double val) {
    double min = -10.0;
    double max = val;
    double dx = 0.001;
    auto f = [](const double &x) { return exp(-p2(x) / 2.); };

    double sum = 0.0f;
    for (double x = min; x < max; x += dx) { sum += dx * f(x); }
    return sum / sqrt(2.0 * pi);
}

/**
 * Full Black Scholes model
 *
 * @param B : Stock price
 *
 * @returns BlackScholesResult
 */
BsOutput black_scholes_no_math(BsInput in) {
    // calculate d1
    double A = ln(in.S / in.K) + in.tau * (in.r + p2(in.sig) / 2.0f);
    double B = in.sig * sqrt(in.tau);
    double d1 = A / B;

    // calculate d2
    double d2 = d1 - in.sig * sqrt(in.tau);

    double call = in.S * N(d1) - N(d2) * in.K * exp(-in.r * in.tau);
    double put = in.K * exp(-in.r * in.tau) - in.S + call;

    return {.call = call, .put = put};
}
