#include "0.common.hpp"
#include "0.fast_math.hpp"

inline double N_naive(double val) {
    double min = -10.0;
    double max = val;
    double dx = 0.001;
    auto f = [](const double &x) { return std::exp(-std::pow(x, 2) / 2.); };

    double sum = 0.0f;
    // no need for abs(f(x)) because e(x) is always positive
    for (double x = min; x < max; x += dx) {
        sum += dx * f(x);
    }
    return sum / sqrt(2.0 * pi);
}

/**
 * Full Black Scholes model
 *
 * @param B : Stock price
 *
 * @returns BlackScholesResult
 */
BsOutput black_scholes(BsInput in) {
    // calculate d1
    double A = std::log(in.S / in.K) + in.tau * (in.r + std::pow(in.sig,2) / 2.0f);
    double B = in.sig * sqrt(in.tau);
    double d1 = A / B;

    // calculate d2
    double d2 = d1 - in.sig * sqrt(in.tau);

    double call = in.S * N_naive(d1) - N_naive(d2) * in.K * exp(-in.r * in.tau);

    double put = in.K * std::exp(-in.r * in.tau) - in.S + call;

    return {.call = call, .put = put};
}
