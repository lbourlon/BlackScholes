#include "0.common.hpp"
#include "0.math.hpp"

/**
 * Full Black Scholes model
 *
 * @param B : Stock price
 *
 * @returns BlackScholesResult
 */
BsOutput black_scholes(BsInput in) {
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
