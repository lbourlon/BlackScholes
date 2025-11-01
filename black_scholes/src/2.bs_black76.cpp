#include "0.common.hpp"
#include "0.math.hpp"

/**
 * Full Black Scholes Alternative Model
 * (Black 76)
 *
 * @returns BlackScholesResult
 */
BsOutput black_scholes_alt(BsInput in) {
    double D = exp(-in.r * in.tau); // discount factor
    double F = in.S / D;          // forward price factor

    // calculate d1
    double A = ln(F / in.K) + (in.tau * p2(in.sig)) / 2.0f;
    double B = in.sig * sqrt(in.tau);
    double d1 = A / B;

    // calculate d2
    double d2 = d1 - in.sig * sqrt(in.tau);

    // call
    double call = D * (N(d1) * F - N(d2) * in.K);

    // put
    double put = in.K * exp(-in.r * in.tau) - in.S + call;

    return {.call = call, .put = put};
}
