#include <cmath>

#include "0.common.hpp"
#include "0.fast_math.hpp"
#include "black_scholes.hpp"

/**
 * Full Black Scholes model
 *
 * @param B : Stock price
 *
 * @returns BlackScholesResult
 */
BsOutput black_scholes_numerical_approx_N(BsInput in) {
    // calculate d1
    const double A = std::log(in.S / in.K) + in.tau * (in.r + std::pow(in.sig, 2) / 2.0f);
    const double B = in.sig * sqrt(in.tau);
    const double d1 = A / B;

    // calculate d2
    double d2 = d1 - in.sig * sqrt(in.tau);

    double call = in.S * fm::N_zelen_severo(d1) -
                  fm::N_zelen_severo(d2) * in.K * exp(-in.r * in.tau);

    double put = in.K * std::exp(-in.r * in.tau) - in.S + call;

    return {call, put};
}
