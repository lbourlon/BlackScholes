#include "0.common.hpp"
#include "0.fast_math.hpp"

/**
 * Full Black Scholes model
 *
 * @param B : Stock price
 *
 * @returns BlackScholesResult
 */
BsOutput black_scholes_pipelined(BsInput in) {
    // calculate d1
    const double A = fm::log(in.S / in.K) + in.tau * (in.r + fm::pow_si(in.sig,2) / 2.0f);
    const double B = in.sig * fm::sqrt(in.tau);
    const double d1 = A / B;

    // calculate d2
    const double d2 = d1 - B;
    const double call = in.S * fm::N_zelen_severo_branchless(d1) - fm::N_zelen_severo_branchless(d2) * in.K * fm::exp(-in.r * in.tau);
    const double put = in.K * fm::exp(-in.r * in.tau) - in.S + call;

    return {call, put};
}
