#pragma once
#include "0.math.hpp"
#include "0.common.hpp"

/**
 * Full Black Scholes model
 *
 * @param S : Stock price
 * @param K : Strike price
 * @param sig : volatility
 * @param r : risk-free interest rate
 * @param tau : time to maturity
 *
 * @returns BlackScholesResult
 */
inline BlackScholesResult black_scholes(ftype S, ftype K, ftype sig, ftype r, ftype tau) {
    // calculate d1
    ftype A = ln(S/K) + tau * (r + p2(sig)/2.0f);
    ftype B = sig * sqrt(tau);
    ftype d1 = A / B;

    // calculate d2
    ftype d2 = d1 - sig * sqrt(tau);

    ftype call = S * N(d1) - N(d2) * K * e(-r * tau);
    ftype put = K * e( -r * tau) - S + call;

    return {.call=call, .put=put};
}
