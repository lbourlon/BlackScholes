#pragma once
#include "0.math.hpp"
#include "0.common.hpp"

/**
 * Full Black Scholes Alternative Model
 * (Black 76)
 *
 * @param S : Stock price
 * @param K : Strike price
 * @param sig : volatility
 * @param r : risk-free interest rate
 * @param tau : time to maturity
 *
 * @returns BlackScholesResult
 */
inline BlackScholesResult black_scholes_alt(ftype S, ftype K, ftype sig, ftype r, ftype tau) {
    ftype D = e(-r * tau); // discount factor
    ftype F = S / D;       // forward price factor

    // calculate d1
    ftype A = ln(F/K) + (tau * p2(sig)) / 2.0f ;
    ftype B = sig * sqrt(tau);
    ftype d1 = A / B;

    // calculate d2
    ftype d2 = d1 - sig * sqrt(tau);

    // call
    ftype call = D * (N(d1) * F - N(d2) * K);

    // put
    ftype put = K * e( -r * tau) - S + call;

    return {.call=call, .put=put};
}
