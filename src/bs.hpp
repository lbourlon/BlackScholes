#pragma once
#include "_math.hpp"

inline ftype calc_d1(ftype S, ftype K, ftype r, ftype sig, ftype tau) {
    ftype A = ln(S/K) + tau * (r + p2(sig)/2.0f);
    ftype B = sig * sqrt(tau);
    return A / B;
}

inline ftype calc_d2(ftype d1, ftype sig, ftype tau) {
    return d1 - sig * sqrt(tau);
}

inline ftype black_scholes_put(ftype S, ftype K, ftype r, ftype tau, ftype call) {
    return K * e( -r * tau) - S + call;
}

/**
 * The full Black Scholes model
 *
 * @param S : Stock price
 * @param K : Strike price
 * @param sig : volatility
 * @param r : risk-free interest rate
 * @param tau : time to maturity
 *
 * @returns the call price for the stock
 */
inline ftype black_scholes_call(ftype S, ftype K, ftype sig, ftype r, ftype tau) {
    ftype d1 = calc_d1(S, K, r, sig, tau);
    ftype d2 = calc_d2(d1, sig, tau);
    return S * N(d1) - N(d2) * K * e(-r * tau);
}
