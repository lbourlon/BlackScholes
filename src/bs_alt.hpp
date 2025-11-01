#pragma once
#include "_math.hpp"

inline ftype calc_d1_alt(ftype F, ftype K, ftype r, ftype sig, ftype tau) {
    ftype A = ln(F/K) + (tau * p2(sig)) / 2.0f ;
    ftype B = sig * sqrt(tau);
    return A / B;
}

inline ftype calc_d2_alt(ftype d1, ftype sig, ftype tau) {
    return d1 - sig * sqrt(tau);
}


inline ftype black_scholes_call_alt(ftype S, ftype K, ftype sig, ftype r, ftype tau) {
    ftype D = e(-r * tau); // discount factor
    ftype F = S / D;       // forward factor factor

    ftype d1 = calc_d1_alt(F, K, r, sig, tau);
    ftype d2 = calc_d2_alt(d1, sig, tau);
    return D * (N(d1) * F - N(d2) * K);
}


inline ftype black_scholes_put_alt(ftype S, ftype K, ftype r, ftype tau, ftype call_alt) {
    float D = e(-r * tau); // discount factor
    return call_alt - S + D * K;
}
