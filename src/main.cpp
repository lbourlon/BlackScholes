#include <cmath>
#include <print>
#include "_math.hpp"
#include "bs.hpp"
#include "bs_alt.hpp"



int main() {
    ftype S = 300.0;
    ftype K = 250.0;
    ftype tau = 1;
    ftype sig = 0.15;
    ftype r = 0.03;

    ftype call = black_scholes_call(S, K, sig, r, tau);
    ftype put = black_scholes_put(S, K, r, tau, call);
    std::println("Call: {}/58.81977, Put: {}/1.43116", call, put);

    call = black_scholes_call_alt(S, K, sig, r, tau);
    put = black_scholes_put_alt(S, K, r, tau, call);
    std::println("Call: {}/58.81977, Put: {}/1.43116 [alt]", call, put);

    S = 821.0;
    K = 781.0;
    tau = 3;
    sig = 0.33;
    r = 0.049;
    call = black_scholes_call(S, K, sig, r, tau);
    put = black_scholes_put(S, K, r, tau, call);
    std::println("Call: {}/251.26628, Put: {}/104.40837", call, put);

    call = black_scholes_call_alt(S, K, sig, r, tau);
    put = black_scholes_put_alt(S, K, r, tau, call);
    std::println("Call: {}/58.81977, Put: {}/1.43116 [alt]", call, put);
}
