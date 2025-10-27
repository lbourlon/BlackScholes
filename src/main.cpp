#include <cmath>
#include <print>

inline float ln(float x) { return std::log10(x);}
inline float p2(float x) { return std::pow(x,2);}
inline float sqrt(float x) { return std::sqrt(x);}
inline float exp(float x) {return std::exp(x);}

constexpr float pi = 3.1451592;

float d1(float S, float K, float r, float sig, float tau) {
    float A = ln(S/K) + tau * (r + (1/2) * p2(sig));
    float B = sig * sqrt(tau);
    return A / B;
}

float d2(float d1, float sig, float tau) {
    return d1 - sig * std::sqrt(tau);
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
float black_sholes_call(float S, float K, float sig, float r, float tau) {
    return pi;
}




int main() {
    std::println("Hello World\n");
    return 1;
}
