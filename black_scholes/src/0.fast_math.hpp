#pragma once

#include <cmath>
#include <cstdint>
#include <numbers>
#include <array>

using std::numbers::e;
using std::numbers::ln2;
using std::numbers::log2e;
using std::numbers::pi;
using std::numbers::sqrt2;
constexpr double sqrt_2pi = 2.5066282746310002;

namespace fm {

constexpr double p2(double x) { return x * x; }
constexpr double p3(double x) { return x * x * x; }

constexpr double pad_exp(double x); //imprecise
constexpr double exp_taylor_expansion(double x);
constexpr double exp_continued_fraction(double x);
constexpr double exp_polynomial(double x);
constexpr double exp_test(double x);
constexpr double exp(double x) { return std::exp(x); }
// constexpr double exp(double x) { return std::exp(x); }


// Newton-Halley method
inline double ln_fracs(double x) {
    constexpr int elements = 2; // converges insanely fast
    double yn = x;
    double yn1 = 0;
    for (int i = 0; i < elements; i += 1) {
        yn1 = yn + 2 *(x-fm::exp(yn))/(x+fm::exp(yn));
        yn = yn1;
    }
    return yn1;
}
// inline double ln(double x) { return std::log(x); }
inline double ln(double x) { return ln_fracs(x); }

constexpr double sqrt_expln(double x) { return fm::exp(fm::ln(x)/2.);  }
constexpr double sqrt_pln(double x) { return std::pow(10, ln(x))/2.; }
inline double sqrt(double x) {return sqrt_expln(x);}

constexpr double N_marsaglia(double x); // slowest
constexpr double N_taylor_expansion(double x); // close second
constexpr double N_zelen_severo(double x); // fastest
constexpr double N(double x) { return N_zelen_severo(x); }


constexpr double pow_i(double val, uint16_t exponent);

/* ----------- Normal Distribution funcs ----------- */

// Standard normal probability density function
// Small phi not to be confused with Big Phi (which I'm notating as N);
constexpr double phi(const double x) { return fm::exp(-0.5 * p2(x)) / sqrt_2pi; }

// Source : George Marsaglia (2004). "Evaluating the Normal Distribution"
// Modified from original implementation to
// Slower than Zelen & Severo
constexpr double N_marsaglia(double x) {
    constexpr double error = 0.000025;
    double s = x, t = 0, b = x, q = x * x, i = 1;
    while (s - t > error) {
        t = s;
        i += 2;
        b *= q / i;
        s = t + b;
    }
    return .5 + s * fm::exp(-.5 * x * x - .91893853320467274178L);
}

// Also based  on the George Marsaglia paper.
// Seems to be faster for the precision level I'm going for
// no dependency on exp function
constexpr double N_taylor_expansion(double x) {
    if (x < 0.0)
        return 1.0 - N_taylor_expansion(-x);
    constexpr int elements = 22;
    double sum = x, factorial = 1., power = x;
    for (int i = 3; i < elements; i += 2) {
        sum += (power *= (x * x)) / (factorial *= i);
    }
    return 0.5 + phi(x) * sum;
}

// Zelen & Severo (1964); error |ε(x)| < 7.5·10−8
// see https://personal.math.ubc.ca/~cbm/aands/page_932.htm
constexpr double N_zelen_severo(double x) {
    using fm::pow_i;
    if (x < 0.0) return 1.0 - N(-x);
    const double t = 1. / (1. + 0.23164'19 * x);
    return 1. - fm::phi(x)*(t * 0.31938'1530 + pow_i(t, 2) * -0.35656'3782 +
                       pow_i(t, 3) * 1.78147'7937 + pow_i(t, 4) * -1.82125'5978 +
                       pow_i(t, 5) * 1.33027'4429);
}


/* ----------- EXP functions ----------- */

// https://mathworld.wolfram.com/PadeApproximant.html
constexpr double pade_exp(double x) {
    using fm::p2;
    using fm::p3;
    return  (120 + 60*x + 12* p2(x) + p3(x))
            /(120 - 60*x + 12* p2(x) - p3(x)) ;// e_3_3
}

// Taylor expansion exponential
// Not considerably slower than std::exp()
// but is actually constexpr ( c++26 not out :/ )
constexpr double exp_taylor_expansion(double x) {
    constexpr int elements = 16;
    double sum = 1., factorial = 1., power = 1;
    for (int i = 1; i < elements; ++i ) { sum += (power *= x) / (factorial *= i); }
    return sum;
}


constexpr double _cnt_frac_recurse(double x, int i) {
    constexpr int step = 4;
    constexpr int elements = step* 5;
    i += step;
    if (i < elements){
        return (i) + p2(x) / _cnt_frac_recurse(x, i);
    } else {
        return x;
    }
}
//
// Continued fraction based on Euler Identity
constexpr double exp_continued_fraction(double x) {
    return 1 + 2*x / ((2-x) + p2(x) /_cnt_frac_recurse(x, 2));
}



// Yaya, D Dia (2023) Approximate incomplete integrals, application to
// complementary error function
// relative error of less than 2^-53
// expects x > 0
constexpr double erfc(double x) {
    constexpr const double b0_0 = 2.92678600515804815402;
    constexpr std::array<double, 5> b1 = {
        8.97280659046817350354, 10.27157061171363078863, 12.72323261907760928036,
        16.88639562007936907786, 24.12333774572479110372};
    constexpr std::array<double, 5> b2 = {
        5.81582518933527390512, 5.70347935898051436684, 5.51862483025707963145,
        5.26184239579604207321, 4.92081346632882032881};
    constexpr std::array<double, 5> c1 = {
        11.61511226260603247078, 18.25323235347346524796, 18.38871225773938486923,
        18.61193318971775795045, 24.14804072812762821134};
    constexpr std::array<double, 5> c2 = {
        3.83362947800146179416, 7.30756258553673541139, 8.42742300458043240405,
        5.66479518878470764762, 4.91396098895240075156};

    double sum = 1;
    for (int i = 0; i < c2.size(); ++i) {
        sum *= ((fm::p2(x) + c2[i] * x + c1[i]) / (fm::p2(x) + b2[i] * x + b1[i]));
    }
    return ( sum / (x + b0_0));
}

/* ----------- Misc ----------- */

// Shore's inverse phi
constexpr double inverse_phi_shore(double x) {
    if (x < 0.5) return -inverse_phi_shore(1 - x);
    return 1 / 5.5556 * (1 - std::pow(((1 - x) / x), 0.1186));
}


constexpr double pow_i(double val, uint16_t exponent) {
    double sum = val;
    for (uint16_t e = 1; e < exponent; ++e)
        sum *= val;
    return sum;
}

} // namespace fm
