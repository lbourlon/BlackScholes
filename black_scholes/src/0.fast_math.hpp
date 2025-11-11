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
constexpr double inv_sqrt_2pi = 1. / sqrt_2pi;
using std::numbers::ln10;
using std::numbers::ln2;

constexpr uint64_t double_sign_bit = 0x8000000000000000;
constexpr uint64_t double_exponent_bits = 0x7ff0000000000000;
constexpr uint64_t double_mantissa_bits = 0x000fffffffffffff;

constexpr double p2(double x) { return x * x; }
constexpr double p3(double x) { return x * x * x; }

constexpr double pad_exp(double x); //imprecise
constexpr double exp_taylor_expansion(double x);
constexpr double exp_continued_fraction(double x);
constexpr double exp_polynomial(double x);
constexpr double exp_test(double x);
constexpr double exp(double x) { return std::exp(x); }
// constexpr double exp(double x) { return std::exp(x); }


constexpr double ln_newton_halley(double x);
constexpr double ln_newton_halley_fast(double x); // not any faster than st:ln
constexpr double ln_power_of_2(double x);         // Fast generalist implem
constexpr double ln_mercator(double x);           // Domain specific optimization
constexpr double ln_sqrt2(double x);              // Just a test, but works lol
constexpr double ln_sqrt2_twice(double x);        // This actually dumb
constexpr double ln(double x) { return ln_mercator(x); }
constexpr double ln_std(double x) { return std::log(x); }
    }
    return yn1;
}
// inline double ln(double x) { return std::log(x); }
inline double ln(double x) { return ln_fracs(x); }

constexpr double sqrt_exp_ln(double x);
constexpr double sqrt_pow_ln(double x);
constexpr double sqrt_newton_raphson(double x);
constexpr double sqrt_kahan_ng(double x);
constexpr double sqrt_std(double x) { return std::sqrt(x); }
constexpr double sqrt(double x) { return sqrt_kahan_ng(x); }

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
constexpr double exp_taylor_expansion(double x) {
    constexpr int elements = 8080;
    double sum = 1., factorial = 1., power = 1;
    for (int i = 1; i < elements; ++i) {
        sum += (power *= x) / (factorial *= i);
    }
    return sum;
}

/* ----------- EXP functions ----------- */

// https://justinwillmert.com/articles/2020/numerically-computing-the-exponential-function-with-polynomial-approximations/#range-reduction
constexpr double exp_reduced(double x) {
    double k = floor(x * log2e + 0.5);
    double r = x - k * ln2;

    const double exp_reduced_taylor =
        1.0 + r + 0.5 * pow_si(r, 2) + 0.16666666666666666 * pow_si(r, 3) +
        0.041666666666666664 * pow_si(r, 4) + 0.008333333333333333 * pow_si(r, 5) +
        0.001388888888888889 * pow_si(r, 6);

    return std::pow(2, k) * exp_reduced_taylor;
}

constexpr double exp_reduced_fast(double x) {
    const double k = static_cast<int>(x * log2e + std::copysign(0.5, x));
    const double r = x - k * ln2;

    const double exp_reduced_taylor =
        1.0 + r + 0.5 * pow_si(r, 2) + 0.16666666666666666 * pow_si(r, 3) +
        0.041666666666666664 * pow_si(r, 4) + 0.008333333333333333 * pow_si(r, 5) +
        0.001388888888888889 * pow_si(r, 6);

    return fm::ldexp(exp_reduced_taylor, k);
}

constexpr double exp_reduced_const_array(double x) {
    const double k = static_cast<int>(x * log2e + std::copysign(0.5, x));
    const double r = x - k * ln2;

    static constexpr std::array<double, 6> coeffs = {
        0.5, 0.16666666666666666, 0.041666666666666664, 0.008333333333333333,
        0.001388888888888889};
    double exp_reduced_taylor = 1.0 + r;
    for (size_t i = 0; i < coeffs.size(); ++i) {
        exp_reduced_taylor += coeffs[i] * pow_si(r, i + 2);
    }

    return fm::ldexp(exp_reduced_taylor, k);
}

/* ----------- Misc ----------- */
constexpr double pow_i(double val, int exponent) {
    if (exponent == 0) {
        return 1.0;
    }
    double sum = 1;
    if (exponent > 0) {
        for (int e = 0; e < exponent; ++e)
            sum *= val;
        return sum;
    } else {
        for (int e = 0; e < -exponent; ++e)
            sum *= val;
        return 1 / sum;
    }
}

constexpr double pow_si(double val, uint64_t expo) {
    double sum = 1;
    for (int e = 0; e < expo; ++e)
        sum *= val;
    return sum;
}
// Compute X * 2^n, without std::ldexp
constexpr double fm_ldexp(double x, int64_t n) {
    // if positive this is straightforward
    if (n > 0) {
        return x * static_cast<double>(1ull << n);
    }

    // if negative we use the property a^{-n} = 1/(a^{n})
    return x / static_cast<double>(1ull << (-n));
};
/* ----------- LN functions ----------- */

// Converges very fast, but is too slow and depends on exp
// https://en.wikipedia.org/wiki/Natural_logarithm#High_precision
constexpr double ln_newton_halley(double x) {
    constexpr int elements = 2; // converges insanely fast
    double yn = x;
    double yn1 = 0;
    for (int i = 0; i < elements; ++i) {
        yn1 = yn + 2 * (x - fm::exp(yn)) / (x + fm::exp(yn));
        yn = yn1;
    }
    return yn1;
}

constexpr double ln_newton_halley_fast(double x) {
    double yn1 = x + 2 * (x - fm::exp(x)) / (x + fm::exp(x));
    return yn1 + 2 * (x - fm::exp(yn1)) / (x + fm::exp(yn1));
}

// This implementation plays on te fact that S/K shouldn't be very big
// And as such reducing the rage of ln may not even be needed
// Precision is yet to be measured, however, its clearly faster than its bigger
// brother right after
// https://en.wikipedia.org/wiki/Mercator_series
constexpr double ln_mercator(double a) {
    a -= 1;
    return a - 0.5 * pow_si(a, 2) + 0.3333333333333333 * pow_si(a, 3) -
           0.25 * pow_si(a, 4) + 0.2 * pow_si(a, 5);
}

// This should be a lot more precise for large numbers from the previous
// approximation however, it doesn't seem needed (CS major eyes lol™)
// https://en.wikipedia.org/wiki/Natural_logarithm#Natural_logarithm_of_10
// https://math.stackexchange.com/questions/3381629/what-is-the-fastest-algorithm-for-finding-the-natural-logarithm-of-a-big-number
// ln(a * 2n) = ln(a) + n * ln(2)
constexpr double ln_power_of_2(double x) {
    // const auto guess = std::bit_cast<double, uint64_t>(x);

    static constexpr int64_t _10bits = 0x3ffuLL; // (2^p)-1

    // This operation shifts away the mantissa, and to get the base 10
    // exponent by subbing the exponent bias from the result (1023 for doubles)
    // https://en.wikipedia.org/wiki/Floating-point_arithmetic#Internal_representation
    // Also not out the sign bit
    const int64_t exponent =
        ((std::bit_cast<uint64_t, double>(x) & double_exponent_bits) >> 52) - 1023ull;

    // Then to factor the mantissa out we zero-out non-mantissa bits
    // And we set the mantissa to be 1023, such that we have 2^0
    // https://en.wikipedia.org/wiki/Double-precision_floating-point_format#Exponent_encoding
    const uint64_t mantissa =
        ((std::bit_cast<uint64_t, double>(x) & double_mantissa_bits) | (_10bits << 52));

    const double a = std::bit_cast<double, uint64_t>(mantissa) - 1.0;

    // Mercator series for ln(a), with only one factor no more factors not really required
    // for the precision I want
    const double ln_a = a - 0.5 * pow_si(a, 2);

    // Solve the new equation
    return ln_a + exponent * ln2;
}
// Note: factoring for power of 10 `ln(a*10n)=ln(a)+n*ln(10)` actually makes no
// sense in my mind, given I already expect |x| < 10

// Because x is so small and we have a pretty fast sqrt approximation
// using the property ln(x) = 2 * ln(sqrt(x)
// We can compute ln_a with a shorter mercator series
// Stems from the same reason that ln_mercator() works as is, S/K should be small
constexpr double ln_sqrt2(double x) {
    const double a = fm::sqrt(x) - 1.0;

    // Mercator series for ln(a), more factors not really required for the precision I
    // want
    const double ln_a = a - 0.5 * pow_si(a, 2) + 0.3333333333333333 * pow_si(a, 3);

    // Solve the new equation
    return 2 * ln_a;
}

// This is so dumb (kinda works™)
// Same thing as the last one, but no mercator
// Stems from the same reason that ln_mercator() works as is, S/K should be small
// https://math.stackexchange.com/questions/3381629/what-is-the-fastest-algorithm-for-finding-the-natural-logarithm-of-a-big-number
constexpr double ln_sqrt2_twice(double x) { return 4 * (fm::sqrt(fm::sqrt(x)) - 1.0); }

/* ----------- Sqrt functions ----------- */

constexpr double sqrt_exp_ln(double x) { return fm::exp(fm::ln_std(x) / 2.); }
constexpr double sqrt_pow_ln(double x) { return std::pow(10, ln(x)) / 2.; }


// Newton(/Raphson)  Square Method
// S.G. Johnson (2015) Square Roots via Newton’s Method
// https://math.mit.edu/~stevenj/18.335/newton-sqrt.pdf
constexpr double sqrt_newton_raphson(double x) {
    double guess = x;
    static constexpr size_t guess_count = 5;

    for (size_t i = 0; i < guess_count; ++i) {
        double g = (guess + x / guess) * 0.5;
        guess = g;
    }
    return guess;
}

// From an unpublished "paper by Prof W. Kahan and K.C. Ng, written in May, 1986"
// I could not find an official source for the paper beyond
// https://www.netlib.org/fdlibm/e_sqrt.c
// The extract also mentions the use of the following correction table
// Though from my needs this seems unneeded, this approach seems was used
// This allows me to only do one newton-raphson step
constexpr double sqrt_kahan_ng(double x) {
    const auto guess_int = (std::bit_cast<uint64_t, double>(x) >> 1) + (1023ull << 51);
    const auto guess = std::bit_cast<double, uint64_t>(guess_int);
    return (guess + x / guess) * 0.5;
}

} // namespace fm
