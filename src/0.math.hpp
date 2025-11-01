#pragma once
#include <cmath>
#include <functional>

using ftype = double;

inline auto ln(auto x) { return std::log(x);}
inline auto p2(auto x) { return std::pow(x,2);}
inline auto sqrt(auto x) { return std::sqrt(x);}
inline auto e(auto x) {return std::exp(x);}
inline auto abs(auto x) {return std::abs(x);}

inline auto pi = 3.141'592'653'589'793;

inline ftype integrate(std::function<ftype(ftype)> f, ftype min, ftype max) {
    ftype intervals = 200'000;
    ftype dx = (abs(max - min)) / intervals;
    ftype sum = 0;
    for (ftype x = min; x < max; x += dx) { sum += abs(dx * f(x)); }
    return sum;
}

/**
 * Standard normal cumulative distribution function
 * (aka Loi normale)
 */
inline ftype N(ftype val) {
    ftype min = -10.0;
    ftype max = val;
    ftype dx = 0.001;
    auto f = [](const ftype& x){return e(-p2(x)/2.);};

    ftype sum = 0.0f;
    // no need for abs(f(x)) because e(x) is always positive
    for (ftype x = min; x < max; x += dx) { sum += dx * f(x); }
    return sum / sqrt(2.0 * pi);
}

/**
 * Standard normal cumulative distribution function
 * Using trapezium integral calculation.
 * (aka Loi normale)
 *
 */
inline ftype N_trap(ftype val) {
    ftype min = -10.0;
    ftype max = val;
    ftype dx = 0.001;
    auto f = [](const ftype& x){return e(-p2(x)/2.);};

    ftype sum = 0.0f;
    for (ftype x = min; x < max; x += dx) {
        sum += 0.5 * dx * (f(x) + f(x+dx));
    }
    return sum / sqrt(2.0 * pi);
}
