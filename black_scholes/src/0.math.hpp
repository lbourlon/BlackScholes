#pragma once
#include <cmath>
#include <functional>

using std::exp;

inline auto ln(auto x) { return std::log(x); }
inline auto p2(auto x) { return std::pow(x, 2); }
inline auto sqrt(auto x) { return std::sqrt(x); }
inline auto abs(auto x) { return std::abs(x); }

using std::numbers::pi;

inline double integrate(std::function<double(double)> f, double min, double max) {
    double intervals = 200'000;
    double dx = (abs(max - min)) / intervals;
    double sum = 0;
    for (double x = min; x < max; x += dx) {
        sum += abs(dx * f(x));
    }
    return sum;
}

/**
 * Standard normal cumulative distribution function
 * (aka Loi normale)
 */
inline double N(double val) {
    double min = -10.0;
    double max = val;
    double dx = 0.001;
    auto f = [](const double &x) { return exp(-p2(x) / 2.); };

    double sum = 0.0f;
    // no need for abs(f(x)) because e(x) is always positive
    for (double x = min; x < max; x += dx) {
        sum += dx * f(x);
    }
    return sum / sqrt(2.0 * pi);
}

/**
 * Standard normal cumulative distribution function
 * Using trapezium integral calculation.
 * (aka Loi normale)
 *
 */
inline double N_trap(double val) {
    double min = -10.0;
    double max = val;
    double dx = 0.001;
    auto f = [](const double &x) { return exp(-p2(x) / 2.); };

    double sum = 0.0f;
    for (double x = min; x < max; x += dx) {
        sum += 0.5 * dx * (f(x) + f(x + dx));
    }
    return sum / sqrt(2.0 * pi);
}
