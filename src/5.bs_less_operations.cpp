#include "0.common.hpp"
#include "0.fast_math.hpp"

std::pair<double,double> N_zelen_severo_branchless_simd(double x1, double x2) {
    const double x1_abs = std::fabs(x1);
    const double x2_abs = std::fabs(x2);

    double t1 = 0.23164'19;
    double t2 = 0.23164'19;

    t1 *= x1_abs;
    t2 *= x2_abs;

    t1 += 1.;
    t2 += 1.;

    t1 = 1. / t1;
    t2 = 1. / t2;

    double taylor_1 = t1 * 0.31938'1530;
    double taylor_2 = t2 * 0.31938'1530;

    const double t1_p2 = t1 * t1;
    const double t2_p2 = t2 * t2;
    taylor_1 -= t1_p2 * 0.35656'3782;
    taylor_2 -= t2_p2 * 0.35656'3782;

    const double t1_p3 = t1_p2 * t1;
    const double t2_p3 = t2_p2 * t2;
    taylor_1 += t1_p3 * 1.78147'7937;
    taylor_2 += t2_p3 * 1.78147'7937;

    const double t1_p4 = t1_p3 * t1;
    const double t2_p4 = t2_p3 * t2;
    taylor_1 -= t1_p4 * 1.82125'5978;
    taylor_2 -= t2_p4 * 1.82125'5978;

    const double t1_p5 = t1_p4 * t1;
    const double t2_p5 = t2_p4 * t2;
    taylor_1 += t1_p5 * 1.33027'4429;
    taylor_2 += t2_p5 * 1.33027'4429;


    double phi_1 = -fm::p2(x1_abs) * 0.5;
    double phi_2 = -fm::p2(x2_abs) * 0.5;

    phi_1 = fm::exp(phi_1);
    phi_2 = fm::exp(phi_2);

    phi_1 *= inv_sqrt_2pi;
    phi_2 *= inv_sqrt_2pi;

    double y1 = 1.;
    double y2 = 1.;

    y1 -= phi_1 * taylor_1;
    y2 -= phi_2 * taylor_2;


    const double mask1 = -static_cast<double>(std::signbit(x1));
    const double mask2 = -static_cast<double>(std::signbit(x2));

    const double res_1 = y1 + mask1 * (1.0 - 2.0 * y1);
    const double res_2 = y2 + mask2 * (1.0 - 2.0 * y2);
    return std::make_pair(res_1, res_2);
}

/**
 * Full Black Scholes model
 *
 * @param B : Stock price
 *
 * @returns BlackScholesResult
 */
BsOutput black_scholes_simd_oriented(BsInput in) {
    // calculate d1
    const double A = fm::log(in.S / in.K) + in.tau * (in.r + fm::pow_si(in.sig,2) / 2.0f);
    const double B = in.sig * fm::sqrt(in.tau);
    const double d1 = A / B;

    // calculate d2
    const double d2 = d1 - B;
    const auto N_res =  N_zelen_severo_branchless_simd(d1, d2);

    const double Ke = in.K * fm::exp(-in.r * in.tau);

    const double call = in.S * N_res.first - N_res.second * Ke;
    const double put = Ke - in.S + call;

    return {call, put};
}
