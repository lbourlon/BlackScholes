#pragma once

#include <format>

struct BsInput {

    double S{};   // Stock price
    double K{};   // Strike price
    double tau{}; // time to maturity
    double sig{}; // volatility
    double r{};   // risk-free interest rate

    BsInput(double _S, double _K, double _tau, double _sig, double _r)
        : S(_S), K(_K), tau(_tau), sig(_sig), r(_r) {};
};

struct BsOutput {
    double call{};
    double put{};
    BsOutput(double _call, double _put) : call(_call), put(_put) {};
};

template <>
struct std::formatter<BsInput> {
    constexpr auto parse(std::format_parse_context &ctx){return ctx.begin();}

    auto format(const BsInput &d, std::format_context &ctx) const {
        return std::format_to(ctx.out(), "BsInput [S:{} | K:{} | tau:{} | sig:{} | r:{}]", d.S, d.K, d.tau, d.sig, d.r);
    }
};

template <>
struct std::formatter<BsOutput> {
    constexpr auto parse(std::format_parse_context &ctx){return ctx.begin();}

    auto format(const BsOutput &d, std::format_context &ctx) const {
        return std::format_to(ctx.out(), "BsOutput [Call:{} | Put:{}]", d.call, d.put);
    }
};
