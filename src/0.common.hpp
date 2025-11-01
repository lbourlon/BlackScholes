#pragma once

struct BsInput {
    double S;   // Stock price
    double K;   // Strike price
    double tau; // time to maturity
    double sig; // volatility
    double r;   // risk-free interest rate
};

struct BsOutput {
    double call;
    double put;
};
