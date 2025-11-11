#pragma once
#include "0.common.hpp"

BsOutput black_scholes(BsInput in);
BsOutput black_scholes_alt(BsInput in); // side-quest

BsOutput black_scholes_numerical_approx_N(BsInput in);
BsOutput black_scholes_fast_math(BsInput in);
BsOutput black_scholes_pipelined(BsInput in);
