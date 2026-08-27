/*
 * sacsma.hpp
 *
 * Public native interface for the SAC-SMA rainfall-runoff model. The
 * interface contains no MATLAB types so it can be used by both the standalone
 * crr_sacsma gateway and the unified crr_model_mex dispatcher.
 *
 * Written by Jasper A. Vrugt.
 * Native-core interface introduced for SAGEhydrology, 2026.
 */

#pragma once
#include <cstddef>

namespace sage_sacsma {
struct Params {
    double uzfwm;
    double uztwm;
    double lzfpm;
    double lzfsm;
    double lztwm;
    double zperc;
    double rexp;
    double uzk;
    double pfree;
    double lzpk;
    double lzsk;
    double acm;
    double kf;
    double T_tr;
    double f_dd;
    double T_sm;
    double eps_m;
    double eps_s;
    double eps;
    double rho;
};

struct Options {
    double InitStep, MaxStep, MinStep, RelTol, AbsTol, Order;
    int maxiter;
};

struct Forcing {
    const double* P = nullptr;
    const double* Ep = nullptr;
    const double* T = nullptr;
    std::size_t n = 0;
};

struct OutputView {
    double* Z = nullptr;
    double* q = nullptr;
    double* J = nullptr;
    std::size_t zrows = 0, zcols = 0, nq = 0, nj = 0;
};

bool run_into(int ns,
              const double* z0,
              std::size_t nz0,
              const Forcing&,
              const Params&,
              const Options&,
              bool mem,
              int ipr,
              bool needJ,
              const OutputView&);
} // namespace sage_sacsma
