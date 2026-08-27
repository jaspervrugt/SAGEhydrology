/*
 * cfe_nwm.hpp
 *
 * Public native interface for the CFE-NWM rainfall-runoff model. The
 * interface contains no MATLAB types so it can be used by both the standalone
 * crr_cfe_nwm gateway and the unified crr_model_mex dispatcher.
 *
 * Written by Jasper A. Vrugt, January 2025.
 * Native-core interface introduced for SAGEhydrology, 2026.
 */

#pragma once
#include <cstddef>

namespace sage_cfe_nwm {
struct Params {
    double s_max;
    double s_fc;
    double s_wp;
    double k_sch;
    double a1;
    double k_perc;
    double lf_thr;
    double a2;
    double k_lf;
    double g_max;
    double c_gw;
    double mm;
    double k_nsh;
    double T_tr;
    double f_dd;
    double dT;
    double T_sm;
    double eps_m;
    double eps;
    double rho;
    int K;
    int L;
    const double* giuh = nullptr;
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
} // namespace sage_cfe_nwm
