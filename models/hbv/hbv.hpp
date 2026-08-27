/*
 * hbv.hpp
 *
 * Public native interface for the HBV rainfall-runoff model. The interface
 * contains no MATLAB types so it can be used by both the standalone crr_hbv
 * gateway and the unified crr_model_mex dispatcher.
 *
 * Written by Jasper A. Vrugt, February 2022.
 * Native-core interface introduced for SAGEhydrology, 2026.
 */

#pragma once

#include <cstddef>
#include <vector>

namespace sage_hbv {

struct Params {
    double f_c;
    double beta;
    double lp;
    double k_0;
    double uzl;
    double k_1;
    double k_2;
    double perc;
    double T_tr;
    double f_dd;
    double sfcf;
    double cfr;
    double eps_s = 1e-3;
    double eps_t = 1e-1;
    double eps_x = 1e-12;
};

struct Options {
    double InitStep;
    double MaxStep;
    double MinStep;
    double RelTol;
    double AbsTol;
    double Order;
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
    std::size_t zrows = 0;
    std::size_t zcols = 0;
    std::size_t nq = 0;
    std::size_t nj = 0;
};

struct Result {
    std::vector<double> Z;
    std::vector<double> q;
    std::vector<double> J;
    std::size_t zrows = 0;
    std::size_t zcols = 0;
    bool fail = false;
};

// High-performance native HBV core. It writes directly into caller-owned
// output buffers, avoiding extra result allocations/copies in MEX gateways.
// Z must be zrows x 65 in MATLAB column-major layout.
bool run_into(int ns,
              const double* z0,
              std::size_t nz0,
              const Forcing& forcing,
              const Params& pars,
              const Options& options,
              bool mem,
              int ipr,
              bool needJ,
              const OutputView& out);

// Owning-result convenience wrapper for callers that do not provide buffers.
Result run(int ns,
           const std::vector<double>& z0,
           const Forcing& forcing,
           const Params& pars,
           const Options& options,
           bool mem,
           int ipr,
           bool needJ);

} // namespace sage_hbv
