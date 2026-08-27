#pragma once

#include <cstddef>

namespace sage_gr4jb {

struct Params {
    double x1;
    double x2;
    double x3;
    double x4;
    double x5;
    double f_p;
    double T_tr;
    double f_dd;
    double eta;
    double tau;
    double kappa;
    double bg;
    double bR;
    double mts;
    double T_sm;
    double eps_m;
    double eps_s;
    double rho;
    int L;
    const double* U = nullptr;
    const double* dUdx4 = nullptr;
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

bool run_into(
    int ns,
    const double* z0,
    std::size_t nz0,
    const Forcing& forcing,
    const Params& parameters,
    const Options& options,
    bool mem,
    int ipr,
    bool needJ,
    const OutputView& output);

} // namespace sage_gr4jb
