// Thin C ABI around SAGEhydrology's unchanged MATLAB-independent HYMOD core.
#include "hymod.hpp"
#include <climits>
#include <exception>

#if defined(_WIN32)
#define HYMOD_API __declspec(dllexport)
#else
#define HYMOD_API __attribute__((visibility("default")))
#endif

extern "C" HYMOD_API int site_hymod_run(
    int days, const double* precipitation, const double* pet,
    const double* temperature, const double* theta, const double* config,
    const double* settings, const double* initial, double* discharge,
    double* jacobian) noexcept {
    if (days < 1 || !precipitation || !pet || !temperature || !theta ||
        !config || !settings || !initial || !discharge || !jacobian ||
        settings[6] < 1 || settings[6] > INT_MAX) return 2;
    try {
        sage_hymod::Params p{theta[0], theta[1], theta[2], theta[3],
                            theta[4], theta[5], theta[6],
                            config[0], config[1], config[2]};
        sage_hymod::Options opt{settings[0], settings[1], settings[2],
                               settings[3], settings[4], settings[5],
                               static_cast<int>(settings[6])};
        sage_hymod::Forcing forcing;
        forcing.P = precipitation; forcing.Ep = pet; forcing.T = temperature;
        forcing.n = static_cast<std::size_t>(days);
        double z0[56] = {}, zfinal[56] = {};
        for (int i = 0; i < 7; ++i) z0[i] = initial[i];
        sage_hymod::OutputView output;
        output.Z = zfinal; output.zrows = 1; output.zcols = 56;
        output.q = discharge; output.nq = forcing.n;
        output.J = jacobian; output.nj = 7;
        return sage_hymod::run_into(days, z0, 56, forcing, p, opt,
                                   false, 1, true, output) ? 1 : 0;
    } catch (...) {
        return 3; // Never allow C++ exceptions to cross ctypes' C boundary.
    }
}
