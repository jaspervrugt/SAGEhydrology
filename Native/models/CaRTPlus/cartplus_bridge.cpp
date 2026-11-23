// ========================================================================
// SAGEhydrology | Model source
// File: Native/models/CaRTPlus/cartplus_bridge.cpp
//
// CaRT / CaRT+ numerical model and supporting routines.
//
// SAGE adaptations: Jasper A. Vrugt, University of California, Irvine.
// Upstream model: Xia (2026). Preserve the adjacent license notices.
// ========================================================================

#ifdef CARTPLUS_UNIFORM_INPUT
#include "cartplus_uniform_catchment.hpp"
#define CARTPLUS_SIMULATE simulate_uniform
#else
#include "cartplus_catchment.hpp"
#define CARTPLUS_SIMULATE simulate_us
#endif

#ifdef CARTPLUS_CALIBRATED_ONLY
#define cartplus_solver_jacobian cartplus_solver_jacobian35
#define cartplus_interval_jacobian cartplus_interval_jacobian35
#define cartplus_regional_jacobian cartplus_regional_jacobian35
#define cartplus_us_jacobian cartplus_us_jacobian35
#endif

namespace {
/**
 * read_parameters: Unpack the native CaRT+ physical parameter order.
 *
 * INPUT ARGUMENTS:
 * p : Model parameter vector; CaRT calibration or CaRT+ physical units.
 *
 * OUTPUT ARGUMENTS:
 * Numerical result or flux structure in the argument time base; non-const
 * state references are advanced in place.
 */
template<class Scalar>
cartplus::CatchmentParametersT<Scalar> read_parameters(const Scalar* p) {
        return {
            {p[0], p[1], p[2], p[3], p[4], p[5]},
            p[6], p[7], p[8], p[9], p[10], p[11], p[12], p[13],
            p[14], p[15], p[16], p[17], p[18], p[19], p[20], p[21],
            p[22], p[23], p[24], p[25], p[26], p[27], p[28], p[29],
            p[30], p[31], p[32], p[33], p[34], p[35], p[36], p[37],
            p[38], p[39], p[40], p[41]
        };
}
}

/**
 * cartplus_solver_forward: Native model entry point for interval fluxes and
 * requested derivatives.
 *
 * INPUT ARGUMENTS:
 * count : Number of intervals or state entries in this helper.
 * forcing : N-by-4 forcing [P, T, PET, training]; depths mm/interval, T deg
 * C.
 * elevation : Five elevation-band heights (m).
 * area : Five nonnegative band fractions, summing to one.
 * parameters : 42 physical CaRT+ parameters in native order.
 * output : Caller-owned numerical output buffer.
 * fixed_spinup : Optional detached baseline stores for derivative replay.
 * captured_spinup : Optional buffer receiving detached spin-up stores.
 * infiltration_floor : Minimum infiltration capacity (mm/day).
 * dt : Integration or observation interval duration (days).
 * integrate_zone : Near-stream forward/sensitivity integration callback.
 *
 * OUTPUT ARGUMENTS:
 * Status code (0 success, nonzero failure); supplied output buffers are
 * filled in observation order.
 */
extern "C" int cartplus_solver_forward(
    std::size_t count, const double* forcing, const double* elevation,
    const double* area, const double* parameters, double* output,
    const cartplus::SpinupStates* fixed_spinup,
    cartplus::SpinupStates* captured_spinup, double infiltration_floor,
    double dt, cartplus::ZoneIntegrator integrate_zone
) {
    try {
        const auto inputs = read_parameters(parameters);
        std::array<double, cartplus::band_count> heights{}, fractions{};
        std::copy_n(elevation, heights.size(), heights.begin());
        std::copy_n(area, fractions.size(), fractions.begin());
        std::vector<cartplus::ForcingDay> days;
        days.reserve(count);
        for (std::size_t i = 0; i < count; ++i) {
            days.push_back({forcing[4*i], forcing[4*i+1], forcing[4*i+2],
                            forcing[4*i+3] != 0.0});
        }
        const auto results = cartplus::CARTPLUS_SIMULATE(
            days, inputs, heights, fractions, fixed_spinup, captured_spinup,
                infiltration_floor,
            dt, integrate_zone
        );
        for (std::size_t i = 0; i < count; ++i) {
            output[4*i] = results[i].discharge;
            output[4*i+1] = results[i].evaporation;
            output[4*i+2] = results[i].external_flux;
            output[4*i+3] = results[i].storage_change;
        }
        return 0;
    } catch (...) {
        return 1;
    }
}

/**
 * cartplus_solver_jacobian: Native model entry point for interval fluxes
 * and requested derivatives.
 *
 * INPUT ARGUMENTS:
 * count : Number of intervals or state entries in this helper.
 * forcing : N-by-4 forcing [P, T, PET, training]; depths mm/interval, T deg
 * C.
 * elevation : Five elevation-band heights (m).
 * area : Five nonnegative band fractions, summing to one.
 * parameters : 42 physical CaRT+ parameters in native order.
 * output : Caller-owned numerical output buffer.
 * fixed_spinup : Optional detached baseline stores for derivative replay.
 * infiltration_floor : Minimum infiltration capacity (mm/day).
 * dt : Integration or observation interval duration (days).
 * integrate_zone : Near-stream forward/sensitivity integration callback.
 *
 * OUTPUT ARGUMENTS:
 * Status code (0 success, nonzero failure); supplied output buffers are
 * filled in observation order.
 */
extern "C" int cartplus_solver_jacobian(
    std::size_t count, const double* forcing, const double* elevation,
    const double* area, const double* parameters, double* output,
    const cartplus::SpinupStates* fixed_spinup, double infiltration_floor,
    double dt, cartplus::ZoneIntegrator integrate_zone
) {
    if (!fixed_spinup) return 1;
    try {
        std::array<double, cartplus::band_count> heights{}, fractions{};
        std::copy_n(elevation, heights.size(), heights.begin());
        std::copy_n(area, fractions.size(), fractions.begin());
        std::vector<cartplus::ForcingDay> days;
        days.reserve(count);
        for (std::size_t i = 0; i < count; ++i) {
            days.push_back({forcing[4*i], forcing[4*i+1], forcing[4*i+2],
                            forcing[4*i+3] != 0.0});
        }
        std::array<cartplus::Tangent,
                   cartplus::physical_parameter_count> seeded{};
        std::size_t column = 0;
        for (std::size_t j = 0; j < seeded.size(); ++j) {
            seeded[j].value = parameters[j];
#ifdef CARTPLUS_CALIBRATED_ONLY
            // Shared vegetation, connectivity and frost constants.
            if ((j >= 22 && j <= 24) || j == 28 || j == 29
                || j == 37 || j == 38) continue;
#endif
            seeded[j].slope[column++] = 1.0;
        }
        // The forward pass already computed the water-balance diagnostics.
        // The tangent pass only needs to return the discharge derivatives.
        cartplus::CARTPLUS_SIMULATE<cartplus::Tangent, true>(
            days, read_parameters(seeded.data()), heights, fractions,
            fixed_spinup, nullptr, infiltration_floor, dt, integrate_zone,
            output);
        return 0;
    } catch (...) {
        return 1;
    }
}

/**
 * cartplus_interval_forward: Native model entry point for interval fluxes
 * and requested derivatives.
 *
 * INPUT ARGUMENTS:
 * count : Number of intervals or state entries in this helper.
 * forcing : N-by-4 forcing [P, T, PET, training]; depths mm/interval, T deg
 * C.
 * elevation : Five elevation-band heights (m).
 * area : Five nonnegative band fractions, summing to one.
 * parameters : 42 physical CaRT+ parameters in native order.
 * output : Caller-owned numerical output buffer.
 * fixed : Detached spin-up stores, or None/empty to determine them.
 * captured : Optional collector receiving baseline spin-up stores.
 * floor : Minimum infiltration capacity in daily-rate units (mm/day).
 * dt : Integration or observation interval duration (days).
 *
 * OUTPUT ARGUMENTS:
 * Status code (0 success, nonzero failure); supplied output buffers are
 * filled in observation order.
 */
extern "C" int cartplus_interval_forward(
    std::size_t count, const double* forcing, const double* elevation,
    const double* area, const double* parameters, double* output,
    const cartplus::SpinupStates* fixed, cartplus::SpinupStates* captured,
    double floor, double dt
) {
    return cartplus_solver_forward(count, forcing, elevation, area,
        parameters, output, fixed, captured, floor, dt, nullptr);
}

/**
 * cartplus_interval_jacobian: Native model entry point for interval fluxes
 * and requested derivatives.
 *
 * INPUT ARGUMENTS:
 * count : Number of intervals or state entries in this helper.
 * forcing : N-by-4 forcing [P, T, PET, training]; depths mm/interval, T deg
 * C.
 * elevation : Five elevation-band heights (m).
 * area : Five nonnegative band fractions, summing to one.
 * parameters : 42 physical CaRT+ parameters in native order.
 * output : Caller-owned numerical output buffer.
 * fixed : Detached spin-up stores, or None/empty to determine them.
 * floor : Minimum infiltration capacity in daily-rate units (mm/day).
 * dt : Integration or observation interval duration (days).
 *
 * OUTPUT ARGUMENTS:
 * Status code (0 success, nonzero failure); supplied output buffers are
 * filled in observation order.
 */
extern "C" int cartplus_interval_jacobian(
    std::size_t count, const double* forcing, const double* elevation,
    const double* area, const double* parameters, double* output,
    const cartplus::SpinupStates* fixed, double floor, double dt
) {
    return cartplus_solver_jacobian(count, forcing, elevation, area,
        parameters, output, fixed, floor, dt, nullptr);
}

/**
 * cartplus_regional_forward: Native model entry point for interval fluxes
 * and requested derivatives.
 *
 * INPUT ARGUMENTS:
 * count : Number of intervals or state entries in this helper.
 * forcing : N-by-4 forcing [P, T, PET, training]; depths mm/interval, T deg
 * C.
 * elevation : Five elevation-band heights (m).
 * area : Five nonnegative band fractions, summing to one.
 * parameters : 42 physical CaRT+ parameters in native order.
 * output : Caller-owned numerical output buffer.
 * fixed : Detached spin-up stores, or None/empty to determine them.
 * captured : Optional collector receiving baseline spin-up stores.
 * infiltration_floor : Minimum infiltration capacity (mm/day).
 *
 * OUTPUT ARGUMENTS:
 * Status code (0 success, nonzero failure); supplied output buffers are
 * filled in observation order.
 */
extern "C" int cartplus_regional_forward(
    std::size_t count, const double* forcing, const double* elevation,
    const double* area, const double* parameters, double* output,
    const cartplus::SpinupStates* fixed, cartplus::SpinupStates* captured,
    double infiltration_floor
) {
    return cartplus_interval_forward(count, forcing, elevation, area,
        parameters, output, fixed, captured, infiltration_floor, 1.0);
}

/**
 * cartplus_regional_jacobian: Native model entry point for interval fluxes
 * and requested derivatives.
 *
 * INPUT ARGUMENTS:
 * count : Number of intervals or state entries in this helper.
 * forcing : N-by-4 forcing [P, T, PET, training]; depths mm/interval, T deg
 * C.
 * elevation : Five elevation-band heights (m).
 * area : Five nonnegative band fractions, summing to one.
 * parameters : 42 physical CaRT+ parameters in native order.
 * output : Caller-owned numerical output buffer.
 * fixed : Detached spin-up stores, or None/empty to determine them.
 * infiltration_floor : Minimum infiltration capacity (mm/day).
 *
 * OUTPUT ARGUMENTS:
 * Status code (0 success, nonzero failure); supplied output buffers are
 * filled in observation order.
 */
extern "C" int cartplus_regional_jacobian(
    std::size_t count, const double* forcing, const double* elevation,
    const double* area, const double* parameters, double* output,
    const cartplus::SpinupStates* fixed, double infiltration_floor
) {
    return cartplus_interval_jacobian(count, forcing, elevation, area,
        parameters, output, fixed, infiltration_floor, 1.0);
}

/**
 * cartplus_us_forward: Native model entry point for interval fluxes and
 * requested derivatives.
 *
 * INPUT ARGUMENTS:
 * count : Number of intervals or state entries in this helper.
 * forcing : N-by-4 forcing [P, T, PET, training]; depths mm/interval, T deg
 * C.
 * elevation : Five elevation-band heights (m).
 * area : Five nonnegative band fractions, summing to one.
 * parameters : 42 physical CaRT+ parameters in native order.
 * output : Caller-owned numerical output buffer.
 * fixed : Detached spin-up stores, or None/empty to determine them.
 * captured : Optional collector receiving baseline spin-up stores.
 *
 * OUTPUT ARGUMENTS:
 * Status code (0 success, nonzero failure); supplied output buffers are
 * filled in observation order.
 */
extern "C" int cartplus_us_forward(
    std::size_t count, const double* forcing, const double* elevation,
    const double* area, const double* parameters, double* output,
    const cartplus::SpinupStates* fixed, cartplus::SpinupStates* captured
) {
    return cartplus_regional_forward(count, forcing, elevation, area,
                                    parameters, output, fixed, captured, 0.0
                                );
}

/**
 * cartplus_us_jacobian: Native model entry point for interval fluxes and
 * requested derivatives.
 *
 * INPUT ARGUMENTS:
 * count : Number of intervals or state entries in this helper.
 * forcing : N-by-4 forcing [P, T, PET, training]; depths mm/interval, T deg
 * C.
 * elevation : Five elevation-band heights (m).
 * area : Five nonnegative band fractions, summing to one.
 * parameters : 42 physical CaRT+ parameters in native order.
 * output : Caller-owned numerical output buffer.
 * fixed : Detached spin-up stores, or None/empty to determine them.
 *
 * OUTPUT ARGUMENTS:
 * Status code (0 success, nonzero failure); supplied output buffers are
 * filled in observation order.
 */
extern "C" int cartplus_us_jacobian(
    std::size_t count, const double* forcing, const double* elevation,
    const double* area, const double* parameters, double* output,
    const cartplus::SpinupStates* fixed
) {
    return cartplus_regional_jacobian(count, forcing, elevation, area,
                                     parameters, output, fixed, 0.0);
}

