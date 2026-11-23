// ========================================================================
// SAGEhydrology | Model source
// File: Native/models/CaRTPlus/cartplus_ode.cpp
//
// CaRT / CaRT+ numerical model and supporting routines.
//
// SAGE adaptations: Jasper A. Vrugt, University of California, Irvine.
// Upstream model: Xia (2026). Preserve the adjacent license notices.
// ========================================================================

/* Adaptive CaRT+ODE with forward sensitivities in physical parameters. */
#include "cartplus_catchment.hpp"
#include <limits>

namespace cartplus::ode {

/* Dormand-Prince coefficients for the accepted and embedded estimates. */
constexpr double a[7][7] = {
    {}, {1./5}, {3./40, 9./40}, {44./45, -56./15, 32./9},
    {19372./6561, -25360./2187, 64448./6561, -212./729},
    {9017./3168, -355./33, 46732./5247, 49./176, -5103./18656},
    {35./384, 0, 500./1113, 125./192, -2187./6784, 11./84}
};
constexpr double low[] = {
    5179./57600, 0, 7571./16695, 393./640, -92097./339200,
    187./2100, 1./40
};

/* Advance to the forcing boundary with state and tangent error control. */
/**
 * integrate: Advance the augmented state with the configured integrator.
 *
 * INPUT ARGUMENTS:
 * state : State vector (including tangent lanes when enabled).
 * rhs : Callback evaluating state derivatives on the model time base.
 * dt : Integration or observation interval duration (days).
 * options : Optional solver controls and regularization settings.
 *
 * OUTPUT ARGUMENTS:
 * No return value; updates referenced states and output buffers.
 */
template<std::size_t N, class S, class Rhs>
void integrate(std::array<S, N>& state, const Rhs& rhs, double dt,
               const double* options) {
    double time = 0;
    double step = std::min({dt, options[0], options[1]});
    auto first = rhs(state);
    for (int trial = 0; time < dt; ++trial) {
        if (trial >= 100000 || step < 1e-14)
            throw std::runtime_error("CaRT+ODE step limit");
        std::array<std::array<S, N>, 7> slopes;
        slopes[0] = first;
        std::array<S, N> next{};
        for (int stage = 1; stage < 7; ++stage) {
            next = state;
            for (int k = 0; k < stage; ++k) {
                if (a[stage][k] == 0) continue;
                for (std::size_t j = 0; j < N; ++j)
                    next[j] += (step*a[stage][k])*slopes[k][j];
            }
            slopes[stage] = rhs(next);
        }
        double error = 0;
        for (std::size_t j = 0; j < N; ++j) {
            double delta = 0;
            for (int k = 0; k < 7; ++k)
                delta += (a[6][k]-low[k])*primal(slopes[k][j]);
            const double scale = options[3] + options[2]*std::max(
                std::abs(primal(state[j])), std::abs(primal(next[j])));
            const double ratio = step*delta/scale;
            error += ratio*ratio;
        }
        error = std::sqrt(error/N);
        if constexpr (std::is_same_v<S, Tangent>) {
            // A steady primal store can still have stiff sensitivities.
            double tangent_error = 0;
            for (std::size_t j = 0; j < N; ++j) {
                for (std::size_t d = 0; d < parameter_count; ++d) {
                    double delta = 0;
                    for (int k = 0; k < 7; ++k)
                        delta += (a[6][k]-low[k])*slopes[k][j].slope[d];
                    const double magnitude = std::max(
                        std::abs(state[j].slope[d]),
                        std::abs(next[j].slope[d]));
                    const double scale = 1e-6 + 1e-4*magnitude;
                    const double ratio = step*delta/scale;
                    tangent_error += ratio*ratio;
                }
            }
            if (std::isfinite(tangent_error))
                error = std::max(error, std::sqrt(
                    tangent_error/(N*parameter_count)));
            else
                error = std::numeric_limits<double>::infinity();
        }
        if (error <= 1) {
            state = next;
            first = slopes[6];
            time += step;
            if (dt-time <= 1e-14*dt) return;
        }
        const double factor = std::isfinite(error) ? std::clamp(
            .9*std::pow(std::max(error, 1e-16), -.2), .2, 5.) : .2;
        step = std::min({step*factor, options[1], dt-time});
    }
}

/* Continuous Pareto partition, evaporation and percolation rates. */
/**
 * production: Evaluate production-store runoff, evaporation, and recharge.
 *
 * INPUT ARGUMENTS:
 * storage : Water store (mm); reference arguments are updated in place.
 * rain : Liquid input depth/rate on this operator time base.
 * pet : Potential evaporation depth/rate on this operator time base.
 * cap : Store capacity (mm).
 * shape : Dimensionless storage/infiltration shape parameter.
 * percolation : Percolation coefficient in the model parameterization.
 * net : Whether to remove interception before the production update.
 *
 * OUTPUT ARGUMENTS:
 * Numerical result or flux structure in the argument time base; non-const
 * state references are advanced in place.
 */
template<class S>
std::array<S, 3> production(S storage, S rain, S pet, S cap, S shape,
                            S percolation, bool net) {
    const S wet = bounded(storage/cap, 0., 1.);
    const S deficit = 1.-wet;
    const S intercept = net ? minimum(rain, pet) : S(0);
    const S exponent = shape/(1.+shape);
    // Keep the retention sensitivity finite as the store fills.
    constexpr double threshold = .01;
    const S retain = deficit < threshold
        ? (deficit/threshold)*exp(exponent*std::log(threshold))
        : exp(exponent*log(deficit));
    return {(rain-intercept)*(1.-retain),
            intercept+(pet-intercept)*(1.-deficit*deficit),
            maximum(storage, 0.)*log1p(percolation*pow(wet, 4.))/4.};
}

/* Continuous near-stream runoff, evaporation and recession. */
/**
 * zone: Evaluate near-stream runoff, evaporation, and drainage.
 *
 * INPUT ARGUMENTS:
 * store : Current water store (mm).
 * rain : Liquid input depth/rate on this operator time base.
 * pet : Potential evaporation depth/rate on this operator time base.
 * capacity : Store or infiltration capacity in this operator.
 * days : Residence time (days).
 *
 * OUTPUT ARGUMENTS:
 * Numerical result or flux structure in the argument time base; non-const
 * state references are advanced in place.
 */
template<class S>
std::array<S, 3> zone(S store, S rain, S pet, S capacity, S days) {
    const S s = bounded(store/capacity, 0., 1.);
    const S intercepted = minimum(rain, pet);
    return {(rain-intercepted)*s*s,
            intercepted+(pet-intercepted)*s*(2.-s),
            maximum(store, 0.)/days};
}

/* Run the coupled catchment, discrete routing and quick-flow ODE. */
/**
 * simulate: Simulate model fluxes and optional discharge sensitivities.
 *
 * INPUT ARGUMENTS:
 * n : Number of intervals or state entries in this helper.
 * forcing : N-by-4 forcing [P, T, PET, training]; depths mm/interval, T deg
 * C.
 * height : Height used by this numerical operator.
 * area : Five nonnegative band fractions, summing to one.
 * p : 42 physical CaRT+ parameters in native order.
 * dt : Integration or observation interval duration (days).
 * floor : Minimum infiltration capacity in daily-rate units (mm/day).
 * options : Optional solver controls and regularization settings.
 *
 * OUTPUT ARGUMENTS:
 * Numerical result or flux structure in the argument time base; non-const
 * state references are advanced in place.
 */
template<class S>
std::vector<std::array<S, 4>> simulate(
    std::size_t n, const double* forcing, const double* height,
    const double* area, const S* p, double dt, double floor,
    const double* options
) {
    std::size_t training = 0, training_end = 0;
    double mean_rain = 0, elevation = 0;
    for (int b = 0; b < 5; ++b) elevation += area[b]*height[b];
    const auto window = std::min(n, static_cast<std::size_t>(365/dt));
    S initial_veg = 0;
    std::vector<S> vegetation(n);
    for (std::size_t i = 0; i < n; ++i) {
        vegetation[i] = logistic((forcing[4*i+1]-p[22])/p[23]);
        if (i < window) initial_veg += vegetation[i]/window;
        if (forcing[4*i+3]) {
            ++training;
            training_end = i+1;
            mean_rain += forcing[4*i];
        }
    }
    if (!training) throw std::invalid_argument("Training mask is empty");
    mean_rain = std::max(mean_rain/(training*dt), .1);
    const S upland = 1.-p[13];
    std::array<S, 5> offset{};
    for (int b = 0; b < 5; ++b)
        offset[b] = p[5]*((height[b]-elevation)/1000.);

    const auto initial = [&](S reference) {
        return std::array<S, 21>{
            0, 0, 0, 0, 0, .5*p[6], .7*p[8], .3*p[14],
            .3*p[25], .5*p[16], mean_rain, initial_veg, reference
        };
    };
    const auto storage = [&](const auto& x) {
        S result = x[5]+upland*x[6]+p[13]*x[7]+x[9]+x[12];
        for (int b = 0; b < 5; ++b) result += area[b]*x[b];
        return result;
    };
    const auto advance = [&](auto& state, std::size_t i, S reference,
                             bool pilot) {
        const double precip = forcing[4*i]/dt;
        const double temp = forcing[4*i+1];
        const double pet = forcing[4*i+2]/dt;
        std::array<S, 5> melt_demand{}, snowfall{};
        S liquid_precip = 0, adjustment = 0;
        for (int b = 0; b < 5; ++b) {
            const S warmth = temp-offset[b]-p[0];
            if (warmth < 0.) {
                snowfall[b] = p[2]*precip;
                adjustment += area[b]*precip*(p[2]-1.);
            } else {
                melt_demand[b] = p[1]*warmth+p[4]*pet;
                liquid_precip += area[b]*precip;
            }
        }
        const S sublimation = p[3]*pet;
        const auto rhs = [&](const auto& x) {
            std::array<S, 21> z{};
            S liquid = liquid_precip, snow_et = 0;
            for (int b = 0; b < 5; ++b) {
                const S snow_rate = maximum(x[b], 0.)/.02;
                const S sub = minimum(sublimation, snow_rate);
                const S melt = minimum(
                    maximum(snow_rate-sub, 0.), melt_demand[b]);
                z[b] = snowfall[b]-sub-melt;
                liquid += area[b]*melt;
                snow_et += area[b]*sub;
            }
            const S demand = pet*(1.-p[21]*(1.-bounded(x[11], 0., 1.)));
            const S surface_et = minimum(demand, maximum(x[5], 0.)/.02);
            const S drain = maximum(x[5], 0.)/p[7];
            const S overflow = maximum(liquid-surface_et-drain, 0.)
                *bounded((x[5]/p[6]-.99)/.01, 0., 1.);
            const S input = drain+overflow;
            const S potential = maximum(demand-surface_et, 0.);
            const S cap = maximum(p[11]*exp(
                -p[12]*maximum(x[10], 0.)/mean_rain),
                std::max(floor, 1e-300));
            const S available = maximum(input-potential, 0.);
            const S excess = available > 40.*cap ? available-cap
                : available+cap*expm1(-available/cap);
            const S rain = input-excess;
            const auto soil = production(x[6], rain, potential, p[8],
                                         p[10], p[9], true);
            const auto near = zone(x[7], rain, potential, p[14], p[15]);
            const auto wet = zone(x[8], rain, potential, p[25], p[26]);
            const S et = upland*soil[1]+p[13]*near[1];
            const S fast = upland*soil[0]+p[13]*near[0]+excess;
            const auto dep = production(x[9], fast,
                maximum(potential-et, 0.), p[16], p[17], p[18], false);
            const S recharge = upland*soil[2]+dep[2];
            z[5] = liquid-surface_et-input;
            z[6] = rain-soil[0]-soil[1]-soil[2];
            z[7] = rain-near[0]-near[1]-near[2];
            z[8] = rain-wet[0]-wet[1]-wet[2];
            z[9] = fast-dep[0]-dep[1]-dep[2];
            z[10] = (input-x[10])/3.;
            z[11] = (vegetation[i]-x[11])/p[24];
            z[20] = recharge;
            if (pilot) return z;
            const S gw = maximum(x[12], 0.);
            const S gw_et = gw*p[20]*maximum(potential-et-dep[1], 0.)
                /(gw+reference);
            const S base = gw*gw/(reference*p[19]);
            const S signal = x[8]/p[25]-p[28]
                +p[30]*bounded(x[6]/p[8], 0., 1.)
                +p[31]*(gw/reference-1.);
            const S connected = p[27]+(1.-p[27])*logistic(signal/p[29]);
            const S bypass = p[32]*connected*dep[0];
            z[12] = recharge-base-gw_et;
            z[13] = dep[0]-bypass;
            z[14] = bypass;
            z[15] = p[13]*near[2];
            z[16] = base;
            z[17] = et+dep[1]+surface_et+snow_et+gw_et;
            z[18] = adjustment;
            z[19] = connected;
            return z;
        };
        std::fill(state.begin()+13, state.end(), S(0));
        integrate(state, rhs, dt, options);
    };

    auto state = initial(0);
    S recharge = 0;
    for (std::size_t i = 0; i < training_end; ++i) {
        advance(state, i, 0, true);
        if (forcing[4*i+3]) recharge += state[20];
    }
    const S reference = maximum(recharge/(training*dt), 1e-6)*p[19];
    state = initial(reference);
    UnitHydrograph<6, S> uh1(p[33], dt);
    UnitHydrograph<13, S> uh2(p[33], dt);
    S quick = .5*p[34];
    std::vector<std::array<S, 4>> output(n);
    for (std::size_t i = 0; i < n; ++i) {
        const S before = storage(state);
        advance(state, i, reference, false);
        const S routed = uh1.advance(state[13]);
        const S bypass = uh2.advance(state[14]);
        const S in_quick = routed/dt+.9*p[36];
        const S in_bypass = bypass/dt+.1*p[36];
        const double temp = forcing[4*i+1];
        const double kelvin = std::max(temp, 0.)+273.15;
        const double viscosity = -std::log(10.)*247.8
            *(1./(kelvin-140.)-1./153.15);
        const S signal = log(maximum(state[19]/dt, 1e-12))+viscosity
            -p[39]*logistic((p[37]-temp)/p[38]);
        const auto quick_rhs = [&](const std::array<S, 3>& x) {
            const S store = maximum(x[0], 0.);
            const S exchange = p[35]*pow(store/p[34], 3.5);
            const S diverted = maximum(in_bypass+.1*exchange, 0.);
            S incoming = in_quick+.9*exchange;
            if (store <= 0.) incoming = maximum(incoming, 0.);
            const S multiplier = signal+p[40]*log(
                maximum(store, 1e-3)/p[34]);
            const S release = .25*store*pow(store/p[34], 4.)
                *exp(6.*tanh(multiplier/6.));
            return std::array<S, 3>{incoming-release, release+diverted,
                incoming-in_quick+diverted-in_bypass};
        };
        std::array<S, 3> x{quick, 0, 0};
        integrate(x, quick_rhs, dt, options);
        const S q = x[1]+state[15]+state[16];
        const S cap = dt*p[41];
        const S loss = -cap*expm1(-q/cap);
        const S change = storage(state)-before + state[13]+state[14]
            -routed-bypass+x[0]-quick;
        output[i] = {q-loss, state[17],
                     state[18]+dt*p[36]+x[2]-loss, change};
        quick = x[0];
    }
    return output;
}

}  // namespace cartplus::ode

/* Return interval fluxes and an optional N-by-42 discharge Jacobian. */
/**
 * cartplus_ode_run: Native model entry point for interval fluxes and
 * requested derivatives.
 *
 * INPUT ARGUMENTS:
 * n : Number of intervals or state entries in this helper.
 * f : Forcing [P, T, PET, training] in the caller layout.
 * elevation : Five elevation-band heights (m).
 * area : Five nonnegative band fractions, summing to one.
 * p : 42 physical CaRT+ parameters in native order.
 * dt : Integration or observation interval duration (days).
 * floor : Minimum infiltration capacity in daily-rate units (mm/day).
 * options : Optional solver controls and regularization settings.
 * output : Caller-owned numerical output buffer.
 * jacobian : Whether to evaluate discharge sensitivities.
 *
 * OUTPUT ARGUMENTS:
 * Status code (0 success, nonzero failure); supplied output buffers are
 * filled in observation order.
 */
extern "C" int cartplus_ode_run(
    std::size_t n, const double* f, const double* elevation,
    const double* area, const double* p, double dt, double floor,
    const double* options, double* output, double* jacobian
) {
    try {
        if (jacobian) {
            std::array<cartplus::Tangent, 42> seeded{};
            for (std::size_t j = 0; j < seeded.size(); ++j) {
                seeded[j].value = p[j];
                seeded[j].slope[j] = 1.;
            }
            const auto rows = cartplus::ode::simulate(
                n, f, elevation, area, seeded.data(), dt, floor, options);
            for (std::size_t i = 0; i < n; ++i) {
                for (int j = 0; j < 4; ++j)
                    output[4*i+j] = rows[i][j].value;
                std::copy(rows[i][0].slope.begin(), rows[i][0].slope.end(),
                          jacobian+42*i);
            }
        } else {
            const auto rows = cartplus::ode::simulate(
                n, f, elevation, area, p, dt, floor, options);
            for (std::size_t i = 0; i < n; ++i)
                std::copy(rows[i].begin(), rows[i].end(), output+4*i);
        }
        return 0;
    } catch (...) {
        return 1;
    }
}

