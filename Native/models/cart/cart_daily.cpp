// ========================================================================
// SAGEhydrology | Model source
// File: Native/models/cart/cart_daily.cpp
//
// CaRT / CaRT+ numerical model and supporting routines.
//
// SAGE adaptations: Jasper A. Vrugt, University of California, Irvine.
// Upstream model: Xia (2026). Preserve the adjacent license notices.
// ========================================================================

// CaRT daily baseline after Xia (2026), adapted under the MIT license.
// Uses the same forward-mode scalar arithmetic as CaRT+.
#define CART_BASELINE_TANGENTS
#include "../CaRTPlus/cartplus_tangent.hpp"
#include <algorithm>
#include <array>
#include <vector>
#include "../diagnostic_hooks.hpp"

namespace cart {
using namespace cartplus;
constexpr int count = parameter_count;
static_assert(parameter_count == count, "CaRT tangent width must match");

// Uniform-input operators: Xia hourly
// f986af265633a7a980836c2e0f2124077d632420 (MIT).
// Legacy cart_daily remains separate; cart_uniform also applies at daily
// resolution.
/**
 * integrated_uh: Evaluate the integral of a unit-hydrograph S-curve.
 *
 * INPUT ARGUMENTS:
 * time : Current model time (days).
 * x4 : First unit-hydrograph time base (days).
 * symmetric : Whether to use the symmetric two-time-base hydrograph.
 *
 * OUTPUT ARGUMENTS:
 * Numerical result or flux structure in the argument time base; non-const
 * state references are advanced in place.
 */
template<class T>
T integrated_uh(double time, T x4, bool symmetric)
{
    if (time <= 0) return T(0);
    const T u = T(time)/x4;
    if (!symmetric) {
        if (u > T(1)) return T(time) - x4*T(1-1./3.5);
        return x4*pow(u,3.5)/T(3.5);
    }
    if (u <= T(1)) return x4*pow(u,3.5)/T(7);
    if (u < T(2)) return T(time)-x4+x4*pow(T(2)-u,3.5)/T(7);
    return T(time)-x4;
}

/**
 * power_release: Compute stable fourth-power store release.
 *
 * INPUT ARGUMENTS:
 * state : State vector (including tangent lanes when enabled).
 * coefficient : Release-law coefficient in the current interval time base.
 *
 * OUTPUT ARGUMENTS:
 * Numerical result or flux structure in the argument time base; non-const
 * state references are advanced in place.
 */
template<class T>
T power_release(T state, T coefficient)
{
    return -state*expm1(T(-.25)*log1p(coefficient*pow(state,4)));
}

/**
 * quick_chain: Apply split quick-store exchange and release substeps.
 *
 * INPUT ARGUMENTS:
 * store : Current water store (mm).
 * x9 : Input to the 90-percent quick-routing branch.
 * x1 : Input to the 10-percent direct-routing branch.
 * x2 : Exchange coefficient in the CaRT parameterization.
 * x3 : Quick-store capacity scale (mm).
 * half_release : Nonlinear release coefficient for each half step.
 * steps : Number of split substeps or integration stages.
 * fraction : Dimensionless partition/drainage fraction.
 *
 * OUTPUT ARGUMENTS:
 * Numerical result or flux structure in the argument time base; non-const
 * state references are advanced in place.
 */
template<class T>
std::array<T,3> quick_chain(T store, T x9, T x1, T x2, T x3,
                            T half_release, int steps, T fraction = T(.9))
{
    T outflow = 0, direct_total = 0;
    for (int step = 0; step < steps; ++step) {
        const T before = power_release(store,half_release);
        store -= before;
        const T exchange = x2*pow(store/x3,3.5);
        const T midpoint = std::max(store+T(.5)*(x9+fraction*exchange),T(0))
            ;
        const T mid_exchange = x2*pow(midpoint/x3,3.5);
        store = std::max(store+x9+fraction*mid_exchange,T(0));
        const T after = power_release(store,half_release);
        store -= after;
        const T direct = std::max(x1+(T(1)-fraction)*mid_exchange,T(0));
        direct_total += direct;
        outflow += before+after+direct;
    }
    return {store,outflow,direct_total};
}

/**
 * route: Convolve interval runoff with unit-hydrograph weights.
 *
 * INPUT ARGUMENTS:
 * runoff : Liquid runoff entering the routing operator.
 * x4 : First unit-hydrograph time base (days).
 * symmetric : Whether to use the symmetric two-time-base hydrograph.
 * resolution : Forcing intervals per day: 1 daily, 24 hourly, 96
 * quarter-hourly.
 * uniform : Whether to use the uniform-input operator.
 *
 * OUTPUT ARGUMENTS:
 * Numerical result or flux structure in the argument time base; non-const
 * state references are advanced in place.
 */
template<class T>
std::vector<T> route(const std::vector<T>& runoff, T x4, bool symmetric,
                     int resolution, bool uniform)
{
    const int length = uniform
        ? static_cast<int>(std::ceil(primal(x4)*(symmetric ? 2 : 1)*
            resolution))+1
        : (symmetric ? 13 : 6)*resolution;
    std::vector<T> weights(length);
    T previous = 0;
    for (int i = 0; i < length; ++i) {
        if (uniform) {
            weights[i] = T(resolution)*(integrated_uh(double(i+1)/resolution
                ,x4,symmetric)
                         - T(2)*integrated_uh(double(i)/resolution,x4,
                             symmetric)
                         + integrated_uh(double(i-1)/resolution,x4,symmetric
                             ));
            continue;
        }
        const T u = T(double(i + 1)/resolution) / x4;
        const T s = symmetric
            ? (u <= T(1) ? T(.5) * pow(u, 2.5)
                          : T(1) - T(.5) * pow(std::max(T(2)-u,T(0)),2.5))
            : pow(std::min(u,T(1)),2.5);
        weights[i] = s - previous;
        previous = s;
    }
    std::vector<T> result(runoff.size());
    for (std::size_t i = 0; i < runoff.size(); ++i)
        for (int j = 0; j < length && std::size_t(j) <= i; ++j)
            result[i] += weights[j] * runoff[i-j];
    return result;
}

/**
 * simulate: Simulate model fluxes and optional discharge sensitivities.
 *
 * INPUT ARGUMENTS:
 * f : Forcing [P, T, PET, training] in the caller layout.
 * n : Number of intervals or state entries in this helper.
 * z : Five elevation-band heights (m).
 * a : Five nonnegative band fractions, summing to one.
 * p : 11 CaRT calibration coordinates, including logarithmic entries.
 * initial : Initial stores (mm), or None/empty to use the spin-up policy.
 * fixed : Detached spin-up stores, or None/empty to determine them.
 * resolution : Forcing intervals per day: 1 daily, 24 hourly, 96
 * quarter-hourly.
 * uniform : Whether to use the uniform-input operator.
 * initialize_only : Whether to return the pilot initialization alone.
 *
 * OUTPUT ARGUMENTS:
 * Numerical result or flux structure in the argument time base; non-const
 * state references are advanced in place.
 */
template<class T>
std::vector<T> simulate(const double* f, int n, const double* z,
                        const double* a, const std::array<T,count>& p,
                        std::array<double,3>& initial, bool fixed, int
                            resolution,
                        bool uniform, bool initialize_only = false)
{
    // Forcing is row-major [precipitation, temperature, PET, training].
    const T dt = 1.0/resolution;
#ifdef CART_SPLIT_MODEL
    const T fraction = p[11];
#else
    const T fraction = .9;
#endif
    const T c = exp(p[0]), x3 = exp(p[2]), tau = exp(p[3]);
    const T kp = T(std::pow(4./9,4)) * exp(p[5]);
    const T cf = exp(p[7]), gamma = exp(p[9]), ks = exp(p[10]);
    std::array<T,5> snow{};
    double mean_z = 0;
    for (int k = 0; k < 5; ++k) mean_z += a[k]*z[k];
    std::vector<T> liquid(n), direct(n), recharge(n), fast(n), q(n);
    for (int i = 0; i < n; ++i) {
        for (int k = 0; k < 5; ++k) {
            const T warm = T(f[4*i+1]) - gamma*T((z[k]-mean_z)/1000)
                         - p[6];
            const bool cold = warm < T(0);
            if (cold) snow[k] += p[8]*T(f[4*i]);
            const T sublimation = ks*T(f[4*i+2]);
            const T potential_melt = std::max(dt*cf*warm,T(0));
            T melt;
            if (uniform) {
                const T demand = sublimation+potential_melt;
                const T fraction = demand > snow[k] ? snow[k]/demand : T(1);
                melt = potential_melt*fraction;
                snow[k] = std::max(snow[k]-sublimation*fraction-melt,T(0));
            } else {
                snow[k] -= std::min(sublimation,snow[k]);
                melt = std::min(potential_melt,snow[k]);
                snow[k] -= melt;
            }
            liquid[i] += T(a[k])*(melt + T(cold ? 0 : f[4*i]));
        }
    }
    // Snow starts empty; the initialized subsurface stores remain
    // explicitly initial water.
    const T half_percolation = dt*kp*T(.5)/pow(c,4);
    const auto production = [&](T store, int steps, bool spinup) {
        for (int i = 0; i < steps; ++i) {
            if (uniform && spinup) store -= power_release(store,
                half_percolation);
            const T pet = f[4*i+2];
            const T pn = std::max(liquid[i]-pet,T(0));
            const T en = std::max(pet-liquid[i],T(0));
            const T s = store/c, tp = tanh(pn/c), te = tanh(en/c);
            const T add = c*(T(1)-s*s)*tp/(T(1)+s*tp);
            const T evap = store*(T(2)-s)*te/(T(1)+(T(1)-s)*te);
            const T v = store - evap + add;
            recharge[i] = uniform
                ? power_release(v,(spinup ? half_percolation : T(2)*
                    half_percolation))
                : v*(T(1)-pow(T(1)+dt*kp*pow(v/c,4),-.25));
            direct[i] = pn-add;
            store = v-recharge[i];
            if (!spinup) {
                CART_TRACE(i, SOIL_END, store); CART_TRACE(i, SOIL_RUNOFF,
                    direct[i]);
                CART_TRACE(i, SOIL_RECHARGE, recharge[i]);
                CART_TRACE(i, SOIL_ET, evap+std::min(liquid[i],pet));
                CART_TRACE(i, SOIL_INPUT, liquid[i]);
                CART_TRACE(i, SOIL_INTERCEPTION, std::min(liquid[i],pet));
            }
        }
        return store;
    };
    if (!fixed) initial[0] = primal(production(T(.7)*c,std::min(n,1460*
        resolution),true));
    CART_INITIAL(1,initial[0]);
    production(T(initial[0]),n,false);
    const int quick_spinup_steps = std::min(n,730*resolution);
    // Causal routing needs only this prefix to initialize the quick store.
    if (initialize_only) direct.resize(quick_spinup_steps);
    const auto main = route(direct,p[1],false,resolution,uniform);
    const auto bypass = route(direct,p[1],true,resolution,uniform);
    const T half_release = T(.5)*dt/pow(x3,4);
    const auto quick = [&](T store, int steps, bool capture=false) {
        for (int i = 0; i < steps; ++i) {
            if (uniform) {
                const T x9 = fraction*main[i], x1 = (T(1)-fraction)*bypass[i
                    ], x2 = dt*p[4];
                const auto fine = quick_chain(store,x9/T(4),x1/T(4),x2/T(4),
                                              x3,half_release/T(4),4,
                                fraction);
                const auto coarse = quick_chain(store,x9/T(2),x1/T(2),x2/T(2
                    ),
                                                x3,half_release/T(2),2,
                                fraction);
                store = std::max((T(4)*fine[0]-coarse[0])/T(3),T(0));
                fast[i] = (T(4)*fine[1]-coarse[1])/T(3);
                if (capture) {
                    CART_TRACE(i, DIRECT_Q, (T(4)*fine[2]-coarse[2])/T(3));
                    CART_TRACE(i, QUICK_Q, fast[i]-(T(4)*fine[2]-coarse[2])/
                        T(3));
                    CART_TRACE(i, QUICK_END, store);
                }
                continue;
            }
            const T loss = dt*p[4]*pow(store/x3,3.5);
            const T v = std::max(store+fraction*(main[i]+loss),T(0));
            store = v*pow(T(1)+dt*pow(v/x3,4),-.25);
            const T qb=std::max((T(1)-fraction)*(bypass[i]+loss),T(0));
            fast[i] = v-store+qb;
            if (capture) { CART_TRACE(i, DIRECT_Q, qb); CART_TRACE(i,
                QUICK_Q, v-store);
                CART_TRACE(i, QUICK_END, store); }
        }
        return store;
    };
    if (!fixed) initial[1] = primal(quick(T(.5)*x3,quick_spinup_steps));
    CART_INITIAL(6,initial[1]);
    if (!initialize_only) quick(T(initial[1]),n,true);
    for (int i=0;i<n && !initialize_only;++i) {
        CART_TRACE(i, ROUTE_IN, fraction*direct[i]);
        CART_TRACE(i, BYPASS_IN, (T(1)-fraction)*direct[i]);
        CART_TRACE(i, ROUTE_OUT, fraction*main[i]);
        CART_TRACE(i, BYPASS_OUT, (T(1)-fraction)*bypass[i]);
    }
    T gbar = 0;
    int days = 0;
    for (int i = 0; i < n; ++i) {
        if (f[4*i+3] == 0) continue;
        gbar += recharge[i];
        ++days;
    }
    gbar = std::max(gbar/(T(days)*dt),T(1e-6));
    const T coefficient = T(1)/(gbar*tau*tau);
    if (!fixed) initial[2] = primal(gbar*tau);
    if (initialize_only) return {};
    T water = initial[2];
    CART_INITIAL(5,water);
    for (int i = 0; i < n; ++i) {
        if (uniform) {
            const T g = std::max(recharge[i],T(0)), a_step = dt*coefficient;
            const T y = g*a_step;
            T ratio, complement;
            if (y < T(.04)) {
                complement = y*(T(1./3)-y*(T(2./15)-y*(T(17./315)-y*T(62./
                    2835))));
                ratio = T(1)-complement;
            } else {
                const T root = sqrt(y);
                ratio = tanh(root)/root;
                complement = T(1)-ratio;
            }
            const T denominator = T(1)+a_step*water*ratio;
            q[i] = fast[i]+(a_step*water*water*ratio+y*water*ratio+g*
                complement)/denominator;
            water = (water+g*ratio)/denominator;
            CART_TRACE(i, GW_END, water); CART_TRACE(i, GW_Q, q[i]-fast[i]);
            continue;
        }
        const T v = water+recharge[i];
        water = v/(T(1)+dt*coefficient*v);
        q[i] = fast[i]+v-water;
        CART_TRACE(i, GW_END, water); CART_TRACE(i, GW_Q, v-water);
    }
    return q;
}
} // namespace cart

#ifdef _WIN32
#define CART_API __declspec(dllexport)
#else
#define CART_API __attribute__((visibility("default")))
#endif

/**
 * run_cart: Integrate the selected CaRT equations and output buffers.
 *
 * INPUT ARGUMENTS:
 * forcing : N-by-4 forcing [P, T, PET, training]; depths mm/interval, T deg
 * C.
 * n : Number of intervals or state entries in this helper.
 * elevation : Five elevation-band heights (m).
 * area : Five nonnegative band fractions, summing to one.
 * parameters : 11 CaRT calibration coordinates, including logarithmic
 * entries.
 * discharge : Caller-owned output buffer for interval discharge (mm).
 * jacobian : Whether to evaluate discharge sensitivities.
 * stores : Store count, state vector, or caller-owned state buffer.
 * fixed : Detached spin-up stores, or None/empty to determine them.
 * resolution : Forcing intervals per day: 1 daily, 24 hourly, 96
 * quarter-hourly.
 * uniform : Whether to use the uniform-input operator.
 *
 * OUTPUT ARGUMENTS:
 * Numerical result or flux structure in the argument time base; non-const
 * state references are advanced in place.
 */
static int run_cart(
    const double* forcing, int n, const double* elevation,
    const double* area, const double* parameters, double* discharge,
    double* jacobian, double* stores, int fixed, int resolution, bool
        uniform)
{
    // Stores are detached for sensitivity calculations, as in Xia's code.
    if (n < 1 || (resolution != 1 && resolution != 24 && resolution != 96))
        return 1;
    bool training = false;
    for (int i = 0; i < n; ++i) training |= forcing[4*i+3] != 0;
    if (!training) return 1;
    std::array<double,cart::count> p{};
    std::copy_n(parameters,cart::count,p.begin());
    std::array<double,3> initial{};
    if (fixed) std::copy_n(stores,3,initial.begin());
    // The hourly tangent replay produces the same discharge with these
    // stores.
    const bool reuse_discharge = jacobian && uniform && resolution == 24;
    if (!reuse_discharge || !fixed) {
        const auto q = cart::simulate(
            forcing,n,elevation,area,p,initial,fixed,resolution,uniform,
                reuse_discharge);
        if (!reuse_discharge) std::copy(q.begin(),q.end(),discharge);
    }
    std::copy(initial.begin(),initial.end(),stores);
    if (jacobian) {
        std::array<cartplus::Tangent,cart::count> seeded{};
        for (int j = 0; j < cart::count; ++j) {
            seeded[j] = p[j];
            seeded[j].slope[j] = 1;
        }
        const auto result = cart::simulate(
            forcing,n,elevation,area,seeded,initial,true,resolution,uniform)
                ;
        for (int i = 0; i < n; ++i) {
            if (reuse_discharge) discharge[i] = result[i].value;
            std::copy_n(result[i].slope.begin(),cart::count,
                        jacobian+i*cart::count);
        }
    }
    return 0;
}

/**
 * cart_daily: Native model entry point for interval fluxes and requested
 * derivatives.
 *
 * INPUT ARGUMENTS:
 * forcing : N-by-4 forcing [P, T, PET, training]; depths mm/interval, T deg
 * C.
 * n : Number of intervals or state entries in this helper.
 * elevation : Five elevation-band heights (m).
 * area : Five nonnegative band fractions, summing to one.
 * parameters : 11 CaRT calibration coordinates, including logarithmic
 * entries.
 * discharge : Caller-owned output buffer for interval discharge (mm).
 * jacobian : Whether to evaluate discharge sensitivities.
 * stores : Store count, state vector, or caller-owned state buffer.
 * fixed : Detached spin-up stores, or None/empty to determine them.
 * resolution : Forcing intervals per day: 1 daily, 24 hourly, 96
 * quarter-hourly.
 *
 * OUTPUT ARGUMENTS:
 * Status code (0 success, nonzero failure); supplied output buffers are
 * filled in observation order.
 */
extern "C" CART_API int cart_daily(
    const double* forcing, int n, const double* elevation,
    const double* area, const double* parameters, double* discharge,
    double* jacobian, double* stores, int fixed, int resolution)
{
    return run_cart(forcing,n,elevation,area,parameters,discharge,jacobian,
                    stores,fixed,resolution,resolution == 24);
}

/**
 * cart_uniform: Native model entry point for interval fluxes and requested
 * derivatives.
 *
 * INPUT ARGUMENTS:
 * forcing : N-by-4 forcing [P, T, PET, training]; depths mm/interval, T deg
 * C.
 * n : Number of intervals or state entries in this helper.
 * elevation : Five elevation-band heights (m).
 * area : Five nonnegative band fractions, summing to one.
 * parameters : 11 CaRT calibration coordinates, including logarithmic
 * entries.
 * discharge : Caller-owned output buffer for interval discharge (mm).
 * jacobian : Whether to evaluate discharge sensitivities.
 * stores : Store count, state vector, or caller-owned state buffer.
 * fixed : Detached spin-up stores, or None/empty to determine them.
 * resolution : Forcing intervals per day: 1 daily, 24 hourly, 96
 * quarter-hourly.
 *
 * OUTPUT ARGUMENTS:
 * Status code (0 success, nonzero failure); supplied output buffers are
 * filled in observation order.
 */
extern "C" CART_API int cart_uniform(
    const double* forcing, int n, const double* elevation,
    const double* area, const double* parameters, double* discharge,
    double* jacobian, double* stores, int fixed, int resolution)
{
    return run_cart(forcing,n,elevation,area,parameters,discharge,jacobian,
                    stores,fixed,resolution,true);
}

