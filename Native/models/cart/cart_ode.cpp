// ========================================================================
// SAGEhydrology | Model source
// File: Native/models/cart/cart_ode.cpp
//
// CaRT / CaRT+ numerical model and supporting routines.
//
// SAGE adaptations: Jasper A. Vrugt, University of California, Irvine.
// Upstream model: Xia (2026). Preserve the adjacent license notices.
// ========================================================================

// Native coupled CaRT with event-aware forward sensitivities.
#define CART_ODE_TANGENTS
#include "../CaRTPlus/cartplus_runge_kutta.hpp"


namespace cartplus::continuous {

/**
 * run_cart: Integrate the selected CaRT equations and output buffers.
 *
 * INPUT ARGUMENTS:
 * n : Number of intervals or state entries in this helper.
 * f : Forcing [P, T, PET, training] in the caller layout.
 * z : Five elevation-band heights (m).
 * a : Five nonnegative band fractions, summing to one.
 * p : 11 CaRT calibration coordinates, including logarithmic entries.
 * dt : Integration or observation interval duration (days).
 * initial : Initial stores (mm), or None/empty to use the spin-up policy.
 * m : Positive routing refinement: stages per X4 time base.
 * rtol : Positive relative integration-error tolerance (-).
 * atol : Positive absolute integration-error tolerance for augmented
 * states.
 * solver : Integration method or solver-control object.
 * q : Discharge or flux values (mm per observation interval).
 * jac : Discharge Jacobian in the selected parameter coordinates.
 * balance : Balance used by this numerical operator.
 *
 * OUTPUT ARGUMENTS:
 * No return value; updates referenced states and output buffers.
 */
template<class S>
void run_cart(std::size_t n, const double* f, const double* z, const double*
    a,
              const S* p, const double* dt, const double* initial, int m,
              double rtol, double atol, int solver, double* q, double* jac,
              double* balance) {
    const int ns = 8+3*m;
    int main_size=ns+5;
#ifdef CART_SPLIT_MODEL
    const S fraction = p[11];
#else
    const S fraction = .9;
#endif
    const S c = exp(p[0]), x3 = exp(p[2]), tau = exp(p[3]);
    const S kp = std::pow(4./9,4)*exp(p[5]), cf = exp(p[7]);
    const S gamma = exp(p[9]), ks = exp(p[10]), lambda = m/p[1];
#ifdef CART_CACHE_CONSTANTS
    const S inverse_c = 1./c;
#endif
    const auto eta9 = exits(m,false), eta1 = exits(m,true);
    double elevation = 0, training_days = 0;
    for (int b = 0; b < 5; ++b) elevation += a[b]*z[b];
    const auto initialize = [&](bool pilot, S gbar) {
        std::vector<S> y(pilot ? 7 : ns+5);
        if (initial) std::copy(initial,initial+(pilot ? 6 : ns),y.begin());
        else {
            y[5] = .7*c;
            if (!pilot) {
                y[6] = .5*x3;
                y[7] = gbar*tau;
            }
        }
        return y;
    };
    Rk2Workspace<S> pilot_steps(solver == 1 ? 7 : 0);
    Rk2Workspace<S> main_steps(solver == 1 ? main_size : 0);
    const auto advance = [&](std::vector<S>& y, std::size_t i, bool pilot,
                             const S& aquifer_scale, double& next_step) {
        const double rain = f[4*i]/dt[i], pet = f[4*i+2]/dt[i];
        std::array<S,5> snowfall{}, melt{}, snow_rate{};
        std::array<bool,5> active{};
        const S sub_demand = ks*pet;
        double liquid_rain = 0, cold_area = 0;
        for (int b = 0; b < 5; ++b) {
            const S warm = f[4*i+1]-gamma*((z[b]-elevation)/1000.)-p[6];
            if (warm < 0.) {
                snowfall[b] = p[8]*rain;
                cold_area += a[b];
            } else {
                melt[b] = cf*warm;
                liquid_rain += a[b]*rain;
            }
            active[b] = y[b] > 0. || snowfall[b]-sub_demand-melt[b] > 0.;
        }
        const S input = rain+(p[8]-1.)*rain*cold_area;
        S pn, en, atmospheric;
        const auto snow_rates = [&] {
            S liquid = liquid_rain, sub_total = 0;
            for (int b = 0; b < 5; ++b) {
                const S sub = active[b] ? sub_demand : minimum(sub_demand,
                    snowfall[b]);
                const S thaw = active[b] ? melt[b]
                    : minimum(melt[b],maximum(snowfall[b]-sub,0.));
                snow_rate[b] = snowfall[b]-sub-thaw;
                liquid += a[b]*thaw;
                sub_total += a[b]*sub;
            }
            const S e0 = minimum(liquid,pet);
            pn = liquid-e0;
            en = pet-e0;
            atmospheric = e0+sub_total;
        };
        snow_rates();
        const auto rhs = [&](const std::vector<S>& state, std::vector<S>&
            out) {
            std::copy(snow_rate.begin(),snow_rate.end(),out.begin());
#ifdef CART_CACHE_CONSTANTS
            const S ratio = state[5]*inverse_c, ratio2 = ratio*ratio;
#else
            const S ratio = state[5]/c, ratio2 = ratio*ratio;
#endif
            const S infiltration = (1.-ratio2)*pn;
            const S evap = ratio*(2.-ratio)*en;
            const S recharge = .25*kp*state[5]*ratio2*ratio2;
            out[5] = infiltration-evap-recharge;
            if (pilot) {
                out[6] = recharge;
                return;
            }
            const auto cascade = [&](int offset, const auto& eta, S incoming
                ) {
                S collected = 0.;
                for (std::size_t j = 0; j < eta.size(); ++j) {
                    const S release = lambda*state[offset+j];
                    out[offset+j] = incoming-release;
                    collected += eta[j]*release;
                    incoming = (1.-eta[j])*release;
                }
                return collected;
            };
            const S direct = pn-infiltration;
            const S u9 = cascade(8,eta9,fraction*direct);
            const S u1 = cascade(8+m,eta1,(1.-fraction)*direct);
            const S quick = maximum(state[6],0.), groundwater = maximum(
                state[7],0.);
            const S loss = -p[4]*pow(quick/x3,3.5);
            const S qr = .25*quick*pow(quick/x3,4.);
            const S qb = maximum(u1-(1.-fraction)*loss,0.);
            const S qg = groundwater*groundwater/aquifer_scale;
            out[6] = u9-qr-fraction*loss;
            out[7] = recharge-qg;
            out[ns] = qr+qb+qg;
            out[ns+1] = atmospheric+evap;
            out[ns+2] = fraction*loss+minimum(u1,(1.-fraction)*loss);
            out[ns+3] = input;
            out[ns+4] = recharge;
        };
        double time = 0;
        S elapsed = 0.;
        const bool fixed_method = solver == 10 || solver == 11;
        while (time < dt[i]) {
            std::array<double,5> candidates{};
            double segment = dt[i]-time;
            S duration = S(dt[i])-elapsed;
            for (int b = 0; b < 5; ++b) {
                candidates[b] = active[b] && snow_rate[b] < 0.
                    ? primal(y[b])/(-primal(snow_rate[b]))
                    : std::numeric_limits<double>::infinity();
                if (candidates[b] < segment) {
                    segment = candidates[b];
                    if (fixed_method) duration = -y[b]/snow_rate[b];
                }
            }
            if (segment <= 0.) throw std::runtime_error(
                "Invalid CaRT snow event");
            if (fixed_method) {
                integrate_fixed_duration(y,rhs,duration,solver == 10);
                elapsed = elapsed+duration;
            } else if (solver == 1)
                integrate_rk2(y,rhs,segment,rtol,atol,next_step,pilot ?
                    pilot_steps : main_steps);
            else
                integrate_segment(y,rhs,segment,rtol,atol,next_step,solver);
            time += segment;
            for (int b = 0; b < 5; ++b) {
                if (std::abs(candidates[b]-segment) >= 1e-11*std::max(1.,dt[
                    i]))
                    continue;
                // Snow is linear within a forcing segment. Correct each
                // tangent for the change of the snow-depletion time.
                if constexpr (std::is_same_v<S,Tangent>) {
                    if (fixed_method) {
                        active[b] = false;
                        y[b] = 0.;
                        snow_rates();
                        continue;
                    }
                    std::vector<S> before(y.size()), after(y.size());
                    rhs(y,before);
                    const auto slopes = y[b].slope;
                    active[b] = false;
                    y[b] = 0.;
                    snow_rates();
                    rhs(y,after);
                    for (std::size_t j = 0; j < y.size(); ++j) {
                        const double jump = primal(before[j])-primal(after[j
                            ]);
                        for (std::size_t k = 0; k < parameter_count; ++k)
                            y[j].slope[k] -= jump*slopes[k]/primal(before[b]
                                );
                    }
                    y[b] = 0.;
                } else {
                    active[b] = false;
                    y[b] = 0.;
                    snow_rates();
                }
            }
            if (dt[i]-time < 1e-14) break;
        }
        const double limit = std::max(1e-7,100*atol);
        for (int j = 0; j < (pilot ? 6 : ns); ++j)
            if (!std::isfinite(primal(y[j])) || primal(y[j]) < -limit)
                throw std::runtime_error("Invalid CaRT storage");
        if (primal(y[5]) > primal(c)+limit)
            throw std::runtime_error("CaRT soil exceeds capacity");
    };
    auto pilot = initialize(true,0.);
    S recharge = 0;
    double next_step = .001;
    // The pilot only supplies recharge over the training mask.
    std::size_t pilot_end = n;
    while (pilot_end > 0 && !f[4*(pilot_end-1)+3]) --pilot_end;
    for (std::size_t i = 0; i < pilot_end; ++i) {
        const S before = pilot[6];
        advance(pilot,i,true,0.,next_step);
        if (f[4*i+3]) {
            training_days += dt[i];
            recharge += pilot[6]-before;
        }
    }
    if (training_days <= 0.) throw std::invalid_argument(
        "Empty CaRT training mask");
    const S gbar = maximum(recharge/training_days,1e-6);
    auto y = initialize(false,gbar);
    const S aquifer_scale = gbar*tau*tau;
    next_step = .001;
    const auto ledger = [&] {
        double storage = 0.;
        for (int b = 0; b < 5; ++b) storage += a[b]*primal(y[b]);
        for (int k = 5; k < ns; ++k) storage += primal(y[k]);
        return storage-primal(y[ns+3]-y[ns+1]-y[ns+2]-y[ns]);
    };
    for (std::size_t i = 0; i < n; ++i) {
        const S before = y[ns];
        const double old_ledger = balance ? ledger() : 0.;
        advance(y,i,false,aquifer_scale,next_step);
        const S flow = y[ns]-before;
        q[i] = primal(flow);
        if constexpr (std::is_same_v<S,Tangent>)
            std::copy(flow.slope.begin(),flow.slope.end(),jac+
                parameter_count*i);
        if (balance) balance[i] = ledger()-old_ledger;
    }
}
} // namespace cartplus::continuous

/**
 * cart_continuous_run: Native model entry point for interval fluxes and
 * requested derivatives.
 *
 * INPUT ARGUMENTS:
 * n : Number of intervals or state entries in this helper.
 * f : Forcing [P, T, PET, training] in the caller layout.
 * z : Five elevation-band heights (m).
 * a : Five nonnegative band fractions, summing to one.
 * p : 11 CaRT calibration coordinates, including logarithmic entries.
 * dt : Integration or observation interval duration (days).
 * initial : Initial stores (mm), or None/empty to use the spin-up policy.
 * m : Positive routing refinement: stages per X4 time base.
 * rtol : Positive relative integration-error tolerance (-).
 * atol : Positive absolute integration-error tolerance for augmented
 * states.
 * solver : Integration method or solver-control object.
 * q : Discharge or flux values (mm per observation interval).
 * jac : Discharge Jacobian in the selected parameter coordinates.
 * balance : Balance used by this numerical operator.
 *
 * OUTPUT ARGUMENTS:
 * Status code (0 success, nonzero failure); supplied output buffers are
 * filled in observation order.
 */
extern "C" __declspec(dllexport) int cart_continuous_run(
    std::size_t n, const double* f, const double* z, const double* a,
    const double* p, const double* dt, const double* initial, int m,
    double rtol, double atol, int solver, double* q, double* jac,
    double* balance) {
    try {
        using namespace cartplus;
        if (jac) {
            std::array<Tangent,parameter_count> params;
            for (std::size_t k = 0; k < params.size(); ++k) {
                params[k] = p[k];
                params[k].slope[k] = 1.;
            }
            continuous::run_cart(n,f,z,a,params.data(),dt,initial,m,rtol,
                atol,
                                 solver,q,jac,balance);
        } else {
            continuous::run_cart(n,f,z,a,p,dt,initial,m,rtol,atol,
                                 solver,q,jac,balance);
        }
        return 0;
    } catch (...) { return 1; }
}

/**
 * cart_configured_run: Native model entry point for interval fluxes and
 * requested derivatives.
 *
 * INPUT ARGUMENTS:
 * n : Number of intervals or state entries in this helper.
 * f : Forcing [P, T, PET, training] in the caller layout.
 * z : Five elevation-band heights (m).
 * a : Five nonnegative band fractions, summing to one.
 * p : 11 CaRT calibration coordinates, including logarithmic entries.
 * dt : Integration or observation interval duration (days).
 * initial : Initial stores (mm), or None/empty to use the spin-up policy.
 * m : Positive routing refinement: stages per X4 time base.
 * rtol : Positive relative integration-error tolerance (-).
 * atol : Positive absolute integration-error tolerance for augmented
 * states.
 * solver : Integration method or solver-control object.
 * q : Discharge or flux values (mm per observation interval).
 * jac : Discharge Jacobian in the selected parameter coordinates.
 * balance : Balance used by this numerical operator.
 * max_step : Maximum/fixed solver step (days), or None for backend default.
 * min_step : Nonnegative lower solver-step bound (days).
 *
 * OUTPUT ARGUMENTS:
 * Status code (0 success, nonzero failure); supplied output buffers are
 * filled in observation order.
 */
extern "C" __declspec(dllexport) int cart_configured_run(
    std::size_t n, const double* f, const double* z, const double* a,
    const double* p, const double* dt, const double* initial, int m,
    double rtol, double atol, int solver, double* q, double* jac,
    double* balance, double max_step, double min_step) {
    try {
        cartplus::continuous::StepLimits limits(max_step,min_step);
        return cart_continuous_run(n,f,z,a,p,dt,initial,m,rtol,atol,solver,q
            ,jac,balance);
    } catch (...) { return 1; }
}

