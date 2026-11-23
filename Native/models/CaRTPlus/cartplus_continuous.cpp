// ========================================================================
// SAGEhydrology | Model source
// File: Native/models/CaRTPlus/cartplus_continuous.cpp
//
// CaRT / CaRT+ numerical model and supporting routines.
//
// SAGE adaptations: Jasper A. Vrugt, University of California, Irvine.
// Upstream model: Xia (2026). Preserve the adjacent license notices.
// ========================================================================

// Fully coupled CaRT+ reference, with the same rate laws as the Python
// port.
#include "cartplus_catchment.hpp"
#include "cartplus_runge_kutta.hpp"


namespace cartplus::continuous {

/**
 * excess_rate: Evaluate smooth infiltration-excess runoff.
 *
 * INPUT ARGUMENTS:
 * water : Available water or store depth (mm).
 * capacity : Store or infiltration capacity in this operator.
 *
 * OUTPUT ARGUMENTS:
 * Numerical result or flux structure in the argument time base; non-const
 * state references are advanced in place.
 */
template<class S>
S excess_rate(S water, S capacity) {
    const S x = water/capacity;
    if (x < .001)
        return capacity*(x*x/2.-pow(x,3.)/6.+pow(x,4.)/24.-pow(x,5.)/120.);
    return water+capacity*expm1(-x);
}

/**
 * production: Evaluate production-store runoff, evaporation, and recharge.
 *
 * INPUT ARGUMENTS:
 * store : Current water store (mm).
 * rain : Liquid input depth/rate on this operator time base.
 * pet : Potential evaporation depth/rate on this operator time base.
 * cap : Store capacity (mm).
 * shape : Dimensionless storage/infiltration shape parameter.
 * perc : Percolation coefficient in the model parameterization.
 * net : Whether to remove interception before the production update.
 * transition : Transition used by this numerical operator.
 *
 * OUTPUT ARGUMENTS:
 * Numerical result or flux structure in the argument time base; non-const
 * state references are advanced in place.
 */
template<class S>
std::array<S,3> production(S store, S rain, S pet, S cap, S shape, S perc,
                           bool net, double transition) {
    const S wet = bounded(store/cap,0.,1.);
    const S deficit = 1.-wet;
    const S intercept = net ? minimum(rain,pet) : S(0);
    const S nu = shape/(1.+shape);
    const S retain = deficit < transition
        ? deficit/transition*exp(nu*std::log(transition))
        : exp(nu*log(deficit));
    return {(rain-intercept)*(1.-retain),
        intercept+(pet-intercept)*(1.-deficit*deficit),
        .25*maximum(store,0.)*log1p(perc*pow(wet,4.))};
}

/**
 * zone: Evaluate near-stream runoff, evaporation, and drainage.
 *
 * INPUT ARGUMENTS:
 * store : Current water store (mm).
 * rain : Liquid input depth/rate on this operator time base.
 * pet : Potential evaporation depth/rate on this operator time base.
 * cap : Store capacity (mm).
 * days : Residence time (days).
 *
 * OUTPUT ARGUMENTS:
 * Numerical result or flux structure in the argument time base; non-const
 * state references are advanced in place.
 */
template<class S>
std::array<S,3> zone(S store, S rain, S pet, S cap, S days) {
    const S wet = bounded(store/cap,0.,1.);
    const S intercept = minimum(rain,pet);
    return {(rain-intercept)*wet*wet,
        intercept+(pet-intercept)*wet*(2.-wet), maximum(store,0.)/days};
}

/**
 * sigmoid: Evaluate the dimensionless logistic response.
 *
 * INPUT ARGUMENTS:
 * x : Numerical operand; scalar, vector, or forward-mode scalar.
 *
 * OUTPUT ARGUMENTS:
 * Numerical result or flux structure in the argument time base; non-const
 * state references are advanced in place.
 */
template<class S>
S sigmoid(S x) {
    if (!(x < S(0))) return 1./(1.+exp(-x));
    const S e = exp(x);
    return e/(1.+e);
}

struct Rk2Mesh {
    std::vector<double> steps;
    std::vector<std::size_t> interval_ends;
    double rejected[2]{};
};

/**
 * simulate: Simulate model fluxes and optional discharge sensitivities.
 *
 * INPUT ARGUMENTS:
 * n : Number of intervals or state entries in this helper.
 * f : Forcing [P, T, PET, training] in the caller layout.
 * z : Five elevation-band heights (m).
 * area : Five nonnegative band fractions, summing to one.
 * p : 42 physical CaRT+ parameters in native order.
 * durations : Positive observation interval durations (days).
 * initial : Initial stores (mm), or None/empty to use the spin-up policy.
 * m : Positive routing refinement: stages per X4 time base.
 * opt : Native numerical-option array in the wrapper order.
 * q : Discharge or flux values (mm per observation interval).
 * jac : Discharge Jacobian in the selected parameter coordinates.
 * balance : Balance used by this numerical operator.
 * solver : Integration method or solver-control object.
 * record : Optional collector of accepted integration steps.
 * replay : Optional accepted-step sequence to replay.
 * stats : Optional integration-count and diagnostic output buffer.
 *
 * OUTPUT ARGUMENTS:
 * No return value; updates referenced states and output buffers.
 */
template<class S>
void simulate(std::size_t n, const double* f, const double* z,
              const double* area, const S* p, const double* durations,
              const double* initial, int m, const double* opt,
              double* q, double* jac, double* balance, int solver,
              Rk2Mesh* record = nullptr, const Rk2Mesh* replay = nullptr,
                  double* stats = nullptr) {
    const int ns = 14+3*m;
    int main_size=ns+7;
    const auto eta9 = exits(m,false), eta1 = exits(m,true);
    const S lambda = m/p[33], upland = 1.-p[13];
    const double availability = opt[3], transition = opt[4];
#ifdef CART_CACHE_CONSTANTS
    const S inverse_surface_time = 1./p[7], inverse_vegetation_time = 1./p[
        24];
#endif
    double train_days = 0, rainfall = 0, elevation = 0;
    double time = 0, year = 0;
    S veg0 = 0;
    std::vector<S> vegetation(n);
    for (int b = 0; b < 5; ++b) elevation += area[b]*z[b];
    for (std::size_t i = 0; i < n; ++i) {
        const double dt = durations[i];
        if (f[4*i+3]) {
            train_days += dt;
            rainfall += f[4*i];
        }
        vegetation[i] = sigmoid((f[4*i+1]-p[22])/p[23]);
        const double w = std::min(dt,std::max(365.-time,0.));
        veg0 += w*vegetation[i];
        year += w;
        time += dt;
    }
    if (train_days <= 0) throw std::invalid_argument("Empty training mask");
    const double mean_rain = std::max(rainfall/train_days,.1);
    const auto initialize = [&](bool pilot, S reference) {
        std::vector<S> y(pilot ? 13 : ns+7);
        if (initial) {
            std::copy(initial,initial+(pilot ? 12 : ns),y.begin());
        } else {
            y[5]=.5*p[6]; y[6]=.7*p[8]; y[7]=.3*p[14]; y[8]=.3*p[25];
            y[9]=.5*p[16]; y[10]=mean_rain; y[11]=veg0/year;
            if (!pilot) { y[12]=reference; y[13]=.5*p[34]; }
        }
        return y;
    };
    const auto storage = [&](const std::vector<S>& y) {
        double value = primal(y[5]+upland*y[6]+p[13]*y[7]+y[9]+y[12]+y[13]);
        for (int b = 0; b < 5; ++b) value += area[b]*primal(y[b]);
        for (int k = 14; k < ns; ++k) value += primal(y[k]);
        return value;
    };
    double next_step = .001;
    const bool exact_vegetation = solver == 4;
    const bool analytic_routing = solver == 5;
    const bool analytic_quick = solver == 8;
    std::vector<S> routing_terms(analytic_routing ? 2*m : 0);
    struct QuickStage { S coefficient, release; };
    struct RoutingStage {
        S inlet, bypass, connectivity, groundwater, near_release;
    };
    Rk2Workspace<S> pilot_work(solver == 1 || solver >= 3 ? 13 : 0);
    Rk2Workspace<S> full_work(solver == 1 || solver >= 3 ? main_size : 0);
    std::size_t mesh_interval = 0, mesh_step = 0;
    const auto advance = [&](std::vector<S>& y, std::size_t i,
                             bool pilot, S reference) {
        const double dt = durations[i], precip = f[4*i]/dt;
        const double temp = f[4*i+1], pet = f[4*i+2]/dt;
        std::array<S,5> snowfall{}, melt_demand{};
        double liquid_rain = 0, cold_area = 0;
        for (int b = 0; b < 5; ++b) {
            const S warm = temp-p[5]*((z[b]-elevation)/1000.)-p[0];
            if (warm < 0.) {
                snowfall[b] = p[2]*precip;
                cold_area += area[b];
            } else {
                melt_demand[b] = p[1]*warm+p[4]*pet;
                liquid_rain += area[b]*precip;
            }
        }
        const S effective = precip+(p[2]-1.)*precip*cold_area;
        const S sub_demand = p[3]*pet;
        const double tk = std::max(temp,0.)+273.15;
        const S thermal = -std::log(10.)*247.8*(1./(tk-140.)-1./153.15)
            -p[39]*sigmoid((p[37]-temp)/p[38]);
        const auto downstream = [&](const S& storage, const RoutingStage&
            stage,
                                     const S& u9, const S& u1, std::vector<S
                                >& out,
                                     QuickStage* quick_stage = nullptr) {
            const S quick = maximum(storage,0.);
            const S exchange = p[35]*pow(quick/p[34],3.5);
            S ir = u9+.9*(p[36]+exchange);
            if (quick <= 0.) ir = maximum(ir,0.);
            const S qb = maximum(u1+.1*(p[36]+exchange),0.);
            const S modulation = log(maximum(stage.connectivity,1e-12))+
                thermal
                +p[40]*log(maximum(quick,1e-3)/p[34]);
            const S qr = .25*quick*pow(quick/p[34],4.)*exp(6.*tanh(
                modulation/6.));
            if (quick_stage)
                *quick_stage = {.25*exp(6.*tanh(modulation/6.))/pow(p[34],4.
                    ),qr};
            const S preloss = qr+qb+stage.near_release+stage.groundwater;
            const S reported = excess_rate(preloss,p[41]);
            out[13]=ir-qr;
            out[ns]=reported;
            out[ns+2]=ir+qb-u9-u1;
            out[ns+5]=preloss-reported;
            out[ns+6]=preloss;
        };
        const auto rhs = [&](const std::vector<S>& x, std::vector<S>& out,
                             RoutingStage* routing = nullptr, QuickStage*
                                quick_stage = nullptr) {
            S liquid = liquid_rain, snow_et = 0;
            for (int b = 0; b < 5; ++b) {
                const S available = maximum(x[b],0.)/availability;
                const S sub = minimum(sub_demand,available);
                const S melt = minimum(maximum(available-sub,0.),melt_demand
                    [b]);
                out[b] = snowfall[b]-sub-melt;
                liquid += area[b]*melt;
                snow_et += area[b]*sub;
            }
            const S demand = pet*(1.-p[21]*(1.-bounded(x[11],0.,1.)));
            const S eh = minimum(demand,maximum(x[5],0.)/availability);
#ifdef CART_CACHE_CONSTANTS
            const S qh = maximum(x[5],0.)*inverse_surface_time;
#else
            const S qh = maximum(x[5],0.)/p[7];
#endif
            const S fh = qh+maximum(liquid-eh-qh,0.)
                *bounded((x[5]/p[6]-(1.-transition))/transition,0.,1.);
            const S ep = maximum(demand-eh,0.);
            const S capacity = maximum(
                p[11]*exp(-p[12]*maximum(x[10],0.)/mean_rain),
                std::max(opt[2],1e-12));
            const S qx = excess_rate(maximum(fh-ep,0.),capacity);
            const S ps = fh-qx;
            const auto soil = production(x[6],ps,ep,p[8],p[10],p[9],true,
                transition);
            const auto near = zone(x[7],ps,ep,p[14],p[15]);
            const auto wet = zone(x[8],ps,ep,p[25],p[26]);
            const S ez = upland*soil[1]+p[13]*near[1];
            const S fast = qx+upland*soil[0]+p[13]*near[0];
            const auto dep = production(x[9],fast,maximum(ep-ez,0.),
                                        p[16],p[17],p[18],false,transition);
            const S recharge = upland*soil[2]+dep[2];
            out[5]=liquid-eh-fh; out[6]=ps-soil[0]-soil[1]-soil[2];
            out[7]=ps-near[0]-near[1]-near[2]; out[8]=ps-wet[0]-wet[1]-wet[2
                ];
            out[9]=fast-dep[0]-dep[1]-dep[2]; out[10]=(fh-x[10])/3.;
#ifdef CART_CACHE_CONSTANTS
            out[11]=exact_vegetation ? S(0) : (vegetation[i]-x[11])*
                inverse_vegetation_time;
#else
            out[11]=exact_vegetation ? S(0) : (vegetation[i]-x[11])/p[24];
#endif
            if (pilot) { out[12]=recharge; return; }
            const S gw = maximum(x[12],0.);
            const S qg = gw*gw/(reference*p[19]);
            const S eg = p[20]*gw*maximum(ep-ez-dep[1],0.)/(gw+reference);
            const S c = p[27]+(1.-p[27])*sigmoid((x[8]/p[25]-p[28]
                +p[30]*bounded(x[6]/p[8],0.,1.)+p[31]*(gw/reference-1.))/p[
                    29]);
            const S bypass = p[32]*c*dep[0];
            const RoutingStage stage{dep[0]-bypass,bypass,c,qg,p[13]*near[2]
                };
            out[12]=recharge-qg-eg;
            out[ns+1]=eh+ez+dep[1]+eg+snow_et;
            out[ns+3]=effective;
            out[ns+4]=recharge;
            if (routing) {
                *routing = stage;
                return;
            }
            const auto cascade = [&](int offset, const auto& eta, S incoming
                ) {
                S collected = 0;
                for (std::size_t j = 0; j < eta.size(); ++j) {
                    const S release = lambda*x[offset+j];
                    out[offset+j] = incoming-release;
                    collected += eta[j]*release;
                    incoming = (1.-eta[j])*release;
                }
                return collected;
            };
            const S u9 = cascade(14,eta9,stage.inlet);
            const S u1 = cascade(14+m,eta1,stage.bypass);
            downstream(x[13],stage,u9,u1,out,quick_stage);
        };
        if (replay) {
            const auto stop = replay->interval_ends.at(mesh_interval++);
            replay_rk2(y,rhs,pilot ? pilot_work : full_work,
                       replay->steps.data()+mesh_step,replay->steps.data()+
                           stop);
            mesh_step = stop;
        } else if (record) {
            integrate_rk2(y,rhs,dt,opt[0],opt[1],next_step,
                          pilot ? pilot_work : full_work,0,Rk2NoExactUpdate{
                              },
                          Rk2DefaultTrial{},&record->steps,&record->rejected
                              [!pilot]);
            record->interval_ends.push_back(record->steps.size());
        } else if (analytic_routing && !pilot) {
            const auto trial = [&](const std::vector<S>& start, auto& work,
                double h) {
                auto& [first, second, euler, heun] = work;
                RoutingStage at_start, at_end;
                rhs(start,first,&at_start);
                for (int j = 0; j < 13; ++j)
                    euler[j] = start[j]+h*first[j];
                rhs(euler,second,&at_end);

                const auto route = [&](int offset, const auto& eta, const S&
                    u0,
                                       const S& u1, std::vector<S>& target)
                                {
                    const S volume = .5*h*(u0+u1);
                    S released = volume;
                    const double pieces = std::ceil(h*primal(lambda)/.5);
                    if (!std::isfinite(pieces) || pieces > 100000)
                        throw std::runtime_error(
                            "Analytic CaRT+ routing step too large");
                    const int parts = std::max(1,static_cast<int>(pieces));
                    const double substep = h/parts;
                    const S rate = substep*lambda;
                    const S change = (u1-u0)/parts;
                    for (std::size_t j = 0; j < eta.size(); ++j)
                        target[offset+j] = start[offset+j];
                    for (int part = 0; part < parts; ++part) {
                        const S input = u0+part*change;
                        for (std::size_t j = 0; j < eta.size(); ++j) {
                            const S incoming = j == 0 ? substep*input
                                : (1.-eta[j-1])*rate*target[offset+j-1];
                            routing_terms[j] = incoming-rate*target[offset+j
                                ];
                        }
                        for (std::size_t j = 0; j < eta.size(); ++j) {
                            target[offset+j] += routing_terms[j];
                            released -= routing_terms[j];
                        }
                        // The scaled cascade generator has row norm at most
                        // one.
                        // Sum its exponential action with the linear inlet
                        // term.
                        for (int degree = 2; degree <= 18; ++degree) {
                            const S coefficient = rate/degree;
                            for (std::size_t j = eta.size()-1; j > 0; --j)
                                routing_terms[j] = coefficient*((1.-eta[j-1]
                                )*routing_terms[j-1]-routing_terms[j]);
                            routing_terms[0] = -coefficient*routing_terms[0]
                                ;
                            if (degree == 2) routing_terms[0] += .5*substep*
                                change;
                            for (std::size_t j = 0; j < eta.size(); ++j) {
                                target[offset+j] += routing_terms[j];
                                released -= routing_terms[j];
                            }
                        }
                    }
                    return released/h;
                };

                const S e9 = route(14,eta9,at_start.inlet,at_start.inlet,
                    euler);
                const S e1 = route(14+m,eta1,at_start.bypass,at_start.bypass
                    ,euler);
                downstream(start[13],at_start,e9,e1,first);
                euler[13] = start[13]+h*first[13];
                for (int j = ns; j < ns+7; ++j)
                    euler[j] = start[j]+h*first[j];

                const S r9 = route(14,eta9,at_start.inlet,at_end.inlet,heun)
                    ;
                const S r1 = route(14+m,eta1,at_start.bypass,at_end.bypass,
                    heun);
                // Use the exact mean cascade releases in both downstream
                // stages.
                // Their trapezoidal integrals then match the cascade water
                // loss.
                downstream(start[13],at_start,r9,r1,first);
                downstream(start[13]+h*first[13],at_end,r9,r1,second);
                for (int j = 0; j < ns+7; ++j) {
                    if (j >= 14 && j < ns) continue;
                    heun[j] = start[j];
                    weighted_add(heun[j],.5*h,first[j]);
                    weighted_add(heun[j],.5*h,second[j]);
                }
            };
            integrate_rk2(y,rhs,dt,opt[0],opt[1],next_step,full_work,0,
                          Rk2NoExactUpdate{},trial);
        } else if (analytic_quick && !pilot) {
            const auto trial = [&](const std::vector<S>& start, auto& work,
                double h) {
                auto& [first, second, euler, heun] = work;
                QuickStage left, right;
                rhs(start,first,nullptr,&left);
                for (std::size_t j = 0; j < start.size(); ++j) {
                    euler[j] = start[j];
                    weighted_add(euler[j],h,first[j]);
                }
                rhs(euler,second,nullptr,&right);
                for (std::size_t j = 0; j < start.size(); ++j) {
                    heun[j] = start[j];
                    weighted_add(heun[j],.5*h,first[j]);
                    weighted_add(heun[j],.5*h,second[j]);
                }
                // Split the mean non-recession flow around the exact
                // power-law recession.
                const S half_input = .25*h*(first[13]+left.release+second[13
                    ]+right.release);
                const S available = start[13]+half_input;
                if (available < 0.) return;
                const S coefficient = .5*(left.coefficient+right.coefficient
                    );
                const S square = available*available;
                const S release = -available*expm1(-.25*log1p(4.*h*
                    coefficient*square*square));
                if (available-release+half_input < 0.) return;
                const S p0 = first[ns+6]-left.release+release/h;
                const S p1 = second[ns+6]-right.release+release/h;
                const S q0 = excess_rate(p0,p[41]), q1 = excess_rate(p1,p[41
                    ]);
                heun[13] = available-release+half_input;
                heun[ns] = start[ns]+.5*h*(q0+q1);
                heun[ns+5] = start[ns+5]+.5*h*(p0-q0+p1-q1);
                heun[ns+6] = start[ns+6]+.5*h*(p0+p1);
            };
            integrate_rk2(y,rhs,dt,opt[0],opt[1],next_step,full_work,0,
                          Rk2NoExactUpdate{},trial,nullptr,stats ? stats+3 :
                              nullptr,
                          stats ? stats+2 : nullptr);
        } else if (exact_vegetation) {
            const auto exact_vegetation_step = [&](const S& start, double h)
                {
                return start + (start-vegetation[i])*expm1(-h/p[24]);
            };
            integrate_rk2(y,rhs,dt,opt[0],opt[1],next_step,
                          pilot ? pilot_work : full_work,11,
                              exact_vegetation_step);
        } else if (solver == 1 || (solver >= 3 && solver <= 8)) {
            double* counts = stats ? stats+2*(!pilot) : nullptr;
            integrate_rk2(y,rhs,dt,opt[0],opt[1],next_step,pilot ?
                pilot_work : full_work,
                          0,Rk2NoExactUpdate{},Rk2DefaultTrial{},nullptr,
                          counts ? counts+1 : nullptr,counts);
        } else {
            integrate_segment(y,rhs,dt,opt[0],opt[1],next_step,solver,pilot
                ? 12 : ns);
        }
        for (int j = 0; j < (pilot ? 12 : ns); ++j)
            if (!std::isfinite(primal(y[j])) || primal(y[j]) < -std::max(
                1e-7,100*opt[1]))
                throw std::runtime_error("Invalid coupled CaRT+ state");
    };
    auto pilot = initialize(true,0.);
    S recharge = 0;
    std::size_t pilot_end = n;
    if (solver == 1 || solver >= 7)
        while (pilot_end > 0 && !f[4*(pilot_end-1)+3]) --pilot_end;
    for (std::size_t i = 0; i < pilot_end; ++i) {
        const S before = pilot[12];
        advance(pilot,i,true,0.);
        if (f[4*i+3]) recharge += pilot[12]-before;
    }
    const S reference = maximum(recharge/train_days,1e-6)*p[19];
    auto y = initialize(false,reference);
    next_step = .001;
    for (std::size_t i = 0; i < n; ++i) {
        const S before = y[ns];
        const double volume = storage(y);
        const auto ledger = [&] { return primal(y[ns+3]+y[ns+2]-y[ns+1]-y[ns
            +5]-y[ns]); };
        const double old_ledger = ledger();
        advance(y,i,false,reference);
        if (q) q[i] = primal(y[ns])-primal(before);
        if constexpr (std::is_same_v<S,Tangent>)
            if (jac)
                for (std::size_t k = 0; k < parameter_count; ++k)
                    jac[parameter_count*i+k] = y[ns].slope[k]-before.slope[k];
        if (balance) balance[i] = storage(y)-volume-(ledger()-old_ledger);
    }
}
} // namespace cartplus::continuous

/**
 * run_continuous: Dispatch the coupled CaRT+ primal/tangent integrator.
 *
 * INPUT ARGUMENTS:
 * n : Number of intervals or state entries in this helper.
 * f : Forcing [P, T, PET, training] in the caller layout.
 * z : Five elevation-band heights (m).
 * a : Five nonnegative band fractions, summing to one.
 * p : 42 physical CaRT+ parameters in native order.
 * durations : Positive observation interval durations (days).
 * initial : Initial stores (mm), or None/empty to use the spin-up policy.
 * m : Positive routing refinement: stages per X4 time base.
 * opt : Native numerical-option array in the wrapper order.
 * seed : Primal-to-parameter derivative seed vector.
 * q : Discharge or flux values (mm per observation interval).
 * jac : Discharge Jacobian in the selected parameter coordinates.
 * balance : Balance used by this numerical operator.
 * solver : Integration method or solver-control object.
 * stats : Optional integration-count and diagnostic output buffer.
 *
 * OUTPUT ARGUMENTS:
 * Numerical result or flux structure in the argument time base; non-const
 * state references are advanced in place.
 */
static int run_continuous(
    std::size_t n, const double* f, const double* z, const double* a,
    const double* p, const double* durations, const double* initial,
    int m, const double* opt, const double* seed, double* q, double* jac,
    double* balance, int solver, double* stats = nullptr,
    bool keep_tangents = false) {
    try {
        if (stats && solver >= 7) std::fill_n(stats,4,0.);
        cartplus::continuous::Rk2Mesh mesh;
        const bool record = solver == 6 && (jac || stats || keep_tangents);
        if (record) {
            mesh.interval_ends.reserve(2*n);
            cartplus::continuous::simulate(n,f,z,a,p,durations,initial,m,opt
                ,
                                            q,nullptr,balance,solver,&mesh);
            if (stats) {
                stats[0] = static_cast<double>(mesh.interval_ends[n-1]);
                stats[1] = mesh.rejected[0];
                stats[2] = static_cast<double>(mesh.steps.size()-mesh.
                    interval_ends[n-1]);
                stats[3] = mesh.rejected[1];
            }
        }
        if (jac || keep_tangents) {
            std::array<cartplus::Tangent,42> params;
            std::size_t column = 0;
            for (int k = 0; k < 42; ++k) {
                params[k] = p[k];
#ifdef CARTPLUS_CALIBRATED_ONLY
                if (k == 22 || k == 23 || k == 24 || k == 28 ||
                    k == 29 || k == 37 || k == 38) continue;
#endif
                params[k].slope[column++] = seed[k];
            }
            cartplus::continuous::simulate(n,f,z,a,params.data(),durations,
                                            initial,m,opt,record ? nullptr :
                                q,jac,
                                            record ? nullptr : balance,
                                solver,nullptr,
                                            record ? &mesh : nullptr,solver
                                >= 7 ? stats : nullptr);
        } else if (!record) {
            cartplus::continuous::simulate(n,f,z,a,p,durations,initial,m,opt
                ,
                                            q,jac,balance,solver,nullptr,
                                nullptr,
                                            solver >= 7 ? stats : nullptr);
        }
        return 0;
    } catch (...) { return 1; }
}

/**
 * cartplus_continuous_run: Native model entry point for interval fluxes and
 * requested derivatives.
 *
 * INPUT ARGUMENTS:
 * n : Number of intervals or state entries in this helper.
 * f : Forcing [P, T, PET, training] in the caller layout.
 * z : Five elevation-band heights (m).
 * a : Five nonnegative band fractions, summing to one.
 * p : 42 physical CaRT+ parameters in native order.
 * durations : Positive observation interval durations (days).
 * initial : Initial stores (mm), or None/empty to use the spin-up policy.
 * m : Positive routing refinement: stages per X4 time base.
 * opt : Native numerical-option array in the wrapper order.
 * seed : Primal-to-parameter derivative seed vector.
 * q : Discharge or flux values (mm per observation interval).
 * jac : Discharge Jacobian in the selected parameter coordinates.
 * balance : Balance used by this numerical operator.
 *
 * OUTPUT ARGUMENTS:
 * Status code (0 success, nonzero failure); supplied output buffers are
 * filled in observation order.
 */
extern "C" __declspec(dllexport) int cartplus_continuous_run(
    std::size_t n, const double* f, const double* z, const double* a,
    const double* p, const double* durations, const double* initial,
    int m, const double* opt, const double* seed, double* q, double* jac,
    double* balance) {
    return run_continuous(n,f,z,a,p,durations,initial,m,opt,seed,q,jac,
        balance,0);
}

/**
 * cartplus_continuous_rk2_run: Native model entry point for interval fluxes
 * and requested derivatives.
 *
 * INPUT ARGUMENTS:
 * n : Number of intervals or state entries in this helper.
 * f : Forcing [P, T, PET, training] in the caller layout.
 * z : Five elevation-band heights (m).
 * a : Five nonnegative band fractions, summing to one.
 * p : 42 physical CaRT+ parameters in native order.
 * durations : Positive observation interval durations (days).
 * initial : Initial stores (mm), or None/empty to use the spin-up policy.
 * m : Positive routing refinement: stages per X4 time base.
 * opt : Native numerical-option array in the wrapper order.
 * seed : Primal-to-parameter derivative seed vector.
 * q : Discharge or flux values (mm per observation interval).
 * jac : Discharge Jacobian in the selected parameter coordinates.
 * balance : Balance used by this numerical operator.
 *
 * OUTPUT ARGUMENTS:
 * Status code (0 success, nonzero failure); supplied output buffers are
 * filled in observation order.
 */
extern "C" __declspec(dllexport) int cartplus_continuous_rk2_run(
    std::size_t n, const double* f, const double* z, const double* a,
    const double* p, const double* durations, const double* initial,
    int m, const double* opt, const double* seed, double* q, double* jac,
    double* balance) {
    return run_continuous(n,f,z,a,p,durations,initial,m,opt,seed,q,jac,
        balance,1);
}

/**
 * cartplus_continuous_rk23_run: Native model entry point for interval
 * fluxes and requested derivatives.
 *
 * INPUT ARGUMENTS:
 * n : Number of intervals or state entries in this helper.
 * f : Forcing [P, T, PET, training] in the caller layout.
 * z : Five elevation-band heights (m).
 * a : Five nonnegative band fractions, summing to one.
 * p : 42 physical CaRT+ parameters in native order.
 * durations : Positive observation interval durations (days).
 * initial : Initial stores (mm), or None/empty to use the spin-up policy.
 * m : Positive routing refinement: stages per X4 time base.
 * opt : Native numerical-option array in the wrapper order.
 * seed : Primal-to-parameter derivative seed vector.
 * q : Discharge or flux values (mm per observation interval).
 * jac : Discharge Jacobian in the selected parameter coordinates.
 * balance : Balance used by this numerical operator.
 *
 * OUTPUT ARGUMENTS:
 * Status code (0 success, nonzero failure); supplied output buffers are
 * filled in observation order.
 */
extern "C" __declspec(dllexport) int cartplus_continuous_rk23_run(
    std::size_t n, const double* f, const double* z, const double* a,
    const double* p, const double* durations, const double* initial,
    int m, const double* opt, const double* seed, double* q, double* jac,
    double* balance) {
    return run_continuous(n,f,z,a,p,durations,initial,m,opt,seed,q,jac,
        balance,2);
}

/**
 * cartplus_continuous_rk2_workspace_run: Native model entry point for
 * interval fluxes and requested derivatives.
 *
 * INPUT ARGUMENTS:
 * n : Number of intervals or state entries in this helper.
 * f : Forcing [P, T, PET, training] in the caller layout.
 * z : Five elevation-band heights (m).
 * a : Five nonnegative band fractions, summing to one.
 * p : 42 physical CaRT+ parameters in native order.
 * durations : Positive observation interval durations (days).
 * initial : Initial stores (mm), or None/empty to use the spin-up policy.
 * m : Positive routing refinement: stages per X4 time base.
 * opt : Native numerical-option array in the wrapper order.
 * seed : Primal-to-parameter derivative seed vector.
 * q : Discharge or flux values (mm per observation interval).
 * jac : Discharge Jacobian in the selected parameter coordinates.
 * balance : Balance used by this numerical operator.
 *
 * OUTPUT ARGUMENTS:
 * Status code (0 success, nonzero failure); supplied output buffers are
 * filled in observation order.
 */
extern "C" __declspec(dllexport) int cartplus_continuous_rk2_workspace_run(
    std::size_t n, const double* f, const double* z, const double* a,
    const double* p, const double* durations, const double* initial,
    int m, const double* opt, const double* seed, double* q, double* jac,
    double* balance) {
    return run_continuous(n,f,z,a,p,durations,initial,m,opt,seed,q,jac,
        balance,3);
}

/**
 * cartplus_continuous_rk2_exact_vegetation_run: Native model entry point
 * for interval fluxes and requested derivatives.
 *
 * INPUT ARGUMENTS:
 * n : Number of intervals or state entries in this helper.
 * f : Forcing [P, T, PET, training] in the caller layout.
 * z : Five elevation-band heights (m).
 * a : Five nonnegative band fractions, summing to one.
 * p : 42 physical CaRT+ parameters in native order.
 * durations : Positive observation interval durations (days).
 * initial : Initial stores (mm), or None/empty to use the spin-up policy.
 * m : Positive routing refinement: stages per X4 time base.
 * opt : Native numerical-option array in the wrapper order.
 * seed : Primal-to-parameter derivative seed vector.
 * q : Discharge or flux values (mm per observation interval).
 * jac : Discharge Jacobian in the selected parameter coordinates.
 * balance : Balance used by this numerical operator.
 *
 * OUTPUT ARGUMENTS:
 * Status code (0 success, nonzero failure); supplied output buffers are
 * filled in observation order.
 */
extern "C" __declspec(dllexport) int
    cartplus_continuous_rk2_exact_vegetation_run(
    std::size_t n, const double* f, const double* z, const double* a,
    const double* p, const double* durations, const double* initial,
    int m, const double* opt, const double* seed, double* q, double* jac,
    double* balance) {
    return run_continuous(n,f,z,a,p,durations,initial,m,opt,seed,q,jac,
        balance,4);
}

/**
 * cartplus_continuous_rk2_analytic_routing_run: Native model entry point
 * for interval fluxes and requested derivatives.
 *
 * INPUT ARGUMENTS:
 * n : Number of intervals or state entries in this helper.
 * f : Forcing [P, T, PET, training] in the caller layout.
 * z : Five elevation-band heights (m).
 * a : Five nonnegative band fractions, summing to one.
 * p : 42 physical CaRT+ parameters in native order.
 * durations : Positive observation interval durations (days).
 * initial : Initial stores (mm), or None/empty to use the spin-up policy.
 * m : Positive routing refinement: stages per X4 time base.
 * opt : Native numerical-option array in the wrapper order.
 * seed : Primal-to-parameter derivative seed vector.
 * q : Discharge or flux values (mm per observation interval).
 * jac : Discharge Jacobian in the selected parameter coordinates.
 * balance : Balance used by this numerical operator.
 *
 * OUTPUT ARGUMENTS:
 * Status code (0 success, nonzero failure); supplied output buffers are
 * filled in observation order.
 */
extern "C" __declspec(dllexport) int
    cartplus_continuous_rk2_analytic_routing_run(
    std::size_t n, const double* f, const double* z, const double* a,
    const double* p, const double* durations, const double* initial,
    int m, const double* opt, const double* seed, double* q, double* jac,
    double* balance) {
    return run_continuous(n,f,z,a,p,durations,initial,m,opt,seed,q,jac,
        balance,5);
}

/**
 * cartplus_continuous_rk2_replay_run: Native model entry point for interval
 * fluxes and requested derivatives.
 *
 * INPUT ARGUMENTS:
 * n : Number of intervals or state entries in this helper.
 * f : Forcing [P, T, PET, training] in the caller layout.
 * z : Five elevation-band heights (m).
 * a : Five nonnegative band fractions, summing to one.
 * p : 42 physical CaRT+ parameters in native order.
 * durations : Positive observation interval durations (days).
 * initial : Initial stores (mm), or None/empty to use the spin-up policy.
 * m : Positive routing refinement: stages per X4 time base.
 * opt : Native numerical-option array in the wrapper order.
 * seed : Primal-to-parameter derivative seed vector.
 * q : Discharge or flux values (mm per observation interval).
 * jac : Discharge Jacobian in the selected parameter coordinates.
 * balance : Balance used by this numerical operator.
 * stats : Optional integration-count and diagnostic output buffer.
 *
 * OUTPUT ARGUMENTS:
 * Status code (0 success, nonzero failure); supplied output buffers are
 * filled in observation order.
 */
extern "C" __declspec(dllexport) int cartplus_continuous_rk2_replay_run(
    std::size_t n, const double* f, const double* z, const double* a,
    const double* p, const double* durations, const double* initial,
    int m, const double* opt, const double* seed, double* q, double* jac,
    double* balance, double* stats) {
    if (!n) return 1;
    return run_continuous(n,f,z,a,p,durations,initial,m,opt,seed,q,jac,
        balance,6,stats);
}

/**
 * cartplus_continuous_rk2_optimized_run: Native model entry point for
 * interval fluxes and requested derivatives.
 *
 * INPUT ARGUMENTS:
 * n : Number of intervals or state entries in this helper.
 * f : Forcing [P, T, PET, training] in the caller layout.
 * z : Five elevation-band heights (m).
 * a : Five nonnegative band fractions, summing to one.
 * p : 42 physical CaRT+ parameters in native order.
 * durations : Positive observation interval durations (days).
 * initial : Initial stores (mm), or None/empty to use the spin-up policy.
 * m : Positive routing refinement: stages per X4 time base.
 * opt : Native numerical-option array in the wrapper order.
 * seed : Primal-to-parameter derivative seed vector.
 * q : Discharge or flux values (mm per observation interval).
 * jac : Discharge Jacobian in the selected parameter coordinates.
 * balance : Balance used by this numerical operator.
 * stats : Optional integration-count and diagnostic output buffer.
 *
 * OUTPUT ARGUMENTS:
 * Status code (0 success, nonzero failure); supplied output buffers are
 * filled in observation order.
 */
extern "C" __declspec(dllexport) int cartplus_continuous_rk2_optimized_run(
    std::size_t n, const double* f, const double* z, const double* a,
    const double* p, const double* durations, const double* initial,
    int m, const double* opt, const double* seed, double* q, double* jac,
    double* balance, double* stats) {
    return run_continuous(n,f,z,a,p,durations,initial,m,opt,seed,q,jac,
        balance,7,stats);
}

/**
 * cartplus_continuous_rk2_analytic_quick_run: Native model entry point for
 * interval fluxes and requested derivatives.
 *
 * INPUT ARGUMENTS:
 * n : Number of intervals or state entries in this helper.
 * f : Forcing [P, T, PET, training] in the caller layout.
 * z : Five elevation-band heights (m).
 * a : Five nonnegative band fractions, summing to one.
 * p : 42 physical CaRT+ parameters in native order.
 * durations : Positive observation interval durations (days).
 * initial : Initial stores (mm), or None/empty to use the spin-up policy.
 * m : Positive routing refinement: stages per X4 time base.
 * opt : Native numerical-option array in the wrapper order.
 * seed : Primal-to-parameter derivative seed vector.
 * q : Discharge or flux values (mm per observation interval).
 * jac : Discharge Jacobian in the selected parameter coordinates.
 * balance : Balance used by this numerical operator.
 * stats : Optional integration-count and diagnostic output buffer.
 *
 * OUTPUT ARGUMENTS:
 * Status code (0 success, nonzero failure); supplied output buffers are
 * filled in observation order.
 */
extern "C" __declspec(dllexport) int
    cartplus_continuous_rk2_analytic_quick_run(
    std::size_t n, const double* f, const double* z, const double* a,
    const double* p, const double* durations, const double* initial,
    int m, const double* opt, const double* seed, double* q, double* jac,
    double* balance, double* stats) {
    return run_continuous(n,f,z,a,p,durations,initial,m,opt,seed,q,jac,
        balance,8,stats);
}

/**
 * cartplus_configured_run: Native model entry point for interval fluxes and
 * requested derivatives.
 *
 * INPUT ARGUMENTS:
 * n : Number of intervals or state entries in this helper.
 * f : Forcing [P, T, PET, training] in the caller layout.
 * z : Five elevation-band heights (m).
 * a : Five nonnegative band fractions, summing to one.
 * p : 42 physical CaRT+ parameters in native order.
 * dt : Integration or observation interval duration (days).
 * initial : Initial stores (mm), or None/empty to use the spin-up policy.
 * m : Positive routing refinement: stages per X4 time base.
 * opt : Native numerical-option array in the wrapper order.
 * seed : Primal-to-parameter derivative seed vector.
 * q : Discharge or flux values (mm per observation interval).
 * jac : Discharge Jacobian in the selected parameter coordinates.
 * balance : Balance used by this numerical operator.
 * solver : Integration method or solver-control object.
 * max_step : Maximum/fixed solver step (days), or None for backend default.
 * min_step : Nonnegative lower solver-step bound (days).
 *
 * OUTPUT ARGUMENTS:
 * Status code (0 success, nonzero failure); supplied output buffers are
 * filled in observation order.
 */
extern "C" __declspec(dllexport) int cartplus_configured_run(
    std::size_t n, const double* f, const double* z, const double* a,
    const double* p, const double* dt, const double* initial, int m,
    const double* opt, const double* seed, double* q, double* jac,
    double* balance, int solver, double max_step, double min_step) {
    try {
        cartplus::continuous::StepLimits limits(max_step,min_step);
        return run_continuous(n,f,z,a,p,dt,initial,m,opt,seed,q,jac,balance,
            solver);
    } catch (...) { return 1; }
}

// Advance the same augmented state used by training, with discharge as
// the only required output. The adaptive error norm still sees tangents.
extern "C" __declspec(dllexport) int cartplus_configured_predict(
    std::size_t n, const double* f, const double* z, const double* a,
    const double* p, const double* dt, const double* initial, int m,
    const double* opt, const double* seed, double* q, double* jac,
    double* balance, int solver, double max_step, double min_step) {
    try {
        cartplus::continuous::StepLimits limits(max_step,min_step);
        return run_continuous(n,f,z,a,p,dt,initial,m,opt,seed,q,jac,balance,
            solver,nullptr,true);
    } catch (...) { return 1; }
}

