// ========================================================================
// SAGEhydrology | Model source
// File: Native/models/CaRTPlus/cartplus_uniform.hpp
//
// CaRT / CaRT+ numerical model and supporting routines.
//
// SAGE adaptations: Jasper A. Vrugt, University of California, Irvine.
// Upstream model: Xia (2026). Preserve the adjacent license notices.
// ========================================================================

// Uniform-input operators adapted from xiaxilin/cartplus (MIT).
// Revision f986af265633a7a980836c2e0f2124077d632420, hourly branch.
// Defaults: prod_lie=true, nsub=4, nsub_slow=1. See LICENSE.
#pragma once
#include "cartplus_catchment.hpp"

namespace cartplus {
/**
 * uniform_release: Compute stable release for uniform interval inputs.
 *
 * INPUT ARGUMENTS:
 * water : Available water or store depth (mm).
 * coefficient : Release-law coefficient in the current interval time base.
 *
 * OUTPUT ARGUMENTS:
 * Numerical result or flux structure in the argument time base; non-const
 * state references are advanced in place.
 */
template<class T> T uniform_release(T water, T coefficient) {
    return -water * expm1(-0.25 * log1p(coefficient * pow(water, 4.0)));
}
/**
 * uniform_psi: Evaluate the stable exponential infiltration-excess factor.
 *
 * INPUT ARGUMENTS:
 * x : Numerical operand; scalar, vector, or forward-mode scalar.
 *
 * OUTPUT ARGUMENTS:
 * Numerical result or flux structure in the argument time base; non-const
 * state references are advanced in place.
 */
template<class T> T uniform_psi(T x) {
    return x < 0.05
        ? x*(0.5-x*(1.0/6-x*(1.0/24-x/120)))
        : 1.0 + expm1(-x)/x;
}
/**
 * uniform_surface: Advance uniform-input surface storage and drainage.
 *
 * INPUT ARGUMENTS:
 * water : Available water or store depth (mm).
 * rain : Liquid input depth/rate on this operator time base.
 * demand : Evaporation demand on this operator time base.
 * capacity : Store or infiltration capacity in this operator.
 * b : Second numerical operand or rate coefficient.
 *
 * OUTPUT ARGUMENTS:
 * Numerical result or flux structure in the argument time base; non-const
 * state references are advanced in place.
 */
template<class T> SurfaceFluxes<T> uniform_surface(
    T& water, T rain, T demand, T capacity, T b) {
    const T net = rain-demand, eb = -expm1(-b);
    const T free = water*(1.0-eb)+net*eb/b;
    const bool over = !(net < 0.0) && free > capacity;
    const T spill = maximum(net-b*capacity, 1e-12);
    const T full_time = minimum(log1p(b*maximum(capacity-water,0.0)/spill)/b
        ,1.0);
    const T empty_time = log1p(b*water/maximum(-net,1e-12))/b;
    const bool empty = net < 0.0 && empty_time < 1.0;
    const T time = over ? full_time : empty ? empty_time : T(1.0);
    const T x = b*time;
    const T drainage = -water*expm1(-x)+net*time*uniform_psi(x)
        + (over ? b*capacity*(1.0-full_time) : T(0.0));
    const T overflow = over ? spill*(1.0-full_time) : T(0.0);
    const T evaporation = empty ? demand*empty_time+rain*(1.0-empty_time) :
        demand;
    water = over ? capacity : empty ? T(0.0) : free;
    return {overflow,drainage,evaporation};
}
/**
 * uniform_pareto: Advance split wet/dry Pareto production operators.
 *
 * INPUT ARGUMENTS:
 * water : Available water or store depth (mm).
 * rain : Liquid input depth/rate on this operator time base.
 * demand : Evaporation demand on this operator time base.
 * capacity : Store or infiltration capacity in this operator.
 * percolation : Percolation coefficient in the model parameterization.
 * shape : Dimensionless storage/infiltration shape parameter.
 * argument 7 : Reserved/unnamed positional argument.
 * net : Whether to remove interception before the production update.
 * dt : Integration or observation interval duration (days).
 * symmetric : Whether to use the symmetric two-time-base hydrograph.
 *
 * OUTPUT ARGUMENTS:
 * Numerical result or flux structure in the argument time base; non-const
 * state references are advanced in place.
 */
template<class T> ProductionFluxes<T> uniform_pareto(
    T& water, T rain, T demand, T capacity, T percolation,
    T shape, double, bool net, double dt=1.0, bool symmetric=false) {
    const T input = net ? maximum(rain-demand,0.0) : rain;
    const T potential = net ? maximum(demand-rain,0.0) : demand;
    const T interception = net ? minimum(rain,demand) : T(0.0);
    const T coefficient = dt*percolation/pow(capacity,4.0);
    T recharge = 0.0;
    if (symmetric) {
        const T g = uniform_release(water,coefficient*0.5);
        water -= g; recharge += g;
    }
    const auto wet = [&](T state) {
        const T deficit = maximum(1.0-bounded(state/capacity,0.0,1.0),1e-7);
        const T exponent = shape+1.0;
        const T ratio = input/(exponent*capacity*exp(log(deficit)/exponent))
            ;
        const T retained = !(ratio < 1.0) ? capacity*deficit
            : -capacity*deficit*expm1(exponent*log1p(-ratio));
        return minimum(retained,input);
    };
    const auto dry = [&](T state, T depth) {
        const T deficit = maximum(1.0-bounded(state/capacity,0.0,1.0),1e-7);
        const T d2 = deficit*deficit, response = tanh(depth/capacity);
        return minimum(capacity*(1.0-d2)*response/(1.0+d2/deficit*response),
            state);
    };
    T retained, evaporation;
    if (symmetric && !net) {
        evaporation = dry(water,0.5*potential); water -= evaporation;
        retained = wet(water); water += retained;
        const T e = dry(water,0.5*potential); water -= e; evaporation += e;
    } else {
        retained = wet(water); evaporation = dry(water,potential);
        water += retained-evaporation;
    }
    const T g = uniform_release(water,coefficient*(symmetric ? 0.5 : 1.0));
    water -= g; recharge += g;
    return {input-retained,recharge,evaporation+interception};
}
/**
 * uniform_groundwater: Advance uniform-input groundwater release and
 * uptake.
 *
 * INPUT ARGUMENTS:
 * water : Available water or store depth (mm).
 * recharge : Recharge supplied to the groundwater store.
 * demand : Evaporation demand on this operator time base.
 * rate : Recession or cascade release rate (1/day).
 * uptake_factor : Dimensionless groundwater evaporation-uptake factor.
 * reference : Pilot mean recharge or reference storage used by the caller.
 *
 * OUTPUT ARGUMENTS:
 * Numerical result or flux structure in the argument time base; non-const
 * state references are advanced in place.
 */
template<class T> GroundwaterFluxes<T> uniform_groundwater(
    T& water, T recharge, T demand, T rate, T uptake_factor, T reference) {
    const T half = uptake_factor*0.5;
    const auto uptake = [&](T state) {
        return -state*expm1(-half*demand/
            (state*exp(-0.5*half*demand/(state+reference))+reference));
    };
    const T e1 = uptake(water); water -= e1;
    const T g = maximum(recharge,0.0), y = g*rate;
    const T complement = y < 0.04
        ? y*(1.0/3-y*(2.0/15-y*(17.0/315-y*62.0/2835)))
        : 1.0-tanh(sqrt(y))/sqrt(y);
    const T ratio = 1.0-complement, den = 1.0+rate*water*ratio;
    const T flow = (rate*water*water*ratio+y*water*ratio+g*complement)/den;
    water = (water+g*ratio)/den;
    const T e2 = uptake(water); water -= e2;
    return {flow,e1+e2};
}
/**
 * uniform_lrel: Evaluate the nonlinear quick-release split operator.
 *
 * INPUT ARGUMENTS:
 * water : Available water or store depth (mm).
 * half : Half-step duration or split-release multiplier.
 * capacity : Store or infiltration capacity in this operator.
 * c0 : Connectivity multiplier in the quick-release law.
 * c1 : Storage-dependent quick-release exponent coefficient.
 *
 * OUTPUT ARGUMENTS:
 * Numerical result or flux structure in the argument time base; non-const
 * state references are advanced in place.
 */
template<class T> T uniform_lrel(T water,T half,T capacity,T c0,T c1) {
    const T th = tanh((c0+c1*log(maximum(water,1e-3)/capacity))/6.0);
    const T a = half*exp(6.0*th), nu = 4.0+c1*(1.0-th*th);
    return -water*expm1(-log1p(nu*a*pow(water,4.0))/nu);
}
/**
 * uniform_quick_chain: Evaluate quick-store splitting at the requested
 * refinement.
 *
 * INPUT ARGUMENTS:
 * water : Available water or store depth (mm).
 * u1 : Input depth to the quick-routing branch.
 * u2 : Input depth to the direct/bypass branch.
 * release : Release used by this numerical operator.
 * exchange : Exchange used by this numerical operator.
 * capacity : Store or infiltration capacity in this operator.
 * c0 : Connectivity multiplier in the quick-release law.
 * c1 : Storage-dependent quick-release exponent coefficient.
 * steps : Number of split substeps or integration stages.
 *
 * OUTPUT ARGUMENTS:
 * Numerical result or flux structure in the argument time base; non-const
 * state references are advanced in place.
 */
template<class T> std::array<T,4> uniform_quick_chain(
    T water,T u1,T u2,T release,T exchange,T capacity,T c0,T c1,int steps) {
    T flow=0.0, external=0.0, direct_sum=0.0;
    for(int i=0;i<steps;++i) {
        const T q1=uniform_lrel(water,release,capacity,c0,c1); water-=q1;
        const T e1=exchange*pow(maximum(water,0.0)/30.0,3.5);
        const T ex=exchange*pow(maximum(water+0.5*(u1+0.9*e1),0.0)/30.0,3.5)
            ;
        const T next=maximum(water+u1+0.9*ex,0.0), direct=maximum(u2+0.1*ex,
            0.0);
        external+=(next-water-u1)+(direct-u2); direct_sum+=direct;
        const T q2=uniform_lrel(next,release,capacity,c0,c1);
        water=next-q2; flow+=q1+q2+direct;
    }
    return {water,flow,external,direct_sum};
}
/**
 * uniform_quick: Advance uniform-input quick routing and direct discharge.
 *
 * INPUT ARGUMENTS:
 * water : Available water or store depth (mm).
 * u1 : Input depth to the quick-routing branch.
 * u2 : Input depth to the direct/bypass branch.
 * release : Release used by this numerical operator.
 * exchange : Exchange used by this numerical operator.
 * capacity : Store or infiltration capacity in this operator.
 * c0 : Connectivity multiplier in the quick-release law.
 * c1 : Storage-dependent quick-release exponent coefficient.
 * direct_output : Optional output pointer receiving direct-branch
 * discharge.
 *
 * OUTPUT ARGUMENTS:
 * Numerical result or flux structure in the argument time base; non-const
 * state references are advanced in place.
 */
template<class T> QuickFluxes<T> uniform_quick(
    T& water,T u1,T u2,T release,T exchange,T capacity,T c0,T c1,T*
        direct_output=nullptr) {
    const auto a=uniform_quick_chain(water,u1/4.0,u2/4.0,release/8.0,
        exchange/4.0,capacity,c0,c1,4);
    const auto b=uniform_quick_chain(water,u1/2.0,u2/2.0,release/4.0,
        exchange/2.0,capacity,c0,c1,2);
    const T next=(4.0*a[0]-b[0])/3.0, correction=maximum(-next,0.0);
    water=next+correction;
    if(direct_output) *direct_output=(4.0*a[3]-b[3])/3.0;
    return {(4.0*a[1]-b[1])/3.0,(4.0*a[2]-b[2])/3.0+correction};
}
/**
 * uniform_snow: Advance uniform-input snow accumulation and melt.
 *
 * INPUT ARGUMENTS:
 * snow : Snow-store depths for the five elevation bands (mm).
 * z : Numerical state or elevation vector, as used by this helper.
 * area : Five nonnegative band fractions, summing to one.
 * p : Model parameter vector; CaRT calibration or CaRT+ physical units.
 * rain : Liquid input depth/rate on this operator time base.
 * temperature : Temperature used by this numerical operator.
 * pet : Potential evaporation depth/rate on this operator time base.
 * dt : Integration or observation interval duration (days).
 *
 * OUTPUT ARGUMENTS:
 * Numerical result or flux structure in the argument time base; non-const
 * state references are advanced in place.
 */
template<class T> SnowFluxes<T> uniform_snow(
    std::array<T,band_count>& snow,const std::array<double,band_count>& z,
    const std::array<double,band_count>& area,const SnowParameters<T>& p,
    double rain,double temperature,double pet,double dt=1.0) {
    double mean=0.0; for(std::size_t k=0;k<band_count;++k) mean+=area[k]*z[k
        ];
    SnowFluxes<T> out;
    for(std::size_t k=0;k<band_count;++k) {
        const T warm=T(temperature)-p.lapse_rate_c_per_km*((z[k]-mean)/
            1000.0)-p.threshold_c;
        const bool cold=warm<0.0;
        if(cold) {snow[k]+=p.snowfall_correction*rain;
            out.snowfall_adjustment+=area[k]*rain*(p.snowfall_correction-1.0
                );}
        const T sub=p.sublimation_fraction*pet;
        const T melt=cold ? T(0.0) : p.degree_day_factor*dt*warm+p.
            pet_melt_factor*pet;
        const T demand=sub+melt;
        const T fraction=demand>snow[k] ? snow[k]/demand : T(1.0);
        snow[k]=maximum(snow[k]-fraction*demand,0.0);
        out.liquid+=area[k]*(fraction*melt+(cold ? 0.0 : rain));
        out.sublimation+=area[k]*fraction*sub;
    }
    return out;
}
/**
 * uniform_uh_integral: Evaluate the uniform-input hydrograph integral.
 *
 * INPUT ARGUMENTS:
 * time : Current model time (days).
 * base : Unit-hydrograph time base (days).
 * symmetric : Whether to use the symmetric two-time-base hydrograph.
 *
 * OUTPUT ARGUMENTS:
 * Numerical result or flux structure in the argument time base; non-const
 * state references are advanced in place.
 */
template<class T> T uniform_uh_integral(double time,T base,bool symmetric) {
    if(time<=0) return T(0.0);
    const T u=time/base;
    if(!symmetric) return u>1.0 ? T(time)-base*(1.0-1.0/3.5) : base*pow(u,
        3.5)/3.5;
    if(u<=1.0) return base*pow(u,3.5)/7.0;
    if(u<2.0) return T(time)-base+base*pow(2.0-u,3.5)/7.0;
    return T(time)-base;
}
template<std::size_t Length,class T> class UniformHydrograph {
    std::vector<T> weights, pending;
    std::vector<std::size_t> seeds;
    std::size_t head=0;
    bool circular;
public:
    /**
     * UniformHydrograph: Initialize the numerical scalar, routing object,
     * or solver workspace.
     *
     * INPUT ARGUMENTS:
     * base : Unit-hydrograph time base (days).
     * dt : Integration or observation interval duration (days).
     * skip_zero : Whether a zero input can use the sparse routing update.
     *
     * OUTPUT ARGUMENTS:
     * Numerical result or flux structure in the argument time base;
     * non-const state references are advanced in place.
     */
    UniformHydrograph(T base,double dt,bool skip_zero=false):circular(
        skip_zero) {
        const auto n=static_cast<std::size_t>(std::ceil(primal(base)*(Length
            ==6?1:2)/dt))+1;
        weights.resize(n); pending.resize(n);
        T previous=0.0;
        for(std::size_t i=0;i<n;++i) {
            const double t=i*dt;
            const T cumulative=(uniform_uh_integral(t+dt,base,Length!=6)
                -uniform_uh_integral(t,base,Length!=6))/dt;
            weights[i]=cumulative-previous; previous=cumulative;
        }
        if constexpr(std::is_same_v<T,Tangent>) {
            for(std::size_t j=0;j<parameter_count;++j)
                if(base.slope[j]!=0.0) seeds.push_back(j);
        }
    }
    /**
     * advance: Advance routing memory and return its outlet depth.
     *
     * INPUT ARGUMENTS:
     * input : Forcing/storage input for the current operator.
     *
     * OUTPUT ARGUMENTS:
     * Numerical result or flux structure in the argument time base;
     * non-const state references are advanced in place.
     */
    T advance(T input) {
        if constexpr(std::is_same_v<T,Tangent>) {
            if(circular) {
                const T result=pending[head]+weights[0]*input;
                pending[head]=T{};
                if(++head==pending.size()) head=0;
                if(input.value==0.0 && std::all_of(input.slope.begin(),input
                    .slope.end(),
                        [](double x){return x==0.0;})) return result;
                std::size_t slot=head;
                for(std::size_t i=1;i<weights.size();++i) {
                    auto& next=pending[slot]; const auto& weight=weights[i];
                    next.value+=weight.value*input.value;
                    for(std::size_t j=0;j<parameter_count;++j)
                        next.slope[j]+=weight.value*input.slope[j];
                    for(auto j:seeds) next.slope[j]+=weight.slope[j]*input.
                        value;
                    if(++slot==pending.size()) slot=0;
                }
                return result;
            }
        }
        const T result=pending[0]+weights[0]*input;
        for(std::size_t i=1;i<weights.size();++i) pending[i-1]=pending[i]+
            weights[i]*input;
        pending.back()=0.0; return result;
    }
    /**
     * storage_change: Compute retained routing water from input minus
     * output.
     *
     * INPUT ARGUMENTS:
     * input : Forcing/storage input for the current operator.
     * output : Caller-owned numerical output buffer.
     *
     * OUTPUT ARGUMENTS:
     * Numerical result or flux structure in the argument time base;
     * non-const state references are advanced in place.
     */
    T storage_change(T input,T output) {return input-output;}
};
} // namespace cartplus
