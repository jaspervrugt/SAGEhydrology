// ========================================================================
// SAGEhydrology | Model source
// File: Native/models/CaRTPlus/cartplus_runge_kutta.hpp
//
// CaRT / CaRT+ numerical model and supporting routines.
//
// SAGE adaptations: Jasper A. Vrugt, University of California, Irvine.
// Upstream model: Xia (2026). Preserve the adjacent license notices.
// ========================================================================

#pragma once

#include "cartplus_tangent.hpp"
#include "cartplus_dop853.hpp"
#include <algorithm>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <vector>

namespace cartplus::continuous {

// Configured entry points set these per call; historical entry points keep
// their original unrestricted adaptive stepping. Units are days.
inline thread_local double maximum_step = std::numeric_limits<double>::
    infinity();
inline thread_local double minimum_step = 0.;
struct StepLimits {
    double old_max = maximum_step, old_min = minimum_step;
    /**
     * StepLimits: Initialize the numerical scalar, routing object, or
     * solver workspace.
     *
     * INPUT ARGUMENTS:
     * maximum : Upper step bound (days).
     * minimum : Lower step bound (days).
     *
     * OUTPUT ARGUMENTS:
     * Numerical result or flux structure in the argument time base;
     * non-const state references are advanced in place.
     */
    StepLimits(double maximum, double minimum) {
        if (!(maximum > 0.) || !std::isfinite(maximum) || minimum < 0.
            || !std::isfinite(minimum) || minimum > maximum)
            throw std::runtime_error("Invalid integration step limits");
        maximum_step = maximum; minimum_step = minimum;
    }
    /**
     * StepLimits: Initialize the numerical scalar, routing object, or
     * solver workspace.
     *
     * INPUT ARGUMENTS:
     * None.
     *
     * OUTPUT ARGUMENTS:
     * Numerical result or flux structure in the argument time base;
     * non-const state references are advanced in place.
     */
    ~StepLimits() { maximum_step = old_max; minimum_step = old_min; }
};

/**
 * bounded_step: Limit a requested step to configured bounds and time left.
 *
 * INPUT ARGUMENTS:
 * requested : Requested numerical step duration (days).
 * remaining : Time remaining in the current interval (days).
 *
 * OUTPUT ARGUMENTS:
 * Numerical result or flux structure in the argument time base; non-const
 * state references are advanced in place.
 */
inline double bounded_step(double requested, double remaining) {
    const double h = std::min({requested, maximum_step, remaining});
    if (h < minimum_step && remaining > minimum_step)
        throw std::runtime_error("Integration requires a step below minimum"
            );
    return h;
}

/**
 * weighted_add: Accumulate a weighted stage value and its derivatives.
 *
 * INPUT ARGUMENTS:
 * target : Target used by this numerical operator.
 * weight : Stage, band, or quadrature weight (-).
 * source : Source stage value or input record.
 *
 * OUTPUT ARGUMENTS:
 * No return value; updates referenced states and output buffers.
 */
inline void weighted_add(double& target, double weight, double source) {
    target += weight*source;
}

/**
 * weighted_add: Accumulate a weighted stage value and its derivatives.
 *
 * INPUT ARGUMENTS:
 * target : Target used by this numerical operator.
 * weight : Stage, band, or quadrature weight (-).
 * source : Source stage value or input record.
 *
 * OUTPUT ARGUMENTS:
 * No return value; updates referenced states and output buffers.
 */
inline void weighted_add(Tangent& target, double weight, const Tangent&
    source) {
    target.value += weight*source.value;
    for (std::size_t i = 0; i < parameter_count; ++i)
        target.slope[i] += weight*source.slope[i];
}

/**
 * integrate: Advance the augmented state with the configured integrator.
 *
 * INPUT ARGUMENTS:
 * y : Numerical state/operand; augmented with sensitivities when enabled.
 * rhs : Callback evaluating state derivatives on the model time base.
 * dt : Integration or observation interval duration (days).
 * rtol : Positive relative integration-error tolerance (-).
 * atol : Positive absolute integration-error tolerance for augmented
 * states.
 *
 * OUTPUT ARGUMENTS:
 * No return value; updates referenced states and output buffers.
 */
template<class S, class F>
void integrate(std::vector<S>& y, const F& rhs, double dt,
               double rtol, double atol) {
    const auto n = y.size();
    std::array<std::vector<S>, 13> k;
    for (auto& stage : k) stage.resize(n);
    std::vector<S> next(n);
    rhs(y, k[0]);
    double t = 0, h = std::min(dt, std::max(.001, minimum_step));
    bool rejected = false;
    for (int trial = 0; t < dt; ++trial) {
        if (trial == 100000 || h < 10*(std::nextafter(t, dt+1)-t))
            throw std::runtime_error("Continuous CaRT+ step limit");
        const double end = t + bounded_step(h, dt-t);
        h = end-t;
        for (int s = 1; s <= 12; ++s) {
            const double* weights = s == 12 ? dop853::b : dop853::a[s];
            for (std::size_t j = 0; j < n; ++j) {
                S sum = 0;
                for (int l = 0; l < s; ++l)
                    if (weights[l] != 0.) weighted_add(sum,weights[l],k[l][j
                        ]);
                next[j] = y[j];
                weighted_add(next[j],h,sum);
            }
            rhs(next, k[s]);
        }
        double norm5 = 0, norm3 = 0;
        const auto accumulate = [&](double d5, double d3,
                                    double old, double candidate) {
            const double scale = atol+rtol*std::max(std::abs(old),
                                                    std::abs(candidate));
            norm5 += (d5/scale)*(d5/scale);
            norm3 += (d3/scale)*(d3/scale);
        };
        for (std::size_t j = 0; j < n; ++j) {
            S err5 = 0, err3 = 0;
            for (int s = 0; s <= 12; ++s) {
                if (dop853::e5[s] != 0.) weighted_add(err5,dop853::e5[s],k[s
                    ][j]);
                if (dop853::e3[s] != 0.) weighted_add(err3,dop853::e3[s],k[s
                    ][j]);
            }
            accumulate(primal(err5), primal(err3), primal(y[j]), primal(next
                [j]));
            if constexpr (std::is_same_v<S, Tangent>)
                for (std::size_t c = 0; c < parameter_count; ++c)
                    accumulate(err5.slope[c], err3.slope[c], y[j].slope[c],
                               next[j].slope[c]);
        }
        constexpr int width = std::is_same_v<S, Tangent> ? parameter_count+1
            : 1;
        const double error = norm5 == 0 && norm3 == 0 ? 0 :
            h*norm5/std::sqrt((norm5+.01*norm3)*n*width);
        const double factor = std::isfinite(error)
            ? std::clamp(.9*std::pow(error, -1./8), .2, 10.) : .2;
        if (error < 1) {
            y.swap(next);
            k[0].swap(k[12]);
            t = end;
            h *= rejected ? std::min(factor, 1.) : factor;
            rejected = false;
        } else {
            h *= factor;
            rejected = true;
        }
    }
}

// Embedded Euler/Heun pair used by the other SAGE C++ RK2 solvers.
// The Euler difference controls the error in states and sensitivities.
template<class S>
struct Rk2Workspace {
    /**
     * Rk2Workspace: Initialize the numerical scalar, routing object, or
     * solver workspace.
     *
     * INPUT ARGUMENTS:
     * n : Number of intervals or state entries in this helper.
     *
     * OUTPUT ARGUMENTS:
     * Numerical result or flux structure in the argument time base;
     * non-const state references are advanced in place.
     */
    explicit Rk2Workspace(std::size_t n)
        : first(n), second(n), euler(n), heun(n) {}
    std::vector<S> first, second, euler, heun;
};

struct Rk2NoExactUpdate {};
struct Rk2DefaultTrial {};

template<class S, class F, class ExactUpdate = Rk2NoExactUpdate,
         class Trial = Rk2DefaultTrial>
/**
 * integrate_rk2: Advance adaptive Heun RK2 with an embedded Euler estimate.
 *
 * INPUT ARGUMENTS:
 * y : Numerical state/operand; augmented with sensitivities when enabled.
 * rhs : Callback evaluating state derivatives on the model time base.
 * dt : Integration or observation interval duration (days).
 * rtol : Positive relative integration-error tolerance (-).
 * atol : Positive absolute integration-error tolerance for augmented
 * states.
 * next_step : Suggested next integration step (days), updated by the
 * solver.
 * work : Reusable Runge-Kutta stage/state workspace.
 * exact_index : Index of the component advanced by an analytic update.
 * exact_update : Optional analytic substep-update callback.
 * custom_trial : Optional trial-step callback.
 * accepted_steps : Optional output of accepted step durations.
 * rejected_steps : Counter of rejected integration trials.
 * accepted_count : Number of accepted integration steps.
 *
 * OUTPUT ARGUMENTS:
 * No return value; updates referenced states and output buffers.
 */
void integrate_rk2(std::vector<S>& y, const F& rhs, double dt,
                   double rtol, double atol, double& next_step,
                   Rk2Workspace<S>& work, std::size_t exact_index = 0,
                   const ExactUpdate& exact_update = {}, const Trial&
                       custom_trial = {},
                   std::vector<double>* accepted_steps = nullptr,
                   double* rejected_steps = nullptr, double* accepted_count
                       = nullptr) {
    const auto n = y.size();
    auto& [first, second, euler, heun] = work;
    if constexpr (std::is_same_v<Trial,Rk2DefaultTrial>) rhs(y,first);
    double t = 0;
    next_step = std::max(next_step, minimum_step);
    for (int trial = 0; t < dt; ++trial) {
        const double end = t + bounded_step(next_step,dt-t);
        const double h = end-t;
        if (trial == 100000 || h <= 0)
            throw std::runtime_error("Continuous CaRT+ RK2 step limit");
        if constexpr (std::is_same_v<Trial,Rk2DefaultTrial>) {
            for (std::size_t j = 0; j < n; ++j) {
                euler[j] = y[j];
                weighted_add(euler[j],h,first[j]);
            }
            if constexpr (!std::is_same_v<ExactUpdate,Rk2NoExactUpdate>)
                euler[exact_index] = exact_update(y[exact_index],h);
            rhs(euler,second);
        } else {
            custom_trial(y,work,h);
        }
        double squared_error = 0;
        const auto accumulate = [&](double candidate, double estimate) {
            const double scale = atol+rtol*std::abs(candidate);
            const double error = (candidate-estimate)/scale;
            squared_error += error*error;
        };
        for (std::size_t j = 0; j < n; ++j) {
            if constexpr (std::is_same_v<Trial,Rk2DefaultTrial>) {
                heun[j] = y[j];
                weighted_add(heun[j],.5*h,first[j]);
                weighted_add(heun[j],.5*h,second[j]);
                if constexpr (!std::is_same_v<ExactUpdate,Rk2NoExactUpdate>)
                    {
                    if (j == exact_index) heun[j] = euler[j];
                }
            }
            accumulate(primal(heun[j]),primal(euler[j]));
            if constexpr (std::is_same_v<S,Tangent>)
                for (std::size_t k = 0; k < parameter_count; ++k)
                    accumulate(heun[j].slope[k],euler[j].slope[k]);
        }
        constexpr int width = std::is_same_v<S,Tangent> ? parameter_count+1
            : 1;
        const double error = std::sqrt(squared_error/(n*width));
        const double factor = std::isfinite(error)
            ? std::clamp(.9/std::sqrt(error),.2,5.) : .2;
        next_step = h*factor;
        if (error <= 1.) {
            if (accepted_steps) accepted_steps->push_back(h);
            if (accepted_count) ++*accepted_count;
            y.swap(heun);
            t = end;
            if (dt-t <= 10*std::numeric_limits<double>::epsilon()*dt) return
                ;
            if constexpr (std::is_same_v<Trial,Rk2DefaultTrial>) rhs(y,first
                );
        } else if (rejected_steps) {
            ++*rejected_steps;
        }
    }
}

/**
 * integrate_rk2: Advance adaptive Heun RK2 with an embedded Euler estimate.
 *
 * INPUT ARGUMENTS:
 * y : Numerical state/operand; augmented with sensitivities when enabled.
 * rhs : Callback evaluating state derivatives on the model time base.
 * dt : Integration or observation interval duration (days).
 * rtol : Positive relative integration-error tolerance (-).
 * atol : Positive absolute integration-error tolerance for augmented
 * states.
 * next_step : Suggested next integration step (days), updated by the
 * solver.
 *
 * OUTPUT ARGUMENTS:
 * No return value; updates referenced states and output buffers.
 */
template<class S, class F>
void integrate_rk2(std::vector<S>& y, const F& rhs, double dt,
                   double rtol, double atol, double& next_step) {
    Rk2Workspace<S> work(y.size());
    integrate_rk2(y,rhs,dt,rtol,atol,next_step,work);
}

// Differentiate the accepted Heun updates on a fixed time mesh.
/**
 * replay_rk2: Replay accepted RK2 steps for forward-sensitivity evaluation.
 *
 * INPUT ARGUMENTS:
 * y : Numerical state/operand; augmented with sensitivities when enabled.
 * rhs : Callback evaluating state derivatives on the model time base.
 * work : Reusable Runge-Kutta stage/state workspace.
 * begin : Begin used by this numerical operator.
 * end : End used by this numerical operator.
 *
 * OUTPUT ARGUMENTS:
 * No return value; updates referenced states and output buffers.
 */
template<class S, class F>
void replay_rk2(std::vector<S>& y, const F& rhs, Rk2Workspace<S>& work,
                const double* begin, const double* end) {
    auto& [first, second, euler, heun] = work;
    for (auto step = begin; step != end; ++step) {
        const double h = *step;
        rhs(y,first);
        for (std::size_t j = 0; j < y.size(); ++j) {
            euler[j] = y[j];
            weighted_add(euler[j],h,first[j]);
        }
        rhs(euler,second);
        for (std::size_t j = 0; j < y.size(); ++j) {
            weighted_add(y[j],.5*h,first[j]);
            weighted_add(y[j],.5*h,second[j]);
        }
    }
}

// Bogacki-Shampine 3(2). The last stage is reused after acceptance.
/**
 * integrate_rk23: Advance the embedded RK2/3 method with admissibility
 * checks.
 *
 * INPUT ARGUMENTS:
 * y : Numerical state/operand; augmented with sensitivities when enabled.
 * rhs : Callback evaluating state derivatives on the model time base.
 * dt : Integration or observation interval duration (days).
 * rtol : Positive relative integration-error tolerance (-).
 * atol : Positive absolute integration-error tolerance for augmented
 * states.
 * next_step : Suggested next integration step (days), updated by the
 * solver.
 * stores : Store count, state vector, or caller-owned state buffer.
 *
 * OUTPUT ARGUMENTS:
 * No return value; updates referenced states and output buffers.
 */
template<class S, class F>
void integrate_rk23(std::vector<S>& y, const F& rhs, double dt,
                    double rtol, double atol, double& next_step, int stores)
                        {
    const auto n = y.size();
    std::vector<S> k1(n), k2(n), k3(n), k4(n), stage(n), candidate(n);
    rhs(y,k1);
    double t = 0;
    next_step = std::max(next_step, minimum_step);
    for (int trial = 0; t < dt; ++trial) {
        const double end = t + bounded_step(next_step,dt-t), h = end-t;
        if (trial == 100000 || h <= 0.)
            throw std::runtime_error("Coupled RK23 step limit");
        for (std::size_t j = 0; j < n; ++j) {
            stage[j] = y[j];
            weighted_add(stage[j],h*.5,k1[j]);
        }
        rhs(stage,k2);
        for (std::size_t j = 0; j < n; ++j) {
            stage[j] = y[j];
            weighted_add(stage[j],h*.75,k2[j]);
            candidate[j] = y[j];
            weighted_add(candidate[j],h*(2./9),k1[j]);
            weighted_add(candidate[j],h/3,k2[j]);
        }
        rhs(stage,k3);
        for (std::size_t j = 0; j < n; ++j)
            weighted_add(candidate[j],h*(4./9),k3[j]);
        // Reject a trial that crosses the admissible storage boundary.
        const double floor = -std::max(1e-7,100*atol);
        bool admissible = true;
        for (int j = 0; j < stores; ++j)
            admissible &= primal(candidate[j]) >= floor;
        if (!admissible) {
            next_step = h*.5;
            continue;
        }
        rhs(candidate,k4);
        double squared_error = 0;
        const auto accumulate = [&](double error, double old, double
            proposed) {
            const double scale = atol+rtol*std::max(std::abs(old),std::abs(
                proposed));
            squared_error += (error/scale)*(error/scale);
        };
        for (std::size_t j = 0; j < n; ++j) {
            S error = 0;
            weighted_add(error,h*(5./72),k1[j]);
            weighted_add(error,-h/12,k2[j]);
            weighted_add(error,-h/9,k3[j]);
            weighted_add(error,h/8,k4[j]);
            accumulate(primal(error),primal(y[j]),primal(candidate[j]));
            if constexpr (std::is_same_v<S,Tangent>)
                for (std::size_t k = 0; k < parameter_count; ++k)
                    accumulate(error.slope[k],y[j].slope[k],candidate[j].
                        slope[k]);
        }
        constexpr int width = std::is_same_v<S,Tangent> ? parameter_count+1
            : 1;
        const double error = std::sqrt(squared_error/(n*width));
        const double factor = std::isfinite(error)
            ? std::clamp(.9*std::pow(error,-1./3),.2,5.) : .2;
        next_step = h*factor;
        if (error <= 1.) {
            y.swap(candidate);
            k1.swap(k4);
            t = end;
            if (dt-t <= 10*std::numeric_limits<double>::epsilon()*dt) return
                ;
        }
    }
}

// Fixed methods deliberately do not adapt or clip stores. An unstable
// user-selected step is reported by the model's state/finite checks.
/**
 * integrate_fixed_duration: Advance a fixed step whose duration may carry
 * derivatives.
 *
 * INPUT ARGUMENTS:
 * y : Numerical state/operand; augmented with sensitivities when enabled.
 * rhs : Callback evaluating state derivatives on the model time base.
 * duration : Duration used by this numerical operator.
 * heun : Whether to use Heun rather than explicit Euler.
 *
 * OUTPUT ARGUMENTS:
 * No return value; updates referenced states and output buffers.
 */
template<class S, class F>
void integrate_fixed_duration(std::vector<S>& y, const F& rhs,
                              const S& duration, bool heun) {
    const double steps = std::ceil(primal(duration)/maximum_step);
    if (steps < 1 || steps > 1000000)
        throw std::runtime_error("Invalid fixed-step count");
    const int count = static_cast<int>(steps);
    const S h = duration/count;
    std::vector<S> k1(y.size()), k2(y.size()), stage(y.size());
    for (int step = 0; step < count; ++step) {
        rhs(y,k1);
        stage = y;
        for (std::size_t j = 0; j < y.size(); ++j)
            stage[j] = stage[j]+h*k1[j];
        if (heun) {
            rhs(stage,k2);
            for (std::size_t j = 0; j < y.size(); ++j) {
                y[j] = y[j]+(h*.5)*(k1[j]+k2[j]);
            }
        } else y.swap(stage);
    }
}

/**
 * integrate_fixed: Advance fixed-step Euler or Heun over one interval.
 *
 * INPUT ARGUMENTS:
 * y : Numerical state/operand; augmented with sensitivities when enabled.
 * rhs : Callback evaluating state derivatives on the model time base.
 * dt : Integration or observation interval duration (days).
 * heun : Whether to use Heun rather than explicit Euler.
 *
 * OUTPUT ARGUMENTS:
 * No return value; updates referenced states and output buffers.
 */
template<class S, class F>
void integrate_fixed(std::vector<S>& y, const F& rhs, double dt, bool heun)
    {
    integrate_fixed_duration(y,rhs,S(dt),heun);
}

// Dormand-Prince 5(4), including the forward sensitivity lanes in the
// error norm, as in the existing adaptive methods.
/**
 * integrate_rk45: Advance the embedded RK4/5 method with admissibility
 * checks.
 *
 * INPUT ARGUMENTS:
 * y : Numerical state/operand; augmented with sensitivities when enabled.
 * rhs : Callback evaluating state derivatives on the model time base.
 * dt : Integration or observation interval duration (days).
 * rtol : Positive relative integration-error tolerance (-).
 * atol : Positive absolute integration-error tolerance for augmented
 * states.
 * next_step : Suggested next integration step (days), updated by the
 * solver.
 *
 * OUTPUT ARGUMENTS:
 * No return value; updates referenced states and output buffers.
 */
template<class S, class F>
void integrate_rk45(std::vector<S>& y, const F& rhs, double dt,
                    double rtol, double atol, double& next_step) {
    constexpr double a[7][7] = {{}, {1./5}, {3./40,9./40},
        {44./45,-56./15,32./9},
        {19372./6561,-25360./2187,64448./6561,-212./729},
        {9017./3168,-355./33,46732./5247,49./176,-5103./18656},
        {35./384,0,500./1113,125./192,-2187./6784,11./84}};
    constexpr double b4[7] = {5179./57600,0,7571./16695,393./640,
                              -92097./339200,187./2100,1./40};
    const auto n = y.size();
    std::array<std::vector<S>,7> k;
    for (auto& v : k) v.resize(n);
    std::vector<S> candidate(n), stage(n);
    double t = 0;
    next_step = std::max(next_step, minimum_step);
    rhs(y,k[0]);
    for (int trial = 0; t < dt; ++trial) {
        const double h = bounded_step(next_step,dt-t);
        if (trial >= 100000 || t+h <= t)
            throw std::runtime_error("RK45 step limit");
        for (int s = 1; s < 7; ++s) {
            stage = y;
            for (std::size_t j = 0; j < n; ++j)
                for (int l = 0; l < s; ++l)
                    if (a[s][l]) weighted_add(stage[j],h*a[s][l],k[l][j]);
            rhs(stage,k[s]);
        }
        candidate = stage;
        double sum = 0;
        const auto accumulate = [&](double error, double old, double next) {
            const double scale = atol+rtol*std::max(std::abs(old),std::abs(
                next));
            sum += (error/scale)*(error/scale);
        };
        for (std::size_t j = 0; j < n; ++j) {
            S error = 0.;
            for (int s = 0; s < 7; ++s)
                weighted_add(error,h*(a[6][s]-b4[s]),k[s][j]);
            accumulate(primal(error),primal(y[j]),primal(candidate[j]));
            if constexpr (std::is_same_v<S,Tangent>)
                for (std::size_t c = 0; c < parameter_count; ++c)
                    accumulate(error.slope[c],y[j].slope[c],candidate[j].
                        slope[c]);
        }
        constexpr int width = std::is_same_v<S,Tangent> ? parameter_count+1
            : 1;
        const double error = std::sqrt(sum/(n*width));
        next_step = h*(std::isfinite(error)
            ? std::clamp(.9*std::pow(error,-.2),.2,5.) : .2);
        if (error <= 1.) {
            y.swap(candidate); k[0].swap(k[6]); t += h;
        }
    }
}

/**
 * integrate_segment: Dispatch the solver for a finite forcing/event
 * segment.
 *
 * INPUT ARGUMENTS:
 * y : Numerical state/operand; augmented with sensitivities when enabled.
 * rhs : Callback evaluating state derivatives on the model time base.
 * dt : Integration or observation interval duration (days).
 * rtol : Positive relative integration-error tolerance (-).
 * atol : Positive absolute integration-error tolerance for augmented
 * states.
 * next_step : Suggested next integration step (days), updated by the
 * solver.
 * solver : Integration method or solver-control object.
 * stores : Store count, state vector, or caller-owned state buffer.
 *
 * OUTPUT ARGUMENTS:
 * No return value; updates referenced states and output buffers.
 */
template<class S, class F>
void integrate_segment(std::vector<S>& y, const F& rhs, double dt,
                       double rtol, double atol, double& next_step, int
                           solver, int stores = 0) {
    if (solver == 1) integrate_rk2(y,rhs,dt,rtol,atol,next_step);
    else if (solver == 2) integrate_rk23(y,rhs,dt,rtol,atol,next_step,stores
        );
    else if (solver == 9) integrate_rk45(y,rhs,dt,rtol,atol,next_step);
    else if (solver == 10 || solver == 11) integrate_fixed(y,rhs,dt,solver
        == 10);
    else integrate(y,rhs,dt,rtol,atol);
}

/**
 * exits: Construct the cascade exit fractions for a unit hydrograph.
 *
 * INPUT ARGUMENTS:
 * m : Positive routing refinement: stages per X4 time base.
 * symmetric : Whether to use the symmetric two-time-base hydrograph.
 *
 * OUTPUT ARGUMENTS:
 * Numerical result or flux structure in the argument time base; non-const
 * state references are advanced in place.
 */
inline std::vector<double> exits(int m, bool symmetric) {
    std::vector<double> eta(symmetric ? 2*m : m);
    double previous = 0;
    for (std::size_t j = 0; j < eta.size(); ++j) {
        const double x = double(j+1)/m;
        const double h = symmetric
            ? (x <= 1 ? .5*std::pow(x,2.5) : 1-.5*std::pow(2-x,2.5))
            : std::pow(x,2.5);
        eta[j] = h-previous;
        previous = h;
    }
    double survival = 0;
    for (std::size_t j = eta.size(); j-- > 0;) {
        survival += eta[j];
        eta[j] /= survival;
    }
    return eta;
}

} // namespace cartplus::continuous
