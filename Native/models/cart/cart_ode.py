# ========================================================================
# SAGEhydrology | Model source
# File: Native/models/cart/cart_ode.py
#
# CaRT / CaRT+ numerical model and supporting routines.
#
# SAGE adaptations: Jasper A. Vrugt, University of California, Irvine.
# Upstream model: Xia (2026). Preserve the adjacent license notices.
# ========================================================================

"""
Coupled continuous-time CaRT with conservative Erlang routing cascades.

This is a continuous counterpart of Xia's interval model, not a replacement
for its finite update maps. Units are mm and days; parameter coordinates are
the same eleven coordinates as cart_python. No GUI solver selection is
changed.
"""

import numpy as np
from scipy.integrate import solve_ivp

from cart_python import Tangent, primal


def routing_fractions(refinement):
    """
    Construct routing-cascade exit fractions.

    SYNOPSIS:
        routing_fractions(...)

    INPUT ARGUMENTS:
    refinement : Positive number of routing stages per X4 time base.

    OUTPUT ARGUMENTS:
    Dimensionless stage-exit fractions for both hydrographs.
    """
    result = []
    for symmetric, count in ((False, refinement), (True, 2*refinement)):
        t = np.arange(count+1, dtype=float)/refinement
        h = (np.where(t <= 1, .5*t**2.5,
                      1-.5*np.maximum(2-t, 0)**2.5) if symmetric
             else np.minimum(t, 1)**2.5)
        weights = np.diff(h)
        survive = np.cumsum(weights[::-1])[::-1]
        eta = weights/survive
        eta[-1] = 1.
        result.append(eta)
    return result


def _inputs(forcing, elevation, area, parameters, initial, resolution,
            interval_days, refinement, rtol, atol):
    """
    Validate model dimensions, units, and integration settings.

    SYNOPSIS:
        _inputs(...)

    INPUT ARGUMENTS:
    forcing : N-by-4 forcing [P, T, PET, training]; depths mm/interval, T
        deg C.
    elevation : Five elevation-band heights (m).
    area : Five nonnegative band fractions, summing to one.
    parameters : 11 CaRT calibration coordinates, including logarithmic
        entries.
    initial : Initial stores (mm), or None/empty to use the spin-up policy.
    resolution : Forcing intervals per day: 1 daily, 24 hourly, 96
        quarter-hourly.
    interval_days : Positive interval duration (days), scalar or length N.
    refinement : Positive number of routing stages per X4 time base.
    rtol : Positive relative integration-error tolerance (-).
    atol : Positive absolute integration-error tolerance for augmented
        states.

    OUTPUT ARGUMENTS:
    Validated forcing, bands, parameters, durations, and initial stores;
    CaRT+ also returns parameter derivative scales.
    """
    f, z, a, p = map(lambda v: np.asarray(v, dtype=float),
                     (forcing, elevation, area, parameters))
    if (f.ndim != 2 or f.shape[1] != 4 or not len(f)
            or z.shape != (5,) or a.shape != (5,) or p.shape != (11,)
            or not all(np.isfinite(v).all() for v in (f, z, a, p))):
        raise ValueError('Expected finite (n,4) forcing,'
            ' five bands, 11 parameters.')
    if (np.any(f[:, [0, 2]] < 0) or np.any(a < 0)
            or not np.isclose(a.sum(), 1, rtol=0, atol=1e-10)
            or not np.isin(f[:, 3], [0, 1]).all() or not f[:, 3].any()):
        raise ValueError('Depths/areas must be nonnegati'
            've; training mask is 0/1.')
    if (resolution not in (1, 24, 96) or isinstance(refinement, bool)
            or not isinstance(refinement, (int, np.integer)) or refinement
                < 1
            or p[1] <= 0 or p[4] > 0 or p[8] < 0
            or not np.isfinite([rtol, atol]).all() or rtol <= 0 or atol <=
                0):
        raise ValueError('Require positive timebase/tole'
            'rances, X2<=0, sf>=0 and integ' 'er refinement.')
    m = int(refinement)
    try:
        durations = np.broadcast_to(1/resolution if interval_days is None
                                    else np.asarray(interval_days, dtype=
                                        float),
                                    (len(f),)).copy()
    except ValueError as exc:
        raise ValueError('interval_days must be scalar or length n.') from \
            exc
    if not np.isfinite(durations).all() or np.any(durations <= 0):
        raise ValueError('Observation interval durations'
            ' must be finite and positive.')
    ns = 8+3*m
    if initial is not None:
        initial = np.asarray(initial, dtype=float)
        if (initial.shape != (ns,) or not np.isfinite(initial).all()
                or np.any(initial < 0) or initial[5] > np.exp(p[0])):
            raise ValueError('Expected ' f'{ns}'
                ' nonnegative initial stores with S<=' 'C.')
    return f, z, a, p, durations, initial


def simulate(forcing, elevation, area, parameters, *, jacobian=False,
             initial=None, resolution=1, refinement=16, rtol=1e-8,
             atol=1e-10, interval_days=None, return_diagnostics=False):
    """
    Simulate model fluxes and optional discharge sensitivities.

    SYNOPSIS:
        simulate(...)

    INPUT ARGUMENTS:
    forcing : N-by-4 forcing [P, T, PET, training]; depths mm/interval, T
        deg C.
    elevation : Five elevation-band heights (m).
    area : Five nonnegative band fractions, summing to one.
    parameters : 11 CaRT calibration coordinates, including logarithmic
        entries.
    jacobian : Whether to evaluate discharge sensitivities.
    initial : Initial stores (mm), or None/empty to use the spin-up policy.
    resolution : Forcing intervals per day: 1 daily, 24 hourly, 96
        quarter-hourly.
    refinement : Positive number of routing stages per X4 time base.
    rtol : Positive relative integration-error tolerance (-).
    atol : Positive absolute integration-error tolerance for augmented
        states.
    interval_days : Positive interval duration (days), scalar or length N.
    return_diagnostics : Whether to append solver/water-balance diagnostics.

    OUTPUT ARGUMENTS:
    q : N interval discharge depths (mm); J : N-by-11 Jacobian or None.
    Optional stores/diagnostics are appended only when requested.

    NOTES:
    Forcing columns: precipitation depth, temperature, PET depth, training
    indicator. Physical initial stores are [N1..N5,S,R,W,Z9...,Z1...].
    Optional interval_days is a positive scalar or length-n vector, allowing
    arbitrary and irregular observation intervals instead of 1/resolution. A
    supplied initial vector is detached (zero initial sensitivities).
    Default initialization is differentiable: zero snow/cascades, S=.7*C,
    R=.5*X3 and W=gbar*tau. There is deliberately no hidden spin-up replay.
    Recharge gbar and its gradient are computed in an independent upstream
    pilot with the same initial snow/soil, forcing, and snow-event
    treatment. Diagnostics include interval water balance residuals and
    event information. Derivatives are branchwise; threshold coincidences
    are flagged explicitly.
    """
    f, z, a, p, durations, initial = _inputs(
        forcing, elevation, area, parameters, initial, resolution,
        interval_days, refinement, rtol, atol)
    m = int(refinement)
    ns = 8+3*m
    use_ad = bool(jacobian)
    zero = np.zeros(11)
    pars = [Tangent(v, seed) for v, seed in zip(p, np.eye(11))] if use_ad \
        else p
    c, x3, tau = (np.exp(pars[k]) for k in (0, 2, 3))
    kp, cf, gamma, ks = ((4/9)**4*np.exp(pars[5]), np.exp(pars[7]),
                          np.exp(pars[9]), np.exp(pars[10]))
    lam = m/pars[1]
    fractions = routing_fractions(m)
    counters = dict(rhs_evaluations=0, snow_events=0,
                    threshold_coincidences=[], boundary_event_coincidences
                        =[])

    def unpack(value, size):
        """
        Separate primal states from forward-sensitivity lanes.

        SYNOPSIS:
            unpack(...)

        INPUT ARGUMENTS:
        value : Primal scalar or numerical operand.
        size : Requested state/array size.

        OUTPUT ARGUMENTS:
        Primal state values and a state-by-direction sensitivity array.
        """
        base = value[:size]
        slopes = value[size:].reshape(size, 11) if use_ad else None
        return base, slopes

    def pack(values, slopes=None):
        """
        Pack primal states and optional forward-sensitivity lanes.

        SYNOPSIS:
            pack(...)

        INPUT ARGUMENTS:
        values : Numerical values or vector of state entries.
        slopes : Forward-sensitivity values for the associated primal
            entries.

        OUTPUT ARGUMENTS:
        Packed primal state vector with optional flattened sensitivity
        lanes.
        """
        return np.r_[values, slopes.ravel()] if use_ad else np.asarray(
            values)

    def value_slope(v):
        """
        Extract a primal value and its forward derivative lanes.

        SYNOPSIS:
            value_slope(...)

        INPUT ARGUMENTS:
        v : Primal numerical value, optionally carrying tangent lanes.

        OUTPUT ARGUMENTS:
        Primal value and its derivative vector in the seeded directions.
        """
        return primal(v), np.broadcast_to(v.slope, (11,)) if isinstance(v,
            Tangent) else zero

    def initialize(pilot, gbar=None):
        """
        Initialize water stores and detach the prescribed spin-up.

        SYNOPSIS:
            initialize(...)

        INPUT ARGUMENTS:
        pilot : Whether to run the reference-recharge initialization pass.
        gbar : Pilot mean recharge defining the groundwater reference.

        OUTPUT ARGUMENTS:
        pack(y, sy) : computed result in the operand units.
        """
        size = 7 if pilot else ns+5
        y, sy = np.zeros(size), np.zeros((size, 11))
        if initial is not None:
            y[:6 if pilot else ns] = initial[:6 if pilot else ns]
        else:
            y[5], sy[5] = value_slope(.7*c)
            if not pilot:
                y[6], sy[6] = value_slope(.5*x3)
                y[7], sy[7] = value_slope(gbar*tau)
        return pack(y, sy)

    def run(pilot, gbar=None):
        """
        Advance the pilot or full model over the forcing sequence.

        SYNOPSIS:
            run(...)

        INPUT ARGUMENTS:
        pilot : Whether to run the reference-recharge initialization pass.
        gbar : Pilot mean recharge defining the groundwater reference.

        OUTPUT ARGUMENTS:
        (start_state, np.asarray(ends), np.asarray(gradients)) : computed
        result in the operand units.
        """
        size = 7 if pilot else ns+5
        state = initialize(pilot, gbar)
        start_state = state[:size].copy()
        ends, gradients = [], []
        for interval, (rain, temperature, pet, _) in enumerate(f):
            dt = durations[interval]
            rain, pet = rain/dt, pet/dt
            warm = [temperature-gamma*((v-a@z)/1000)-pars[6] for v in z]
            cold = np.array([primal(v) < 0 for v in warm])
            if not pilot and any(abs(primal(v)) < 1e-9 for v in warm):
                counters['threshold_coincidences'].append(interval)
            snowfall = [pars[8]*rain if cold[k] else 0. for k in range(5)]
            melt_potential = [cf*max(v, 0.) for v in warm]
            snow0 = state[:5]
            active = np.array([snow0[k] > 0 or primal(snowfall[k]-ks*pet-
                melt_potential[k]) > 0
                               for k in range(5)])

            def rates(flags):
                """
                Evaluate state and cumulative-flux derivatives.

                SYNOPSIS:
                    rates(...)

                INPUT ARGUMENTS:
                flags : Active snow-store event flags.

                OUTPUT ARGUMENTS:
                (snow_rate, liquid - e0, pet - e0, e0 + sum((a[k] * sub[k]
                for k in range(5)))) : computed result in the operand units.
                """
                sub = [ks*pet if flags[k] else min(ks*pet, snowfall[k])
                    for k in range(5)]
                melt = [melt_potential[k] if flags[k] else
                        min(melt_potential[k], max(snowfall[k]-sub[k], 0.)
                            ) for k in range(5)]
                snow_rate = [snowfall[k]-sub[k]-melt[k] for k in range(5)]
                liquid = sum(a[k]*((0. if cold[k] else rain)+melt[k]) for
                    k in range(5))
                e0 = min(liquid, pet)
                return snow_rate, liquid-e0, pet-e0, e0+sum(a[k]*sub[k]
                    for k in range(5))

            def rhs(_, augmented, flags=active):
                """
                Evaluate augmented water-state and sensitivity derivatives.

                SYNOPSIS:
                    rhs(...)

                INPUT ARGUMENTS:
                _ : Reserved/unused callback argument.
                augmented : Packed state vector with optional
                    forward-sensitivity lanes.
                flags : Active snow-store event flags.

                OUTPUT ARGUMENTS:
                Computed fluxes, stores, or numerical values described in
                the notes.
                """
                counters['rhs_evaluations'] += 1
                y, sy = unpack(augmented, size)
                dy, dsy = np.zeros(size), np.zeros((size, 11))

                def assign(k, v):
                    """
                    Write a rate and its sensitivity lanes to the output
                    buffers.

                    SYNOPSIS:
                        assign(...)

                    INPUT ARGUMENTS:
                    k : State index or release coefficient in this helper.
                    v : Primal numerical value, optionally carrying tangent
                        lanes.

                    OUTPUT ARGUMENTS:
                    None. Writes primal and derivative rates to the
                    enclosing output.
                    """
                    dy[k], dsy[k] = value_slope(v)

                def dual(k):
                    """
                    Read a state as a primal or forward-mode scalar.

                    SYNOPSIS:
                        dual(...)

                    INPUT ARGUMENTS:
                    k : State index or release coefficient in this helper.

                    OUTPUT ARGUMENTS:
                    State entry as a plain or forward-mode scalar.
                    """
                    return Tangent(y[k], sy[k]) if use_ad else y[k]

                sn, pn, en, atmospheric = rates(flags)
                for k in range(5):
                    assign(k, sn[k])
                soil = dual(5)
                ratio = soil/c
                infiltrate = (1-ratio**2)*pn
                evap = ratio*(2-ratio)*en
                recharge = .25*kp*soil*ratio**4
                direct = pn-infiltrate
                assign(5, infiltrate-evap-recharge)
                if pilot:
                    assign(6, recharge)
                    return pack(dy, dsy)
                outlets = []
                for offset, count, alpha, eta in ((8, m, .9, fractions[0]),
                                                  (8+m, 2*m, .1, fractions
                                        [1])):
                    stores = y[offset:offset+count]
                    release = primal(lam)*stores
                    through = (1-eta)*release
                    dy[offset:offset+count] = np.r_[primal(alpha*direct),
                        through[:-1]]-release
                    if use_ad:
                        drelease = primal(lam)*sy[offset:offset+count]+ \
                            stores[:, None]*lam.slope
                        dthrough = (1-eta)[:, None]*drelease
                        dsy[offset:offset+count] = np.vstack((value_slope(
                            alpha*direct)[1], dthrough[:-1]))-drelease
                        outlets.append(Tangent(eta@release, eta@drelease))
                    else:
                        outlets.append(eta@release)
                u9, u1 = outlets
                # Negative internal solver stages are extended by positive
                # part;
                # accepted trajectories are checked rather than clipped
                # silently.
                quick, groundwater = max(dual(6), 0.), max(dual(7), 0.)
                loss = -pars[4]*(quick/x3)**3.5
                qr = .25*quick*(quick/x3)**4
                qb = max(u1-.1*loss, 0.)
                qg = groundwater**2/(gbar*tau**2)
                assign(6, u9-qr-.9*loss)
                assign(7, recharge-qg)
                assign(ns, qr+qb+qg)
                assign(ns+1, atmospheric+evap)
                assign(ns+2, .9*loss+min(u1, .1*loss))
                assign(ns+3, rain+(pars[8]-1)*rain*sum(a[cold]))
                assign(ns+4, recharge)
                return pack(dy, dsy)

            time = 0.
            while time < dt:
                sn = np.array([primal(v) for v in rates(active)[0]])
                candidates = np.full(5, np.inf)
                draining = active & (sn < 0)
                candidates[draining] = state[:5][draining]/(-sn[draining])
                duration = min(dt-time, float(candidates.min()))
                if duration <= 0:
                    raise RuntimeError("Nonpositive snow-event interval.")
                solution = solve_ivp(rhs, (time, time+duration), state,
                                     method='DOP853', rtol=rtol, atol=atol)
                if not solution.success or not np.isfinite(solution.y[:, -
                    1]).all():
                    raise RuntimeError(
                        f"CaRT coupled ODE failed: {solution.message}")
                state = solution.y[:, -1]
                time += duration
                hits = np.flatnonzero(np.isfinite(candidates) &
                                      (np.abs(candidates-duration) < 1e-11
                                        *max(1., dt)))
                for k in hits:
                    before = rhs(time, state).copy()
                    # Event-time derivative at fixed absolute time. The
                    # reset
                    # is identity; only the vector field changes
                    # (saltation).
                    te = -unpack(state, size)[1][k]/before[k] if use_ad \
                        else zero
                    active[k] = False
                    state[k] = 0.
                    after = rhs(time, state)
                    if use_ad:
                        state[size:] += ((before[:size]-after[:size])[:,
                            None]*te).ravel()
                        state[size:].reshape(size, 11)[k] = 0.
                    if not pilot:
                        counters['snow_events'] += 1
                        if abs(time-dt) < 1e-10:
                            counters['boundary_event_coincidences'].append \
                                ((interval, int(k)))
                if dt-time < 1e-14:
                    time = dt
            y, sy = unpack(state, size)
            if np.min(y[:6 if pilot else ns]) < -max(1e-7, 100*atol):
                raise RuntimeError('Accepted trajectory has a nega'
                    'tive storage; tighten toleranc' 'es.')
            if y[5] > primal(c)+max(1e-7, 100*atol):
                raise RuntimeError('Accepted production storage ex'
                    'ceeds capacity.')
            ends.append(y.copy())
            if use_ad:
                gradients.append(sy.copy())
        return start_state, np.asarray(ends), np.asarray(gradients)

    _, upstream, upstream_grad = run(True)
    train = f[:, 3].astype(bool)
    recharge = np.diff(np.r_[0., upstream[:, 6]])
    mean = recharge[train].sum()/durations[train].sum()
    if use_ad:
        dg = np.diff(np.vstack((zero, upstream_grad[:, 6])), axis=0)
        mean = Tangent(mean, dg[train].sum(axis=0)/durations[train].sum())
    gbar = max(mean, 1e-6)
    start, ends, gradients = run(False, gbar)
    depths = np.diff(np.r_[0., ends[:, ns]])
    jac = np.diff(np.vstack((zero, gradients[:, ns])), axis=0) if use_ad \
        else None
    if not return_diagnostics:
        return depths, jac
    volume = ends[:, :5]@a+ends[:, 5:ns].sum(axis=1)
    volume0 = start[:5]@a+start[5:ns].sum()
    ledger = ends[:, ns+3]-ends[:, ns+1]-ends[:, ns+2]-ends[:, ns]
    residual = np.diff(np.r_[0., volume-volume0-ledger])
    diagnostics = dict(counters, initial=start[:ns], states=ends[:, :ns],
                       state_jacobian=gradients[:, :ns] if use_ad else None,
                       interval_days=durations,
                       gbar=primal(gbar), gbar_gradient=value_slope(gbar)[
                           1],
                       balance_residual_mm=residual, refinement=m,
                       rtol=rtol, atol=atol,
                       initialization='detached' if initial is not None
                           else 'differentiated',
                       state_order='N1..N5,S,R,W,Z9_1..Z9_m,Z1_1..Z1_2m')
    return depths, jac, diagnostics
