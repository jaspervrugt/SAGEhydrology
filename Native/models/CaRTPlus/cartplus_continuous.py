# ========================================================================
# SAGEhydrology | Model source
# File: Native/models/CaRTPlus/cartplus_continuous.py
#
# CaRT / CaRT+ numerical model and supporting routines.
#
# SAGE adaptations: Jasper A. Vrugt, University of California, Irvine.
# Upstream model: Xia (2026). Preserve the adjacent license notices.
# ========================================================================

"""
Fully coupled continuous-time CaRT+ research reference.

All stores, routing cascades, flux integrals and forward sensitivities are
advanced in one ODE. The older cartplus_ode module remains the hybrid model.
Physical parameters follow cartplus_inputs.PARAMETER_MAP; all times are
days.
"""

from pathlib import Path
import sys

import numpy as np
from scipy.integrate import solve_ivp

from cartplus_inputs import PARAMETER_MAP
from cartplus_python import Tangent, primal

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'cart'))
try:
    from cart_ode import routing_fractions
finally:
    sys.path.pop(0)


def logistic(x):
    """
    Evaluate the logistic response with forward derivatives.

    SYNOPSIS:
        logistic(...)

    INPUT ARGUMENTS:
    x : Numerical operand; scalar, vector, or forward-mode scalar.

    OUTPUT ARGUMENTS:
    Dimensionless logistic response, preserving derivative lanes.
    """
    if primal(x) >= 0:
        return 1/(1+np.exp(-x))
    e = np.exp(x)
    return e/(1+e)


def bounded(x, low=0., high=1.):
    """
    Clip a scalar while retaining the selected tangent branch.

    SYNOPSIS:
        bounded(...)

    INPUT ARGUMENTS:
    x : Numerical operand; scalar, vector, or forward-mode scalar.
    low : Lower bound for the clipping operator.
    high : Upper bound for the clipping operator.

    OUTPUT ARGUMENTS:
    Bounded value and derivatives of the selected branch.
    """
    return min(max(x, low), high)


def excess_rate(water, capacity):
    """
    Evaluate smooth infiltration-excess runoff.

    SYNOPSIS:
        excess_rate(...)

    INPUT ARGUMENTS:
    water : Available water or store depth (mm).
    capacity : Store or infiltration capacity in this operator.

    OUTPUT ARGUMENTS:
    Infiltration-excess rate on the supplied forcing time base.
    """
    x = water/capacity
    if primal(x) < 1e-3:
        return capacity*(x*x/2-x**3/6+x**4/24-x**5/120)
    return water+capacity*np.expm1(-x)


def _inputs(forcing, elevation, area, parameters, initial, resolution,
            interval_days, refinement, rtol, atol, infiltration_floor,
            availability_days, transition_fraction,
                sensitivity_coordinates):
    """
    Validate model dimensions, units, and integration settings.

    SYNOPSIS:
        _inputs(...)

    INPUT ARGUMENTS:
    forcing : N-by-4 forcing [P, T, PET, training]; depths mm/interval, T
        deg C.
    elevation : Five elevation-band heights (m).
    area : Five nonnegative band fractions, summing to one.
    parameters : 42 physical CaRT+ parameters in native order.
    initial : Initial stores (mm), or None/empty to use the spin-up policy.
    resolution : Forcing intervals per day: 1 daily, 24 hourly, 96
        quarter-hourly.
    interval_days : Positive interval duration (days), scalar or length N.
    refinement : Positive number of routing stages per X4 time base.
    rtol : Positive relative integration-error tolerance (-).
    atol : Positive absolute integration-error tolerance for augmented
        states.
    infiltration_floor : Minimum infiltration capacity (mm/day).
    availability_days : Water-availability regularization time (days).
    transition_fraction : Dimensionless transition smoothing width.
    sensitivity_coordinates : Sensitivity directions: physical or
        transformed coordinates.

    OUTPUT ARGUMENTS:
    Validated forcing, bands, parameters, durations, and initial stores;
    CaRT+ also returns parameter derivative scales.
    """
    f, z, a, p = [np.asarray(v, dtype=float) for v in (forcing, elevation,
        area, parameters)]
    if (f.ndim != 2 or f.shape[1] != 4 or not len(f) or z.shape != (5,)
            or a.shape != (5,) or p.shape != (42,)
            or not all(np.isfinite(v).all() for v in (f, z, a, p))):
        raise ValueError('Expected finite (n,4) forcing,'
            ' five bands and 42 physical pa' 'rameters.')
    if (np.any(f[:, [0, 2]] < 0) or np.any(a < 0)
            or not np.isclose(a.sum(), 1., rtol=0, atol=1e-10)
            or not np.isin(f[:, 3], [0, 1]).all() or not f[:, 3].any()):
        raise ValueError('Require nonnegative depths/are'
            'as, normalized bands, and a 0/' '1 training mask.')
    positive = [1, 5, 6, 7, 8, 9, 10, 11, 14, 15, 16, 17, 18, 19,
                23, 24, 25, 26, 29, 33, 34, 38, 41]
    if (np.any(p[positive] <= 0) or np.any(p[[2, 3, 4, 12, 20, 30, 31, 36,
        39, 40]] < 0)
            or np.any((p[[13, 21, 27, 32]] < 0) | (p[[13, 21, 27, 32]] > 1))
            or p[35] > 0):
        raise ValueError('Physical capacities/times/widt'
            'hs must be positive and fracti' 'ons in [0,1]; exchange<=0.')
    if (resolution not in (1, 24, 96) or isinstance(refinement, bool)
            or not isinstance(refinement, (int, np.integer)) or refinement
                < 1
            or sensitivity_coordinates not in ('physical', 'calibration')
            or not np.isfinite([rtol, atol, infiltration_floor,
                availability_days, transition_fraction]).all()
            or min(rtol, atol, availability_days) <= 0 or
                infiltration_floor < 0
            or not 0 < transition_fraction <= 1):
        raise ValueError('Invalid resolution, refinement'
            ', coordinates, solver or regul' 'arization controls.')
    try:
        durations = np.broadcast_to(1/resolution if interval_days is None
                                    else np.asarray(interval_days, dtype=
                                        float), (len(f),)).copy()
    except ValueError as exc:
        raise ValueError('interval_days must be scalar or length n.') from \
            exc
    if not np.isfinite(durations).all() or np.any(durations <= 0):
        raise ValueError('Observation durations must be '
            'finite and positive.')
    ns = 14+3*int(refinement)
    if initial is not None:
        initial = np.asarray(initial, dtype=float)
        if (initial.shape != (ns,) or not np.isfinite(initial).all() or np
            .any(initial < 0)
                or np.any(initial[[5, 6, 7, 8, 9]] > p[[6, 8, 14, 25, 16]])
                or initial[11] > 1):
            raise ValueError('Expected ' f'{ns}'
                ' nonnegative initial states within t' 'heir capacities.')
    seed_scale = np.array([v if mult is not None else 1.
                           for v, (_, _, mult) in zip(p, PARAMETER_MAP)])
    if sensitivity_coordinates == 'physical':
        seed_scale[:] = 1.
    return f, z, a, p, durations, initial, seed_scale


def simulate(forcing, elevation, area, parameters, *, jacobian=False,
             initial=None, resolution=1, interval_days=None, refinement=16,
             rtol=1e-8, atol=1e-10, infiltration_floor=0.,
             availability_days=.02, transition_fraction=.01,
             sensitivity_coordinates='physical', return_diagnostics=False,
             first_step=None, sensitivity_columns=None):
    """
    Simulate model fluxes and optional discharge sensitivities.

    SYNOPSIS:
        simulate(...)

    INPUT ARGUMENTS:
    forcing : N-by-4 forcing [P, T, PET, training]; depths mm/interval, T
        deg C.
    elevation : Five elevation-band heights (m).
    area : Five nonnegative band fractions, summing to one.
    parameters : 42 physical CaRT+ parameters in native order.
    jacobian : Whether to evaluate discharge sensitivities.
    initial : Initial stores (mm), or None/empty to use the spin-up policy.
    resolution : Forcing intervals per day: 1 daily, 24 hourly, 96
        quarter-hourly.
    interval_days : Positive interval duration (days), scalar or length N.
    refinement : Positive number of routing stages per X4 time base.
    rtol : Positive relative integration-error tolerance (-).
    atol : Positive absolute integration-error tolerance for augmented
        states.
    infiltration_floor : Minimum infiltration capacity (mm/day).
    availability_days : Water-availability regularization time (days).
    transition_fraction : Dimensionless transition smoothing width.
    sensitivity_coordinates : Sensitivity directions: physical or
        transformed coordinates.
    return_diagnostics : Whether to append solver/water-balance diagnostics.
    first_step : Optional initial solver step (days).
    sensitivity_columns : Selected parameter directions, or None for all
        directions.

    OUTPUT ARGUMENTS:
    q : N discharge depths (mm/interval); J : discharge Jacobian or None.
    Diagnostics are appended only when requested.

    NOTES:
    Forcing columns are P depth, temperature, PET depth, training mask
    (0/1). interval_days optionally supplies arbitrary scalar or irregular
    durations. State order: N1..N5,H,S,N,B,A,J,V,W,R,Z9_1..m,Z1_1..2m. B,J,V
    are diagnostic/memory states, excluded from physical water storage.
    Supplied initial states are detached. Default capacity fractions, V's
    first-year forcing target mean and Wref are differentiated. No hidden
    spin-up is performed. A differentiated upstream pilot computes Wref.
    Returned derivatives default to physical parameter coordinates;
    selecting 'calibration' applies the log-coordinate chain rule in
    PARAMETER_MAP, retaining all 42 columns (35 fitted and seven shared
    constants).
    The fixed availability time and transition fraction regularize snow,
    surface availability/overflow and Pareto retention, as documented in the
    paper. They never scale with the observation interval or ODE step size.
    Piecewise rates are continuous at internal clips, requiring no
    saltation; rain/snow threshold coincidences at forcing boundaries are
    flagged.
    """
    f, z, a, p, durations, initial, seed_scale = _inputs(
        forcing, elevation, area, parameters, initial, resolution,
            interval_days,
        refinement, rtol, atol, infiltration_floor, availability_days,
        transition_fraction, sensitivity_coordinates)
    columns = np.arange(42) if sensitivity_columns is None else np.asarray \
        (sensitivity_columns)
    if (columns.ndim != 1 or columns.dtype.kind not in "iu" or not len(
        columns)
            or np.any((columns < 0) | (columns >= 42))
            or len(np.unique(columns)) != len(columns)):
        raise ValueError('Expected distinct CaRT+ sensit'
            'ivity column indices in [0, 42' ').')
    m, ns, npar = int(refinement), 14+3*int(refinement), len(columns)
    params = [Tangent(v, seed) for v, seed in zip(p, np.diag(seed_scale)[:
        , columns])] if jacobian else p
    (threshold, cf, sf, ks, km, lapse, ch, th, cs, kp, bs, ki, gi, an,
     cn, tn, ca, ba, ka, tg, eg_factor, veg_reduction, tv, wv, tveg,
     cb, tb, cmin, ccenter, cw, soil_connection, groundwater_connection,
     beta, routing_days, cr, exchange, exterior, tf, wf, frost, exponent,
         cch) = params
    upland, lam = 1-an, m/routing_days
    etas = routing_fractions(m)
    train = f[:, 3].astype(bool)
    mean_rain = max(f[train, 0].sum()/durations[train].sum(), .1)
    veg_targets = [logistic((temp-tv)/wv) for temp in f[:, 1]]
    starts = np.r_[0., np.cumsum(durations)[:-1]]
    year_weights = np.minimum(durations, np.maximum(365.-starts, 0.))
    veg0 = sum(v*w for v, w in zip(veg_targets, year_weights))/ \
        year_weights.sum()
    zero = np.zeros(npar)
    counters = dict(rhs_evaluations=0, threshold_coincidences=[])

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
        return primal(v), np.broadcast_to(v.slope, (npar,)) if isinstance(
            v, Tangent) else zero

    def pack(y, sy):
        """
        Pack primal states and optional forward-sensitivity lanes.

        SYNOPSIS:
            pack(...)

        INPUT ARGUMENTS:
        y : Numerical state/operand; augmented with sensitivities when
            enabled.
        sy : State sensitivities in selected parameter directions.

        OUTPUT ARGUMENTS:
        Packed primal state vector with optional flattened sensitivity
        lanes.
        """
        return np.r_[y, sy.ravel()] if jacobian else y

    def unpack(y, count):
        """
        Separate primal states from forward-sensitivity lanes.

        SYNOPSIS:
            unpack(...)

        INPUT ARGUMENTS:
        y : Numerical state/operand; augmented with sensitivities when
            enabled.
        count : Number of intervals or state entries in this helper.

        OUTPUT ARGUMENTS:
        Primal state values and a state-by-direction sensitivity array.
        """
        return y[:count], y[count:].reshape(count, npar) if jacobian else \
            None

    def pareto(store, rain, pet, capacity, shape, perc, net):
        """
        Advance the Pareto production-store water balance.

        SYNOPSIS:
            pareto(...)

        INPUT ARGUMENTS:
        store : Current water store (mm).
        rain : Liquid input depth/rate on this operator time base.
        pet : Potential evaporation depth/rate on this operator time base.
        capacity : Store or infiltration capacity in this operator.
        shape : Dimensionless storage/infiltration shape parameter.
        perc : Percolation coefficient in the model parameterization.
        net : Whether to remove interception before the production update.

        OUTPUT ARGUMENTS:
        Updated storage, runoff, recharge, and evaporation for the interval;
        the continuous helper returns flux rates.
        """
        wet = bounded(store/capacity)
        deficit = 1-wet
        nu = shape/(1+shape)
        intercept = min(rain, pet) if net else 0.
        retain = (deficit/transition_fraction*np.exp(nu*np.log(
            transition_fraction))
                  if primal(deficit) < transition_fraction else deficit**nu)
        return ((rain-intercept)*(1-retain),
                intercept+(pet-intercept)*(1-deficit**2),
                .25*max(store, 0.)*np.log1p(perc*wet**4))

    def zone(store, rain, pet, capacity, days):
        """
        Evaluate near-stream runoff, evaporation, and drainage.

        SYNOPSIS:
            zone(...)

        INPUT ARGUMENTS:
        store : Current water store (mm).
        rain : Liquid input depth/rate on this operator time base.
        pet : Potential evaporation depth/rate on this operator time base.
        capacity : Store or infiltration capacity in this operator.
        days : Residence time (days).

        OUTPUT ARGUMENTS:
        Near-stream runoff, evaporation, and drainage rates.
        """
        wet = bounded(store/capacity)
        intercept = min(rain, pet)
        return ((rain-intercept)*wet**2,
                intercept+(pet-intercept)*wet*(2-wet), max(store, 0.)/days)

    def initialize(pilot, reference=None):
        """
        Initialize water stores and detach the prescribed spin-up.

        SYNOPSIS:
            initialize(...)

        INPUT ARGUMENTS:
        pilot : Whether to run the reference-recharge initialization pass.
        reference : Pilot mean recharge or reference storage used by the
            caller.

        OUTPUT ARGUMENTS:
        pack(y, sy) : computed result in the operand units.
        """
        count = 13 if pilot else ns+7
        y, sy = np.zeros(count), np.zeros((count, npar))
        if initial is not None:
            y[:12 if pilot else ns] = initial[:12 if pilot else ns]
        else:
            for k, v in ((5, .5*ch), (6, .7*cs), (7, .3*cn), (8, .3*cb),
                         (9, .5*ca), (10, mean_rain), (11, veg0)):
                y[k], sy[k] = value_slope(v)
            if not pilot:
                y[12], sy[12] = value_slope(reference)
                y[13], sy[13] = value_slope(.5*cr)
        return pack(y, sy)

    def run(pilot, reference=None):
        """
        Advance the pilot or full model over the forcing sequence.

        SYNOPSIS:
            run(...)

        INPUT ARGUMENTS:
        pilot : Whether to run the reference-recharge initialization pass.
        reference : Pilot mean recharge or reference storage used by the
            caller.

        OUTPUT ARGUMENTS:
        (first, np.asarray(ends), np.asarray(gradients)) : computed result
        in the operand units.
        """
        count = 13 if pilot else ns+7
        state = initialize(pilot, reference)
        first = state[:count].copy()
        ends, gradients = [], []
        pilot_end = np.flatnonzero(train)[-1]+1 if pilot else len(f)
        for i, (rain, temp, pet, _) in enumerate(f[:pilot_end]):
            dt = durations[i]
            rain, pet = rain/dt, pet/dt
            warm = [temp-lapse*((height-a@z)/1000)-threshold for height in
                z]
            cold = [primal(t) < 0 for t in warm]
            snowfall = [sf*rain if is_cold else 0. for is_cold in cold]
            melt_demand = [0. if cold[b] else cf*warm[b]+km*pet for b in
                range(5)]
            liquid_rain = rain*sum(a[b] for b in range(5) if not cold[b])
            snow_adjustment = (sf-1)*rain*sum(a[b] for b in range(5) if
                cold[b])
            if not pilot and any(abs(primal(t)) < 1e-9 for t in warm):
                counters['threshold_coincidences'].append(i)
            kelvin = max(temp, 0.)+273.15
            viscosity = -np.log(10)*247.8*(1/(kelvin-140)-1/153.15)
            frozen = frost*logistic((tf-temp)/wf)

            def rhs(_, augmented):
                """
                Evaluate augmented water-state and sensitivity derivatives.

                SYNOPSIS:
                    rhs(...)

                INPUT ARGUMENTS:
                _ : Reserved/unused callback argument.
                augmented : Packed state vector with optional
                    forward-sensitivity lanes.

                OUTPUT ARGUMENTS:
                Computed fluxes, stores, or numerical values described in
                the notes.
                """
                counters['rhs_evaluations'] += 1
                y, sy = unpack(augmented, count)
                dy, dsy = np.zeros(count), np.zeros((count, npar))

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
                    return Tangent(y[k], sy[k]) if jacobian else y[k]

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

                subs, melts = [], []
                for b in range(5):
                    available = max(dual(b), 0.)/availability_days
                    sub = min(ks*pet, available)
                    melt = min(max(available-sub, 0.), melt_demand[b])
                    subs.append(sub)
                    melts.append(melt)
                    assign(b, snowfall[b]-sub-melt)
                liquid = liquid_rain+sum(a[b]*melts[b] for b in range(5))
                h, s, n, wet, depression, memory, veg = [dual(k) for k in
                    range(5, 12)]
                demand = pet*(1-veg_reduction*(1-bounded(veg)))
                eh = min(demand, max(h, 0.)/availability_days)
                qh = max(h, 0.)/th
                fh = qh+max(liquid-eh-qh, 0.)*bounded((h/ch-(1-
                    transition_fraction))/transition_fraction)
                ep = max(demand-eh, 0.)
                capacity = max(ki*np.exp(-gi*max(memory, 0.)/mean_rain),
                    infiltration_floor, 1e-12)
                qx = excess_rate(max(fh-ep, 0.), capacity)
                ps = fh-qx
                qs, es, ds = pareto(s, ps, ep, cs, bs, kp, True)
                qn, en, fn = zone(n, ps, ep, cn, tn)
                qw, ew, fw = zone(wet, ps, ep, cb, tb)
                ez = upland*es+an*en
                fast = qx+upland*qs+an*qn
                qa, ea, da = pareto(depression, fast, max(ep-ez, 0.), ca,
                    ba, ka, False)
                recharge = upland*ds+da
                assign(5, liquid-eh-fh)
                assign(6, ps-qs-es-ds)
                assign(7, ps-qn-en-fn)
                assign(8, ps-qw-ew-fw)
                assign(9, fast-qa-ea-da)
                assign(10, (fh-memory)/3)
                assign(11, (veg_targets[i]-veg)/tveg)
                if pilot:
                    assign(12, recharge)
                    return pack(dy, dsy)
                groundwater, quick = max(dual(12), 0.), max(dual(13), 0.)
                qg = groundwater**2/(reference*tg)
                eg = eg_factor*groundwater*max(ep-ez-ea, 0.)/(groundwater+
                    reference)
                connected = cmin+(1-cmin)*logistic((wet/cb-ccenter+
                    soil_connection*bounded(s/cs)+groundwater_connection*(
                        groundwater/reference-1))/cw)
                bypass_input = beta*connected*qa
                main_input = qa-bypass_input
                outlets = []
                for offset, length, incoming, eta in ((14, m, main_input,
                    etas[0]),
                                                     (14+m, 2*m,
                                        bypass_input, etas[1])):
                    stores = y[offset:offset+length]
                    release = primal(lam)*stores
                    through = (1-eta)*release
                    dy[offset:offset+length] = np.r_[primal(incoming),
                        through[:-1]]-release
                    if jacobian:
                        drelease = primal(lam)*sy[offset:offset+length]+ \
                            stores[:, None]*lam.slope
                        dsy[offset:offset+length] = np.vstack((value_slope
                            (incoming)[1],
                            ((1-eta)[:, None]*drelease)[:-1]))-drelease
                        outlets.append(Tangent(eta@release, eta@drelease))
                    else:
                        outlets.append(eta@release)
                u9, u1 = outlets
                exchange_rate = exchange*(quick/cr)**3.5
                ir = u9+.9*(exterior+exchange_rate)
                if primal(quick) <= 0:
                    ir = max(ir, 0.)
                qb = max(u1+.1*(exterior+exchange_rate), 0.)
                signal = np.log(max(connected, 1e-12))+viscosity-frozen
                modulation = signal+exponent*np.log(max(quick, 1e-3)/cr)
                qr = .25*quick*(quick/cr)**4*np.exp(6*np.tanh(modulation/6))
                preloss = qr+qb+an*fn+qg
                reported = excess_rate(preloss, cch)
                channel_loss = preloss-reported
                net_external = ir+qb-u9-u1
                assign(12, recharge-qg-eg)
                assign(13, ir-qr)
                assign(ns, reported)
                assign(ns+1, eh+ez+ea+eg+sum(a[b]*subs[b] for b in range(5
                    )))
                assign(ns+2, net_external)
                assign(ns+3, rain+snow_adjustment)
                assign(ns+4, recharge)
                assign(ns+5, channel_loss)
                assign(ns+6, preloss)
                return pack(dy, dsy)

            result = solve_ivp(rhs, (0., dt), state, method='DOP853', rtol
                =rtol, atol=atol,
                               first_step=None if first_step is None else
                                   min(dt, first_step))
            if not result.success or not np.isfinite(result.y[:, -1]).all():
                raise RuntimeError('CaRT+ coupled integration failed: '
                    f'{result.message}')
            state = result.y[:, -1]
            y, sy = unpack(state, count)
            if np.min(y[:12 if pilot else ns]) < -max(1e-7, 100*atol):
                raise RuntimeError('Negative accepted CaRT+ state;'
                    ' tighten solver tolerances.')
            if np.any(y[[5, 6, 7, 8, 9]] > p[[6, 8, 14, 25, 16]]+max(1e-7,
                100*atol)):
                raise RuntimeError(
                    'A CaRT+ accepted state exceeds capacity.')
            ends.append(y.copy())
            if jacobian:
                gradients.append(sy.copy())
        return first, np.asarray(ends), np.asarray(gradients)

    _, pilot, spilot = run(True)
    recharge = np.diff(np.r_[0., pilot[:, 12]])
    mean_recharge = recharge[train[:len(recharge)]].sum()/durations[train] \
        .sum()
    if jacobian:
        ds = np.diff(np.vstack((zero, spilot[:, 12])), axis=0)
        mean_recharge = Tangent(mean_recharge, ds[train[:len(ds)]].sum(
            axis=0)/durations[train].sum())
    reference = max(mean_recharge, 1e-6)*tg
    first, ends, gradients = run(False, reference)
    q = np.diff(np.r_[0., ends[:, ns]])
    j = np.diff(np.vstack((zero, gradients[:, ns])), axis=0) if jacobian \
        else None
    if not return_diagnostics:
        return q, j
    # Diagnostic B and memories J,V carry no basin water.
    weights = np.r_[a, 1., primal(upland), primal(an), 0., 1., 0., 0., 1.,
        1., np.ones(3*m)]
    volume = ends[:, :ns]@weights
    volume0 = first[:ns]@weights
    ledger = ends[:, ns+3]+ends[:, ns+2]-ends[:, ns+1]-ends[:, ns+5]-ends[
        :, ns]
    residual = np.diff(np.r_[0., volume-volume0-ledger])
    diagnostics = dict(counters, initial=first[:ns], states=ends[:, :ns],
                       state_jacobian=gradients[:, :ns] if jacobian else
                           None,
                       balance_residual_mm=residual, reference_storage=
                           primal(reference),
                       reference_gradient=value_slope(reference)[1],
                       channel_loss_mm=np.diff(np.r_[0., ends[:, ns+5]]),
                       preloss_depth_mm=np.diff(np.r_[0., ends[:, ns+6]]),
                       interval_days=durations, refinement=m, rtol=rtol,
                           atol=atol,
                       sensitivity_coordinates=sensitivity_coordinates,
                       initialization='detached' if initial is not None
                           else 'differentiated',
                       availability_days=availability_days,
                           transition_fraction=transition_fraction,
                       state_order=
                           'N1..N5,H,S,N,B,A,J,V,W,R,Z9_1..m,Z1_1..2m')
    return q, j, diagnostics
