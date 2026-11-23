# ========================================================================
# SAGEhydrology | Model source
# File: Native/models/CaRTPlus/cartplus_python.py
#
# CaRT / CaRT+ numerical model and supporting routines.
#
# SAGE adaptations: Jasper A. Vrugt, University of California, Irvine.
# Upstream model: Xia (2026). Preserve the adjacent license notices.
# ========================================================================

"""
Standalone CaRT+ water balance and routing for Python cases 1-3.
"""


from math import ceil

import numpy as np


def simulate(forcing, elevation, area, p, resolution, floor, zone,
             fixed=None):
    """
    Simulate model fluxes and optional discharge sensitivities.

    SYNOPSIS:
        simulate(...)

    INPUT ARGUMENTS:
    forcing : N-by-4 forcing [P, T, PET, training]; depths mm/interval, T
        deg C.
    elevation : Five elevation-band heights (m).
    area : Five nonnegative band fractions, summing to one.
    p : 42 physical CaRT+ parameters in native order.
    resolution : Forcing intervals per day: 1 daily, 24 hourly, 96
        quarter-hourly.
    floor : Minimum infiltration capacity in daily-rate units (mm/day).
    zone : Near-stream zone integration callback.
    fixed : Detached spin-up stores, or None/empty to determine them.

    OUTPUT ARGUMENTS:
    output : N-by-4 [Q, ET, external flux, storage change] (mm/interval);
    captured : detached spin-up stores.
    """
    dtype = object if isinstance(p[0], Tangent) else float
    dt = 1 / resolution
    n = len(forcing)
    spin = min(n, 1460 * resolution)
    quick_spin = min(n, 730 * resolution)
    train = forcing[:, 3].astype(bool)
    if not train.any():
        raise ValueError("CaRT+ requires training intervals.")
    (threshold, melt_rate, snowfall_factor, sublimation_factor, pet_melt,
     lapse, surface_cap, surface_days, soil_cap, soil_perc, soil_shape,
     infiltration_cap, infiltration_shape, zone_fraction, zone_cap,
     zone_days, depression_cap, depression_shape, depression_perc,
     groundwater_days, groundwater_uptake, vegetation_reduction,
     vegetation_temperature, vegetation_width, vegetation_days,
     wet_cap, wet_days, connected_min, connected_center, connected_width,
     soil_connection, groundwater_connection, bypass_fraction, routing_days,
     quick_cap, exchange_rate, exterior_rate, frost_temperature,
     frost_width, frost_reduction, quick_exponent, channel_cap) = p
    precip, temperature, pet = forcing[:, :3].T
    vegetation = logistic(
        (temperature - vegetation_temperature) / vegetation_width
    )
    baseline = vegetation[:365 * resolution].mean()
    mean_rain = max(precip[train].mean(), .1 * dt)
    vegetation_memory = memory(
        vegetation - baseline, np.exp(-dt / vegetation_days),
        2048 * resolution,
    )
    demand = pet * (1 - vegetation_reduction * (
        1 - np.clip(baseline + vegetation_memory, 0, 1)
    ))
    snow = np.full(5, Tangent(0) if dtype is object else 0., dtype=dtype)
    snow_fluxes = np.zeros((4, n), dtype=dtype)
    rain, sublimation, snow_extra, snow_change = snow_fluxes
    band_offset = lapse * (elevation - area @ elevation) / 1000
    for i in range(n):
        before = area @ snow
        warmth = temperature[i] - band_offset - threshold
        cold = warmth < 0
        snow += cold * snowfall_factor * precip[i]
        sub = np.minimum(snow, sublimation_factor * pet[i])
        snow -= sub
        melt = np.minimum(
            snow, np.where(cold, 0, melt_rate*dt*warmth + pet_melt*pet[i]),
        )
        snow -= melt
        rain[i] = area @ (melt + (~cold)*precip[i])
        sublimation[i] = area @ sub
        snow_extra[i] = (area @ cold) * precip[i] * (snowfall_factor - 1)
        snow_change[i] = area @ snow - before

    drainage = -np.expm1(-dt / surface_days)
    surface, warm_rain, warm_pet, _, _ = surface_series(
        .5*surface_cap, rain[:spin], demand[:spin], surface_cap, drainage,
    )
    surface = fixed[0] if fixed is not None else surface
    surface_start = surface
    _, rain, demand, surface_et, change = surface_series(
        surface, rain, demand, surface_cap, drainage,
    )
    warm_rain, _ = infiltrate(
        warm_rain, warm_pet, mean_rain, infiltration_cap,
        infiltration_shape, floor, resolution,
    )
    rain, excess = infiltrate(
        rain, demand, mean_rain, infiltration_cap, infiltration_shape,
        floor, resolution,
    )
    soil, near, wet = .7*soil_cap, .3*zone_cap, .3*wet_cap
    if fixed is None:
        for water, potential in zip(warm_rain, warm_pet):
            soil = pareto(
                soil, water, potential, soil_cap, soil_perc,
                soil_shape, dt, True,
            )[0]
            near = zone_step(near, water, potential, zone_cap,
                             zone_days, dt, zone)[0]
            wet = zone_step(wet, water, potential, wet_cap,
                            wet_days, dt, zone)[0]
    else:
        soil, near, wet = fixed[1:4]
    soil_start, near_start, wet_start = soil, near, wet
    fluxes = np.zeros((6, n), dtype=dtype)
    fast, recharge, et, near_q, saturation, wetness = fluxes
    upland = 1 - zone_fraction
    for i, (water, potential) in enumerate(zip(rain, demand)):
        before = upland*soil + zone_fraction*near
        soil, runoff, percolation, evap = pareto(
            soil, water, potential, soil_cap, soil_perc, soil_shape,
            dt, True,
        )
        near, zone_runoff, zone_q, zone_et = zone_step(
            near, water, potential, zone_cap, zone_days, dt, zone,
        )
        wet = zone_step(wet, water, potential, wet_cap,
                        wet_days, dt, zone)[0]
        fast[i] = upland*runoff + zone_fraction*zone_runoff + excess[i]
        recharge[i] = upland*percolation
        et[i] = upland*evap + zone_fraction*zone_et
        near_q[i] = zone_fraction*zone_q
        saturation[i] = np.clip(soil/soil_cap, 0, 1)
        wetness[i] = wet/wet_cap
        change[i] += upland*soil + zone_fraction*near - before

    depression = .5*depression_cap
    if fixed is None:
        for i in range(spin):
            depression = pareto(
                depression, fast[i], max(demand[i]-et[i], 0),
                depression_cap, depression_perc, depression_shape,
                dt, False,
            )[0]
    else:
        depression = fixed[4]
    depression_start = depression
    for i in range(n):
        before = depression
        depression, fast[i], add, evap = pareto(
            depression, fast[i], max(demand[i]-et[i], 0),
            depression_cap, depression_perc, depression_shape, dt, False,
        )
        recharge[i] += add
        et[i] += evap
        change[i] += depression - before

    equilibrium = max(recharge[train].mean()/dt, 1e-6)*groundwater_days
    groundwater_start = equilibrium if fixed is None else fixed[5]
    groundwater = groundwater_start
    recession = dt/(equilibrium*groundwater_days)
    base_q, connected = np.zeros((2, n), dtype=dtype)
    for i in range(n):
        before = groundwater
        groundwater += recharge[i]
        x = groundwater_uptake*max(demand[i]-et[i], 0)
        x /= groundwater + groundwater_start
        uptake = groundwater * (
            x - .5*x*x + x*x*x/6 if x < 1e-3 else -np.expm1(-x)
        )
        groundwater -= uptake
        next_store = groundwater/(1+recession*groundwater)
        base_q[i] = groundwater-next_store
        groundwater = next_store
        change[i] += groundwater-before
        et[i] += uptake
        signal = wetness[i] - connected_center
        signal += soil_connection*saturation[i]
        signal += groundwater_connection*(groundwater/groundwater_start-1)
        connected[i] = connected_min + (1-connected_min)*logistic(
            signal/connected_width
        )

    direct = bypass_fraction*connected*fast
    routed = hydrograph(fast-direct, routing_days, dt, False)
    bypass = hydrograph(direct, routing_days, dt, True)
    change += fast-routed-bypass
    exterior = exterior_rate*dt
    routed += .9*exterior
    bypass += .1*exterior
    kelvin = np.maximum(temperature, 0)+273.15
    viscosity = -np.log(10)*247.8*(1/(kelvin-140)-1/153.15)
    connectivity = np.log(connected)+viscosity-frost_reduction*logistic(
        (frost_temperature-temperature)/frost_width
    )
    quick = .5*quick_cap
    release = .25*dt*quick_cap**-4
    exchange = dt*exchange_rate*(30/quick_cap)**3.5

    def advance_quick(store, i):
        """
        Advance the quick store for the selected observation.

        SYNOPSIS:
            advance_quick(...)

        INPUT ARGUMENTS:
        store : Current water store (mm).
        i : Observation/state index.

        OUTPUT ARGUMENTS:
        (next_store, available - next_store + direct_q, exchanged) :
        computed result in the operand units.
        """
        flux = exchange*(max(store, 0)/30)**3.5
        available = max(store+routed[i]+.9*flux, 0)
        direct_q = max(bypass[i]+.1*flux, 0)
        log_factor = connectivity[i]+quick_exponent*np.log(
            max(available, 1e-3)/quick_cap
        )
        rate = release*np.exp(6*np.tanh(log_factor/6))
        next_store = available*np.exp(-.25*np.log1p(4*rate*available**4))
        exchanged = available-store-routed[i]+direct_q-bypass[i]
        return next_store, available-next_store+direct_q, exchanged

    if fixed is None:
        for i in range(quick_spin):
            quick = advance_quick(quick, i)[0]
    else:
        quick = fixed[6]
    captured = (surface_start, soil_start, near_start, wet_start,
                depression_start, groundwater_start, quick)
    output = np.empty((n, 4), dtype=dtype)
    for i in range(n):
        before = quick
        quick, flow, exchanged = advance_quick(quick, i)
        flow += base_q[i]+near_q[i]
        cap = dt*channel_cap
        x = max(flow, 0)/cap
        loss = cap*(x-.5*x*x if x < 1e-3 else -np.expm1(-x))
        output[i] = (flow-loss, et[i]+surface_et[i]+sublimation[i],
                     snow_extra[i]+exterior+exchanged-loss,
                     snow_change[i]+change[i]+quick-before)
    return output, captured


def logistic(value):
    """
    Evaluate the logistic response with forward derivatives.

    SYNOPSIS:
        logistic(...)

    INPUT ARGUMENTS:
    value : Primal scalar or numerical operand.

    OUTPUT ARGUMENTS:
    Dimensionless logistic response, preserving derivative lanes.
    """
    return 1/(1+np.exp(-value))


def memory(values, decay, window):
    """
    Apply the finite exponentially weighted memory kernel.

    SYNOPSIS:
        memory(...)

    INPUT ARGUMENTS:
    values : Numerical values or vector of state entries.
    decay : Memory-retention coefficient (-).
    window : Finite memory-kernel length in observation intervals.

    OUTPUT ARGUMENTS:
    Exponentially weighted memory series in the input value units.
    """
    out = np.empty_like(values)
    previous = 0.
    tail = (1-decay)*decay**window
    for i, value in enumerate(values):
        previous = decay*previous+(1-decay)*value
        if i >= window:
            previous -= tail*values[i-window]
        out[i] = previous
    return out


def surface_series(store, rain, pet, capacity, fraction):
    """
    Advance surface storage, overflow, drainage, and evaporation.

    SYNOPSIS:
        surface_series(...)

    INPUT ARGUMENTS:
    store : Current water store (mm).
    rain : Liquid input depth/rate on this operator time base.
    pet : Potential evaporation depth/rate on this operator time base.
    capacity : Store or infiltration capacity in this operator.
    fraction : Dimensionless partition/drainage fraction.

    OUTPUT ARGUMENTS:
    Final store and interval runoff, residual demand, evaporation, and
    storage-change series.
    """
    runoff, demand, et, change = np.zeros((4, len(rain)), dtype=rain.dtype)
    for i, (water, potential) in enumerate(zip(rain, pet)):
        before = store
        store += water
        et[i] = min(potential, store)
        store -= et[i]
        overflow = max(store-capacity, 0)
        store -= overflow
        drainage = fraction*store
        store -= drainage
        runoff[i] = overflow+drainage
        demand[i] = max(potential-et[i], 0)
        change[i] = store-before
    return store, runoff, demand, et, change


def infiltrate(rain, pet, mean_rain, capacity, shape, floor, resolution):
    """
    Partition liquid input into infiltration and excess runoff.

    SYNOPSIS:
        infiltrate(...)

    INPUT ARGUMENTS:
    rain : Liquid input depth/rate on this operator time base.
    pet : Potential evaporation depth/rate on this operator time base.
    mean_rain : Mean training rainfall on the forcing time base.
    capacity : Store or infiltration capacity in this operator.
    shape : Dimensionless storage/infiltration shape parameter.
    floor : Minimum infiltration capacity in daily-rate units (mm/day).
    resolution : Forcing intervals per day: 1 daily, 24 hourly, 96
        quarter-hourly.

    OUTPUT ARGUMENTS:
    Infiltrating liquid input and infiltration-excess runoff.
    """
    dt = 1/resolution
    decay = np.exp(-dt/3)
    history = memory(rain, decay, 36*resolution)
    previous = np.r_[0., history[:-1]]
    cap = dt*np.maximum(capacity*np.exp(-shape*(
        previous/mean_rain+decay**np.arange(len(rain))
    )), floor)
    available = np.maximum(rain-pet, 0)
    x = available/cap
    excess = available-cap*np.where(x < 1e-3, x-.5*x*x, -np.expm1(-x))
    return rain-excess, excess


def pareto(store, rain, demand, cap, percolation, shape, dt, net):
    """
    Advance the Pareto production-store water balance.

    SYNOPSIS:
        pareto(...)

    INPUT ARGUMENTS:
    store : Current water store (mm).
    rain : Liquid input depth/rate on this operator time base.
    demand : Evaporation demand on this operator time base.
    cap : Store capacity (mm).
    percolation : Percolation coefficient in the model parameterization.
    shape : Dimensionless storage/infiltration shape parameter.
    dt : Integration or observation interval duration (days).
    net : Whether to remove interception before the production update.

    OUTPUT ARGUMENTS:
    Updated storage, runoff, recharge, and evaporation for the interval; the
    continuous helper returns flux rates.
    """
    water = max(rain-demand, 0) if net else rain
    potential = max(demand-rain, 0) if net else demand
    interception = min(rain, demand) if net else 0
    deficit = max(1-np.clip(store/cap, 0, 1), 1e-7)
    exponent = shape+1
    remaining = deficit**(1/exponent)-water/(exponent*cap)
    retained = np.clip(cap*(deficit-max(remaining, 0)**exponent), 0, water)
    response = np.tanh(potential/cap)
    dry = deficit**2
    evap = min(cap*(1-dry)*response/(1+dry/deficit*response), store)
    store += retained-evap
    recharge = -store*np.expm1(-.25*dt*np.log1p(
        percolation*(store/cap)**4
    ))
    return store-recharge, water-retained, recharge, evap+interception


def zone_step(store, rain, demand, cap, residence, dt, integrate):
    """
    Advance near-stream storage and its interval fluxes.

    SYNOPSIS:
        zone_step(...)

    INPUT ARGUMENTS:
    store : Current water store (mm).
    rain : Liquid input depth/rate on this operator time base.
    demand : Evaporation demand on this operator time base.
    cap : Store capacity (mm).
    residence : Residence time (days).
    dt : Integration or observation interval duration (days).
    integrate : Zone-integration callback supplied by the backend.

    OUTPUT ARGUMENTS:
    Updated store, runoff, drainage, and evaporation for the interval.
    """
    net = (rain-demand)/cap
    rate = dt/residence
    next_fraction, integral = integrate(store/cap, net, rate, dt)
    next_store = cap*next_fraction
    discharge = cap*rate*integral
    change = next_store-store
    runoff = rain-demand-discharge-change if net > 0 else 0
    return next_store, runoff, discharge, rain-runoff-discharge-change


def hydrograph(inflow, days, dt, symmetric):
    """
    Route input depths through the selected unit hydrograph.

    SYNOPSIS:
        hydrograph(...)

    INPUT ARGUMENTS:
    inflow : Input series to be routed, in interval depths (mm).
    days : Residence time (days).
    dt : Integration or observation interval duration (days).
    symmetric : Whether to use the symmetric two-time-base hydrograph.

    OUTPUT ARGUMENTS:
    Routed interval depths (mm), with the input series length.
    """
    length = (13 if symmetric else 6) if dt == 1 else ceil(
        primal(days)*(2 if symmetric else 1)/dt
    )
    lag = np.arange(length+1)*dt/days
    if symmetric:
        cumulative = np.where(lag <= 1, .5*lag**2.5,
                              1-.5*np.maximum(2-lag, 0)**2.5)
    else:
        cumulative = np.minimum(lag, 1)**2.5
    return np.convolve(inflow, np.diff(cumulative))[:len(inflow)]


def primal(value):
    """
    Extract the real/primal value without derivative lanes.

    SYNOPSIS:
        primal(...)

    INPUT ARGUMENTS:
    value : Primal scalar or numerical operand.

    OUTPUT ARGUMENTS:
    Primal scalar value without forward derivative lanes.
    """
    return value.value if isinstance(value, Tangent) else value


class Tangent:
    'Scalar forward-mode value and ' 'its parameter gradient.'

    def __init__(self, value, slope=0.):
        """
        Construct a scalar with primal value and derivative lanes.

        SYNOPSIS:
            __init__(...)

        INPUT ARGUMENTS:
        self : Current forward-mode scalar, including its derivative lanes.
        value : Primal scalar or numerical operand.
        slope : Forward derivative or vector of parameter-direction
            derivatives.

        OUTPUT ARGUMENTS:
        None. Updates the supplied object or checks the stated conditions.
        """
        self.value = value
        self.slope = slope

    @staticmethod
    def lift(value):
        """
        Promote an ordinary value to a zero-derivative tangent.

        SYNOPSIS:
            lift(...)

        INPUT ARGUMENTS:
        value : Primal scalar or numerical operand.

        OUTPUT ARGUMENTS:
        value if isinstance(value, Tangent) else Tangent(value) : computed
        result in the operand units.
        """
        return value if isinstance(value, Tangent) else Tangent(value)

    def __add__(self, other):
        """
        Add operands and propagate forward derivatives.

        SYNOPSIS:
            __add__(...)

        INPUT ARGUMENTS:
        self : Current forward-mode scalar, including its derivative lanes.
        other : Second scalar/array operand, promoted when required.

        OUTPUT ARGUMENTS:
        Numerical result with propagated derivative lanes; comparisons
        return bool.
        """
        if isinstance(other, np.ndarray):
            return NotImplemented
        other = self.lift(other)
        return Tangent(self.value+other.value, self.slope+other.slope)

    __radd__ = __add__

    def __neg__(self):
        """
        Negate the primal value and its forward derivatives.

        SYNOPSIS:
            __neg__(...)

        INPUT ARGUMENTS:
        self : Current forward-mode scalar, including its derivative lanes.

        OUTPUT ARGUMENTS:
        Numerical result with propagated derivative lanes; comparisons
        return bool.
        """
        return Tangent(-self.value, -self.slope)

    def __sub__(self, other):
        """
        Subtract operands and propagate forward derivatives.

        SYNOPSIS:
            __sub__(...)

        INPUT ARGUMENTS:
        self : Current forward-mode scalar, including its derivative lanes.
        other : Second scalar/array operand, promoted when required.

        OUTPUT ARGUMENTS:
        Numerical result with propagated derivative lanes; comparisons
        return bool.
        """
        if isinstance(other, np.ndarray):
            return NotImplemented
        return self + -self.lift(other)

    def __rsub__(self, other):
        """
        Evaluate the reversed arithmetic operation with forward derivatives.

        SYNOPSIS:
            __rsub__(...)

        INPUT ARGUMENTS:
        self : Current forward-mode scalar, including its derivative lanes.
        other : Second scalar/array operand, promoted when required.

        OUTPUT ARGUMENTS:
        Numerical result with propagated derivative lanes; comparisons
        return bool.
        """
        return -self + other

    def __mul__(self, other):
        """
        Multiply operands using the forward product rule.

        SYNOPSIS:
            __mul__(...)

        INPUT ARGUMENTS:
        self : Current forward-mode scalar, including its derivative lanes.
        other : Second scalar/array operand, promoted when required.

        OUTPUT ARGUMENTS:
        Numerical result with propagated derivative lanes; comparisons
        return bool.
        """
        if isinstance(other, np.ndarray):
            return NotImplemented
        other = self.lift(other)
        return Tangent(self.value*other.value,
                       self.slope*other.value+self.value*other.slope)

    __rmul__ = __mul__

    def __truediv__(self, other):
        """
        Divide operands using the forward quotient rule.

        SYNOPSIS:
            __truediv__(...)

        INPUT ARGUMENTS:
        self : Current forward-mode scalar, including its derivative lanes.
        other : Second scalar/array operand, promoted when required.

        OUTPUT ARGUMENTS:
        Numerical result with propagated derivative lanes; comparisons
        return bool.
        """
        if isinstance(other, np.ndarray):
            return NotImplemented
        other = self.lift(other)
        value = self.value/other.value
        return Tangent(value, (self.slope-value*other.slope)/other.value)

    def __rtruediv__(self, other):
        """
        Evaluate the reversed arithmetic operation with forward derivatives.

        SYNOPSIS:
            __rtruediv__(...)

        INPUT ARGUMENTS:
        self : Current forward-mode scalar, including its derivative lanes.
        other : Second scalar/array operand, promoted when required.

        OUTPUT ARGUMENTS:
        Numerical result with propagated derivative lanes; comparisons
        return bool.
        """
        if isinstance(other, np.ndarray):
            return NotImplemented
        return self.lift(other)/self

    def __pow__(self, other):
        """
        Raise an operand to a power with forward derivatives.

        SYNOPSIS:
            __pow__(...)

        INPUT ARGUMENTS:
        self : Current forward-mode scalar, including its derivative lanes.
        other : Second scalar/array operand, promoted when required.

        OUTPUT ARGUMENTS:
        Numerical result with propagated derivative lanes; comparisons
        return bool.
        """
        if isinstance(other, np.ndarray):
            return NotImplemented
        other = self.lift(other)
        if self.value == 0 and other.value > 1:
            return Tangent(0., self.slope*0.)
        value = self.value**other.value
        slope = other.value*self.value**(other.value-1)*self.slope
        if np.any(other.slope):
            slope += value*np.log(self.value)*other.slope
        return Tangent(value, slope)

    def __rpow__(self, other):
        """
        Evaluate the reversed arithmetic operation with forward derivatives.

        SYNOPSIS:
            __rpow__(...)

        INPUT ARGUMENTS:
        self : Current forward-mode scalar, including its derivative lanes.
        other : Second scalar/array operand, promoted when required.

        OUTPUT ARGUMENTS:
        Numerical result with propagated derivative lanes; comparisons
        return bool.
        """
        return self.lift(other)**self

    def __lt__(self, other):
        """
        Compare primal values without differentiating the branch.

        SYNOPSIS:
            __lt__(...)

        INPUT ARGUMENTS:
        self : Current forward-mode scalar, including its derivative lanes.
        other : Second scalar/array operand, promoted when required.

        OUTPUT ARGUMENTS:
        Numerical result with propagated derivative lanes; comparisons
        return bool.
        """
        return self.value < primal(other)

    def __le__(self, other):
        """
        Compare primal values without differentiating the branch.

        SYNOPSIS:
            __le__(...)

        INPUT ARGUMENTS:
        self : Current forward-mode scalar, including its derivative lanes.
        other : Second scalar/array operand, promoted when required.

        OUTPUT ARGUMENTS:
        Numerical result with propagated derivative lanes; comparisons
        return bool.
        """
        return self.value <= primal(other)

    def __gt__(self, other):
        """
        Compare primal values without differentiating the branch.

        SYNOPSIS:
            __gt__(...)

        INPUT ARGUMENTS:
        self : Current forward-mode scalar, including its derivative lanes.
        other : Second scalar/array operand, promoted when required.

        OUTPUT ARGUMENTS:
        Numerical result with propagated derivative lanes; comparisons
        return bool.
        """
        return self.value > primal(other)

    def __ge__(self, other):
        """
        Compare primal values without differentiating the branch.

        SYNOPSIS:
            __ge__(...)

        INPUT ARGUMENTS:
        self : Current forward-mode scalar, including its derivative lanes.
        other : Second scalar/array operand, promoted when required.

        OUTPUT ARGUMENTS:
        Numerical result with propagated derivative lanes; comparisons
        return bool.
        """
        return self.value >= primal(other)

    def exp(self):
        """
        Evaluate exp and propagate forward derivatives.

        SYNOPSIS:
            exp(...)

        INPUT ARGUMENTS:
        self : Current forward-mode scalar, including its derivative lanes.

        OUTPUT ARGUMENTS:
        Numerical result with propagated derivative lanes; comparisons
        return bool.
        """
        value = np.exp(self.value)
        return Tangent(value, value*self.slope)

    def expm1(self):
        """
        Evaluate expm1 and propagate forward derivatives.

        SYNOPSIS:
            expm1(...)

        INPUT ARGUMENTS:
        self : Current forward-mode scalar, including its derivative lanes.

        OUTPUT ARGUMENTS:
        Numerical result with propagated derivative lanes; comparisons
        return bool.
        """
        return Tangent(np.expm1(self.value), np.exp(self.value)*self.slope)

    def log(self):
        """
        Evaluate log and propagate forward derivatives.

        SYNOPSIS:
            log(...)

        INPUT ARGUMENTS:
        self : Current forward-mode scalar, including its derivative lanes.

        OUTPUT ARGUMENTS:
        Numerical result with propagated derivative lanes; comparisons
        return bool.
        """
        return Tangent(np.log(self.value), self.slope/self.value)

    def log1p(self):
        """
        Evaluate log1p and propagate forward derivatives.

        SYNOPSIS:
            log1p(...)

        INPUT ARGUMENTS:
        self : Current forward-mode scalar, including its derivative lanes.

        OUTPUT ARGUMENTS:
        Numerical result with propagated derivative lanes; comparisons
        return bool.
        """
        return Tangent(np.log1p(self.value), self.slope/(1+self.value))

    def tanh(self):
        """
        Evaluate tanh and propagate forward derivatives.

        SYNOPSIS:
            tanh(...)

        INPUT ARGUMENTS:
        self : Current forward-mode scalar, including its derivative lanes.

        OUTPUT ARGUMENTS:
        Numerical result with propagated derivative lanes; comparisons
        return bool.
        """
        value = np.tanh(self.value)
        return Tangent(value, (1-value*value)*self.slope)
