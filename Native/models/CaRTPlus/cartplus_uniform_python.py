# ========================================================================
# SAGEhydrology | Model source
# File: Native/models/CaRTPlus/cartplus_uniform_python.py
#
# CaRT / CaRT+ numerical model and supporting routines.
#
# SAGE adaptations: Jasper A. Vrugt, University of California, Irvine.
# Upstream model: Xia (2026). Preserve the adjacent license notices.
# ========================================================================

"""
Uniform-input CaRT+ port of Xia f986af2 (MIT), Python reference backend.
"""
from pathlib import Path
from importlib.util import spec_from_file_location, module_from_spec
from math import ceil
import numpy as np
_spec=spec_from_file_location('_cartplus_legacy_scalar',Path(__file__).
    with_name('cartplus_python.py'))
_legacy=module_from_spec(_spec)
_spec.loader.exec_module(_legacy)
Tangent,primal=_legacy.Tangent,_legacy.primal
logistic,memory,zone_step=_legacy.logistic,_legacy.memory,_legacy.zone_step

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
        sub_demand = sublimation_factor*pet[i]
        melt_demand = np.where(cold, 0, melt_rate*dt*warmth+pet_melt*pet[i])
        total = sub_demand+melt_demand
        fraction = np.array([min(1., snow[k]/total[k]) if total[k]>0 else 1.
                             for k in range(5)], dtype=dtype)
        sub = sub_demand*fraction
        melt = melt_demand*fraction
        snow -= sub+melt
        rain[i] = area @ (melt + (~cold)*precip[i])
        sublimation[i] = area @ sub
        snow_extra[i] = (area @ cold) * precip[i] * (snowfall_factor - 1)
        snow_change[i] = area @ snow - before

    drainage = dt / surface_days
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
                soil_shape, dt, True, symmetric=True,
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
        soil_before, wet_before = soil, wet
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
        saturation[i] = np.clip(.5*(soil_before+soil)/soil_cap, 0, 1)
        wetness[i] = .5*(wet_before+wet)/wet_cap
        change[i] += upland*soil + zone_fraction*near - before

    depression = .5*depression_cap
    if fixed is None:
        for i in range(spin):
            depression = pareto(
                depression, fast[i], max(demand[i]-et[i], 0),
                depression_cap, depression_perc, depression_shape,
                dt, False, symmetric=True,
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
        groundwater, base_q[i], uptake = groundwater_step(
            groundwater, recharge[i], max(demand[i]-et[i],0), recession,
            groundwater_uptake, groundwater_start)
        change[i] += groundwater-before
        et[i] += uptake
        signal = wetness[i] - connected_center
        signal += soil_connection*saturation[i]
        signal += groundwater_connection*(.5*(before+groundwater)/
            groundwater_start-1)
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
        quick_step(store, routed[i], bypass[i], release, exchange,
        quick_cap, connectivity[i], quick_exponent) : computed result in the
        operand units.
        """
        return quick_step(store,routed[i],bypass[i],release,exchange,
                          quick_cap,connectivity[i],quick_exponent)

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


def release(water, coefficient):
    """
    Evaluate stable nonlinear release over a split step.

    SYNOPSIS:
        release(...)

    INPUT ARGUMENTS:
    water : Available water or store depth (mm).
    coefficient : Release-law coefficient in the current interval time base.

    OUTPUT ARGUMENTS:
    Water released by the nonlinear split-step operator (mm).
    """
    return -water*np.expm1(-.25*np.log1p(coefficient*water**4))


def surface_series(store, rain, pet, capacity, b):
    """
    Advance surface storage, overflow, drainage, and evaporation.

    SYNOPSIS:
        surface_series(...)

    INPUT ARGUMENTS:
    store : Current water store (mm).
    rain : Liquid input depth/rate on this operator time base.
    pet : Potential evaporation depth/rate on this operator time base.
    capacity : Store or infiltration capacity in this operator.
    b : Second numerical operand or rate coefficient.

    OUTPUT ARGUMENTS:
    Final store and interval runoff, residual demand, evaporation, and
    storage-change series.
    """
    dtype=object if any(isinstance(x,Tangent) for x in (store,capacity,b)) \
        else float
    runoff,demand,et,change=np.zeros((4,len(rain)),dtype=dtype)
    for i,(p,e) in enumerate(zip(rain,pet)):
        before=store; net=p-e; eb=-np.expm1(-b)
        free=store*(1-eb)+net*eb/b
        over=net>=0 and free>capacity
        spill=max(net-b*capacity,1e-12)
        tc=min(np.log1p(b*max(capacity-store,0)/spill)/b,1.)
        t0=np.log1p(b*store/max(-net,1e-12))/b
        empty=net<0 and t0<1
        time=tc if over else t0 if empty else 1.
        x=b*time
        psi=x*(.5-x*(1/6-x*(1/24-x/120))) if x<.05 else 1+np.expm1(-x)/x
        drainage=-store*np.expm1(-x)+net*time*psi+(b*capacity*(1-tc) if
            over else 0)
        overflow=spill*(1-tc) if over else 0
        et[i]=e*t0+p*(1-t0) if empty else e
        store=capacity if over else 0. if empty else free
        runoff[i]=overflow+drainage; demand[i]=max(e-et[i],0)
        change[i]=store-before
    return store,runoff,demand,et,change


def infiltrate(rain,pet,mean,capacity,shape,floor,resolution):
    """
    Partition liquid input into infiltration and excess runoff.

    SYNOPSIS:
        infiltrate(...)

    INPUT ARGUMENTS:
    rain : Liquid input depth/rate on this operator time base.
    pet : Potential evaporation depth/rate on this operator time base.
    mean : Reference rainfall/elevation mean used by this operator.
    capacity : Store or infiltration capacity in this operator.
    shape : Dimensionless storage/infiltration shape parameter.
    floor : Minimum infiltration capacity in daily-rate units (mm/day).
    resolution : Forcing intervals per day: 1 daily, 24 hourly, 96
        quarter-hourly.

    OUTPUT ARGUMENTS:
    Infiltrating liquid input and infiltration-excess runoff.
    """
    dt=1/resolution; decay=np.exp(-dt/3)
    history=memory(rain,decay,36*resolution)
    previous=np.r_[0.,history[:-1]]+decay**np.arange(len(rain))*mean
    wet=rain+(previous-rain)*(1-decay)/(dt/3)
    cap=dt*np.maximum(capacity*np.exp(-shape*wet/mean),floor)
    available=np.maximum(rain-pet,0); x=available/cap
    excess=available-cap*np.where(x<1e-3,x-.5*x*x,-np.expm1(-x))
    return rain-excess,excess


def pareto(store,rain,demand,cap,percolation,shape,dt,net,symmetric=False):
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
    symmetric : Whether to use the symmetric two-time-base hydrograph.

    OUTPUT ARGUMENTS:
    Updated storage, runoff, recharge, and evaporation for the interval; the
    continuous helper returns flux rates.
    """
    water=max(rain-demand,0) if net else rain
    potential=max(demand-rain,0) if net else demand
    interception=min(rain,demand) if net else 0
    coefficient=dt*percolation/cap**4
    recharge=0.
    if symmetric:
        recharge=release(store,.5*coefficient); store-=recharge
    def wet(s):
        """
        Retain liquid input according to the production deficit.

        SYNOPSIS:
            wet(...)

        INPUT ARGUMENTS:
        s : S used by this numerical operator.

        OUTPUT ARGUMENTS:
        min(cap * deficit if ratio >= 1 else -cap * deficit *
        np.expm1(exponent * np.log1p(-ratio)), water) : computed result in
        the operand units.
        """
        deficit=max(1-np.clip(s/cap,0,1),1e-7); exponent=shape+1
        ratio=water/(exponent*cap*np.exp(np.log(deficit)/exponent))
        return min(cap*deficit if ratio>=1 else -cap*deficit*np.expm1(
            exponent*np.log1p(-ratio)),water)
    def dry(s,depth):
        """
        Remove soil evaporation for the specified demand.

        SYNOPSIS:
            dry(...)

        INPUT ARGUMENTS:
        s : S used by this numerical operator.
        depth : Depth used by this numerical operator.

        OUTPUT ARGUMENTS:
        min(cap * (1 - d2) * response / (1 + d2 / deficit * response), s) :
        computed result in the operand units.
        """
        deficit=max(1-np.clip(s/cap,0,1),1e-7); d2=deficit**2
        response=np.tanh(depth/cap)
        return min(cap*(1-d2)*response/(1+d2/deficit*response),s)
    if symmetric and not net:
        evap=dry(store,.5*potential); store-=evap
        retained=wet(store); store+=retained
        add=dry(store,.5*potential); store-=add; evap+=add
    else:
        retained=wet(store); evap=dry(store,potential); store+=retained-evap
    g=release(store,coefficient*(.5 if symmetric else 1))
    return store-g,water-retained,recharge+g,evap+interception


def groundwater_step(store,recharge,demand,rate,lam,reference):
    """
    Advance groundwater storage, discharge, and uptake.

    SYNOPSIS:
        groundwater_step(...)

    INPUT ARGUMENTS:
    store : Current water store (mm).
    recharge : Recharge supplied to the groundwater store.
    demand : Evaporation demand on this operator time base.
    rate : Recession or cascade release rate (1/day).
    lam : Groundwater evaporation-uptake coefficient.
    reference : Pilot mean recharge or reference storage used by the caller.

    OUTPUT ARGUMENTS:
    Updated groundwater store, discharge, and uptake depths.
    """
    half=.5*lam
    def uptake(w):
        """
        Evaluate groundwater uptake during a split step.

        SYNOPSIS:
            uptake(...)

        INPUT ARGUMENTS:
        w : Integration workspace or scalar operand in this helper.

        OUTPUT ARGUMENTS:
        -w * np.expm1(-half * demand / (w * np.exp(-0.5 * half * demand / (w
        + reference)) + reference)) : computed result in the operand units.
        """
        return -w*np.expm1(-half*demand/(w*np.exp(-.5*half*demand/(w+
            reference))+reference))
    e1=uptake(store); store-=e1
    g=max(recharge,0); y=g*rate
    om=y*(1/3-y*(2/15-y*(17/315-y*62/2835))) if y<.04 else 1-np.tanh(np.
        sqrt(y))/np.sqrt(y)
    tc=1-om; den=1+rate*store*tc
    q=(rate*store*store*tc+y*store*tc+g*om)/den
    store=(store+g*tc)/den
    e2=uptake(store)
    return store-e2,q,e1+e2


def quick_step(store,u1,u2,p1,p2,cap,c0,c1):
    """
    Advance quick routing, exchange, and bypass discharge.

    SYNOPSIS:
        quick_step(...)

    INPUT ARGUMENTS:
    store : Current water store (mm).
    u1 : Input depth to the quick-routing branch.
    u2 : Input depth to the direct/bypass branch.
    p1 : Quick-store release coefficient for the split operator.
    p2 : Quick-store exchange coefficient for the split operator.
    cap : Store capacity (mm).
    c0 : Connectivity multiplier in the quick-release law.
    c1 : Storage-dependent quick-release exponent coefficient.

    OUTPUT ARGUMENTS:
    Updated quick store, outlet discharge, and exchange depth.
    """
    def lrel(w,half):
        """
        Evaluate nonlinear quick-store release during a split step.

        SYNOPSIS:
            lrel(...)

        INPUT ARGUMENTS:
        w : Integration workspace or scalar operand in this helper.
        half : Half-step duration or split-release multiplier.

        OUTPUT ARGUMENTS:
        -w * np.expm1(-np.log1p(nu * a * w ** 4) / nu) : computed result in
        the operand units.
        """
        th=np.tanh((c0+c1*np.log(max(w,1e-3)/cap))/6)
        a=half*np.exp(6*th); nu=4+c1*(1-th*th)
        return -w*np.expm1(-np.log1p(nu*a*w**4)/nu)
    def chain(w,ns):
        """
        Combine scalar values with chain-rule derivative lanes.

        SYNOPSIS:
            chain(...)

        INPUT ARGUMENTS:
        w : Integration workspace or scalar operand in this helper.
        ns : Number of hydrological state entries.

        OUTPUT ARGUMENTS:
        (w, q, ex) : computed result in the operand units.
        """
        q=ex=0.; a,b,j=u1/ns,u2/ns,p2/ns; half=.5*p1/ns
        for _ in range(ns):
            q1=lrel(w,half); w-=q1
            e1=j*(max(w,0)/30)**3.5
            flux=j*(max(w+.5*(a+.9*e1),0)/30)**3.5
            v=max(w+a+.9*flux,0); direct=max(b+.1*flux,0)
            ex+=v-w-a+direct-b
            q2=lrel(v,half); w=v-q2; q+=q1+q2+direct
        return w,q,ex
    a,b=chain(store,4),chain(store,2)
    w,q,ex=[(4*x-y)/3 for x,y in zip(a,b)]
    correction=max(-w,0)
    return w+correction,q,ex+correction


def hydrograph(inflow,days,dt,symmetric):
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
    def integral(t):
        """
        Evaluate the time integral of a unit-hydrograph S-curve.

        SYNOPSIS:
            integral(...)

        INPUT ARGUMENTS:
        t : T used by this numerical operator.

        OUTPUT ARGUMENTS:
        Computed fluxes, stores, or numerical values described in the notes.
        """
        if t<=0:return 0.
        u=t/days
        if not symmetric:
            return t-days*(1-1/3.5) if u>1 else days*u**3.5/3.5
        if u<=1:return days*u**3.5/7
        return t-days+days*(2-u)**3.5/7 if u<2 else t-days
    length=ceil(primal(days)*(2 if symmetric else 1)/dt)+1
    cumulative=[(integral((i+1)*dt)-integral(i*dt))/dt for i in range(
        length)]
    return np.convolve(inflow,np.diff([0.,*cumulative]))[:len(inflow)]

