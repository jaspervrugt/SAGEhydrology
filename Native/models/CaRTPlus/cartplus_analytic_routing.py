# ========================================================================
# SAGEhydrology | Model source
# File: Native/models/CaRTPlus/cartplus_analytic_routing.py
#
# CaRT / CaRT+ numerical model and supporting routines.
#
# SAGE adaptations: Jasper A. Vrugt, University of California, Irvine.
# Upstream model: Xia (2026). Preserve the adjacent license notices.
# ========================================================================

"""
Analytic propagation of one continuous CaRT+ routing cascade.
"""

import numpy as np
from scipy.linalg import expm, expm_frechet


def routing_step(exit_fractions, refinement, routing_days, step_days, *,
                 differentiate=False):
    """
    Evaluate a cascade transition and its optional derivative.

    SYNOPSIS:
        routing_step(...)

    INPUT ARGUMENTS:
    exit_fractions : Routing-stage exit fractions (-).
    refinement : Positive number of routing stages per X4 time base.
    routing_days : Routing time base (days).
    step_days : Requested propagation duration (days).
    differentiate : Whether to differentiate the cascade transition.

    OUTPUT ARGUMENTS:
    Transition matrix and optional derivative with respect to the routing
    time base.

    NOTES:
    Apply the map to [initial stores, start rate, end rate]. The outputs are
    final stores followed by integrated release. Rates are in mm/day, stores
    and release in mm. The optional second map differentiates the first with
    respect to log(routing_days).
    """
    exits = np.asarray(exit_fractions, dtype=float)
    if (exits.ndim != 1 or not exits.size or not np.isfinite(exits).all()
            or np.any((exits < 0) | (exits > 1)) or exits[-1] != 1):
        raise ValueError('Expected cascade exit fraction'
            's in [0, 1], ending in one.')
    if (not np.isfinite(routing_days) or routing_days <= 0
            or not np.isfinite(step_days) or step_days < 0):
        raise ValueError('Expected positive routing time'
            ' and nonnegative step duration' '.')

    stages = exits.size
    if refinement < 1 or stages not in (refinement, 2*refinement):
        raise ValueError('Expected refinement matching t'
            'he short or long cascade.')
    # Work in normalized step time. The last two coordinates carry
    # the inlet rate and its constant derivative in normalized time.
    rate = step_days * refinement / routing_days
    if not np.isfinite(rate):
        raise ValueError('Routing rate times step durati'
            'on must be finite.')
    generator = np.zeros((stages + 3, stages + 3))
    i = np.arange(stages)
    generator[i, i] = -rate
    generator[i[1:], i[:-1]] = rate * (1 - exits[:-1])
    generator[stages, i] = rate * exits
    log_time_direction = -generator.copy()
    generator[0, stages + 1] = step_days
    generator[stages + 1, stages + 2] = 1

    def endpoint_map(matrix):
        """
        Select the endpoint storage and discharge transition rows.

        SYNOPSIS:
            endpoint_map(...)

        INPUT ARGUMENTS:
        matrix : State-transition matrix for the routing cascade.

        OUTPUT ARGUMENTS:
        Transition rows for endpoint stores and cumulative outflow.
        """
        result = matrix[:stages + 1].copy()
        result[:, stages + 1] -= result[:, stages + 2]
        return np.delete(result, stages, axis=1)

    if differentiate:
        transition, derivative = expm_frechet(generator, log_time_direction)
        return endpoint_map(transition), endpoint_map(derivative)
    return endpoint_map(expm(generator)), None
