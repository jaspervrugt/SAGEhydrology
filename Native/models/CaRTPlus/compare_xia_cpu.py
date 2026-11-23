# ========================================================================
# SAGEhydrology | Model source
# File: Native/models/CaRTPlus/compare_xia_cpu.py
#
# CaRT / CaRT+ numerical model and supporting routines.
#
# SAGE adaptations: Jasper A. Vrugt, University of California, Irvine.
# Upstream model: Xia (2026). Preserve the adjacent license notices.
# ========================================================================

"""
Compare Xia's unmodified daily forward code with SAGE on a CPU.

Run with a local checkout of xiaxilin/cartplus and installed CAMELS-US
Project 6 data. The CPU shim evaluates the upstream Triton kernel bodies
with NumPy. It does not emulate GPU rounding or run training.
"""

import argparse
import json
from pathlib import Path
import sys
from types import ModuleType, SimpleNamespace

import numpy as np
import torch

from cartplus_cpu import simulate as sage_plus
from cartplus_inputs import PARAMETER_MAP, elevation_bands, \
    physical_parameters


class Pointer:
    'Represent a flat kernel buffer' ' and its lane offsets.'

    def __init__(self, array, offset=0):
        """
        Construct a scalar with primal value and derivative lanes.

        SYNOPSIS:
            __init__(...)

        INPUT ARGUMENTS:
        self : Current forward-mode scalar, including its derivative lanes.
        array : Contiguous numerical array, or None for a null native
            pointer.
        offset : Offset of a routing-state block in the full state vector.

        OUTPUT ARGUMENTS:
        None. Updates the supplied object or checks the stated conditions.
        """
        self.array = array
        self.offset = offset

    def __add__(self, offset):
        """
        Add operands and propagate forward derivatives.

        SYNOPSIS:
            __add__(...)

        INPUT ARGUMENTS:
        self : Current forward-mode scalar, including its derivative lanes.
        offset : Offset of a routing-state block in the full state vector.

        OUTPUT ARGUMENTS:
        Numerical result with propagated derivative lanes; comparisons
        return bool.
        """
        return Pointer(self.array, self.offset + offset)

    def __sub__(self, offset):
        """
        Subtract operands and propagate forward derivatives.

        SYNOPSIS:
            __sub__(...)

        INPUT ARGUMENTS:
        self : Current forward-mode scalar, including its derivative lanes.
        offset : Offset of a routing-state block in the full state vector.

        OUTPUT ARGUMENTS:
        Numerical result with propagated derivative lanes; comparisons
        return bool.
        """
        return Pointer(self.array, self.offset - offset)


def load(pointer, mask=True, other=0):
    """
    Load active lanes and fill masked lanes with the upstream default.

    SYNOPSIS:
        load(...)

    INPUT ARGUMENTS:
    pointer : Pointer used by this numerical operator.
    mask : Mask used by this numerical operator.
    other : Second scalar/array operand, promoted when required.

    OUTPUT ARGUMENTS:
    Computed fluxes, stores, or numerical values described in the notes.
    """
    mask = np.broadcast_to(mask, np.shape(pointer.offset))
    return np.where(mask, pointer.array[np.where(mask, pointer.offset, 0)],
                    np.asarray(other, dtype=pointer.array.dtype))


def store(pointer, value, mask=True):
    """
    Write active lanes into the original PyTorch buffer.

    SYNOPSIS:
        store(...)

    INPUT ARGUMENTS:
    pointer : Pointer used by this numerical operator.
    value : Primal scalar or numerical operand.
    mask : Mask used by this numerical operator.

    OUTPUT ARGUMENTS:
    None. Updates the supplied object or checks the stated conditions.
    """
    mask = np.broadcast_to(mask, np.shape(pointer.offset))
    pointer.array[np.asarray(pointer.offset)[mask]] = (
        np.broadcast_to(value, mask.shape)[mask]
    )


def install_cpu_shim():
    """
    Provide only the Triton operations used by the upstream kernels.

    SYNOPSIS:
        install_cpu_shim(...)

    INPUT ARGUMENTS:
    None.

    OUTPUT ARGUMENTS:
    None. Updates the supplied object or checks the stated conditions.
    """
    triton = ModuleType("triton")
    language = ModuleType("triton.language")
    block_id = 0

    class Kernel:
        def __init__(self, function):
            """
            Construct a scalar with primal value and derivative lanes.

            SYNOPSIS:
                __init__(...)

            INPUT ARGUMENTS:
            self : Current forward-mode scalar, including its derivative
                lanes.
            function : Numerical callback or native callable.

            OUTPUT ARGUMENTS:
            None. Updates the supplied object or checks the stated
            conditions.
            """
            self.function = function

        def __call__(self, *args, **kwargs):
            """
            Invoke the comparison shim with forwarded model arguments.

            SYNOPSIS:
                __call__(...)

            INPUT ARGUMENTS:
            self : Current forward-mode scalar, including its derivative
                lanes.
            args : Arguments forwarded to the selected native/model
                callback.
            kwargs : Additional keyword arguments forwarded to the solver.

            OUTPUT ARGUMENTS:
            Numerical result with propagated derivative lanes; comparisons
            return bool.
            """
            return self.function(*args, **kwargs)

        def __getitem__(self, grid):
            """
            Bind the requested launch grid for the comparison shim.

            SYNOPSIS:
                __getitem__(...)

            INPUT ARGUMENTS:
            self : Current forward-mode scalar, including its derivative
                lanes.
            grid : Optional numerical verification grid.

            OUTPUT ARGUMENTS:
            Numerical result with propagated derivative lanes; comparisons
            return bool.
            """
            def launch(*args, **kwargs):
                """
                Launch the selected native comparison or test executable.

                SYNOPSIS:
                    launch(...)

                INPUT ARGUMENTS:
                args : Arguments forwarded to the selected native/model
                    callback.
                kwargs : Additional keyword arguments forwarded to the
                    solver.

                OUTPUT ARGUMENTS:
                None. Updates the supplied object or checks the stated
                conditions.
                """
                nonlocal block_id
                args = [
                    Pointer(x.detach().numpy().reshape(-1))
                    if isinstance(x, torch.Tensor) else x for x in args
                ]
                for block_id in range(grid[0]):
                    self.function(*args, **kwargs)
            return launch

    triton.jit = Kernel
    triton.cdiv = lambda a, b: (a+b-1)//b
    for name in ("abs", "arange", "exp", "log", "maximum", "minimum",
                 "sqrt", "sum", "where", "zeros", "zeros_like", "float32"):
        setattr(language, name, getattr(np, name))
    language.sigmoid = lambda x: 1 / (1 + np.exp(-x))
    language.static_range = range
    language.constexpr = object()
    language.program_id = lambda axis: block_id
    language.load = load
    language.store = store
    triton.language = language
    sys.modules.update({"triton": triton, "triton.language": language})


def compare(upstream, data, gauges, days, jacobian_rows=0):
    """
    Compare day-by-day discharge at upstream initial parameter values.

    SYNOPSIS:
        compare(...)

    INPUT ARGUMENTS:
    upstream : Upstream used by this numerical operator.
    data : Basin forcing, elevation bands, and model settings.
    gauges : Gauges used by this numerical operator.
    days : Residence time (days).
    jacobian_rows : Number of observation rows selected for Jacobian
        comparison.

    OUTPUT ARGUMENTS:
    None. Updates the supplied object or checks the stated conditions.
    """
    sys.path.insert(0, str(upstream / "src"))
    install_cpu_shim()
    from cart.config import PRESETS
    from cart.model import CaRT
    from cart.parameters import static_table
    import importlib.util

    path = Path(__file__).resolve().parents[1] / "cart/cart_python.py"
    spec = importlib.util.spec_from_file_location("sage_cart", path)
    baseline = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(baseline)
    root = Path(__file__).resolve().parent
    params = json.loads((root / "parameters.json").read_text())
    torch.set_num_threads(1)

    with np.load(data / "camels_us_daymet671_v4full.npz") as archive:
        ids = list(archive["ids"].astype(str))
        rain = archive["P"][:days]
        pet = archive["PET"][:days]
        temp = .5 * (archive["Tmax"][:days] + archive["Tmin"][:days])
    mask = np.arange(len(rain)) >= 365
    if not mask.any():
        raise ValueError('Select more than 365 days for ' \
            'the training mask.')
    for gauge in gauges:
        col = ids.index(gauge)
        z, area = elevation_bands(data, gauge)
        f = np.column_stack((rain[:, col], temp[:, col], pet[:, col], mask))
        dataset = SimpleNamespace(
            P=f[:, :1], Tm=f[:, 1:2]+273.15, E=f[:, 2:3],
            static=np.zeros((1, 1), np.float32), band_H=z[None],
            band_A=area[None], mask=lambda _: mask,
        )
        for name in ("cart", "cartplus"):
            cfg = PRESETS[name]
            static = {k: v[2] for k, v in static_table(cfg).items()}
            law = {k: v[2] for k, v in params["law"].items()}
            shared = params["shared"]
            model = CaRT(dataset, cfg, "cpu")
            law_parameters = {
                k: torch.tensor([[v]], dtype=torch.float32,
                                requires_grad=bool(jacobian_rows))
                for k, v in law.items()
            } if name == "cartplus" else {}
            model.net = SimpleNamespace(
                phg={k: torch.tensor(shared[k], dtype=torch.float32)
                     for k in ("Tv", "lsv", "ltv")},
                lawg={k: torch.tensor(shared[k], dtype=torch.float32)
                      for k in ("TF", "lsF", "thc", "lsth")},
                law=lambda _: law_parameters.copy(),
            )
            model.net = FixedNetwork(model.net, static, bool(jacobian_rows))
            with torch.set_grad_enabled(bool(jacobian_rows)), np.errstate( \
                over="ignore", invalid="ignore",
                                               divide="ignore"):
                reference = model.simulate()[0].ravel()
                q = reference.detach().numpy()
            if name == "cart":
                actual, native = baseline.simulate(
                    f, z, area, list(static.values()), mcode=4, jacobian= \
                        bool(jacobian_rows)
                )
            else:
                p = physical_parameters(static, law, shared)
                actual, native = sage_plus(root / "cartplus_cpu.dll", f, z \
                    , area, p,
                                           jacobian=bool(jacobian_rows))
                actual = actual[:, 0]
            if not (np.isfinite(actual).all() and np.isfinite(q).all()):
                raise ValueError(f"Nonfinite discharge for {gauge} {name}.")
            delta = actual - q
            print(f"{gauge} {name}: n={len(q)}, "
                  f"max |difference|={np.max(np.abs(delta)):.6g}, "
                  f"RMSE={np.sqrt(np.mean(delta**2)):.6g} mm/day",
                  flush=True)

            if jacobian_rows:
                leaves = model.net.static | law_parameters
                if name == "cartplus":
                    mapping = {
                        key: (i, p[i] if scale is not None else 1)
                        for i, (group, key, scale) in enumerate( \
                            PARAMETER_MAP)
                        if group != 2
                    }
                    native = np.column_stack([
                        native[:, mapping[key][0]] * mapping[key][1]
                        for key in leaves
                    ])
                sampled_days = np.unique(np.r_[
                    np.linspace(0, len(q)-1, jacobian_rows, dtype=int),
                    np.argmax(q), 365,
                ])
                check_jacobian(reference, native, leaves, sampled_days)


def check_jacobian(reference, native, leaves, days):
    """
    Evaluate upstream adjoint rows and compare in calibration space.

    SYNOPSIS:
        check_jacobian(...)

    INPUT ARGUMENTS:
    reference : Pilot mean recharge or reference storage used by the caller.
    native : Native used by this numerical operator.
    leaves : Leaves used by this numerical operator.
    days : Residence time (days).

    OUTPUT ARGUMENTS:
    None. Updates the supplied object or checks the stated conditions.
    """
    rows = []
    for i in days:
        gradient = torch.autograd.grad(
            reference[i], tuple(leaves.values()), retain_graph=True,
            allow_unused=True,
        )
        rows.append([float(g.item()) if g is not None else 0 for g in \
            gradient])
    upstream = np.asarray(rows)
    actual = native[days]
    if not (np.isfinite(upstream).all() and np.isfinite(actual).all()):
        raise ValueError("Nonfinite discharge Jacobian.")
    error = actual - upstream
    i, j = np.unravel_index(np.abs(error).argmax(), error.shape)
    print(f"  Jacobian {error.shape}: max abs={abs(error[i,j]):.6g}, "
          f"RMSE={np.sqrt(np.mean(error**2)):.6g}, "
          f"""relative L2={np.linalg.norm(error)/
              max(np.linalg.norm(upstream),1e-30):.6g}; """
          f"worst day={days[i]+1}, parameter={list(leaves)[j]}, "
          f"upstream={upstream[i,j]:.6g}, C++={actual[i,j]:.6g}", flush= \
              True)


class FixedNetwork:
    'Supply identical fixed catchme' 'nt parameters to the upstream ' \
        'model.'

    def __init__(self, network, static, gradients):
        """
        Construct a scalar with primal value and derivative lanes.

        SYNOPSIS:
            __init__(...)

        INPUT ARGUMENTS:
        self : Current forward-mode scalar, including its derivative lanes.
        network : Network used by this numerical operator.
        static : Mapping from static parameter names to calibration values.
        gradients : Gradients used by this numerical operator.

        OUTPUT ARGUMENTS:
        None. Updates the supplied object or checks the stated conditions.
        """
        self.__dict__.update(vars(network))
        self.static = {k: torch.tensor([v], dtype=torch.float32, \
            requires_grad=gradients)
                       for k, v in static.items()}

    def __call__(self, attributes):
        """
        Invoke the comparison shim with forwarded model arguments.

        SYNOPSIS:
            __call__(...)

        INPUT ARGUMENTS:
        self : Current forward-mode scalar, including its derivative lanes.
        attributes : Regional basin attribute records.

        OUTPUT ARGUMENTS:
        Numerical result with propagated derivative lanes; comparisons
        return bool.
        """
        return attributes, self.static.copy()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("upstream", type=Path)
    parser.add_argument("data", type=Path)
    parser.add_argument("gauges", nargs="+")
    parser.add_argument("--days", type=int, default=3650)
    parser.add_argument(
        "--jacobian-rows", type=int, default=0,
        help=('Number of evenly spaced Jacobi' \
            'an rows to compare. A positive' ' '
              'value also includes peak flow ' \
                  'and the first training day.'),
    )
    args = parser.parse_args()
    compare(args.upstream, args.data, args.gauges, args.days, args. \
        jacobian_rows)
