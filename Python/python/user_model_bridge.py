# ==========================================================================
# SITEhydrology | Python source
# File: python/user_model_bridge.py
#
# MATLAB-free execution bridge shared with SAGE Python user models.
#
# Copyright (c) 2026 Jasper A. Vrugt
# University of California, Irvine | jasper@uci.edu
# ==========================================================================

"""MATLAB-free execution bridge shared with SAGE Python user models."""

import ctypes as C
import hashlib, json
from pathlib import Path
import numpy as np


class NativeModel:

    def __init__(self, folder):
        self.folder = Path(folder)
        manifest = json.loads(
            (self.folder / "python_library.json").read_text()
        )
        path = self.folder / manifest["file"]

        if (
            path.parent != self.folder
            or hashlib.sha256(path.read_bytes()).hexdigest()
            != manifest["sha256"]
        ):
            raise ValueError("Custom model library checksum differs.")

        self.lib = C.CDLL(str(path))
        p = C.c_void_p
        d = C.POINTER(C.c_double)

        for name, args, ret in [
            ("sage_array", [C.c_int, C.c_int, d], p),
            ("sage_object", [], p),
            ("sage_field", [p, C.c_char_p, p], None),
            ("sage_free", [p], None),
            ("sage_rows", [p], C.c_int),
            ("sage_cols", [p], C.c_int),
            ("sage_values", [p], d),
            ("sage_call", [C.POINTER(p), C.POINTER(p)], C.c_int),
            ("sage_error", [], C.c_char_p),
        ]:

            fn = getattr(self.lib, name)
            fn.argtypes = args
            fn.restype = ret

    def simulate(
        self,
        theta,
        forcing,
        spec,
        *,
        y0=None,
        model_config=None,
        settings=(0.01, 1.0, 1e-4, 1e-3, 1e-3, 2.0, 10000.0),
        extra_data=None,
        all_outputs=False
    ):

        lib = self.lib

        def wrap(value):

            if isinstance(value, dict):

                obj = lib.sage_object()

                for key, item in value.items():
                    lib.sage_field(obj, key.encode(), wrap(item))

                return obj

            a = np.asarray(value, float)

            if a.ndim < 2:
                a = a.reshape(-1, 1)

            a = np.asfortranarray(a)

            return lib.sage_array(
                *a.shape, a.ctypes.data_as(C.POINTER(C.c_double))
            )

        forcing = np.asarray(forcing, float)
        theta = np.asarray(theta, float)

        if (
            forcing.ndim != 2
            or forcing.shape[1] != 3
            or not np.isfinite(forcing).all()
            or not np.isfinite(theta).all()
        ):
            raise ValueError("Invalid custom model inputs.")

        n = len(forcing)
        parameters = spec["parameters"]
        cfg = model_config or {}

        mdl = dict(
            th_min=[p["lower_bound"] for p in parameters],
            th_max=[p["upper_bound"] for p in parameters],
            pspace=0,
            tout=n,
            idx=[1, n + 1],
            L=cfg.get("L", 20),
            config=cfg,
            y0=(
                y0
                if y0 is not None
                else [s["initial_value"] for s in spec["states"]]
            ),
        )

        data = dict(P=forcing[:, 0], Ep=forcing[:, 1], T=forcing[:, 2])
        data.update(extra_data or {})

        ode = dict(
            zip(
                (
                    "InitStep",
                    "MaxStep",
                    "MinStep",
                    "RelTol",
                    "AbsTol",
                    "Order",
                    "maxiter",
                ),
                settings,
            )
        )
        ode["mem"] = 0

        inputs = [wrap(value) for value in (theta, mdl, data, ode)]
        outputs = (C.c_void_p * 9)()

        try:

            if lib.sage_call((C.c_void_p * 4)(*inputs), outputs):
                raise ValueError(lib.sage_error().decode())

            result = []

            for item in outputs:

                m, k = lib.sage_rows(item), lib.sage_cols(item)

                result.append(
                    np.ctypeslib.as_array(
                        lib.sage_values(item), shape=(m * k,)
                    )
                    .copy()
                    .reshape((m, k), order="F")
                    if m * k
                    else np.empty((m, k))
                )

            if result[4].size and result[4].ravel()[0]:
                raise ValueError("Custom model integration failed.")

            q = result[0].ravel()
            jac = result[1]

            if (
                q.shape != (n,)
                or jac.shape != (n, len(theta))
                or not np.isfinite(q).all()
                or not np.isfinite(jac).all()
            ):
                raise ValueError(
                    "Custom model returned invalid discharge or"
                    " sensitivities."
                )

            if all_outputs:
                observations = {"Q": q}
                jacobians = {"Q": jac}
                for name, value_index, jac_index in [
                    ("SWE", 5, 6),
                    ("SM", 7, 8),
                ]:
                    value = result[value_index].ravel()
                    derivative = result[jac_index]
                    if value.size:
                        if (
                            value.shape != (n,)
                            or derivative.shape != (n, len(theta))
                            or not np.isfinite(value).all()
                            or not np.isfinite(derivative).all()
                        ):
                            raise ValueError(
                                "Invalid custom-model " + name + " outputs."
                            )
                        observations[name] = value
                        jacobians[name] = derivative
                return observations, jacobians
            return q, jac

        finally:

            for item in outputs:

                if item:
                    lib.sage_free(item)

            for item in inputs:
                lib.sage_free(item)
