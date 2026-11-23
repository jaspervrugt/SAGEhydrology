# ==========================================================================
# SITEhydrology | Python source
# File: python/user_model_native.py
#
# Discover and execute the installed MATLAB SITE user-model C++ sources.
#
# Copyright (c) 2026 Jasper A. Vrugt
# University of California, Irvine | jasper@uci.edu
# ==========================================================================

"""Discover and execute the installed MATLAB SITE user-model C++ sources."""

import ast
from functools import lru_cache
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys

import numpy as np

FILES = (
    "crr_user_model_mex.cpp",
    "user_model_prepare.cpp",
    "user_model_prepare.hpp",
    "user_model.cpp",
    "user_model.hpp",
    "delta_n.hpp",
    "log_fail_params.hpp",
)


def directory():
    from sage_paths import SOFTWARE_ROOT, PYTHON_ROOT

    candidates = [
        PYTHON_ROOT.parent / "user_model",
        SOFTWARE_ROOT / "user_model",
        PYTHON_ROOT / "user_model",
    ]
    if getattr(sys, "frozen", False):
        executable_dir = Path(sys.executable).parent
        if executable_dir.parent.parent.suffix == ".app":
            candidates.insert(0, executable_dir.parents[2] / "user_model")
        else:
            candidates.insert(0, executable_dir / "user_model")
    required = (*FILES[:5], "user_model_info.mat")
    return next(
        (path for path in candidates if all((path / name).is_file() for name in required)),
        None,
    )


def _number(text):
    def visit(node):
        if isinstance(node, ast.Constant) and isinstance(
            node.value, (int, float)
        ):
            return node.value
        if isinstance(node, ast.Name) and node.id in ("eps", "Inf", "pi"):
            return {
                "eps": np.finfo(float).eps,
                "Inf": float("inf"),
                "pi": np.pi,
            }[node.id]
        if isinstance(node, ast.UnaryOp) and isinstance(
            node.op, (ast.USub, ast.UAdd)
        ):
            return (-1 if isinstance(node.op, ast.USub) else 1) * visit(
                node.operand
            )
        if isinstance(node, ast.BinOp):
            a, b = visit(node.left), visit(node.right)
            if isinstance(node.op, ast.Div):
                return a / b
            if isinstance(node.op, ast.Mult):
                return a * b
            if isinstance(node.op, ast.Add):
                return a + b
            if isinstance(node.op, ast.Sub):
                return a - b
        raise ValueError(
            "Unsupported numeric user-model configuration: " + text
        )

    return float(visit(ast.parse(text.strip(), mode="eval").body))


def configuration(folder):
    path = Path(folder) / "user_model_configuration.m"
    if not path.is_file():
        return []
    text = path.read_text(encoding="utf-8-sig")

    def cells(name):
        match = re.search(r"\b" + name + r"\s*=\s*\{(.*?)\}", text, re.S)
        if not match:
            return []
        return [
            "".join(
                re.findall(r"'([^']*)'", item.replace("char(176)", "'°'"))
            )
            for item in match[1].split(";")
        ]

    def numbers(name):
        match = re.search(r"\b" + name + r"\s*=\s*\[(.*?)\]", text, re.S)
        return (
            [_number(value) for value in match[1].split(";")]
            if match
            else []
        )

    symbols, defaults = cells("symbols"), numbers("defaults")
    names, units, lower, refs = (
        cells("names"),
        cells("units"),
        numbers("lower"),
        cells("references"),
    )
    editable_match = re.search(r"\beditable\s*=\s*\[(.*?)\]", text, re.S)
    editable = []
    if editable_match:
        for term in editable_match[1].split(";"):
            match = re.fullmatch(
                r"\s*(true|false)(?:\((\d+),1\))?\s*", term
            )
            if not match:
                raise ValueError(
                    "Unsupported editable configuration metadata."
                )
            editable.extend([match[1] == "true"] * int(match[2] or 1))
    if len(defaults) != len(symbols):
        raise ValueError(
            "User-model configuration symbols/defaults differ in length."
        )
    return [
        dict(
            Symbol=symbol,
            Name=names[i] if i < len(names) else symbol,
            Unit=units[i] if i < len(units) else "-",
            Default=defaults[i],
            Editable=editable[i] if i < len(editable) else True,
            Lower=lower[i] if i < len(lower) else -float("inf"),
            Reference=refs[i] if i < len(refs) else "",
        )
        for i, symbol in enumerate(symbols)
    ]


def metadata(folder=None, resolution=1):
    from scipy.io import loadmat

    folder = Path(folder or directory() or "user_model")
    info = folder / "user_model_info.mat"
    if not info.is_file():
        raise FileNotFoundError(
            "Missing user_model_info.mat in " + str(folder)
        )
    data = loadmat(info, simplify_cells=True)
    rows = np.asarray(data["par_info"], dtype=object).reshape(-1, 6)
    parameters = []
    for _, symbol, description, unit, lo, hi in rows:
        unit = str(unit).replace("�", "°")
        factor = (
            resolution
            if unit.lower() in ("d", "day", "days")
            else (
                1 / resolution
                if re.search(r"/d(?=$|/)|/day", unit.lower())
                else 1
            )
        )
        if resolution != 1:
            token = "h" if resolution == 24 else "15 min"
            unit = (
                token
                if factor == resolution
                else re.sub(r"/days?(?=$|/)|/d(?=$|/)", "/" + token, unit)
            )
        parameters.append(
            dict(
                symbol=str(symbol),
                description=str(description),
                unit=unit,
                lower=float(lo) * factor,
                upper=float(hi) * factor,
            )
        )
    count = int(data["m"])
    routing = 20 if resolution == 1 else 100
    if data.get("model_template") == "gr4jb":
        count = 4 + 2 * routing
    initial = np.atleast_1d(data["y0"]).astype(float)
    if len(initial) != count:
        initial = np.full(count, initial[0])
    values = configuration(folder)
    config = {item["Symbol"]: item["Default"] for item in values}
    config["L"] = routing
    names = list(np.atleast_1d(data.get("state_name", [])))
    if data.get("model_template") == "gr4jb":
        names = [
            "snow_water_equivalent",
            "production_store",
            "routing_store",
            *[f"unit_hydrograph_1_storage_{i+1}" for i in range(routing)],
            *[f"unit_hydrograph_2_storage_{i+1}" for i in range(routing)],
            "discharge",
        ]
    states = [
        dict(
            Symbol=str(names[i]) if i < len(names) else f"S_{i+1}",
            Description=str(names[i]) if i < len(names) else f"State {i+1}",
            Unit="mm",
            Default=float(initial[i]),
            Editable=True,
            Capacity="",
        )
        for i in range(count)
    ]
    return dict(
        name="user_model",
        parameters=parameters,
        states=count,
        initial_states=initial.tolist(),
        definition=dict(
            states=states, values=values, config=config, y0=initial.tolist()
        ),
        folder=str(folder.resolve()),
    )


def compile_model(folder, cache=None, compiler=None):
    """
    Compile the actual user sources into a content-addressed local cache.
    """
    from sage_paths import OUTPUT_ROOT, PYTHON_ROOT

    folder = Path(folder).resolve()
    bridge = PYTHON_ROOT / "resources/user_model_bridge"
    sources = [folder / name for name in FILES if (folder / name).is_file()]
    for name in FILES[:5]:
        if not (folder / name).is_file():
            raise FileNotFoundError(folder / name)
    inputs = sources + [bridge / "mex.h", bridge / "python_gateway.cpp"]
    digest = hashlib.sha256()
    for file in inputs:
        digest.update(file.name.encode())
        digest.update(file.read_bytes())
    digest.update(sys.platform.encode())
    base = Path(cache) if cache else OUTPUT_ROOT.parent / "user_models"
    target = base / digest.hexdigest()[:24]
    manifest = target / "python_library.json"
    if manifest.is_file():
        return target
    # Accept a precompiled matching local bridge on machines without a
    # compiler.
    local = folder / "python_library.json"
    if local.is_file():
        content = json.loads(local.read_text())
        if content.get("source_sha256") == digest.hexdigest():
            return folder
    executable = compiler or shutil.which("g++") or shutil.which("clang++")
    if not executable and sys.platform == "win32":
        candidates = sorted(
            Path(os.environ.get("PROGRAMDATA", "C:/ProgramData")).glob(
                "MATLAB/SupportPackages/*/3P.instrset/mingw*/bin/g++.exe"
            ),
            reverse=True,
        )
        executable = str(candidates[0]) if candidates else None
    if not executable:
        raise RuntimeError(
            "Compile the user-model Python library with a C++17 compiler"
            " first."
        )
    target.mkdir(parents=True, exist_ok=True)
    for file in inputs:
        shutil.copy2(file, target / file.name)
    library = target / (
        "site_user_model.dll"
        if sys.platform == "win32"
        else (
            "libsite_user_model.dylib"
            if sys.platform == "darwin"
            else "libsite_user_model.so"
        )
    )
    flags = (
        ["-static", "-static-libgcc", "-static-libstdc++"]
        if sys.platform == "win32"
        else ["-fPIC"]
    )
    command = [
        executable,
        *flags,
        "-std=c++17",
        "-O2",
        "-shared",
        "-I",
        str(target),
        *[
            str(target / name)
            for name in (
                "python_gateway.cpp",
                "crr_user_model_mex.cpp",
                "user_model_prepare.cpp",
                "user_model.cpp",
            )
        ],
        "-o",
        str(library),
    ]
    completed = subprocess.run(
        command, capture_output=True, text=True, timeout=180
    )
    (target / "compile_log.txt").write_text(
        completed.stdout + completed.stderr
    )
    if completed.returncode:
        raise RuntimeError(
            "User-model compilation failed. See "
            + str(target / "compile_log.txt")
        )
    manifest.write_text(
        json.dumps(
            dict(
                file=library.name,
                source_sha256=digest.hexdigest(),
                sha256=hashlib.sha256(library.read_bytes()).hexdigest(),
            )
        )
    )
    return target


@lru_cache(maxsize=8)
def _model(library_folder):
    from user_model_bridge import NativeModel

    return NativeModel(library_folder)


@lru_cache(maxsize=8)
def _run_metadata(folder, resolution):
    return metadata(folder, resolution)


def simulate(
    theta,
    forcing,
    *,
    y0=None,
    model_config=None,
    settings=(0.01, 1.0, 1e-4, 1e-3, 1e-3, 2.0, 10000.0),
    routing_start=0,
):
    cfg = dict(model_config or {})
    folder = cfg.pop("_user_model_folder", None) or directory()
    resolution = cfg.pop("_resolution", 1)
    info = cfg.pop("_user_model_metadata", None) or _run_metadata(
        str(folder), resolution
    )
    library = cfg.pop("_user_model_library", None) or compile_model(folder)
    spec = dict(
        parameters=[
            dict(lower_bound=p["lower"], upper_bound=p["upper"])
            for p in info["parameters"]
        ],
        states=[
            dict(initial_value=value) for value in info["initial_states"]
        ],
    )
    defaults = dict(
        info["definition"]["config"],
        **{k: v for k, v in cfg.items() if not k.startswith("_")},
    )
    return _model(str(library)).simulate(
        theta,
        forcing,
        spec,
        y0=y0,
        model_config=defaults,
        settings=settings,
    )
