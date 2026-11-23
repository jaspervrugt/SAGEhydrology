# ==========================================================================
# SAGEhydrology | Python source
# File: python/build_hbv_routing.py
#
# Build the C adapter to the HBV router also used by MATLAB.
#
# Copyright (c) 2026 Jasper A. Vrugt
# University of California, Irvine | jasper@uci.edu
# ==========================================================================

"""Build a small routing library without rebuilding the ODE model cores."""

import argparse
import hashlib
import json
from pathlib import Path
from native_sources import native_root, library_directory
import shutil
import subprocess
import sys


def library_path():
    """Return the platform-specific native library path."""
    name = (
        "hbv_routing.dll"
        if sys.platform == "win32"
        else (
            "libhbv_routing.dylib"
            if sys.platform == "darwin"
            else "libhbv_routing.so"
        )
    )
    return library_directory() / name


def build(sage_root, compiler="g++", output_dir=None):
    """Compile the public native numerical library."""
    source = native_root(sage_root) / "models/hbv"
    files = [source / "hbv_routing_c.cpp", source / "hbv_routing.hpp"]
    hashes = {
        p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in files
    }
    output = library_path().resolve()
    if output_dir is not None:
        output = Path(output_dir).resolve() / output.name
    metadata = output.with_suffix(".json")
    if output.exists() and metadata.exists():
        old = json.loads(metadata.read_text(encoding="utf-8"))
        if old.get("sources") == hashes:
            if (
                old.get("library_sha256")
                == hashlib.sha256(output.read_bytes()).hexdigest()
            ):
                return output
    executable = shutil.which(compiler)
    if not executable:
        raise FileNotFoundError(compiler)
    output.parent.mkdir(parents=True, exist_ok=True)
    flags = [
        "-std=c++17",
        "-O3",
        "-DNDEBUG",
        "-shared",
        "-ffp-contract=off",
    ]
    if sys.platform == "win32":
        flags += ["-static", "-static-libgcc", "-static-libstdc++"]
    else:
        flags += ["-fPIC"]
    command = [executable, *flags, str(files[0]), "-o", str(output)]
    subprocess.run(command, check=True)
    metadata.write_text(
        json.dumps(
            dict(
                sources=hashes,
                command=command,
                library_sha256=hashlib.sha256(
                    output.read_bytes()
                ).hexdigest(),
            ),
            indent=2,
        ),
        encoding="utf-8",
    )
    print(output)
    return output


if __name__ == "__main__":
    resources = Path(__file__).resolve().parents[1] / "resources/sage"

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sage-root", type=Path, default=resources)
    parser.add_argument("--compiler", default="g++")
    parser.add_argument("--output-dir", type=Path)
    args = parser.parse_args()
    build(args.sage_root, args.compiler, args.output_dir)
