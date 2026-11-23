# ==========================================================================
# SAGEhydrology | Python source
# File: python/build_native.py
#
# Build a ctypes library from the existing SAGE HYMOD C++ core (requires
# g++).
#
# Copyright (c) 2026 Jasper A. Vrugt
# University of California, Irvine | jasper@uci.edu
# ==========================================================================

"""
Build a ctypes library from the existing SAGE HYMOD C++ core
(requires g++).
"""

import argparse
import hashlib
import json
from pathlib import Path
from native_sources import native_root, library_directory
import platform
import shutil
import subprocess
import sys


def library_path():
    """Return the platform-specific native library path."""
    name = (
        "site_hymod.dll"
        if sys.platform == "win32"
        else (
            "libsite_hymod.dylib"
            if sys.platform == "darwin"
            else "libsite_hymod.so"
        )
    )
    return library_directory() / name


def build(sage_root, compiler="g++", output_dir=None):
    """Compile the public native numerical library."""
    core = native_root(sage_root) / "models" / "hymod"
    wrapper = native_root(sage_root) / "adapters/python/hymod_c_api.cpp"
    for path in (core / "hymod.cpp", core / "hymod.hpp"):
        if not path.is_file():
            raise FileNotFoundError(path)
    executable = shutil.which(compiler)
    if executable is None:
        raise FileNotFoundError(f"Compiler unavailable: {compiler}")
    output = library_path().resolve()
    if output_dir is not None:
        output = Path(output_dir).resolve() / output.name
    output.parent.mkdir(parents=True, exist_ok=True)
    flags = ["-std=c++17", "-O3", "-DNDEBUG", "-shared"]
    if sys.platform == "win32":
        flags += ["-static", "-static-libgcc", "-static-libstdc++"]
    else:
        flags += ["-fPIC"]
    command = [
        executable,
        *flags,
        "-I",
        str(core),
        str(wrapper),
        str(core / "hymod.cpp"),
        "-o",
        str(output),
    ]
    subprocess.run(command, check=True)
    metadata = dict(
        compiler=subprocess.check_output(
            [executable, "--version"], text=True
        ),
        command=command,
        platform=platform.platform(),
        sources={
            str(p): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in (core / "hymod.cpp", core / "hymod.hpp", wrapper)
        },
        library_sha256=hashlib.sha256(output.read_bytes()).hexdigest(),
    )
    (output.parent / "build.json").write_text(
        json.dumps(metadata, indent=2) + "\n"
    )
    print(output)
    return output


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--sage-root",
        type=Path,
        default=__import__("sage_paths").SAGE_RESOURCES,
    )
    parser.add_argument("--compiler", default="g++")
    parser.add_argument("--output-dir", type=Path)
    args = parser.parse_args()
    build(args.sage_root, args.compiler, args.output_dir)
