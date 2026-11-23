# ========================================================================
# SAGEhydrology | Model source
# File: Native/models/CaRTPlus/build_uniform.py
#
# CaRT / CaRT+ numerical model and supporting routines.
#
# SAGE adaptations: Jasper A. Vrugt, University of California, Irvine.
# Upstream model: Xia (2026). Preserve the adjacent license notices.
# ========================================================================

"""
Build versioned Xia f986af2 uniform-input libraries without replacing legacy
DLLs.
"""
import argparse
import hashlib
import json
from pathlib import Path
import platform
import shutil
import subprocess

REVISION = "f986af265633a7a980836c2e0f2124077d632420"


def build(compiler="g++"):
    """
    Build the isolated native model libraries.

    SYNOPSIS:
        build(...)

    INPUT ARGUMENTS:
    compiler : C++ compiler executable used by this build command.

    OUTPUT ARGUMENTS:
    None. Writes isolated native libraries and their build metadata.
    """
    root = Path(__file__).resolve().parent
    executable = shutil.which(compiler)
    if executable is None:
        raise FileNotFoundError(f"Compiler unavailable: {compiler}")
    flags = ["-std=c++17", "-O3", "-shared", "-DCART_FAST_SCALAR",
             "-DCARTPLUS_UNIFORM_INPUT"]
    flags += ["-static"] if platform.system() == "Windows" else ["-fPIC"]
    extension = ".dll" if platform.system() == "Windows" else ".so"
    outputs = {}
    for calibrated in (False, True):
        name = "cartplus_uniform_f986af2" + ("_35" if calibrated else "") \
            + extension
        command = [executable, *flags]
        if calibrated:
            command += ["-DCARTPLUS_CALIBRATED_ONLY"]
        command += [str(root / "cartplus_bridge.cpp"), "-o", str(root /
            name)]
        subprocess.run(command, check=True)
        outputs[name] = hashlib.sha256((root / name).read_bytes()). \
            hexdigest()
    manifest = dict(upstream="https://github.com/xiaxilin/cartplus",
        branch="hourly", revision=REVISION,
        defaults=dict(prod_lie=True, nsub=4, nsub_slow=1),
        compiler=subprocess.check_output([executable,"--version"],text=
            True).splitlines()[0],
        sources={p.name:hashlib.sha256(p.read_bytes()).hexdigest()
                 for p in [root/"cartplus_bridge.cpp", *sorted(root.glob(
                     "*.hpp"))]},
        libraries=outputs)
    (root / "cartplus_uniform_build.json").write_text(json.dumps(manifest,
        indent=2)+"\n")
    print(json.dumps(outputs,indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--compiler",default="g++")
    build(parser.parse_args().compiler)
