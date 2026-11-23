# ==========================================================================
# SAGEhydrology | Python source
# File: python/build_crr.py
#
# Build a MATLAB-free library from all seven public SAGE C++ model cores.
#
# Copyright (c) 2026 Jasper A. Vrugt
# University of California, Irvine | jasper@uci.edu
# ==========================================================================

"""
Build a MATLAB-free library from all seven public SAGE C++ model
cores.
"""

import argparse
import hashlib
import json
from pathlib import Path
from native_sources import native_root, library_directory
import re
import shutil
import subprocess
import sys

MODELS = {
    1: ("hymod", "sage_hymod", 7, 7),
    2: ("hmodel", "sage_hmodel", 6, 9),
    3: ("sacsma", "sage_sacsma", 10, 15),
    4: ("Xinanjiang", "sage_xinanjiang", 9, 16),
    5: ("gr4jA", "sage_gr4ja", 10, 8),
    6: ("hbv", "sage_hbv", 5, 12),
    7: ("cfe_nwm", "sage_cfe_nwm", 14, 15),
}


def library_path():
    """Return the platform-specific native library path."""
    name = (
        "site_crr.dll"
        if sys.platform == "win32"
        else (
            "libsite_crr.dylib"
            if sys.platform == "darwin"
            else "libsite_crr.so"
        )
    )
    return library_directory() / name


def build(sage_root, compiler="g++", output_dir=None):
    """Compile the public native numerical library."""
    root = native_root(sage_root) / "models"
    output = library_path().resolve()
    if output_dir is not None:
        output = Path(output_dir).resolve() / output.name
    output.parent.mkdir(parents=True, exist_ok=True)
    headers = []
    sources = []
    cases = []
    manifest = {}
    for model, (folder, namespace, m, d) in MODELS.items():
        stem = "xinanjiang" if folder == "Xinanjiang" else folder
        header = root / folder / (stem + ".hpp")
        source = root / folder / (stem + ".cpp")
        if not header.exists() or not source.exists():
            raise FileNotFoundError(header)
        block = re.search(
            r"struct Params\s*\{(.*?)\};", header.read_text(), re.S
        ).group(1)
        fields = re.findall(r"\b(double|int)\s+(\w+)\s*(?:=[^;]+)?;", block)
        headers.append(
            "#include " + json.dumps(str(header).replace("\\", "/"))
        )
        sources.append(source)
        assignments = "\n".join(
            f"p.{name} = "
            + (
                f"static_cast<int>(v[{i}]);"
                if kind == "int"
                else f"v[{i}];"
            )
            for i, (kind, name) in enumerate(fields)
        )
        if model == 7:
            assignments += "\np.giuh=giuh;"
        state_expression = "3+p.L+p.K+1" if model == 7 else str(m)
        cases.append(
            (
                "case "
                f"{model}"
                ": {\n            if(np!="
                f"{len(fields)}"
                " || d!="
                f"{d}"
                ") return 2;\n            "
                f"{namespace}"
                "::Params p{}; "
                f"{assignments}"
                "\n            const int expected="
                f"{state_expression}"
                "; if(m!=expected) return 2;\n            if(model==7 "
                "&& ng!=p.L) return 2;\n            "
                f"{namespace}"
                "::Options opt{s[0],s[1],s[2],s[3],s[4],s[5],static_cas"
                "t<int>(s[6])};\n            "
                f"{namespace}"
                "::Forcing forcing{P,Ep,T,static_cast<std::size_t>(days"
                ")};\n            "
                f"{namespace}"
                "::OutputView out;out.Z=z.data();out.zrows=1;out.zcols="
                "z.size();out.q=q;out.nq=days;out.J=J;out.nj=d;\n      "
                "      return "
                f"{namespace}"
                "::run_into(days,z0.data(),z0.size(),forcing,p,opt,fals"
                "e,1,true,out)?1:0;\n        }"
            ).replace(
                "if(model==7 && ng!=p.L) return 2;",
                "if(ng!=p.L) return 2;" if model == 7 else "",
            )
        )
        manifest[str(model)] = dict(
            folder=folder,
            states=m,
            core_parameters=d,
            fields=[name for _, name in fields],
            source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
        )
    code = (
        "\n".join(headers)
        + (
            "\n#include <vector>\n#include <cmath>\n#include "
            "<climits>\n#ifdef _WIN32\n#define API "
            "__declspec(dllexport)\n#else\n#define API "
            '__attribute__((visibility("default")))\n#endif\nextern '
            '"C" API int site_crr_run(int model,int days,int m,int '
            "d,const double*P,const double*Ep,const double*T,const "
            "double*v,int np,const double*s,const double*initial,const "
            "double*giuh,int ng,double*q,double*J) noexcept {\n    "
            "if(days<1||m<1||d<1||!P||!Ep||!T||!v||!s||!initial||!q||!J"
            "||s[6]<1||s[6]>INT_MAX)return 2;\n    try { "
            "std::vector<double> z0(m*(d+1),0.),z(m*(d+1),0.);for(int "
            "i=0;i<m;++i)z0[i]=initial[i];\n    switch(model) {\n"
        )
        + "\n".join(cases)
        + (
            "default:return 2;}\n    }catch(...){return 3;}\n}\nextern "
            '"C" API void site_crr_rss(int n,int d,const '
            "double*q,const double*J,const double*y,const unsigned "
            "char*mask,const double*scale,double*f,double*g) noexcept "
            "{\n    *f=0.;for(int j=0;j<d;++j)g[j]=0.;\n    for(int "
            "i=0;i<n;++i)if(mask[i]){double "
            "r=q[i]-y[i];*f+=r*r;for(int j=0;j<d;++j)g[j]+=2*r*J[i+n*j]"
            "*scale[j];}\n}\n"
        )
    )
    gateway = output.parent / "crr_generated.cpp"
    gateway.write_text(code)
    executable = shutil.which(compiler)
    if not executable:
        raise FileNotFoundError(compiler)
    flags = ["-std=c++17", "-O2", "-fwrapv", "-DNDEBUG", "-shared"] + (
        ["-static", "-static-libgcc", "-static-libstdc++"]
        if sys.platform == "win32"
        else ["-fPIC"]
    )
    command = [
        executable,
        *flags,
        str(gateway),
        *[str(s) for s in sources],
        "-o",
        str(output),
    ]
    subprocess.run(command, check=True)
    metadata = dict(
        models=manifest,
        command=command,
        compiler=subprocess.check_output(
            [executable, "--version"], text=True
        ),
        library_sha256=hashlib.sha256(output.read_bytes()).hexdigest(),
    )
    (output.parent / "crr_build.json").write_text(
        json.dumps(metadata, indent=2) + "\n"
    )
    from build_hbv_routing import build as build_routing

    build_routing(sage_root, compiler, output_dir=output.parent)
    print(output)


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
