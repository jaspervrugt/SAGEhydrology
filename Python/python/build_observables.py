# ==========================================================================
# SAGEhydrology | Python source
# File: python/build_observables.py
#
# Build a separate state-history gateway to SAGE's unchanged C++ cores.
#
# Copyright (c) 2026 Jasper A. Vrugt
# University of California, Irvine | jasper@uci.edu
# ==========================================================================

"""Build a separate state-history gateway to SAGE's unchanged C++ cores."""

import argparse
import hashlib
import json
from pathlib import Path
from native_sources import native_root, library_directory
import re
import shutil
import subprocess
import sys
from build_crr import MODELS


def library_path():
    """Return the platform-specific native library path."""
    name = (
        "sage_observables.dll"
        if sys.platform == "win32"
        else (
            "libsage_observables.dylib"
            if sys.platform == "darwin"
            else "libsage_observables.so"
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
    soils = {
        1: [],
        2: [2],
        3: [1, 2, 3, 4, 5],
        4: [1, 2],
        5: [1],
        6: [1],
        7: [1],
    }
    for model, (folder, namespace, m, d) in MODELS.items():
        stem = "xinanjiang" if folder == "Xinanjiang" else folder
        header = root / folder / (stem + ".hpp")
        source = root / folder / (stem + ".cpp")
        block = re.search(
            r"struct Params\s*\{(.*?)\};", header.read_text(), re.S
        ).group(1)
        fields = re.findall(r"\b(double|int)\s+(\w+)\s*(?:=[^;]+)?;", block)
        headers.append(
            "#include "
            + json.dumps(str(header.resolve()).replace("\\", "/"))
        )
        sources.append(str(source.resolve()))
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
        expected = "4+p.K+p.L" if model == 7 else str(m)
        soil = ",".join(map(str, soils[model]))
        cases.append(
            (
                "case "
                f"{model}"
                ": {\n          if(np!="
                f"{len(fields)}"
                " || d!="
                f"{d}"
                ")return 2;\n          "
                f"{namespace}"
                "::Params p{};"
                f"{assignments}"
                "\n          if(m!="
                f"{expected}"
                ")return 2;\n          "
                f"{namespace}"
                "::Options opt{s[0],s[1],s[2],s[3],s[4],s[5],static_cast<in"
                "t>(s[6])};\n          "
                f"{namespace}"
                "::Forcing forcing{P,Ep,T,static_cast<std::size_t>(days)};"
                "\n          "
                f"{namespace}"
                "::OutputView out;out.Z=z.data();out.zrows=nt;out.zcols=m*("
                "d+1);\n          bool failed="
                f"{namespace}"
                "::run_into(days,z0.data(),z0.size(),forcing,p,opt,true,1,t"
                "rue,out);\n          if(failed)return 1;\n          "
                "std::vector<int> defaults={"
                f"{soil}"
                "};\n          extract(days,m,d,z.data(),defaults,sm_indice"
                "s,nsm,observed,jac);\n          return 0;\n        }"
            )
        )
        manifest[str(model)] = dict(
            fields=[name for _, name in fields],
            source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
        )
    code = (
        "\n".join(headers)
        + (
            "\n#include <vector>\n#include <cmath>\n#ifdef "
            "_WIN32\n#define API __declspec(dllexport)\n#else\n#define "
            'API __attribute__((visibility("default")))\n#endif\nstatic'
            " void extract(int n,int m,int d,const double*z,const "
            "std::vector<int>&defaults,const int*custom,int "
            "ncustom,double*out,double*J){\n const int "
            "nt=n+1;std::vector<int> soil=defaults;\n "
            "if(ncustom>=0)soil.assign(custom,custom+ncustom);\n "
            "for(int i=0;i<n;++i){\n  out[i]=z[i+1+nt*(m-1)]-z[i+nt*(m-"
            "1)];out[i+n]=z[i+1];out[i+2*n]=0;\n  for(int "
            "state:soil)out[i+2*n]+=z[i+1+nt*state];\n  for(int "
            "j=0;j<d;++j){\n   J[i+n*j]=z[i+1+nt*((j+1)*m+m-1)]-z[i+nt*"
            "((j+1)*m+m-1)];\n   J[i+n*j+n*d]=z[i+1+nt*((j+1)*m)];J[i+n"
            "*j+2*n*d]=0;\n   for(int state:soil)J[i+n*j+2*n*d]+=z[i+1+"
            'nt*((j+1)*m+state)];\n  }\n }\n}\nextern "C" API int '
            "sage_observables_run(int model,int days,int m,int d,const "
            "double*P,const double*Ep,const double*T,const "
            "double*v,int np,const double*s,const double*initial,const "
            "double*giuh,int ng,const int*sm_indices,int "
            "nsm,double*observed,double*jac) noexcept {\n "
            "if(days<1||m<1||d<1||!observed||!jac||!P||!Ep||!T||!v||!s|"
            "|!initial)return 2;\n for(int "
            "i=0;i<nsm;++i)if(sm_indices[i]<0||sm_indices[i]>=m-1)retur"
            "n 2;\n try{const std::size_t nt=static_cast<std::size_t>(d"
            "ays)+1;\n std::vector<double> "
            "z0(m*(d+1),0.),z(nt*m*(d+1),0.);for(int "
            "i=0;i<m;++i)z0[i]=initial[i];\n switch(model){\n"
        )
        + "\n".join(cases)
        + """default:return 2;}
 }catch(...){return 3;}
}
"""
    )
    gateway = output.with_suffix(".cpp")
    gateway.write_text(code, encoding="utf-8")
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
        *sources,
        "-o",
        str(output),
    ]
    subprocess.run(command, check=True)
    output.with_suffix(".json").write_text(
        json.dumps(
            dict(
                models=manifest,
                sha256=hashlib.sha256(output.read_bytes()).hexdigest(),
                command=command,
            ),
            indent=2,
        )
    )
    from build_hbv_routing import build as build_routing

    build_routing(sage_root, compiler, output_dir=output.parent)
    print(output)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--sage-root",
        type=Path,
        default=__import__("sage_paths").SAGE_RESOURCES,
    )
    p.add_argument("--compiler", default="g++")
    p.add_argument("--output-dir", type=Path)
    a = p.parse_args()
    build(a.sage_root, a.compiler, a.output_dir)
