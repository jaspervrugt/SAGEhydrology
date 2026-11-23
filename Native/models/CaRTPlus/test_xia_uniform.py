# ========================================================================
# SAGEhydrology | Model source
# File: Native/models/CaRTPlus/test_xia_uniform.py
#
# CaRT / CaRT+ numerical model and supporting routines.
#
# SAGE adaptations: Jasper A. Vrugt, University of California, Irvine.
# Upstream model: Xia (2026). Preserve the adjacent license notices.
# ========================================================================

"""
Compare SAGE with unmodified Xia f986af2 kernel bodies on a CPU.

The reference executes upstream Triton forward kernels through the existing
NumPy shim, including float32 output rounding. No GPU or training is needed.
"""
import argparse
import importlib.util
import json
from pathlib import Path
import sys
from types import ModuleType, SimpleNamespace
import numpy as np
import torch
from compare_xia_cpu import FixedNetwork, install_cpu_shim
from cartplus_cpu import simulate, _entrypoint
from cartplus_inputs import physical_parameters


def main(upstream):
    """
    Run this command-line build, comparison, or validation task.

    SYNOPSIS:
        main(...)

    INPUT ARGUMENTS:
    upstream : Upstream used by this numerical operator.

    OUTPUT ARGUMENTS:
    results : computed result in the operand units.
    """
    install_cpu_shim()
    ld = ModuleType("triton.language.extra.libdevice")
    for name in ("expm1", "log1p", "tanh", "exp", "log", "pow"):
        setattr(ld, name, np.power if name == "pow" else getattr(np, name))
    extra = ModuleType("triton.language.extra")
    extra.libdevice = ld
    sys.modules.update({"triton.language.extra": extra,
                        "triton.language.extra.libdevice": ld})
    sys.path.insert(0, str(upstream / "src"))
    from cart.config import config_for
    from cart.model import CaRT
    from cart.parameters import static_table
    root = Path(__file__).resolve().parent
    spec = json.loads((root / "parameters.json").read_text())
    shared = spec["shared"]
    law = {k: v[2] for k, v in spec["law"].items()}
    z = np.array([100.,220.,350.,510.,740.])
    area = np.array([.1,.2,.3,.25,.15])
    torch.set_num_threads(1)
    baseline_spec = importlib.util.spec_from_file_location(
        "uniform_cart_reference", root.parent / "cart/cart_python.py")
    baseline = importlib.util.module_from_spec(baseline_spec)
    baseline_spec.loader.exec_module(baseline)
    results = []
    for resolution in (1, 24, 96):
        t = np.arange(2200 if resolution == 1 else 480, dtype=float)
        f = np.column_stack(((2+1.5*np.sin(t/7)+12*(np.sin(t/3)>.8))/
            resolution,
                             4-14*np.cos(t/25),
                             (1+.5*np.cos(t/9))/resolution, t>=20))
        dataset = SimpleNamespace(P=f[:,:1], Tm=f[:,1:2]+273.15,
            E=f[:,2:3], static=np.zeros((1,1),np.float32),
            band_H=z[None], band_A=area[None], mask=lambda _: f[:,3].
                astype(bool))
        for name in ("cart", "cartplus"):
            cfg = config_for(name, "gb", dt=1/resolution)
            static = {k:v[2] for k,v in static_table(cfg).items()}
            model = CaRT(dataset,cfg,"cpu")
            network = SimpleNamespace(
                phg={k:torch.tensor(shared[k],dtype=torch.float32)
                     for k in ("Tv","lsv","ltv")},
                lawg={k:torch.tensor(shared[k],dtype=torch.float32)
                      for k in ("TF","lsF","thc","lsth")},
                law=lambda _: {k:torch.tensor([[v]],dtype=torch.float32)
                    for k,v in law.items()})
            model.net=FixedNetwork(network,static,False)
            with torch.no_grad(), np.errstate(all="ignore"):
                reference=model.simulate()[0].ravel().numpy()
            if name=="cart":
                actual,jac=baseline.simulate(f,z,area,list(static.values()),
                    mcode=4,jacobian=True,resolution=resolution,
                        uniform_input=True)
                balance=None
            else:
                p=physical_parameters(static,law,shared)
                output,jac=simulate(root/"cartplus_uniform_f986af2.dll",f,
                    z,area,p,
                    jacobian=True,resolution=resolution)
                actual=output[:,0]
                balance=float(np.max(np.abs(f[:,0]+output[:,2]-output[:,0]
                    -output[:,1]-output[:,3])))
                assert balance<1e-7, balance
            assert np.isfinite(jac).all()
            error=float(np.max(np.abs(actual-reference)))
            relative=float(np.linalg.norm(actual-reference)/max(np.linalg.
                norm(reference),1e-20))
            row=dict(model=name,resolution=resolution,max_abs=error,
                relative_l2=relative,
                     balance_max_abs=balance)
            print(json.dumps(row),flush=True)
            results.append(row)
            np.testing.assert_allclose(actual,reference,rtol=8e-4,atol=2e-4)
            if name == "cartplus":
                short = f[:12].copy()
                short[:,3] = np.arange(len(short)) >= 2
                a, ja = simulate(root/"cartplus_uniform_f986af2.dll",short
                    ,z,area,p,
                                 jacobian=True,resolution=resolution)
                b, jb = simulate(root/"missing.dll",short,z,area,p,
                                 jacobian=True,resolution=resolution,mcode
                                     =2,
                                 uniform_input=True)
                np.testing.assert_allclose(a,b,rtol=1e-6,atol=3e-6)
                np.testing.assert_allclose(ja,jb,rtol=1e-4,atol=2e-4)
                c, jc = simulate(root/"cartplus_uniform_f986af2_35.dll",
                    short,z,area,p,
                                  jacobian=True,resolution=resolution,
                                      calibrated_only=True)
                columns=[i for i in range(42) if i not in (22,23,24,28,29,
                    37,38)]
                np.testing.assert_allclose(c,a,rtol=1e-12,atol=1e-12)
                np.testing.assert_allclose(jc,ja[:,columns],rtol=1e-12,
                    atol=1e-12)
    from test_cartplus_jacobian import check
    for resolution in (1,24,96):
        print("finite-difference Jacobian",resolution,
              check(resolution,4,"cartplus_uniform_f986af2.dll",1e-8),
                  flush=True)
    return results


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("upstream",type=Path)
    args=parser.parse_args()
    main(args.upstream.resolve())
