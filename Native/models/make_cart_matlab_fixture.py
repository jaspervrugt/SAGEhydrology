"""Write native Python expectations for the MATLAB gateway parity test."""
from pathlib import Path
import sys
import numpy as np
from scipy.io import savemat

ROOT = Path(__file__).resolve().parent
sys.path[:0] = [str(ROOT/'cart'), str(ROOT/'CaRTPlus')]
import cart_ode_cpu
import cart_python
import cartplus_continuous_cpu
import cartplus_cpu
from test_cart_ode import P, Z, A, forcing
from test_cartplus_continuous import parameters, forcing as plus_forcing


def build():
    cases = []
    for model,p,weather in ((19,P,forcing()), (18,parameters(),plus_forcing())):
        for resolution in (1,24,96):
            f = np.repeat(weather,resolution,axis=0)
            f[:,[0,2]] /= resolution
            for method in ('dop853','rk2','rk23','rk45','heun_fixed','euler'):
                options = dict(refinement=2,rtol=1e-7,atol=1e-9,max_step=1/768)
                native = cart_ode_cpu if model==19 else cartplus_continuous_cpu
                q,j = native.simulate(f,Z,A,p,jacobian=True,solver=method,
                                      resolution=resolution,**options)
                cases.append(dict(model=model,p=p,f=f,resolution=resolution,
                    formulation='ode',method=method,uniform=False,q=q,j=j))
            for uniform in (False,True):
                if model==19:
                    q,j = cart_python.simulate(f,Z,A,p,jacobian=True,mcode=4,
                        resolution=resolution,uniform_input=uniform)
                else:
                    library = 'cartplus_uniform_f986af2.dll' if uniform else 'cartplus_cpu.dll'
                    out,j = cartplus_cpu.simulate(ROOT/'CaRTPlus'/library,f,Z,A,p,
                        jacobian=True,resolution=resolution)
                    q = out[:,0]
                cases.append(dict(model=model,p=p,f=f,resolution=resolution,
                    formulation='sequential',method='',uniform=uniform,q=q,j=j))
    destination = ROOT.parents[1]/'Matlab/tests/cart_numerics_expected.mat'
    savemat(destination,dict(cases=np.array(cases,dtype=object),z=Z,a=A),
            do_compression=True,oned_as='column')
    print(f'{len(cases)} cases: {destination}')


if __name__ == '__main__':
    build()
