"""Check new methods, derivatives, convergence and forcing subdivision."""
import sys
from pathlib import Path
import unittest
import numpy as np

ROOT = Path(__file__).resolve().parent
sys.path[:0] = [str(ROOT/'cart'), str(ROOT/'CaRTPlus')]
import cart_ode_cpu as cart
import cartplus_continuous_cpu as plus
from test_cart_ode import P, Z, A, forcing
from test_cartplus_continuous import parameters, forcing as plus_forcing


class NumericsTests(unittest.TestCase):
    def test_cartplus_prediction_preserves_augmented_solver_output(self):
        for solver in ("dop853", "rk2", "rk23", "rk45", "heun_fixed", "euler"):
            with self.subTest(solver=solver):
                options = dict(
                    solver=solver, max_step=.01, refinement=2,
                    calibrated_only=True, rtol=1e-3, atol=1e-3,
                )
                inputs = (plus_forcing(), Z, A, parameters())
                expected, _ = plus.simulate(*inputs, jacobian=True, **options)
                discharge, jacobian = plus.simulate(
                    *inputs, keep_tangents=True, **options
                )
                np.testing.assert_array_equal(discharge, expected)
                self.assertIsNone(jacobian)

    def test_accuracy_convergence_and_directional_gradients(self):
        for model, p, f in ((cart,P,forcing()), (plus,parameters(),plus_forcing())):
            reference = model.simulate(f,Z,A,p,refinement=2,rtol=1e-11,atol=1e-13)[0]
            for method in ('rk45','heun_fixed','euler'):
                step = 1/768 if method == 'euler' else 1/192
                opts = dict(refinement=2,solver=method,max_step=step,
                            rtol=1e-10,atol=1e-12)
                q,j,d = model.simulate(f,Z,A,p,jacobian=True,
                                      return_diagnostics=True,**opts)
                self.assertLess(abs(d['balance_residual_mm']).max(),1e-8)
                self.assertLess(abs(q-reference).max(),.004)
                direction = np.cos(np.arange(len(p))+.3)*np.maximum(abs(p),.01)/10
                h = 1e-4
                fd = (model.simulate(f,Z,A,p+h*direction,**opts)[0]
                      -model.simulate(f,Z,A,p-h*direction,**opts)[0])/(2*h)
                np.testing.assert_allclose(j@direction,fd,atol=3e-5,rtol=3e-4,
                                           err_msg=f'{model.__name__}: {method}')
                if method != 'rk45':
                    coarse = model.simulate(f,Z,A,p,refinement=2,solver=method,max_step=2*step)[0]
                    ratio = np.linalg.norm(coarse-reference)/np.linalg.norm(q-reference)
                    self.assertGreater(ratio,1.7 if method=='euler' else 3.1)

    def test_rk45_daily_hourly_and_quarter_hour(self):
        for model,p,f in ((cart,P,forcing()),(plus,parameters(),plus_forcing())):
            options = dict(refinement=2,solver='rk45',rtol=1e-10,atol=1e-12)
            daily,j = model.simulate(f,Z,A,p,jacobian=True,**options)
            for resolution in (24,96):
                sub = np.repeat(f,resolution,axis=0)
                sub[:,[0,2]] /= resolution
                q,jac = model.simulate(sub,Z,A,p,resolution=resolution,jacobian=True,**options)
                np.testing.assert_allclose(q.reshape(-1,resolution).sum(1),daily,atol=2e-7,rtol=2e-7)
                np.testing.assert_allclose(jac.reshape(-1,resolution,len(p)).sum(1),j,atol=2e-6,rtol=2e-5)

    def test_calibrated_columns_and_step_validation(self):
        dry = np.array([[0.,15.,0.,1.]])
        q,_ = cart.simulate(dry,Z,A,P,solver='rk45',refinement=2,
                            initial=np.zeros(14),max_step=.02,min_step=.01)
        np.testing.assert_array_equal(q,0.)
        p,f = parameters(),plus_forcing()
        for method in ('rk45','heun_fixed','euler'):
            options = dict(solver=method,refinement=2,max_step=1/192,jacobian=True)
            q,full = plus.simulate(f,Z,A,p,**options)
            q35,j35 = plus.simulate(f,Z,A,p,calibrated_only=True,**options)
            np.testing.assert_allclose(q35,q,atol=2e-8,rtol=2e-8)
            np.testing.assert_allclose(j35,full[:,plus.CALIBRATED_COLUMNS],atol=2e-7,rtol=2e-5)
        for model,p,f in ((cart,P,forcing()),(plus,p,f)):
            for maximum,minimum in ((0,0),(.1,.2),(float('nan'),0)):
                with self.assertRaises(ValueError):
                    model.simulate(f,Z,A,p,max_step=maximum,min_step=minimum)
            with self.assertRaises(RuntimeError):
                model.simulate(f,Z,A,p,solver='rk45',max_step=1,min_step=.1,
                               rtol=1e-12,atol=1e-14,refinement=2)


if __name__ == '__main__':
    unittest.main()
