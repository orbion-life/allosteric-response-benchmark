"""Independent numerical checks before the response-operator study."""
import json,unittest
from pathlib import Path
import numpy as np
from scipy.linalg import eigh,qr
from scipy.integrate import quad
from scipy.sparse import csr_matrix
import operator_study as study

def dense_action(H,F,t):
    lam,U=eigh(H.toarray())
    c=study.mm(U.T,F)
    return study.mm(c.T*np.exp(-t*lam),c)-study.mm(F.T,F)

class OperatorTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.protocol=json.loads((study.ROOT/'preanalysis.json').read_text())
    def test_original_grid_response_parity(self):
        geometry=self.protocol['geometries'][0]
        H,F,times,sd,_=study.build(self.protocol,geometry,10,5,4)
        measured,_=study.full_reference(H,F,times)
        p=dict(self.protocol,geometry_A=geometry['coordinates_A'])
        r0,edges,lam,B,oldsd=study.reference.static_model(p,10)
        old,_=study.reference.grid(r0,edges,lam,B,10,3,5,4,'biquadratic',times,oldsd)
        for key in ['C','equilibrium','delayed']:
            np.testing.assert_allclose(measured[key],old[key],atol=1e-12,rtol=0)
        dense=np.array([dense_action(H,F,t) for t in times])
        np.testing.assert_allclose(measured['C'],dense,atol=1e-12,rtol=0)
    def test_residual_bounds_against_independent_dense_exponential(self):
        rng=np.random.default_rng(20260912)
        for scale in [.01,1,10]:
            raw=rng.normal(size=(16,16));H=csr_matrix(scale*study.mm(raw.T,raw)/16)
            F=rng.normal(size=(16,3));V=qr(np.column_stack((F,H@F)),mode='economic')[0]
            times=np.array([.1,1,10]);a,r=study.projected(H,F,V,times)
            truth=np.array([dense_action(H,F,t) for t in times])
            error=abs(truth-a['C'])
            self.assertTrue(np.all(error<=a['bound']+1e-10))
            self.assertLess(r['static_max_error'],1e-11)
            for C in a['C']:self.assertLess(np.linalg.eigvalsh(C).max(),1e-10)
    def test_nonorthogonal_basis_bound_accounts_for_mass_matrix(self):
        H=csr_matrix(np.diag([.1,.3,1.,2.,4.,6.]))
        F=np.eye(6)[:,:3]+.1
        V=qr(F,mode='economic')[0]*np.array([.99,1.01,1.02])
        times=np.array([.1,1,10]);a,r=study.projected(H,F,V,times)
        truth=np.array([dense_action(H,F,t) for t in times])
        self.assertGreater(r['mass_matrix_spectral_error'],.03)
        self.assertTrue(np.all(abs(truth-a['C'])<=a['bound']+1e-10))

    def test_integral_kernels_against_quadrature(self):
        lam=np.array([0,1e-13,.1,3.])
        for t in [.1,1,10]:
            phi,psi=study.kernels(lam,t)
            for i in range(4):
                for j in range(4):
                    rate=lam[i]+lam[j]
                    exact=quad(lambda s:np.exp(-rate*s),0,t,epsabs=1e-12)[0]
                    weighted=quad(lambda s:(t-s)*np.exp(-rate*s),0,t,epsabs=1e-12)[0]
                    self.assertAlmostEqual(phi[i,j],exact,places=10)
                    self.assertAlmostEqual(psi[i,j],weighted,places=10)
    def test_invariant_observable_subspace_stops_at_seed(self):
        H=csr_matrix(np.diag(np.arange(1,13,dtype=float)))
        F=np.eye(12)[:,:3]*np.array([.5,1,2])
        arrays,records,selected,_=study.block_basis(H,F,self.protocol,np.array([.1,1,10]))
        self.assertTrue(selected);self.assertEqual(records[-1]['rank'],3)
        self.assertLess(records[-1]['max_response_bound'],1e-9)
    def test_non_psd_operator_is_rejected(self):
        H=csr_matrix(np.diag([-1.,2,3,4]));F=np.eye(4)[:,:3]
        with self.assertRaises(ArithmeticError):study.projected(H,F,F,np.array([1.]))
    def test_missing_seed_component_is_counted(self):
        H=csr_matrix(np.diag([1.,2.,3.,4.]));F=np.eye(4)[:,:3]
        V=np.eye(4)[:,:2];a,r=study.projected(H,F,V,np.array([.1,1,10]))
        truth=np.array([dense_action(H,F,t) for t in [.1,1,10]])
        self.assertGreater(r['static_max_error'],.9)
        self.assertTrue(np.all(abs(truth-a['C'])<=a['bound']+1e-12))

if __name__=='__main__':unittest.main(verbosity=2)
