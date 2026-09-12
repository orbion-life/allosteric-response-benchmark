"""Focused tests of semantic portability and rejection of changed physics."""
import unittest
import json
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
import numpy as np
from verify import physical_comparison, gauge_consistency, product


class SemanticVerificationTests(unittest.TestCase):
    def test_absolute_physical_tolerance(self):
        a = np.array([0., 1.])
        self.assertTrue(physical_comparison(a, a+1e-12)['passed'])
        self.assertFalse(physical_comparison(a, a+1e-8)['passed'])

    def test_large_values_do_not_relax_absolute_tolerance(self):
        a = np.array([1e6])
        self.assertFalse(physical_comparison(a, a+1e-4)['passed'])

    def test_nonfinite_shape_dtype_and_discrete_changes_fail(self):
        a = np.array([0., 1.])
        for b in [np.array([np.nan, 1.]), np.array([np.inf, 1.]),
                  np.array([[0., 1.]]), a.astype(np.float32)]:
            self.assertFalse(physical_comparison(a, b)['passed'])
        self.assertFalse(physical_comparison(np.array([1, 2]), np.array([2, 1]))['passed'])

    def fixture(self):
        rng = np.random.default_rng(20260912)
        h = np.diag([0., 1., 4.]); b = rng.normal(size=(3,3)); r = np.diag([.2, .7, 1.2])
        times = np.array([.1, 1., 10.]); static = product(b.T,b)
        delayed = np.array([product(b.T*np.exp(-t*np.diag(h)),b) for t in times])
        arrays = {'rank3_Hr':h,'rank3_coefficients':b,'rank3_residual_Gram':r,
                  'rank3_equilibrium':static,'rank3_delayed':delayed,'rank3_C':delayed-static}
        return arrays,times,rng

    def test_rotated_basis_is_physically_equivalent(self):
        arrays,times,rng=self.fixture()
        self.assertTrue(gauge_consistency(arrays,'rank3_',times)['passed'])
        q,_=np.linalg.qr(rng.normal(size=(3,3)))
        rotated={k:v.copy() for k,v in arrays.items()}
        for field in ['Hr','residual_Gram']:
            rotated['rank3_'+field]=product(q.T,product(arrays['rank3_'+field],q))
        rotated['rank3_coefficients']=product(q.T,arrays['rank3_coefficients'])
        self.assertFalse(np.allclose(arrays['rank3_coefficients'],rotated['rank3_coefficients']))
        self.assertTrue(gauge_consistency(rotated,'rank3_',times)['passed'])

    def test_changed_observable_or_response_is_rejected(self):
        for field in ['coefficients','C','equilibrium','delayed']:
            arrays,times,_=self.fixture(); arrays['rank3_'+field].flat[0]+=1e-4
            self.assertFalse(gauge_consistency(arrays,'rank3_',times)['passed'],field)
        a=np.array([.002]); self.assertFalse(physical_comparison(a,a+1e-5)['passed'])

    def test_invalid_reduced_operator_or_residual_is_rejected(self):
        for field,value in [('Hr',-.1),('residual_Gram',-.1),('coefficients',np.nan)]:
            arrays,times,_=self.fixture(); arrays['rank3_'+field][0,0]=value
            self.assertFalse(gauge_consistency(arrays,'rank3_',times)['passed'],field)

    def test_full_cli_writes_diagnostics_before_rejecting_changed_response(self):
        root=Path(__file__).resolve().parent
        with tempfile.TemporaryDirectory(prefix='operator-verifier-test-') as temporary:
            replay=Path(temporary)/'replay';shutil.copytree(root/'replay-independent',replay)
            target=replay/'development-k1-n17-e4.npz'
            with np.load(target,allow_pickle=False) as archive:
                arrays={k:archive[k].copy() for k in archive.files}
            arrays['rank96_C'][0,0,0]+=1e-4
            np.savez_compressed(target,**arrays)
            output=Path(temporary)/'verification.json'
            run=subprocess.run([sys.executable,'-W','error',str(root/'verify.py'),
                                '--replay',str(replay),'--output',str(output)],capture_output=True,text=True)
            self.assertEqual(run.returncode,1,run.stdout+run.stderr)
            receipt=json.loads(output.read_text())
            self.assertEqual(receipt['status'],'FAIL')
            self.assertEqual(receipt['physical_arrays_compared'],786)
            failures=[x for x in receipt['failures'] if x['context']=='semantic physical array' and x.get('field')=='rank96_C']
            self.assertEqual(len(failures),1)
            self.assertGreater(failures[0]['diagnostic']['max_absolute_difference'],9.9e-5)


if __name__=='__main__':
    unittest.main()
