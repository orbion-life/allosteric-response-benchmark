"""Compare operator tensor implementation with independent direct-energy GH.

This imports the implementation under test, not for the expected coefficients.
Expected values use hermite_quadrature_check.py's direct quartic contact energies.
Both sides use the same explicitly supplied covariance eigenvector chart.
"""
from pathlib import Path
import importlib.util
import json
import math
import sys
import numpy as np
import hermite_quadrature_check as gh

HERE=Path(__file__).resolve().parent
BASE=HERE.parents[1]
SOURCE=BASE/'operator/source/core.py'
spec=importlib.util.spec_from_file_location('pulsar_operator_under_test',SOURCE)
operator=importlib.util.module_from_spec(spec)
spec.loader.exec_module(operator)

def main():
    control_path=BASE/'nonlinear/results/controls.npz'
    snapshot_path=BASE/'nonlinear/inputs/frozen-gaussian-result.npz'
    control=gh.load(control_path);snapshot=gh.load(snapshot_path)
    model={'r0':control['r0'],'edges':control['edges'],'B':control['basis']}
    sigma=control['Sigma_Gaussian'];sd=control['sd']
    nodes=np.array([0.,.1,1.])/control['eigenvalues'][0]
    tested=operator.Hermite(model,sigma,sd)
    expected=gh.quadrature_coefficients(model['r0'],model['edges'],model['B'],sigma,sd,
                                        eigenvectors=tested.V)
    alphas=expected['alphas'];occupations=[tuple(np.repeat(np.arange(len(alpha)),alpha)) for alpha in alphas]
    h=np.array([tested.coeff(indices) for indices in occupations])
    Z=np.concatenate([expected['coefficients']*np.exp(-t*expected['rates'])[:,None] for t in nodes],axis=1)
    expected_psi=gh.mm(Z,snapshot['W'])
    expected_factor=np.sqrt(expected['rates'])[:,None]*expected_psi
    rows=[tested.row(indices,snapshot['W'],nodes) for indices in occupations]
    actual_psi=np.array([r[0] for r in rows]);actual_factor=np.array([r[1] for r in rows]);rates=np.array([r[2] for r in rows])
    checks={
        'all_34_coefficients_maxabs':float(abs(h-expected['coefficients']).max()),
        'degree_coefficients_maxabs':{str(k):float(abs(h[np.sum(alphas,axis=1)==k]-expected['coefficients'][np.sum(alphas,axis=1)==k]).max()) for k in range(1,5)},
        'all_34_Psi_rows_maxabs':float(abs(actual_psi-expected_psi).max()),
        'all_34_factor_rows_maxabs':float(abs(actual_factor-expected_factor).max()),
        'all_34_rates_maxabs':float(abs(rates-expected['rates']).max()),
        'factor_Gram_vs_saved_Hr_tau_maxabs':float(abs(gh.mm(actual_factor.T,actual_factor)-snapshot['Hr']).max()/control['eigenvalues'][0]),
        'mass_vs_identity_maxabs':float(abs(gh.mm(actual_psi.T,actual_psi)-np.eye(len(snapshot['Hr']))).max()),
        'observable_embedding_vs_saved_B_maxabs':float(abs(gh.mm(actual_psi.T,h)-snapshot['B']).max()),
    }
    # Check the mu/beta factor in spectral rates separately, with fixed Sigma.
    changed=operator.Hermite(model,sigma,sd,mu=3.,beta=2.)
    checks['nonunit_mu_beta_rate_scaling_maxabs']=float(abs(changed.rates-1.5*tested.rates).max())
    coefficient_tolerance=1e-10;embedding_tolerance=1e-8
    passed=checks['all_34_coefficients_maxabs']<coefficient_tolerance
    passed=passed and all(v<embedding_tolerance for k,v in checks.items() if isinstance(v,float))
    receipt={'status':'PASS' if passed else 'FAIL','source_under_test':str(SOURCE),
       'scope':'All 34 centered degree1-4 coefficients and factor rows for existing three-coordinate triangle only.',
       'independence':'Expected coefficients are order5 GH integrals of direct energy values, without normal-ordering formulas; operator implementation is only the tested side.',
       'coordinate_chart':'Covariance eigenvectors supplied by implementation to independent quadrature; avoids meaningless eigenvector sign differences.',
       'quadrature_order':gh.ORDER,'randomness':'None','basis_size':34,'constants':{'beta':1.,'mu':1.,'kappa':1.},
       'nonunit_rate_check':{'mu':3.,'beta':2.,'Sigma':'held fixed'},
       'tolerances':{'coefficients':coefficient_tolerance,'embedding':embedding_tolerance},'checks':checks,
       'source_hashes':{'script':gh.digest(__file__),'independent_helper':gh.digest(gh.__file__),
                        'operator_source':gh.digest(SOURCE),'controls':gh.digest(control_path),'snapshot_model':gh.digest(snapshot_path)},
       'python':sys.version,'numpy':np.__version__}
    (HERE/'operator-Hermite-comparison.json').write_text(json.dumps(receipt,indent=2,allow_nan=False)+'\n')
    print(json.dumps(receipt,indent=2));assert passed

if __name__=='__main__':main()
