"""Check replay arrays and the confidence-bound logic without executing circuits."""
import argparse
import json
from pathlib import Path
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import expm_multiply
from analyze import coefficients, tail_bound, choose_degree, interval_decision, mm


def main(reference, repeat, receipt):
    a=np.load(reference/'analysis-arrays.npz')
    b=np.load(repeat/'analysis-arrays.npz')
    assert set(a.files)==set(b.files)
    comparisons={}
    for key in a.files:
        x,y=a[key],b[key]
        assert x.shape==y.shape
        if x.dtype.kind in 'iu':
            passed=bool(np.array_equal(x,y))
        else:
            passed=bool(np.all(np.isfinite(x)) and np.all(np.isfinite(y)) and
                        np.allclose(x,y,rtol=1e-11,atol=1e-11))
        comparisons[key]=dict(passed=passed,max_absolute_difference=float(np.max(abs(x-y))))
        assert passed,key
    # The derived normalized physical response has its own strict criterion.
    assert np.max(abs(a['C_exact']-b['C_exact'])) < 1e-10
    altered=a['C_exact'].copy()
    altered[0,0]+=1e-6
    assert np.max(abs(a['C_exact']-altered)) > 1e-10
    # The cheaper nonorthogonal family has its own preparation and reconstruction.
    root=Path(__file__).parent
    primary=np.load(root/'inputs/primary.npz')
    metadata=json.loads((root/'inputs/primary.json').read_text())
    F=np.sqrt(primary['pi'])[:,None]*primary['monomial_centered']
    norms=np.linalg.norm(F,axis=0)
    raw_basis=F/norms
    raw_D=primary['A']*norms[None,:]
    G=len(primary['pi'])
    H=sparse.csr_matrix((primary['H_data'],primary['H_indices'],primary['H_indptr']),shape=(G,G))
    tau=metadata['tau_model_time']
    evolved=expm_multiply(-tau*H,raw_basis,traceA=-tau*H.diagonal().sum())
    raw_delayed=mm(raw_basis.T,evolved)
    raw_covariance=mm(raw_basis.T,raw_basis)
    raw_C=-mm(mm(raw_D,raw_covariance-raw_delayed),raw_D.T)/np.outer(a['scale'],a['scale'])
    raw_error=float(np.max(abs(raw_C-primary['C'])))
    assert raw_error < 1e-10
    r=json.loads((reference/'results.json').read_text())
    z=r['model']['z']
    for K in [0,21,25,34,41,42,47,64,100]:
        actual=max(0.,1-float(coefficients(z,K).sum()))
        assert actual <= tail_bound(z,K)+1e-13
    rho2=r['observables']['rho_max_squared']
    K=choose_degree(z,rho2,.0005)
    assert rho2*tail_bound(z,K) <= .0005
    assert K==0 or rho2*tail_bound(z,K-1) > .0005
    # This test supplies measured means and known reconstruction constants only.
    # It verifies that insufficient precision cannot certify a chosen shortlist.
    means=np.zeros((8,8))
    means[:7,7]=[.08,.07,.06,.05,.04,.01,.005]
    means[7,:7]=means[:7,7]
    common=dict(D=np.eye(8),equilibrium=np.zeros((8,8)),scale=np.ones(8),
                retained_sum=1.,tail=0.,gamma=np.ones(8),rho=np.ones(8),
                receiver=np.array([7]),candidates=np.array([True]*7+[False]),
                canonical=np.arange(1,9))
    tight=interval_decision(means,epsilon=.001,**common)
    loose=interval_decision(means,epsilon=.02,**common)
    noisy=interval_decision(means,epsilon=.001,device_score_bound=.02,**common)
    assert tight['certified'] and tight['top5_canonical']==[1,2,3,4,5]
    assert not loose['certified'] and not noisy['certified']
    receipt.write_text(json.dumps(dict(status='passed',array_comparisons=comparisons,
        independent_normalized_monomial_response_error=raw_error,
        tests=['All replay array shapes/finite entries and values checked',
               'Integer state indices and ranking identities match exactly',
               'Normalized response replay threshold 1e-10 rejects planted 1e-6 change',
               'Cheaper normalized-monomial family independently reproduces the same physical response',
               'Nine numerical Bessel tails lie below the analytic Chernoff bound',
               'Selected fixed degree is minimal under the declared Chernoff criterion',
               'Measured-mean interface certifies separated intervals and rejects overlapping/noisy intervals'],
        tight_interval_test=tight,loose_interval_test=loose,noise_interval_test=noisy),indent=2)+'\n')
    print('Replay and confidence-bound checks passed.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    root=Path(__file__).parent
    p.add_argument('--reference',type=Path,default=root/'results')
    p.add_argument('--repeat',type=Path,default=root/'results-repeat')
    p.add_argument('--receipt',type=Path,default=root/'results/verification.json')
    a=p.parse_args()
    main(a.reference,a.repeat,a.receipt)
