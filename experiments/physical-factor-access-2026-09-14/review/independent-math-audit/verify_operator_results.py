"""Independent bounded audit of saved ABL force and selected Hermite outputs.

No imports from operator production code. Local force blocks are integrated from
the Hessian of the original quartic contact energy. Selected Hermite coefficients
use conditional Gaussian contact moments and low-dimensional GH quadrature, not
normal-ordered coefficient tensors. All writes remain in this audit directory.
"""
from pathlib import Path
import json, math, sys, time
import numpy as np
from scipy.linalg import eigh, solve
import hermite_quadrature_check as gh

HERE=Path(__file__).resolve().parent
BASE=HERE.parents[1]
OP=BASE/'operator'

def main():
    start=time.perf_counter()
    protocol=json.loads((OP/'protocol.json').read_text())
    for filename,spec in protocol['input']['files'].items():
        assert gh.digest(OP/'inputs'/filename)==spec['sha256']
    paths={name:OP/'results'/name for name in ['physical-factor.npz','selected-hermite-rows.npz','matched-responses.npz','small-physical-model.npz']}
    model=gh.load(OP/'inputs/static-model.npz');saved=gh.load(OP/'inputs/operator.npz')
    Sigma=gh.load(OP/'inputs/fitted-covariance.npz')['covariance']
    sd=gh.load(OP/'inputs/all-mode-harmonic-scales.npz')['harmonic_sd']
    force=gh.load(paths['physical-factor.npz']);selected=gh.load(paths['selected-hermite-rows.npz'])
    response=gh.load(paths['matched-responses.npz']);small=gh.load(paths['small-physical-model.npz'])
    reported=json.loads((OP/'results/classical-summary.json').read_text())
    access=json.loads((OP/'results/hermite-access.json').read_text())
    N,d,rank=len(sd),len(Sigma),len(saved['Hr']);r0=model['r0'];edges=model['edges'];B0=model['B']
    i,j=edges.T;T=B0.reshape(N,3,d)[i]-B0.reshape(N,3,d)[j]
    a=r0[i]-r0[j];length2=np.sum(a*a,axis=1);degree=np.bincount(edges.ravel(),minlength=N)
    eig,Q=eigh(Sigma);A=np.einsum('eij,jk->eik',T,Q,optimize=True)*np.sqrt(eig)
    S=np.einsum('eik,ejk->eij',A,A,optimize=False)

    # Integrate d²u/dx² = (rr^T + .5*(r²-l²)I)/l², for beta=mu=kappa=1.
    localvals,localvectors=np.linalg.eigh(S)
    assert localvals.min()>-1e-10
    roots=localvectors*np.sqrt(np.maximum(localvals,0))[:,None,:]
    nodes,weights=gh.standard_nodes(3)
    displacement=np.einsum('ab,ecb->eac',nodes,roots,optimize=False)
    vector=a[:,None,:]+displacement;change=np.sum(vector*vector,axis=2)-length2[:,None]
    blocks=(np.einsum('a,eai,eaj->eij',weights,vector,vector,optimize=False)
            +.5*np.einsum('a,ea->e',weights,change,optimize=False)[:,None,None]*np.eye(3))/length2[:,None,None]
    bval,bvec=np.linalg.eigh(blocks)
    bsignroots=(bvec*np.sqrt(bval)[:,None,:])@bvec.transpose(0,2,1)
    F=np.einsum('eij,ejk->eik',bsignroots,T,optimize=False).reshape(-1,d)
    Gamma=F.T@F;Sinv=solve(Sigma,np.eye(d),assume_a='pos')
    mismatch=Gamma-Sinv
    force_checks={
        'GH_expected_Hessian_vs_saved_blocks_maxabs':float(abs(blocks-force['blocks']).max()),
        'GH_local_factor_Gram_vs_saved_Gamma_maxabs':float(abs(Gamma-force['Gamma']).max()),
        'local_expected_Hessian_min_eigenvalue':float(bval.min()),
        'stationarity_Gamma_minus_inverse_Sigma_maxabs':float(abs(mismatch).max()),
        'stationarity_relative_Frobenius':float(np.linalg.norm(mismatch)/np.linalg.norm(Sinv)),
        'stationarity_spectral_norm':float(abs(eigh((mismatch+mismatch.T)/2,eigvals_only=True)).max()),
        'saved_factor_covariance_identity_maxabs':float(abs(force['Gamma']@force['Sigma_physical']-np.eye(d)).max()),
        'dense_projected_factor_Frobenius_squared':float(np.sum(F*F)),
    }

    # The conditional residual xi is Gaussian in R^3. For x=m+xi, R=Cov(xi),
    # E[(2a.x+x²)² | selected y] = (2a.m+m²+trR)²
    #                              +4(a+m)^T R(a+m)+2 tr(R²).
    # Integrating this polynomial times He_alpha/sqrt(alpha!) gives h_alpha.
    expected_h=[];degrees=[]
    for row in selected['indices']:
        ix=row[row>=0];unique,counts=np.unique(ix,return_counts=True)
        z,w=gh.standard_nodes(len(unique));As=A[:,:,unique]
        m=np.einsum('na,eia->eni',z,As,optimize=False)
        residual=S-np.einsum('eia,eja->eij',As,As,optimize=False)
        trace=np.trace(residual,axis1=1,axis2=2)
        delta=2*np.einsum('ei,eni->en',a,m,optimize=False)+np.sum(m*m,axis=2)
        shifted=a[:,None,:]+m
        conditional=(delta+trace[:,None])**2+4*np.einsum('eni,eij,enj->en',shifted,residual,shifted,optimize=False)+2*np.einsum('eij,eji->e',residual,residual,optimize=False)[:,None]
        basis=gh.normalized_hermites(z,counts[None,:])[:,0]
        coefficient=np.einsum('n,en->e',w*basis,conditional,optimize=False)/(8*length2)
        h=np.bincount(i,weights=coefficient,minlength=N)+np.bincount(j,weights=coefficient,minlength=N)
        expected_h.append(h/(degree*sd));degrees.append(len(ix))
    expected_h=np.array(expected_h);degrees=np.array(degrees)
    omega=np.array([sum(1/eig[row[row>=0]]) for row in selected['indices']])
    nodes=selected['nodes'];Z=np.concatenate([expected_h*np.exp(-t*omega)[:,None] for t in nodes],axis=1)
    expected_psi=Z@saved['W'];expected_factor=np.sqrt(omega)[:,None]*expected_psi
    coefficient_checks={
        'selected64_coefficients_maxabs':float(abs(expected_h-selected['coefficients']).max()),
        'degree_maxabs':{str(k):float(abs(expected_h[degrees==k]-selected['coefficients'][degrees==k]).max()) for k in range(1,5)},
        'selected64_embedding_rows_maxabs':float(abs(expected_psi-selected['embedding']).max()),
        'selected64_factor_rows_maxabs':float(abs(expected_factor-selected['factor_rows']).max()),
        'selected64_rates_maxabs':float(abs(omega-selected['omega']).max()),
    }
    rng=np.random.default_rng(20260914);expected_indices=[]
    for k in range(1,5):
        rows=[tuple([0]*k),tuple([d-1]*k)]
        while len(rows)<16:
            ix=tuple(sorted(rng.integers(0,d,size=k).tolist()))
            if ix not in rows:rows.append(ix)
        expected_indices.extend(rows)
    assert np.array_equal(selected['indices'],np.array([list(ix)+[-1]*(4-len(ix)) for ix in expected_indices]))
    assert np.max(abs(nodes/float(saved['tau'])-np.array([0,.0001,.001,.01,.1,1.])))<1e-14

    D=math.comb(d+4,4)-1;rowpad=2**math.ceil(math.log2(D));colpad=2**math.ceil(math.log2(rank))
    dimensions={'degree_counts':{str(k):math.comb(d+k-1,k) for k in range(1,5)},
                'full_centered_Hermite_dimension':D,'full_observable_coefficient_float64_bytes':8*D*N,
                'full_snapshot_embedding_float64_bytes':8*D*rank,'full_factor_float64_bytes':8*D*rank,
                'padded_factor_rows':rowpad,'padded_factor_columns':colpad,
                'generic_table_loader_qubits':rowpad.bit_length()-1+colpad.bit_length()-1+1,
                'generic_loader_Ry_angles_per_UL':rowpad*colpad,
                'generic_loader_table_float64_bytes':8*rowpad*colpad,
                'full_output_scalar_count_at_four_times':N*N*4,'symmetric_scalar_count_at_four_times':N*(N+1)//2*4}
    assert all(access[k]==v for k,v in dimensions.items())
    assert access['dense_whitening_bytes']==saved['W'].nbytes
    assert np.isclose(access['actual_saved_reduced_factor_Frobenius_squared'],np.trace(saved['Hr']))
    assert access['sequential_full_row_time_extrapolation_seconds']==float(np.median(access['selected_coefficient_and_embedding_seconds'])*D)
    costs={k:dimensions[k] for k in dimensions if not isinstance(dimensions[k],dict)}
    costs['basis_dimension_verified']=D
    costs['generic_loader_scope']='Fixed explicit table implementation only; not an algorithm-independent lower bound or competitive full response cost.'

    # Recompute complete saved reduced responses independently from dense Hr,
    # and verify the reporting arithmetic on every saved output array.
    lam,U=eigh(saved['Hr']);projected=U.T@saved['B']
    Cs=np.array([projected.T@(np.exp(-t*lam)[:,None]*projected)-saved['G0'] for t in response['times']])
    response_checks={
       'recomputed_spectral_vs_saved_spectral_maxabs':float(abs(Cs-response['reduced_spectral']).max()),
       'matrix_free_vs_spectral_maxabs':float(abs(response['reduced_matrix_free']-response['reduced_spectral']).max()),
       'reduced_vs_full_Gaussian_maxabs':float(abs(response['reduced_spectral']-response['full_Gaussian']).max()),
       'factor_consistent_vs_original_Gaussian_maxabs':float(abs(response['factor_consistent_Gaussian']-response['full_Gaussian']).max()),
       'Gram_static_vs_saved_G0_maxabs':float(abs(saved['B'].T@saved['B']-saved['G0']).max()),
    }
    for key,reported_key in [('matrix_free_vs_spectral_maxabs','matrix_free_error_maxabs'),('reduced_vs_full_Gaussian_maxabs','reduction_error_maxabs_vs_full_Gaussian'),('factor_consistent_vs_original_Gaussian_maxabs','factor_covariance_response_change_maxabs')]:
        assert response_checks[key]==reported[reported_key]
    # Small actual-ABL fixture: fully direct GH covers all14 terms, independently
    # of operator tensor formulas and without re-fitting its covariance.
    smallgh=gh.quadrature_coefficients(small['r0'],small['edges'],small['B'],small['Sigma'],small['sd'],eigenvectors=small['covariance_vectors'])
    coefficient_lookup={tuple(np.repeat(np.arange(2),a)):h for a,h in zip(smallgh['alphas'],smallgh['coefficients'])}
    smallh=np.array([coefficient_lookup[tuple(row[row>=0])] for row in small['occupations']])
    small_checks={'all14_direct_GH_coefficients_maxabs':float(abs(smallh-small['h']).max()),
                  'force_Gram_maxabs':float(abs(small['force_factor'].T@small['force_factor']-small['force_Gamma']).max()),
                  'source_coordinates_match_first_three_ABL_residues':bool(np.array_equal(small['r0'],r0[:3])),
                  'source_contacts_match_inherited_ABL':bool(np.array_equal(small['edges'],edges[(edges<3).all(axis=1)]))}
    tolerances={'force_absolute':1e-10,'coefficient_absolute':1e-10,'embedding_absolute':1e-7,'saved_response_recompute':1e-10}
    passed=(force_checks['GH_expected_Hessian_vs_saved_blocks_maxabs']<1e-10
      and force_checks['GH_local_factor_Gram_vs_saved_Gamma_maxabs']<1e-10
      and force_checks['stationarity_Gamma_minus_inverse_Sigma_maxabs']<1e-10
      and coefficient_checks['selected64_coefficients_maxabs']<1e-10
      and coefficient_checks['selected64_factor_rows_maxabs']<1e-7
      and response_checks['recomputed_spectral_vs_saved_spectral_maxabs']<1e-10
      and small_checks['all14_direct_GH_coefficients_maxabs']<1e-10
      and small_checks['source_coordinates_match_first_three_ABL_residues']
      and small_checks['source_contacts_match_inherited_ABL'])
    checked_paths=[__file__,gh.__file__,OP/'protocol.json',OP/'protocol-amendment-1.json',
       OP/'source/core.py',OP/'source/run_classical.py',OP/'results/classical-summary.json',
       OP/'results/hermite-access.json',*paths.values()]
    receipt={'status':'PASS' if passed else 'FAIL','scope':'Independent numeric checks of force blocks, selected64 ABL coefficient rows, saved output arithmetic and complete14-term actual-ABL submodel. Full13.36-billion-term embedding identity and full quantum access remain unverified/unimplemented.',
       'methods':{'force':'Order5 GH of direct local quartic-energy Hessian in3 variables',
                  'coefficients':'Order5 GH over1..4 selected standard-normal modes after exact conditional Gaussian integration of other coordinates; no normal-ordered tensor formula',
                  'small_fixture':'Direct quartic energy GH with all2 coordinates and all14 Hermite coefficients',
                  'response':'Fresh eigendecomposition of saved Hr and complete matrix outputs'},
       'randomness':'No integration randomness; only prescribed seed20260914 used to verify selected row identities.',
       'dimensions':dimensions,'tolerances':tolerances,'force_checks':force_checks,'coefficient_checks':coefficient_checks,
       'response_checks':response_checks,'small_ABL_checks':small_checks,
       'source_hashes':{str(path):gh.digest(path) for path in checked_paths},
       'input_hashes':{filename:gh.digest(OP/'inputs'/filename) for filename in protocol['input']['files']},
       'seconds':time.perf_counter()-start,'python':sys.version,'numpy':np.__version__}
    (HERE/'operator-results-check.json').write_text(json.dumps(receipt,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:v for k,v in receipt.items() if k not in ['source_hashes','input_hashes']},indent=2));assert passed

if __name__=='__main__':main()
