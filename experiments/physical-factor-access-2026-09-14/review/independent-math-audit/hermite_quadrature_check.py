"""Independent degree-1--4 Hermite/gradient factor check on saved triangle.

Coefficients are obtained by deterministic Gaussian quadrature of direct contact
energies, not by normal-ordering tensor formulas. A second, nested quadrature
constructs the conditional snapshot functions and their gradients directly.
No fitting, cloud jobs, or imports from the production/operator implementation.
Writes only into this audit directory. Main execution follows fixture agreement.
"""
from pathlib import Path
import hashlib, itertools, json, math, sys, time, warnings
import numpy as np
from numpy.polynomial.hermite import hermgauss
from scipy.linalg import expm
warnings.filterwarnings('error',category=RuntimeWarning)
HERE=Path(__file__).resolve().parent
BASE=HERE.parents[1]
ORDER=5
RAW_TOL=1e-10
EMBEDDING_TOL=1e-8

def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def mm(a,b):return np.einsum('ij,jk->ik',a,b,optimize=False)
def load(path):
    with np.load(path) as f:return {k:f[k].copy() for k in f.files}

def standard_nodes(d,order=ORDER):
    x,w=hermgauss(order)
    indices=np.array(list(itertools.product(range(order),repeat=d)))
    return np.sqrt(2)*x[indices],np.prod(w[indices]/np.sqrt(np.pi),axis=1)

def normalized_hermites(nodes,alphas):
    # Probabilists' recurrence He_(k+1)=x He_k-k He_(k-1), normalized by sqrt(k!).
    table=np.empty((5,len(nodes),nodes.shape[1]));table[0]=1;table[1]=nodes
    for k in range(1,4):table[k+1]=nodes*table[k]-k*table[k-1]
    for k in range(5):table[k]/=math.sqrt(math.factorial(k))
    return np.stack([np.prod(np.stack([table[a,:,j] for j,a in enumerate(alpha)],axis=1),axis=1) for alpha in alphas],axis=1)

def direct_contacts(q,r0,edges,physical_basis,kappa=1.):
    n=len(r0);d=physical_basis.shape[1];degree=np.bincount(edges.ravel(),minlength=n)
    r=r0[None]+mm(q,physical_basis.T).reshape(-1,n,3)
    values=np.zeros((len(q),n));gradient=np.zeros((len(q),n,d))
    for i,j in edges:
        rest2=np.sum((r0[i]-r0[j])**2);vector=r[:,i]-r[:,j]
        change=np.sum(vector*vector,axis=1)-rest2;c=kappa/(8*rest2)
        value=c*change**2;T=physical_basis[3*i:3*i+3]-physical_basis[3*j:3*j+3]
        grad=4*c*change[:,None]*mm(vector,T)
        for residue in (i,j):values[:,residue]+=value/degree[residue];gradient[:,residue]+=grad/degree[residue]
    return values,gradient

def quadrature_coefficients(r0,edges,physical_basis,Sigma,sd,beta=1.,mu=1.,kappa=1.,eigenvectors=None):
    d=len(Sigma);nodes,weights=standard_nodes(d)
    if eigenvectors is None:
        eigenvalues,eigenvectors=np.linalg.eigh(Sigma)
        # Fix coordinate phases explicitly so coefficients are reusable.
        for column in eigenvectors.T:
            if column[np.argmax(abs(column))]<0:column*=-1
    else:
        eigenvalues=np.diag(mm(eigenvectors.T,mm(Sigma,eigenvectors)))
        assert np.max(abs(mm(Sigma,eigenvectors)-eigenvectors*eigenvalues))<1e-10
    assert eigenvalues.min()>0
    root=eigenvectors*np.sqrt(eigenvalues)
    q=mm(nodes,root.T);values,grad=direct_contacts(q,r0,edges,physical_basis,kappa)
    means=np.einsum('i,ij->j',weights,values,optimize=False)
    f=(values-means)/sd
    alphas=np.array(sorted([a for a in itertools.product(range(5),repeat=d) if 1<=sum(a)<=4],key=lambda a:(sum(a),a)),dtype=int)
    basis_values=normalized_hermites(nodes,alphas)
    h=mm(basis_values.T,weights[:,None]*f)
    rates=np.einsum('ij,j->i',alphas,mu/(beta*eigenvalues),optimize=False)
    return dict(alphas=alphas,coefficients=h,rates=rates,nodes=nodes,weights=weights,
                physical_q=q,eigenvalues=eigenvalues,eigenvectors=eigenvectors,root=root,
                means=means,normalized_values=f,basis_values=basis_values,
                polynomial_reconstruction_error=float(abs(mm(basis_values,h)-f).max()),
                normalized_basis_gram_error=float(abs(mm(basis_values.T,weights[:,None]*basis_values)-np.eye(len(alphas))).max()))

def main():
    started=time.perf_counter()
    control_path=BASE/'nonlinear/results/controls.npz'
    snapshot_path=BASE/'nonlinear/inputs/frozen-gaussian-result.npz'
    control,snapshot=load(control_path),load(snapshot_path)
    r0,edges,B0,Sigma,sd=(control[k] for k in ['r0','edges','basis','Sigma_Gaussian','sd'])
    assert np.array_equal(Sigma,snapshot['Sigma'])
    assert np.max(abs(sd-snapshot['sd']))<1e-12
    beta=mu=kappa=1.;tau=1/control['eigenvalues'][0]
    ratios=np.array([0.,.1,1.]);times=ratios*tau
    assert abs(snapshot['times'][0]-times[1])<1e-12 and abs(snapshot['times'][1]-times[2])<1e-12
    hdata=quadrature_coefficients(r0,edges,B0,Sigma,sd,beta,mu,kappa)
    h,rates=hdata['coefficients'],hdata['rates'];alpha=hdata['alphas'];weights=hdata['weights']
    assert len(alpha)==34
    Z=np.concatenate([h*np.exp(-t*rates)[:,None] for t in times],axis=1)
    W=snapshot['W'];Psi=mm(Z,W);factor=np.sqrt(rates)[:,None]*Psi
    mass=mm(Psi.T,Psi);Hr=mm(factor.T,factor);observables=mm(Psi.T,h)
    # Independent exact conditional Gaussian quadrature of direct contact gradients.
    q0=hdata['physical_q'];nodes=hdata['nodes'];root=hdata['root'];v=hdata['eigenvalues'];Q=hdata['eigenvectors']
    snapshot_values=[];snapshot_gradients=[]
    for t in times:
        decay=np.exp(-mu*t/(beta*v));P=mm(Q*decay,Q.T)
        innovation=mm(nodes,(root*np.sqrt(-np.expm1(-2*mu*t/(beta*v)))).T)
        qfuture=mm(q0,P.T)[:,None,:]+innovation[None,:,:]
        energy,grad=direct_contacts(qfuture.reshape(-1,3),r0,edges,B0,kappa)
        energy=energy.reshape(len(q0),len(nodes),3);grad=grad.reshape(len(q0),len(nodes),3,3)
        expected=np.einsum('j,ijk->ik',weights,energy,optimize=False)
        expected_grad=np.einsum('j,ijab->iab',weights,grad,optimize=False)
        # Chain derivative of qfuture with respect to q0.
        initial_grad=np.einsum('iab,bc->iac',expected_grad,P,optimize=False)
        snapshot_values.append((expected-hdata['means'])/sd)
        snapshot_gradients.append(initial_grad/sd[None,:,None])
    direct_values=np.concatenate(snapshot_values,axis=1)
    direct_gradients=np.concatenate(snapshot_gradients,axis=1)
    M=mm(direct_values.T,weights[:,None]*direct_values)
    A=mu/beta*np.einsum('i,ijk,ilk->jl',weights,direct_gradients,direct_gradients,optimize=False)
    D=mm(direct_values.T,weights[:,None]*hdata['normalized_values'])
    spectral_M=mm(Z.T,Z);spectral_A=mm(Z.T,rates[:,None]*Z);spectral_D=mm(Z.T,h)
    response=[];difference=[]
    static=mm(observables.T,observables)
    for t,Csaved in zip(snapshot['times'],snapshot['C']):
        C=mm(observables.T,mm(expm(-t*Hr)-np.eye(len(Hr)),observables))
        response.append(C);difference.append(float(abs(C-Csaved).max()))
    checks={'normalized_Hermite_Gram':hdata['normalized_basis_gram_error'],
            'pointwise_original_contact_reconstruction':hdata['polynomial_reconstruction_error'],
            'pointwise_snapshot_reconstruction':float(abs(mm(hdata['basis_values'],Z)-direct_values).max()),
            'direct_GH_M_vs_saved':float(abs(M-snapshot['M']).max()),
            'direct_gradient_GH_tauA_vs_saved':float(tau*abs(A-snapshot['A']).max()),
            'direct_GH_D_vs_saved':float(abs(D-snapshot['D']).max()),
            'Hermite_M_vs_direct_GH':float(abs(spectral_M-M).max()),
            'Hermite_tauA_vs_direct_gradient_GH':float(tau*abs(spectral_A-A).max()),
            'Hermite_D_vs_direct_GH':float(abs(spectral_D-D).max()),
            'Psi_Gram_vs_identity':float(abs(mass-np.eye(len(Hr))).max()),
            'factor_tauHr_vs_saved':float(tau*abs(Hr-snapshot['Hr']).max()),
            'PsiT_h_vs_saved_B':float(abs(observables-snapshot['B']).max()),
            'response_vs_saved_three_times':max(difference),
            'positive_local_factor_min_rate':float(rates.min())}
    rawkeys=['normalized_Hermite_Gram','pointwise_original_contact_reconstruction','pointwise_snapshot_reconstruction','direct_GH_M_vs_saved','direct_gradient_GH_tauA_vs_saved','direct_GH_D_vs_saved','Hermite_M_vs_direct_GH','Hermite_tauA_vs_direct_gradient_GH','Hermite_D_vs_direct_GH']
    embedkeys=['Psi_Gram_vs_identity','factor_tauHr_vs_saved','PsiT_h_vs_saved_B','response_vs_saved_three_times']
    passed=all(checks[k]<=RAW_TOL for k in rawkeys) and all(checks[k]<=EMBEDDING_TOL for k in embedkeys)
    output=HERE/'small-Hermite-factor.npz'
    np.savez_compressed(output,alphas=alpha,h=h,rates=rates,Psi=Psi,factor=factor,M_quadrature=M,A_gradient_quadrature=A,D_quadrature=D,
                        snapshot_coefficients=Z,coordinate_eigenvalues=v,coordinate_eigenvectors=Q,Hr_from_factor=Hr,B_from_embedding=observables,C_from_factor=np.array(response),snapshot_times=times)
    receipt={'status':'PASS' if passed else 'FAIL','scope':'Complete34-chaos-term check of existing9D Gaussian triangle snapshot model only; no full ABL factor, nonlinear fidelity, efficient input oracle, or QPU claim',
             'basis':'Normalized probabilists Hermites He_alpha/sqrt(alpha!)','degrees':[1,2,3,4],'basis_size':len(alpha),
             'coefficient_method':'Order5 GH integration of direct quartic contact energies; no normal-ordering tensor formula',
             'gradient_Gram_method':'Nested125x125 GH conditional averages of direct contact gradients; Dirichlet coefficient mu/beta',
             'coordinate_phases':'Largest absolute component of each covariance eigenvector is positive','quadrature_order':ORDER,
             'outer_nodes':125,'conditional_innovation_nodes':125,'randomness':'None','snapshot_times_over_tau':ratios.tolist(),
             'constants':{'beta':beta,'mu':mu,'kappa':kappa,'tau':float(tau)},'raw_tolerance':RAW_TOL,'embedding_tolerance':EMBEDDING_TOL,
             'checks':checks,'degree1_coefficient_norm':float(np.linalg.norm(h[np.sum(alpha,axis=1)==1])),
             'charged_preprocessing':'Saved covariance fit, original physical basis, coordinate covariance diagonalization, and saved snapshot whitening are inputs/preprocessing, not free oracles.',
             'source_hashes':{'script':digest(__file__),'controls':digest(control_path),'snapshot_model':digest(snapshot_path),'output':digest(output)},
             'seconds':time.perf_counter()-started,'python':sys.version,'numpy':np.__version__}
    (HERE/'small-Hermite-factor-check.json').write_text(json.dumps(receipt,indent=2,allow_nan=False)+'\n')
    print(json.dumps(receipt,indent=2));assert passed

if __name__=='__main__':main()
