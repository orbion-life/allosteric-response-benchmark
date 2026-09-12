"""Observable-seeded response-operator canary with full finite-grid references.

All physical coordinates and original quartic node observables are retained.
This is a finite numerical model study, not a protein or hardware validation.
"""
from pathlib import Path
import argparse, datetime, hashlib, json, resource, sys, time
import numpy as np
from scipy.linalg import eigh, svd
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import expm_multiply
from scipy.special import expit, logsumexp
import reference_source as reference
ROOT=Path(__file__).resolve().parent

def digest(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(p,obj): Path(p).write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n')
def mm(a,b): return np.einsum('ij,jk->ik',a,b,optimize=False)
def rss(): return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform=='darwin' else 1024)

def build(protocol,geometry,kappa,n,extent):
    """Reproduce the v0.3.0 full-grid energy/rates, exposing H and observables."""
    started=time.perf_counter()
    local=dict(protocol,geometry_A=geometry['coordinates_A'])
    r0,edges,lam,B,sd=reference.static_model(local,kappa)
    axis=np.linspace(-extent/np.sqrt(lam[0]),extent/np.sqrt(lam[0]),n)
    delta=axis[1]-axis[0];inds=np.indices((n,n,n)).reshape(3,-1).T
    q=axis[inds];G=len(q);rr=r0[None]+mm(q,B.T).reshape(G,3,3)
    Ei=np.zeros((G,3));energy=np.zeros(G);degree=np.bincount(edges.ravel(),minlength=3)
    for i,j in edges:
        rest2=np.sum((r0[i]-r0[j])**2)
        length2=np.sum((rr[:,i]-rr[:,j])**2,axis=1)
        u=kappa*(length2-rest2)**2/(8*rest2);energy+=u
        Ei[:,i]+=u/degree[i];Ei[:,j]+=u/degree[j]
    pi=np.exp(-energy-logsumexp(-energy));sp=np.sqrt(pi)
    means=np.einsum('i,ij->j',pi,Ei,optimize=False)
    F=sp[:,None]*(Ei-means)/sd[None,:]
    row=[];col=[];data=[];diag=np.zeros(G);flat=np.arange(G).reshape((n,n,n))
    for dim in range(3):
        lo=[slice(None)]*3;hi=[slice(None)]*3;lo[dim]=slice(0,-1);hi[dim]=slice(1,None)
        x=flat[tuple(lo)].ravel();y=flat[tuple(hi)].ravel()
        xy=2/delta**2*expit(-(energy[y]-energy[x]));yx=2/delta**2*expit(energy[y]-energy[x])
        h=-np.sqrt(xy*yx);row.extend([x,y]);col.extend([y,x]);data.extend([h,h])
        np.add.at(diag,x,xy);np.add.at(diag,y,yx)
    row.append(np.arange(G));col.append(np.arange(G));data.append(diag)
    H=coo_matrix((np.concatenate(data),(np.concatenate(row),np.concatenate(col))),shape=(G,G)).tocsr()
    times=np.array(protocol['times_over_tau'])/lam[0]
    meta={'common_setup_seconds':time.perf_counter()-started,'states':G,'nnz':H.nnz,
          'stationarity_residual':float(np.max(abs(H@sp))),
          'centered_observable_residual':float(np.max(abs(np.einsum('i,ij->j',sp,F,optimize=False)))),
          'boundary_probability':float(pi[np.any((inds==0)|(inds==n-1),axis=1)].sum())}
    assert meta['stationarity_residual']<1e-9
    assert np.isfinite(F).all() and np.all(sd>0)
    return H,F,times,sd,meta

def kernels(eigenvalues,t):
    z=eigenvalues[:,None]+eigenvalues[None,:];x=t*z;small=abs(x)<1e-4
    a=np.empty_like(x);b=np.empty_like(x)
    y=x[small]
    a[small]=t*(1-y/2+y*y/6-y**3/24+y**4/120-y**5/720)
    b[small]=t*t*(.5-y/6+y*y/24-y**3/120+y**4/720-y**5/5040)
    a[~small]=-np.expm1(-x[~small])/z[~small]
    b[~small]=(x[~small]+np.expm1(-x[~small]))/(z[~small]**2)
    return a,b

def positive_quadratic_upper(matrix,vector):
    """Conservative summation allowance; not directed-rounding certification."""
    terms=np.outer(vector,vector)*matrix
    value=float(np.sum(terms,dtype=np.longdouble));mass=float(np.sum(abs(terms),dtype=np.longdouble))
    allowance=128*np.finfo(float).eps*max(mass,np.finfo(float).tiny)
    if value < -max(allowance,1e-12*mass):
        raise ArithmeticError('Residual integral lost positive semidefiniteness.')
    return max(0.,value)+allowance

def projected(H,F,V,times):
    started=time.perf_counter();HV=H@V;Hr=mm(V.T,HV);Hr=(Hr+Hr.T)/2
    lam,U=eigh(Hr);negative=max(0.,-float(lam[0]));scale=max(1.,float(np.max(abs(lam))))
    if negative>1e-10*scale:raise ArithmeticError('Projected operator is not positive semidefinite.')
    # Tiny negative Ritz values are numerical. Account for the spectral clipping
    # separately; the exact Galerkin operator is positive semidefinite.
    lam_used=np.maximum(lam,0.)
    coeff=mm(V.T,F);rot=mm(U.T,coeff);P=mm(V,coeff);defect=F-P
    norms=np.linalg.norm(F,axis=0);pnorms=np.linalg.norm(P,axis=0);dnorms=np.linalg.norm(defect,axis=0)
    initial=np.outer(dnorms,norms)+np.outer(pnorms,dnorms);initial=np.minimum(initial,initial.T)
    gram=mm(F.T,F);reduced_gram=mm(coeff.T,coeff)
    static_error=abs(gram-reduced_gram)
    residual=HV-mm(V,Hr);Rgram=mm(residual.T,residual);Rgram=(Rgram+Rgram.T)/2
    coupling=mm(mm(U.T,Rgram),U);coupling=(coupling+coupling.T)/2
    normR=np.sqrt(max(0.,float(eigh(Rgram,eigvals_only=True)[-1])))
    # Numerical nonorthogonality is measured. A first-order integral term covers
    # any nonzero V^T R before using the two-sided residual inequality.
    ortho=float(np.max(abs(mm(V.T,V)-np.eye(V.shape[1]))))
    residual_parallel=float(np.linalg.norm(mm(V.T,residual),ord=2))
    delayed=[];one_bounds=[];two_bounds=[];bounds=[];uniform=[]
    for t in times:
        K=mm(rot.T*np.exp(-t*lam_used),rot);delayed.append(K)
        phi,psi=kernels(lam_used,float(t));I=[];J=[]
        for j in range(F.shape[1]):
            I.append(positive_quadratic_upper(coupling*phi,rot[:,j]))
            J.append(positive_quadratic_upper(coupling*psi,rot[:,j]))
        one=np.outer(norms,np.sqrt(t*np.array(I)));one=np.minimum(one,one.T)
        two=np.sqrt(np.outer(J,J))+t*residual_parallel*np.outer(pnorms,pnorms)
        numerical=t*negative*np.outer(pnorms,pnorms)
        # Initial-span error affects delayed covariance. The directly calculated
        # Gram discrepancy accounts for the static term in response reconstruction.
        response_bound=np.minimum(one,two)+initial+static_error+numerical
        response_bound+=256*np.finfo(float).eps*(1+t*scale)*np.outer(norms,norms)
        one_bounds.append(one+initial+static_error+numerical)
        two_bounds.append(two+initial+static_error+numerical)
        uniform.append(t*normR*np.outer(norms,norms)+initial+static_error+numerical)
        bounds.append(response_bound)
    delayed=np.array(delayed);C=delayed-reduced_gram
    row={'rank':V.shape[1],'orthogonality_error':ortho,'seed_projection_norm':float(np.linalg.norm(defect)),
         'static_max_error':float(static_error.max()),'min_projected_eigenvalue':float(lam[0]),
         'clipped_eigenvalue_magnitude':negative,'residual_spectral_norm':float(normR),
         'residual_parallel_norm':residual_parallel,'max_response_bound':float(np.max(bounds)),
         'projection_query_seconds':time.perf_counter()-started}
    arrays={'C':C,'delayed':delayed,'equilibrium':reduced_gram,'bound':np.array(bounds),
            'one_sided_bound':np.array(one_bounds),'two_sided_bound':np.array(two_bounds),
            'uniform_bound':np.array(uniform),'Hr':Hr,'coefficients':coeff,'residual_Gram':Rgram}
    return arrays,row

def block_basis(H,F,protocol,times):
    started=time.perf_counter();checks=protocol['basis']['rank_checkpoints'];maximum=protocol['basis']['maximum_rank']
    cutoff=protocol['basis']['relative_dependency_threshold']
    left,s,_=svd(F,full_matrices=False);keep=s>cutoff*s[0];V=left[:,keep]
    if V.shape[1]==0:raise ArithmeticError('No observable seed.')
    frontier=V.copy();records=[];snapshots=[];selected=False;max_time=protocol['resources']['max_total_wall_seconds']
    while True:
        rank=V.shape[1]
        if rank in checks or rank==maximum:
            arrays,row=projected(H,F,V,times);row['operator_seconds_to_checkpoint']=time.perf_counter()-started
            row['bound_stop_pass']=bool(row['max_response_bound']<=protocol['gates']['calculated_response_error_bound'] and row['static_max_error']<=protocol['gates']['normalized_static_error'] and row['orthogonality_error']<=protocol['gates']['basis_orthogonality'])
            records.append(row);snapshots.append(arrays)
            if row['bound_stop_pass']:selected=True;break
        if rank>=maximum:break
        block=H@frontier;original_scale=float(np.linalg.norm(block,ord=2))
        for _ in range(2):block-=mm(V,mm(V.T,block))
        left,s,_=svd(block,full_matrices=False);keep=s>cutoff*max(original_scale,np.finfo(float).tiny)
        left=left[:,keep]
        if left.shape[1]==0:
            if not records or records[-1]['rank']!=rank:
                arrays,row=projected(H,F,V,times);row['operator_seconds_to_checkpoint']=time.perf_counter()-started;row['bound_stop_pass']=bool(row['max_response_bound']<=protocol['gates']['calculated_response_error_bound'] and row['static_max_error']<=protocol['gates']['normalized_static_error'] and row['orthogonality_error']<=protocol['gates']['basis_orthogonality']);records.append(row);snapshots.append(arrays)
            selected=records[-1]['bound_stop_pass'];break
        # Preserve checkpoint boundaries if a rank-deficient block changes size.
        next_check=min([c for c in checks if c>rank]+[maximum]);left=left[:,:min(left.shape[1],next_check-rank,maximum-rank)]
        for column in left.T:
            if column[np.argmax(abs(column))]<0:column*=-1
        V=np.column_stack((V,left));frontier=left
        if time.perf_counter()-started>max_time:raise TimeoutError('Basis construction exceeded total allocation.')
    return snapshots,records,selected,time.perf_counter()-started

def full_reference(H,F,times):
    started=time.perf_counter();equilibrium=mm(F.T,F);delayed=[]
    for t in times:
        propagated=expm_multiply(-t*H,F,traceA=-t*H.diagonal().sum())
        delayed.append(mm(F.T,propagated))
    delayed=np.array(delayed)
    return {'C':delayed-equilibrium,'delayed':delayed,'equilibrium':equilibrium},time.perf_counter()-started

def case(protocol,geometry,kappa,n,output,deadline,extent=4,reference_only=False):
    name=f"{geometry['id']}-k{kappa}-n{n}-e{extent}";started=time.perf_counter();cpu=time.process_time()
    H,F,times,sd,meta=build(protocol,geometry,kappa,n,extent)
    if not reference_only:
        snapshots,records,selected,operator_seconds=block_basis(H,F,protocol,times)
    exact,reference_seconds=full_reference(H,F,times)
    result=dict(geometry=geometry['id'],kappa=kappa,n=n,extent=extent,times=times.tolist(),reference_only=reference_only,**meta,reference_seconds=reference_seconds)
    values={f'reference_{k}':v for k,v in exact.items()};values.update(times=times,harmonic_sd=sd)
    if not reference_only:
        for index,(approx,row) in enumerate(zip(snapshots,records)):
            err=abs(approx['C']-exact['C']);delayerr=abs(approx['delayed']-exact['delayed'])
            row.update(response_max_error=float(err.max()),delayed_max_error=float(delayerr.max()),minimum_bound_slack=float(np.min(approx['bound']-err)),
                       actual_error_within_calculated_bound=bool(np.all(err<=approx['bound']+1e-9)),dimension_reduction=float(H.shape[0]/row['rank']))
            for k,v in approx.items():values[f'rank{row["rank"]}_{k}']=v
        last=records[-1];total_operator=meta['common_setup_seconds']+operator_seconds;total_reference=meta['common_setup_seconds']+reference_seconds
        # The operational choice uses the bound only. Exact errors are revealed
        # here, after construction, solely for validation.
        finite_pass=bool(selected and last['response_max_error']<=protocol['gates']['normalized_response_error'] and last['delayed_max_error']<=protocol['gates']['normalized_delayed_error'] and last['actual_error_within_calculated_bound'])
        ratio=total_operator/total_reference
        result.update(checkpoints=records,selected_by_bound=selected,selected_rank=last['rank'],operator_seconds=operator_seconds,
                      total_operator_seconds=total_operator,total_reference_seconds=total_reference,total_time_ratio=ratio,
                      finite_fidelity_pass=finite_pass,practical_pass=bool(finite_pass and last['dimension_reduction']>=4 and ratio<=1),
                      scope='Finite-grid response only; no biological inference or continuum accuracy claim.')
    result.update(case_seconds=time.perf_counter()-started,cpu_seconds=time.process_time()-cpu,peak_process_bytes=rss())
    assert result['peak_process_bytes']<=protocol['resources']['max_aggregate_bytes']
    np.savez_compressed(output/(name+'.npz'),**values);dump(output/(name+'.json'),result)
    print(json.dumps({k:result[k] for k in ['geometry','kappa','n','reference_only','case_seconds']},sort_keys=True),flush=True)
    if time.perf_counter()>deadline:raise TimeoutError('Total experiment allocation exhausted.')
    return result

def run(output):
    output.mkdir(parents=True,exist_ok=False)
    protocol=json.loads((ROOT/'preanalysis.json').read_text());assert digest(ROOT/'preanalysis.json')==(ROOT/'preanalysis.sha256').read_text().strip()
    for name,h in json.loads((ROOT/'source-receipts.json').read_text()).items():assert digest(ROOT/name)==h,name
    receipt={'status':'RUNNING','started_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'protocol_sha256':digest(ROOT/'preanalysis.json'),'script_sha256':digest(__file__),'cases':[]}
    started=time.perf_counter();deadline=started+protocol['resources']['max_total_wall_seconds'];dump(output/'run.json',receipt)
    try:
        for geometry in protocol['geometries']:
            for kappa in protocol['kappa']:
                for n in protocol['grid_n']:
                    receipt['cases'].append(case(protocol,geometry,kappa,n,output,deadline))
                    dump(output/'run.json',receipt)
        for geometry in protocol['geometries']:
            for kappa in protocol['kappa']:
                d=protocol['domain_check'];receipt['cases'].append(case(protocol,geometry,kappa,d['n'],output,deadline,extent=d['extent'],reference_only=True));dump(output/'run.json',receipt)
        receipt['status']='COMPLETE'
    except BaseException as exc:
        receipt.update(status='STOPPED',failure_type=type(exc).__name__,failure=str(exc));raise
    finally:
        receipt.update(wall_seconds=time.perf_counter()-started,cpu_seconds=time.process_time(),peak_process_bytes=rss(),finished_utc=datetime.datetime.now(datetime.timezone.utc).isoformat());dump(output/'run.json',receipt)

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,required=True);args=parser.parse_args();run(args.output)
