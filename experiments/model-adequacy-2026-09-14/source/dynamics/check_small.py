"""Independent finite-matrix and field-response implementation checks."""
from core import *
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import expm_multiply
from scipy.special import expit


def csr(grid):
    n=grid.G; flat=np.arange(n).reshape(grid.shape); V=host(grid.V)
    diag=np.zeros(n);rows=[];cols=[];vals=[]
    for k in range(3):
        lo=[slice(None)]*3;hi=[slice(None)]*3;lo[k]=slice(0,-1);hi[k]=slice(1,None)
        aa=flat[tuple(lo)].ravel();bb=flat[tuple(hi)].ravel();delta=V[bb]-V[aa]
        base=2*grid.model['mu']*grid.model['eigenvalues'][k]/grid.h**2
        ab=base*expit(-delta);ba=base*expit(delta);off=-np.sqrt(ab*ba)
        np.add.at(diag,aa,ab);np.add.at(diag,bb,ba)
        rows.extend([aa,bb]);cols.extend([bb,aa]);vals.extend([off,off])
    rows.append(np.arange(n));cols.append(np.arange(n));vals.append(diag)
    return coo_matrix((np.concatenate(vals),(np.concatenate(rows),np.concatenate(cols))),shape=(n,n)).tocsr()


def check(output,backend_name='numpy'):
    out=Path(output);out.mkdir(parents=True,exist_ok=False);start=time.perf_counter();xp=backend(backend_name)
    g=Grid([4,4,18],2,xp);H=csr(g);F=host(g.F)
    action_error=float(np.max(abs(host(g.apply(g.F))- (H@F.T).T)))
    records=[]; arrays={'G0':g.G0}
    for tt in [.0075,.1,1,10]:
        t=tt*g.tau;ref=F@expm_multiply(-t*H,F.T,traceA=-t*H.diagonal().sum())
        c,cm=g.propagate(t,'chebyshev');u,um=g.propagate(t,'uniformization')
        rec={'time_over_tau':tt,'chebyshev_error':float(np.max(abs(c-ref))),
             'uniformization_error':float(np.max(abs(u-ref))),'chebyshev':cm,'uniformization':um}
        records.append(rec);arrays['t'+str(tt)]=np.stack([ref,c,u])
    field=[]
    for source in range(3):
        eta0=min(.01,.05*g.model['beta']*g.model['degree'][source]*g.normal['sd'][source])
        central=[]
        for divisor in [1,2,4]:
            eta=eta0/divisor;pm=[]
            for sign in [1,-1]:
                gg=Grid([4,4,18],2,xp,(source,sign*eta))
                val,_=gg.propagate(.1*g.tau,left=gg.left,tolerance=1e-12)
                pm.append(val[0])
            central.append((pm[0]-pm[1])/(2*eta))
        target=arrays['t0.1'][0][source]-g.G0[source]
        field.append({'source':source,'eta0':eta0,'successive_changes':np.max(abs(np.diff(central,axis=0)),axis=1).tolist(),
                      'FDT_error_final':float(np.max(abs(central[-1]-target))),'central':np.array(central).tolist(),'reference':target.tolist()})
    # An independent finite-chain Duhamel residual check, with all domains valid.
    V,_=np.linalg.qr(F.T); B=V.T@F.T; Hr=V.T@(H@V);R=H@V-V@Hr;G=R.T@R
    lam,Q=np.linalg.eigh(Hr);times=np.array([.0075,.025,.1,1])*g.tau
    from scipy.integrate import quad
    residual=[]
    for t in times:
        ref=F@expm_multiply(-t*H,F.T,traceA=-t*H.diagonal().sum())
        reduced=B.T@(Q*np.exp(-t*lam))@Q.T@B
        bounds=np.empty((3,3))
        for i in range(3):
            def integrand(s):
                v=Q@(np.exp(-s*lam)*(Q.T@B[:,i]));return float(np.linalg.norm(R@v))
            val,err=quad(integrand,0,t,epsabs=1e-11,epsrel=1e-10)
            bounds[:,i]=np.linalg.norm(F,axis=1)*(val+err)
        residual.append({'time_over_tau':t/g.tau,'max_error':float(np.max(abs(ref-reduced))),
                         'max_bound':float(bounds.max()),'all_errors_below_computed_bound':bool(np.all(abs(ref-reduced)<=bounds+1e-10))})
    np.savez_compressed(out/'raw-checks.npz',**arrays,residual_Gram=G,reduced_H=Hr,B=B)
    passed=action_error<=1e-12 and max(max(x['chebyshev_error'],x['uniformization_error']) for x in records)<=1e-8
    passed=passed and max(x['FDT_error_final'] for x in field)<=.0005 and all(x['all_errors_below_computed_bound'] for x in residual)
    rec={'status':'PASS' if passed else 'FAIL','scope':'Finite implementation checks only; coarse grid is not continuum reference',
         'grid':g.receipt(),'CSR_action_max_error':action_error,'propagation':records,'field':field,'finite_residual':residual,
         'seconds':time.perf_counter()-start,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
         'core_sha256':hashlib.sha256((Path(__file__).parent/'core.py').read_bytes()).hexdigest()}
    dump(out/'receipt.json',rec);print(json.dumps({k:rec[k] for k in ['status','CSR_action_max_error','seconds']}),flush=True)
    if not passed: raise RuntimeError('Small implementation checks failed; receipts retained')
    return rec


if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('--output',required=True);p.add_argument('--backend',default='numpy')
    a=p.parse_args();check(a.output,a.backend)
