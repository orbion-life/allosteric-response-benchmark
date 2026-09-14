"""Original finite reversible model; FP64 CPU and GPU implementations.

The new algorithms change propagation and storage, not the original potential,
coordinate measure, diffusion or reflecting finite-volume discretization.
"""
import os
for _key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ[_key]='1'
from pathlib import Path
import sys, json, math, time, resource, hashlib, datetime
import numpy as np
from scipy.special import ive
from scipy.stats import poisson
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'vendor'))
from model import build_model, analytic_gaussian, grid_shape
from chebyshev import required_degree, log_tail_bound


def dump(p,obj):
    Path(p).write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n')


def host(x):
    return x.get() if hasattr(x,'get') else np.asarray(x)


def scalar(x):
    return float(host(x))


def sync(xp):
    if xp.__name__=='cupy': xp.cuda.Stream.null.synchronize()


def backend(name):
    if name=='cupy':
        import cupy
        return cupy
    return np


GPU_CODE=r'''
extern "C" __global__ void apply(
 const double* diag,const double* l0,const double* l1,const double* l2,
 const double* cur,const double* old,double* out,double* result,
 const long long G,const int n0,const int n1,const int n2,const int nvec,
 const double a,const double b,const double c,const double coeff,const int accum){
 long long z=(long long)blockDim.x*blockIdx.x+threadIdx.x;
 if(z>=G*nvec) return;
 long long g=z%G; long long st0=(long long)n1*n2;
 int i0=g/st0, i1=(g/n2)%n1, i2=g%n2;
 double h=diag[g]*cur[z];
 if(i0<n0-1) h+=l0[g]*cur[z+st0];
 if(i0>0) h+=l0[g-st0]*cur[z-st0];
 if(i1<n1-1) h+=l1[g]*cur[z+n2];
 if(i1>0) h+=l1[g-n2]*cur[z-n2];
 if(i2<n2-1) h+=l2[g]*cur[z+1];
 if(i2>0) h+=l2[g-1]*cur[z-1];
 double val=a*cur[z]+b*h;
 if(c!=0) val+=c*old[z];
 out[z]=val;
 if(accum) result[z]+=coeff*val;
}
'''


class Grid:
    def __init__(self,faces,spacing,xp=np,field=None):
        start=time.perf_counter(); self.xp=xp
        self.model=build_model(); self.normal=analytic_gaussian(self.model)
        self.faces=np.asarray(faces,float); self.h=spacing
        self.shape=grid_shape(faces,spacing); self.G=int(np.prod(self.shape))
        self.sd=xp.asarray(self.normal['sd']); self.field=field
        self.U=xp.empty(self.G,dtype=xp.float64)
        self.E=xp.empty((3,self.G),dtype=xp.float64)
        self.Bx=xp.asarray(self.model['Bx']); self.r0=xp.asarray(self.model['r0'])
        self.tau=1/(self.model['mu']*self.model['eigenvalues'][0])
        n0,n1,n2=self.shape
        # Stream geometry to avoid an Nstates x Cartesian-dimension allocation.
        for lo in range(0,self.G,262144):
            hi=min(self.G,lo+262144); ids=xp.arange(lo,hi,dtype=xp.int64)
            ix=xp.stack((ids//(n1*n2),(ids//n2)%n1,ids%n2),axis=1)
            x=-xp.asarray(faces)+(ix+.5)*spacing
            rr=x@self.Bx.T
            rr=rr.reshape(-1,3,3)+self.r0
            U=xp.zeros(hi-lo); E=xp.zeros((3,hi-lo))
            for i,j in self.model['edges']:
                rest2=float(np.sum((self.model['r0'][i]-self.model['r0'][j])**2))
                delta=xp.sum((rr[:,i]-rr[:,j])**2,axis=1)-rest2
                u=self.model['kappa']*delta**2/(8*rest2)
                U+=u; E[i]+=u/self.model['degree'][i]; E[j]+=u/self.model['degree'][j]
            self.U[lo:hi]=U; self.E[:,lo:hi]=E
        V=self.model['beta']*self.U
        vmax=xp.max(-V); lz=vmax+xp.log(xp.sum(xp.exp(-V-vmax)))
        self.logpi0=-V-lz
        self.pi0=xp.exp(self.logpi0)
        self.mean=self.E@self.pi0
        self.f=(self.E-self.mean[:,None])/self.sd[:,None]
        if field is not None:
            source,eta=field
            V=V+eta*self.f[source]
        vmax=xp.max(-V); self.logZ=vmax+xp.log(xp.sum(xp.exp(-V-vmax)))
        self.logpi=-V-self.logZ; self.pi=xp.exp(self.logpi); self.sp=xp.exp(.5*self.logpi)
        self.F=xp.ascontiguousarray(self.sp[None,:]*self.f)
        self.left=xp.ascontiguousarray(xp.exp(self.logpi0-.5*self.logpi))
        self.V=V
        self.diag=xp.zeros(self.G); self.links=[xp.zeros(self.G) for _ in range(3)]
        vv=V.reshape(self.shape); dd=self.diag.reshape(self.shape)
        for k in range(3):
            lo=[slice(None)]*3; hi=[slice(None)]*3
            lo[k]=slice(0,-1); hi[k]=slice(1,None); lo,hi=tuple(lo),tuple(hi)
            delta=vv[hi]-vv[lo]
            # Stable logistic pairs preserve exact rate sum (up to roundoff).
            base=2*self.model['mu']*self.model['eigenvalues'][k]/spacing**2
            z=xp.exp(-xp.abs(delta)); small=z/(1+z)
            pab=xp.where(delta>=0,small,1-small); pba=1-pab
            dd[lo]+=base*pab; dd[hi]+=base*pba
            # Evaluate sqrt(pab*pba) without catastrophic 1-p cancellation.
            off=-base*xp.exp(-xp.abs(delta)/2)/(1+z)
            self.links[k].reshape(self.shape)[lo]=off
        self.gamma=scalar(xp.max(self.diag))
        row=self.diag.copy().reshape(self.shape)
        for k in range(3):
            lo=[slice(None)]*3; hi=[slice(None)]*3
            lo[k]=slice(0,-1); hi[k]=slice(1,None);lo,hi=tuple(lo),tuple(hi)
            off=self.links[k].reshape(self.shape)[lo]
            row[lo]-=off; row[hi]-=off
        self.upper=scalar(xp.max(row)); del row
        self.kernel=xp.RawKernel(GPU_CODE,'apply') if xp.__name__=='cupy' else None
        station=self.apply(self.sp[None,:],a=0,b=1)
        self.stationarity=scalar(xp.max(xp.abs(station)))
        self.G0=host(self.F@self.F.T)
        sync(xp)
        self.build_seconds=time.perf_counter()-start
        if self.stationarity>1e-9: raise RuntimeError('Stationarity gate failed: '+str(self.stationarity))

    def apply(self,cur,a=0.,b=1.,old=None,c=0.,out=None,result=None,coeff=0.):
        xp=self.xp
        cur=xp.ascontiguousarray(cur)
        if out is None: out=xp.empty_like(cur)
        if old is None: old=cur
        if self.kernel is not None:
            if result is None: result=out; accum=0
            else: accum=1
            args=(self.diag,*self.links,cur,old,out,result,np.int64(self.G),
                  *[np.int32(n) for n in self.shape],np.int32(cur.shape[0]),
                  np.float64(a),np.float64(b),np.float64(c),np.float64(coeff),np.int32(accum))
            self.kernel(((cur.size+255)//256,),(256,),args)
        else:
            X=cur.reshape((len(cur),*self.shape))
            Y=out.reshape(X.shape); Y[:]=self.diag.reshape(self.shape)*X
            for k in range(3):
                lo=[slice(None)]*3; hi=[slice(None)]*3
                lo[k]=slice(0,-1);hi[k]=slice(1,None);lo,hi=tuple(lo),tuple(hi)
                off=self.links[k].reshape(self.shape)[lo]
                Y[(slice(None),*lo)]+=off*X[(slice(None),*hi)]
                Y[(slice(None),*hi)]+=off*X[(slice(None),*lo)]
            out*=b; out+=a*cur
            if c: out+=c*old
            if result is not None: result+=coeff*out
        return out

    def propagate(self,t,method='chebyshev',F=None,left=None,tolerance=1e-9,progress=None):
        xp=self.xp; F=self.F if F is None else xp.ascontiguousarray(F)
        left=F if left is None else xp.atleast_2d(left)
        normF=xp.sqrt(xp.sum(F*F,axis=1)); normL=xp.sqrt(xp.sum(left*left,axis=1))
        amp=scalar(xp.max(normF))*scalar(xp.max(normL)); tol=min(.5,tolerance/max(amp,1e-300))
        start=time.perf_counter();old=F.copy(); temp=xp.empty_like(F)
        if method=='chebyshev':
            z=t*self.upper/2; degree=required_degree(z,tol)
            coeff=ive(np.arange(degree+1),z); result=coeff[0]*old
            if degree:
                cur=self.apply(F,a=1,b=-2/self.upper)
                result+=2*coeff[1]*cur
                for k in range(2,degree+1):
                    self.apply(cur,a=2,b=-4/self.upper,old=old,c=-1,out=temp,result=result,coeff=2*coeff[k])
                    old,cur,temp=cur,temp,old
                    if progress is not None and k%200==0: progress(k,degree)
            tail=math.exp(log_tail_bound(z,degree)); calls=degree
        elif method=='uniformization':
            z=t*self.gamma
            degree=int(poisson.isf(tol,z))+1
            while poisson.sf(degree,z)>tol: degree+=1
            weights=poisson.pmf(np.arange(degree+1),z)
            result=weights[0]*old;cur=old
            for k in range(1,degree+1):
                self.apply(cur,a=1,b=-1/self.gamma,out=temp,result=result,coeff=weights[k])
                cur,temp=temp,cur
                if progress is not None and k%1000==0: progress(k,degree)
            tail=float(poisson.sf(degree,z)); calls=degree
        else: raise ValueError(method)
        value=host(left@result.T)
        sync(xp)
        if not np.isfinite(value).all(): raise FloatingPointError('nonfinite response')
        return value,{'method':method,'degree':degree,'time':t,'time_over_tau':t/self.tau,
                      'seconds':time.perf_counter()-start,'operator_actions':calls,
                      'tail_bound_times_norms':tail*amp,'norm_amplification':amp,
                      'arithmetic_scope':'Floating-point computation; no interval arithmetic certificate.'}

    def receipt(self):
        xp=self.xp
        r={'faces':self.faces.tolist(),'spacing':self.h,'shape':list(self.shape),'states':self.G,
           'field':self.field,'backend':xp.__name__,'build_seconds':self.build_seconds,
           'stationarity':self.stationarity,'G0':self.G0.tolist(),'gamma':self.gamma,'spectral_upper':self.upper,
           'eigenvalues':self.model['eigenvalues'].tolist(),'sd':host(self.sd).tolist(),
           'log_partition_cell_measure':scalar(self.logZ)+3*math.log(self.h),
           'numpy_version':np.__version__,'backend_version':xp.__version__}
        if xp.__name__=='cupy':
            r['GPU']=xp.cuda.runtime.getDeviceProperties(0)['name'].decode()
            r['GPU_pool_reserved_bytes']=xp.get_default_memory_pool().total_bytes()
        return r


def run_reference(out,faces,spacing,backend_name='numpy',times=(.1,1,10),independent=True):
    out=Path(out);out.mkdir(parents=True,exist_ok=False)
    started=time.perf_counter();xp=backend(backend_name)
    grid=Grid(faces,spacing,xp)
    rec=grid.receipt(); allK={};work=[]
    methods=['chebyshev','uniformization'] if independent else ['chebyshev']
    for method in methods:
        kk=[]
        for tt in times:
            def progress(k,degree):
                dump(out/'progress.json',{'method':method,'time_over_tau':tt,'action':k,'total_actions':degree,'seconds':time.perf_counter()-started})
            val,meta=grid.propagate(tt*grid.tau,method,progress=progress)
            kk.append(val);work.append(meta)
            print(json.dumps({'method':method,'time_over_tau':tt,'seconds':meta['seconds'],'degree':meta['degree']}),flush=True)
        allK[method]=np.array(kk)
    data={'G0':grid.G0,'times_over_tau':np.array(times),'sd':host(grid.sd)}
    for method,k in allK.items(): data['K_'+method]=k;data['C_'+method]=k-grid.G0
    np.savez_compressed(out/'response.npz',**data)
    error=float(np.max(abs(allK['chebyshev']-allK['uniformization']))) if independent else None
    rec.update(engineering_status='COMPLETE',scientific_scope='Stated finite reflecting model; reference convergence and coverage assessed separately',
               independent_error=error,independent_status='PASS' if error is not None and error<=1e-8 else 'NOT_RUN' if error is None else 'FAIL',
               propagation=work,total_seconds=time.perf_counter()-started,
               process_peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform=='darwin' else 1024),
               source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    dump(out/'receipt.json',rec)
    return rec


if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('--output',required=True);p.add_argument('--faces',type=float,nargs=3,required=True)
    p.add_argument('--spacing',type=float,required=True);p.add_argument('--backend',default='numpy');p.add_argument('--no-independent',action='store_true')
    a=p.parse_args();run_reference(a.output,a.faces,a.spacing,a.backend,independent=not a.no_independent)
