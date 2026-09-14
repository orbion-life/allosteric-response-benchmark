"""Independent 64-state on/off field propagation in the original quartic model."""
from pathlib import Path
import os
for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS']:os.environ[k]='3'
import sys,time,json,hashlib
import numpy as np
from scipy.linalg import eigh,expm
ROOT=Path(__file__).resolve().parent;sys.path.insert(0,str(ROOT/'vendor'))
from finite_core import Grid
from scipy.linalg.blas import dgemm

def mm(a,b):return dgemm(1.,a,b)
def dump(p,x):Path(p).write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')
def run():
 p=json.loads((ROOT/'protocol.json').read_text())['tiny_pulse'];out=ROOT/'tiny-field';out.mkdir(exist_ok=False);start=time.perf_counter()
 base=Grid(p['faces'],p['spacing']);H=base.apply(np.eye(base.G));v,V=eigh(H);d=p['duration_over_tau']*base.tau;times=np.array(p['times_over_tau']);G0=base.G0
 def action(hv,hV,t,x):return mm(hV,np.exp(-t*hv)[:,None]*mm(hV.T,np.atleast_2d(x).reshape(len(hv),-1))).reshape(np.asarray(x).shape)
 def K(t):return mm(mm(base.F,V)*np.exp(-t*v),mm(V.T,base.F.T))
 step=np.array([K(t*base.tau)-G0 for t in times]);pulse=step.copy()
 for i,t in enumerate(times):
  if t>=p['duration_over_tau']:pulse[i]-=K((t-p['duration_over_tau'])*base.tau)-G0
 raw=np.zeros((3,3,2,len(times),3));independent=np.zeros_like(raw);conservation=[]
 for source in p['source_residues']:
  for e,eta in enumerate(p['eta']):
   for q,sign in enumerate([1,-1]):
    g=Grid(p['faces'],p['spacing'],field=(source,sign*eta));Hg=g.apply(np.eye(g.G));a,U=eigh(Hg)
    Qg=-(g.sp[:,None]*Hg)/g.sp[None,:];Q0=-(base.sp[:,None]*H)/base.sp[None,:]
    for i,ratio in enumerate(times):
     t=ratio*base.tau;first=min(t,d)
     density=g.sp*action(a,U,first,base.pi/g.sp)
     direct=mm(expm(first*Qg),base.pi[:,None])[:,0]
     if t>d:
      density=base.sp*action(v,V,t-d,density/base.sp)
      direct=mm(expm((t-d)*Q0),direct[:,None])[:,0]
     raw[source,e,q,i]=mm(base.f,density[:,None])[:,0];independent[source,e,q,i]=mm(base.f,direct[:,None])[:,0]
     conservation.append({'mass_error':float(abs(density.sum()-1)),'minimum_probability':float(density.min()),'forward_generator_max_abs':float(np.max(abs(density-direct)))})
 central=(raw[:,:,0]-raw[:,:,1])/(2*np.array(p['eta'])[None,:,None,None])
 ref=pulse.transpose(1,0,2)[:,None,:,:];errors=np.max(abs(central-ref),axis=(0,2,3));successive=np.max(abs(np.diff(central,axis=1)),axis=(0,2,3))
 independentcentral=(independent[:,:,0]-independent[:,:,1])/(2*np.array(p['eta'])[None,:,None,None])
 np.savez_compressed(out/'raw.npz',step_reference=step,pulse_reference=pulse,raw_on_off_means=raw,independent_forward_means=independent,central_response=central,times_over_tau=times,duration_over_tau=p['duration_over_tau'],eta=p['eta'],base_H=H,base_f=base.f,base_pi=base.pi,base_sp=base.sp)
 rec={'status':'COMPLETE','scientific_status':'PASS' if errors[-1]<=p['acceptance_max_abs'] and max(successive)<=p['acceptance_max_abs'] else 'FAIL','scope':p['scope'],'states':base.G,'max_identity_error_by_amplitude':errors.tolist(),'successive_central_changes':successive.tolist(),'independent_forward_central_max_abs':float(np.max(abs(central-independentcentral))),'max_probability_mass_error':max(x['mass_error'] for x in conservation),'minimum_probability':min(x['minimum_probability'] for x in conservation),'max_forward_probability_discrepancy':max(x['forward_generator_max_abs'] for x in conservation),'seconds':time.perf_counter()-start,'protocol_sha256':hashlib.sha256((ROOT/'protocol.json').read_bytes()).hexdigest(),'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
 dump(out/'receipt.json',rec);print(json.dumps(rec,indent=2))
if __name__=='__main__':run()
