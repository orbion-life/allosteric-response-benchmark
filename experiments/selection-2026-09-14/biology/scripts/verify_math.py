"""Independent Wick polynomial and finite-difference checks on the actual GRB2 fit."""
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS']:os.environ[k]='1'
from pathlib import Path
import sys,json,functools,itertools
import numpy as np
P=Path(__file__).resolve().parents[1];sys.path.insert(0,str(P/'vendor'));import gaussian_protein as gp
m=dict(np.load(P/'results/primary-model/static-model.npz'));fit=np.load(P/'results/primary-model/gaussian/fitted-covariance.npz');B=m['B'];Sigma=fit['covariance'];Y=fit['whitened_covariance'];r0=m['r0'];edges=m['edges'];lam=m['eigenvalues'];tau=1/lam[0];N=len(r0)
values,V=np.linalg.eigh(Sigma);cart=B@Sigma@B.T;lag=(B@V*(values*np.exp(-tau/values)))@(B@V).T
# Expand p(x)=(2*a.x+x.x)^2 by polynomial multiplication, independently of Hermite contraction.
def polynomial(a):
 base={}
 for i in range(3):
  e=[0]*3;e[i]=1;base[tuple(e)]=2*a[i];e[i]=2;base[tuple(e)]=1.
 out={}
 for p,c in base.items():
  for q,d in base.items():
   k=tuple(x+y for x,y in zip(p,q));out[k]=out.get(k,0.)+c*d
 return out

def wick_cov(a,b,S,T,C):
 cov=np.block([[S,C],[C.T,T]])
 @functools.lru_cache(None)
 def moment(ids):
  if not ids:return 1.
  if len(ids)%2:return 0.
  return sum(cov[ids[0],ids[j]]*moment(ids[1:j]+ids[j+1:]) for j in range(1,len(ids)))
 def ids(p,offset):return tuple(offset+i for i,k in enumerate(p) for _ in range(k))
 pa=polynomial(a);pb=polynomial(b);ex=sum(c*moment(ids(p,0)) for p,c in pa.items());ey=sum(c*moment(ids(p,3)) for p,c in pb.items());xy=sum(c*d*moment(ids(p,0)+ids(q,3)) for p,c in pa.items() for q,d in pb.items());return xy-ex*ey

def edgeblock(M,e,f):
 i,j=edges[e];u,v=edges[f]
 def block(i,j):return M[3*i:3*i+3,3*j:3*j+3]
 return block(i,u)-block(i,v)-block(j,u)+block(j,v)
checks=[]
for e,f in [(0,0),(0,len(edges)//2),(len(edges)//3,len(edges)-1)]:
 a=r0[edges[e,0]]-r0[edges[e,1]];b=r0[edges[f,0]]-r0[edges[f,1]];S=edgeblock(cart,e,e);T=edgeblock(cart,f,f)
 for time,M in [('static',cart),('tau',lag)]:
  C=edgeblock(M,e,f);ind=wick_cov(a,b,S,T,C);prod=float(gp.contracted_covariance(a,b,S,T,C));err=abs(ind-prod)/max(1.,abs(ind));assert err<1e-11;checks.append(dict(edge_pair=[e,f],time=time,independent_covariance=ind,production_covariance=prod,scaled_absolute_error=err))
W=B/np.sqrt(lam)
def independent_objective(y):
 V=W@y@W.T;total=.5*np.trace(y)-.5*np.linalg.slogdet(y)[1]
 for i,j in edges:
  s=V[3*i:3*i+3,3*i:3*i+3]+V[3*j:3*j+3,3*j:3*j+3]-V[3*i:3*i+3,3*j:3*j+3]-V[3*j:3*j+3,3*i:3*i+3];l2=np.sum((r0[i]-r0[j])**2);total+=(2*np.trace(s@s)+np.trace(s)**2)/(8*l2)
 return float(total)
obj=gp.CovarianceObjective(r0,edges,B,lam);test=.8*Y;F,G,*_=obj.evaluate(test);assert abs(F-independent_objective(test))<1e-10;rng=np.random.default_rng(9189);grad=[]
for _ in range(6):
 D=rng.normal(size=Y.shape);D=(D+D.T)/2;D/=np.linalg.norm(D);eps=1e-5;finite=(independent_objective(test+eps*D)-independent_objective(test-eps*D))/(2*eps);analytic=float(np.sum(G*D));err=abs(finite-analytic);assert err<2e-7;grad.append(dict(finite_difference=finite,production_gradient=analytic,absolute_error=err))
# A second full harmonic implementation has an independently coded incidence contraction.
from run_pilot import harmonic_covariances
hc=harmonic_covariances(r0,edges,B,lam,tau);gpc,_,_=gp.moments_for_covariance(m,np.diag(1/lam),[tau]);he=max(np.max(abs(hc[0]-gpc[0])),np.max(abs(hc[1]-gpc[1])));assert he<1e-10
print(json.dumps({'status':'PASS','actual_model_N':N,'actual_model_D':len(lam),'independent_polynomial_Wick_checks':checks,'independent_free_energy_gradient_checks':grad,'independent_full_harmonic_covariance_max_error':float(he),'scope':'Actual finite Gaussian model implementation checks. No original-nonlinear physical fidelity or biological advantage inferred.'},indent=2))
