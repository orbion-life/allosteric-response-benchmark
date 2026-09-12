"""Independent harmonic-coordinate verification, 2026-09-12.
Run python3 verify_math.py. NumPy/SciPy required.
Geometry and law match the synthetic benchmark in the repository root; stiffness is 100.
OU moments are an independently derived unbounded harmonic benchmark.
"""
import itertools, json, math, warnings
from pathlib import Path
import numpy as np
from numpy.polynomial.hermite import hermgauss
from scipy.special import logsumexp
warnings.filterwarnings("error", category=RuntimeWarning)
np.seterr(all="raise")
def mm(a,b):
 a,b=np.asarray(a),np.asarray(b)
 return np.einsum({(1,1):"i,i->",(1,2):"i,ij->j",(2,1):"ij,j->i",(2,2):"ij,jk->ik"}[a.ndim,b.ndim],a,b,optimize=False)
ROOT=Path(__file__).resolve().parent
R0=np.array([[0,0,0],[1.1,0,0],[.45,.95,0],[.35,.3,.95],[.65,.28,-.85]],float)
EDGES=[(i,j) for i in range(5) for j in range(i+1,5) if (i,j)!=(3,4)]
DEG=np.bincount(np.array(EDGES).ravel(),minlength=5)
KAPPA=100.; BETA=MU=1.
H=np.zeros((15,15))
for i,j in EDGES:
 v=np.zeros(15);a=R0[i]-R0[j];v[3*i:3*i+3]=a/np.linalg.norm(a);v[3*j:3*j+3]=-v[3*i:3*i+3]
 H+=KAPPA*np.outer(v,v)
vals,vecs=np.linalg.eigh(H);ix=np.flatnonzero(vals>1e-8)
lam=vals[ix];B=vecs[:,ix]
for j in range(B.shape[1]):
 if B[np.argmax(abs(B[:,j])),j]<0:B[:,j]*=-1
TAU=1/(MU*lam[0])

def fixture(d):
 powers=[]
 for degree in (2,3,4):
  for indices in itertools.combinations_with_replacement(range(d),degree):
   p=np.bincount(indices,minlength=d);powers.append(tuple(p))
 where={p:i for i,p in enumerate(powers)};A=np.zeros((5,len(powers)))
 for i,j in EDGES:
  a=R0[i]-R0[j];T=B[3*i:3*i+3,:d]-B[3*j:3*j+3,:d];l2=mm(a,a)
  coefficients={}
  for u in range(d):
   p=[0]*d;p[u]=1;coefficients[tuple(p)]=2*mm(a,T[:,u])
   for v in range(u,d):
    p=[0]*d;p[u]+=1;p[v]+=1;coefficients[tuple(p)]=mm(T[:,u],T[:,v])*(1 if u==v else 2)
  for p,c in coefficients.items():
   for q,z in coefficients.items():
    k=where[tuple(x+y for x,y in zip(p,q))];coeff=KAPPA*c*z/(8*l2)
    A[i,k]+=coeff/DEG[i];A[j,k]+=coeff/DEG[j]
 return np.array(powers),A

def joint_moment(a,b,variance,rho):
 # Wick/Isserlis pairings: k cross-time pairs, all other pairs same-time.
 if (a+b)%2:return 0.
 ans=0.
 for k in range(min(a,b)+1):
  if (a-k)%2 or (b-k)%2:continue
  left,right=(a-k)//2,(b-k)//2
  count=math.factorial(a)*math.factorial(b)/(math.factorial(k)*math.factorial(left)*math.factorial(right)*2**(left+right))
  ans+=count*variance**((a+b)/2)*rho**k
 return ans

def ou(d,t):
 powers,A=fixture(d);count=len(powers);equal=np.ones((count,count));delayed=np.ones((count,count));means=np.ones(count)
 for u in range(d):
  variance=1/(BETA*lam[u]);rho=np.exp(-MU*lam[u]*t) if MU*lam[u]*t<700 else 0.
  now=np.array([[joint_moment(a,b,variance,1.) for b in range(5)] for a in range(5)])
  later=np.array([[joint_moment(a,b,variance,rho) for b in range(5)] for a in range(5)])
  equal*=now[powers[:,u,None],powers[None,:,u]]
  delayed*=later[powers[:,u,None],powers[None,:,u]]
  means*=np.array([joint_moment(a,0,variance,1.) for a in powers[:,u]])
 cov=mm(mm(A,equal-np.outer(means,means)),A.T)
 lag=mm(mm(A,delayed-np.outer(means,means)),A.T)
 R=-BETA*(cov-lag);sh=np.sqrt(np.diag(cov));C=R/(BETA*np.outer(sh,sh))
 return {'d':d,'monomials':count,'R':R,'C':C,'sh':sh,'cov':cov,'delayed':lag,'A':A,'powers':powers}

def energy(q):
 d=q.shape[1];r=R0[None,:,:]+mm(q,B[:,:d].T).reshape(-1,5,3);E=np.zeros((len(q),5))
 for i,j in EDGES:
  l2=np.sum((R0[i]-R0[j])**2);u=KAPPA/(8*l2)*(np.sum((r[:,i]-r[:,j])**2,axis=1)-l2)**2
  E[:,i]+=u/DEG[i];E[:,j]+=u/DEG[j]
 return E

def quadrature_check(d=2):
 # Independent direct contact evaluation; exact degree-eight Gaussian quadrature.
 x,w=hermgauss(5);indices=np.array(list(itertools.product(range(5),repeat=2*d)))
 z=np.sqrt(2)*x[indices];weights=np.prod(w[indices]/np.sqrt(np.pi),axis=1)
 rho=np.exp(-MU*lam[:d]*TAU);sigma=np.sqrt(1/(BETA*lam[:d]))
 q0=z[:,:d]*sigma;qt=(rho*z[:,:d]+np.sqrt(1-rho**2)*z[:,d:])*sigma
 e0,et=energy(q0),energy(qt);mu0=mm(weights,e0);mut=mm(weights,et)
 cov=mm(e0.T,weights[:,None]*e0)-np.outer(mu0,mu0)
 lag=mm(e0.T,weights[:,None]*et)-np.outer(mu0,mut)
 R=-BETA*(cov-lag)
 return {'points':len(weights),'R':R.tolist(),'maximum_difference_from_moment_route':float(abs(R-ou(d,TAU)['R']).max())}

def hermite_ou(d,t):
 # Independent spectral route: expand monomials in orthogonal Hermite modes.
 powers,A=fixture(d);coeff={}
 for col,p in enumerate(powers):
  options=[]
  for u,a in enumerate(p):
   scale=(1/(BETA*lam[u]))**(a/2)
   options.append([(a-2*j,scale*math.factorial(int(a))/(2**j*math.factorial(j)*math.factorial(int(a-2*j)))) for j in range(a//2+1)])
  for term in itertools.product(*options):
   degree=tuple(v[0] for v in term)
   if not any(degree):continue
   value=np.prod([v[1] for v in term])
   coeff.setdefault(degree,np.zeros(5))
   coeff[degree]+=A[:,col]*value
 R=np.zeros((5,5))
 for degree,c in coeff.items():
  norm=math.prod(math.factorial(int(a)) for a in degree)
  decay=1-np.exp(-MU*sum(a*l for a,l in zip(degree,lam))*t)
  R-=BETA*norm*decay*np.outer(c,c)
 return R

models=[ou(d,TAU) for d in (2,3,4,9)];reference=models[-1]['sh']
results=[]
for o in models:
 c9=o['R']/(BETA*np.outer(reference,reference))
 results.append({'d':o['d'],'monomials':o['monomials'],'R54':float(o['R'][4,3]),'C54_own_normalization':float(o['C'][4,3]),'C54_common_d9_normalization':float(c9[4,3]),'C_common_d9':c9.tolist(),'max_R_asymmetry':float(abs(o['R']-o['R'].T).max()),'hermite_spectral_max_error':float(abs(o['R']-hermite_ou(o['d'],TAU)).max()),'largest_R_eigenvalue':float(np.linalg.eigvalsh(o['R']).max())})
# Independent recreation of uniform-monomial measurement amplification.
amplification=[]
for d,n in ((2,32),(4,16)):
 powers,A=fixture(d);o=ou(d,TAU);half=4/np.sqrt(BETA*lam[0]);axis=np.linspace(-half,half,n)
 q=np.stack(np.meshgrid(*[axis]*d,indexing='ij'),axis=-1).reshape(-1,d)
 mon=np.stack([np.prod(q**p,axis=1) for p in powers],axis=1);E=energy(q)
 for model in ('harmonic','biquadratic'):
  U=.5*np.sum(q*q*lam[:d],axis=1) if model=='harmonic' else mm(E,DEG)/2
  pi=np.exp(-BETA*U-logsumexp(-BETA*U));means=mm(pi,mon)
  norms=np.sqrt(mm(pi,(mon-means)**2));gamma=mm(abs(A),norms)/o['sh']
  amplification.append({'d':d,'n':n,'model':model,'gamma':gamma.tolist(),'maximum_gamma_product':float(gamma.max()**2),'apex_gamma_product':float(gamma[3]*gamma[4]),'energy_identity_max_error':float(abs(mm(mon,A.T)-E).max())})
counts=[]
for d in (2,4,6,8,10):
 b=math.comb(d+4,4)-1-d;m=b*(b+1)//2;shots=math.ceil(2*math.log(2*m/.05)/.005**2)
 counts.append({'d':d,'grid_side_16_G':16**d,'grid_side_32_G':32**d,'basis':b,'unique_symmetric_overlaps':m,'illustrative_overlap_error':.005,'total_failure':.05,'shots_per_overlap':shots,'total_shots':m*shots})
zero=ou(2,0);eq=ou(2,1e4*TAU)
out={'date':'2026-09-12','scope':'New audit calculations: exact Gaussian-moment harmonic continuum comparator and shot-bound arithmetic; no quantum-device execution or biological validation','tau':TAU,'harmonic_results':results,'quadrature_check':quadrature_check(),'R0_max_abs':float(abs(zero['R']).max()),'equilibrium_delayed_covariance_max_abs':float(abs(eq['delayed']).max()),'amplification':amplification,'shot_counts':counts}
(ROOT/'math-calculations-independent.json').write_text(json.dumps(out,indent=2))
print(json.dumps({k:v for k,v in out.items() if k not in ('amplification','shot_counts')},indent=2))
