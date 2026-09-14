"""Classical loading-coefficient rounding check for the named dense fallback."""
from core import *
from scipy.linalg import cholesky
from threadpoolctl import threadpool_limits

def main():
    s=dict(np.load(ROOT/'inputs/operator.npz'));raw=dict(np.load(ROOT/'results/matched-responses.npz'));L=cholesky(s['Hr'],lower=False);B=s['B'];bmax=float(np.max(np.linalg.norm(B,axis=0)));rows=[]
    for dtype in [np.float32,np.float16]:
        Lq=L.astype(dtype).astype(float);Bq=B.astype(dtype).astype(float);Hq=mm(Lq.T,Lq);dH=Hq-s['Hr'];dHnorm=float(np.max(abs(eigh((dH+dH.T)/2,eigvals_only=True))));bdelta=float(np.max(np.linalg.norm(Bq-B,axis=0)));Gq=mm(Bq.T,Bq);gdelta=float(np.max(abs(Gq-s['G0'])));lam,Q=eigh(Hq);weights=mm(Q.T,Bq);details=[]
        for t,ref in zip(raw['times'],raw['reduced_spectral']):
            val=mm(weights.T,np.exp(-t*lam)[:,None]*weights)-Gq
            bound=t*dHnorm*bmax*bmax+bdelta*(2*bmax+bdelta)+gdelta
            details.append({'time':float(t),'response_maxabs_error':float(np.max(abs(val-ref))),'conservative_fixed_response_bound':float(bound)})
        rows.append({'loading_coefficient_dtype':np.dtype(dtype).name,'factor_stored_bytes':L.astype(dtype).nbytes,'observable_stored_bytes':B.astype(dtype).nbytes,'factor_Gram_perturbation_spectral_norm':dHnorm,'max_observable_vector_L2_error':bdelta,'static_Gram_maxabs_error':gdelta,'by_time':details})
    dump(ROOT/'results/precision-checks.json',{'scope':'Numerical coefficient-rounding sensitivity for explicitly dense Cholesky fallback. Does not validate physical gate synthesis precision, noise, compact factor oracle or full ABL circuit.','reference':'Frozen rank1261 ABL Gaussian reduced response, fixed harmonic normalization, same4times','results':rows})

if __name__=='__main__':
    with threadpool_limits(limits=4):main()
