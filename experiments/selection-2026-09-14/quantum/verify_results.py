"""Independent emitted-gate accounting, high-precision rotation audit and saved response replay."""
import os
for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS']:os.environ[key]='1'
import hashlib,json,math,time
from pathlib import Path
import mpmath as mp
import numpy as np
from scipy.special import ive
from qiskit import qpy
ROOT=Path(__file__).resolve().parent
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()

def circuit_counts(p):
    with Path(p).open('rb') as f:q=qpy.load(f)[0]
    ops=q.count_ops();return dict(T=ops.get('t',0)+ops.get('tdg',0),CX=ops.get('cx',0),operations=sum(ops.values()),qubits=q.num_qubits)

def mp_error(record):
    theta=mp.mpf(str(record['angle']));gates=record['gates'];a=mp.mpc(1);b=mp.mpc(0);c=mp.mpc(0);d=mp.mpc(1)
    w=mp.exp(mp.j*mp.pi/4);root2=mp.sqrt(2)
    if gates is None:
        k=round(record['angle']/(math.pi/4));a=mp.exp(-mp.j*k*mp.pi/8);d=mp.exp(mp.j*k*mp.pi/8)
    else:
        for g in gates:
            if g=='H':a,b,c,d=(a+b)/root2,(a-b)/root2,(c+d)/root2,(c-d)/root2
            elif g=='T':b*=w;d*=w
            elif g=='S':b*=mp.j;d*=mp.j
            elif g=='X':a,b=b,a;c,d=d,c
            elif g=='W':a*=w;b*=w;c*=w;d*=w
            else:raise ValueError(g)
    a-=mp.exp(-mp.j*theta/2);d-=mp.exp(mp.j*theta/2)
    x=abs(a)**2+abs(c)**2;z=abs(b)**2+abs(d)**2;y=mp.conj(a)*b+mp.conj(c)*d
    return mp.sqrt((x+z+mp.sqrt((x-z)**2+4*abs(y)**2))/2)

def main():
    start=time.perf_counter();mp.mp.dps=70;checks={};cases=[]
    for record in json.loads((ROOT/'inputs/provenance.json').read_text())['files']:
        assert sha(ROOT/record['path'])==record['sha256']
    for path,h in json.loads((ROOT/'source-seal.json').read_text())['files'].items():assert sha(ROOT/path)==h
    inp=np.load(ROOT/'inputs/triangle.npz',allow_pickle=False);H,B=inp['Hr'],inp['B'];ev,V=np.linalg.eigh(H);QB=V.T@B
    exact=np.array([(QB.T*np.expm1(-t*ev))@QB for t in inp['times']]);scale=np.outer(np.linalg.norm(B,axis=0),np.linalg.norm(B,axis=0));projection=float(np.max(abs(exact-inp['reference_C'])))
    receipts={}
    for name in ['select_swap','dense_norm_optimal','dense_factor_normalization']:
        folder=ROOT/'results'/name;r=json.loads((folder/'receipt.json').read_text());z=np.load(folder/'mathematical-results.npz');receipts[name]=r
        coeff=[]
        for ti,t in enumerate(inp['times']):
            k=np.arange(r['polynomial_degrees'][ti]+1);c=ive(k,r['alpha']*t);c[1:]*=2;coeff.append(c)
        C=np.array([scale*np.einsum('k,kij->ij',c,z['terms'][:len(c)]).real-B.T@B for c in coeff])
        err=float(np.max(abs(C-z['C'])));assert err<1e-12
        assert np.max(abs(C-inp['reference_C']))<.002
        if (folder/'controlled-walk.qpy').exists():
            wc=circuit_counts(folder/'controlled-walk.qpy')
            for key in wc:assert wc[key]==r['cost']['controlled_walk'][key]
        prep=[circuit_counts(folder/f'prepare-{i}.qpy') for i in range(3)]
        for p,q in zip(prep,r['cost']['preparations']):
            for key in p:assert p[key]==q[key]
        shots=math.ceil(2*float(scale.max())**2*math.log(36/.05)/(.002-projection-2e-6-2e-8)**2)
        assert shots==r['sampling']['shots_per_entry']
        totalT=0.
        for ti,c in enumerate(coeff):
            order=float(np.dot(np.arange(len(c)),c/c.sum()))
            for i in range(3):
                for j in range(i,3):totalT+=shots*(prep[i]['T']+prep[j]['T']+order*r['cost']['controlled_walk']['T'])
        relative=abs(totalT/r['cost']['expected_total_T']-1);assert relative<1e-12
        cases.append(dict(name=name,response_reconstruction_error=err,expected_T_reconstruction_relative_error=relative,shots_per_entry=shots,
                          source_sha256=r['source_sha256'],source_hash_matches=any(sha(ROOT/p)==r['source_sha256'] for p in ['experiment.py','candidate_experiment.py'])))
    rotation=[]
    for p in sorted((ROOT/'results').glob('**/rotation-records.json')):
        records=json.loads(p.read_text());worst=mp.mpf(0);maximum=mp.mpf(0)
        for rec in records:
            error=mp_error(rec);ratio=error/mp.mpf(str(rec['epsilon']));worst=max(worst,ratio);maximum=max(maximum,error)
            assert ratio<=1+mp.mpf('1e-10'),(str(p),rec['angle'],str(error),rec['epsilon'])
            if rec['gates'] is not None:assert rec['gates'].count('T')==rec['T']
        rotation.append(dict(path=str(p.relative_to(ROOT)),rotations=len(records),max_operator_error=float(maximum),max_error_over_declared_tolerance=float(worst),sha256=sha(p)))
    protein=[]
    for target in ['kras','abl','myc','myh7']:
        folder=ROOT/'results/proteins'/target;r=json.loads((folder/'receipt.json').read_text());arrays=np.load(ROOT/'inputs'/(target+'-measured-arrays.npz'));d=r['padded_rank'];rank=r['rank'];n=(d-1).bit_length();N=r['residues']
        assert np.max(abs(arrays['L'].T@arrays['L']-arrays['H']))<1e-10
        assert sha(folder/'factor-angle-tables.npz')==r['factor_tables_sha256']
        assert sha(folder/'observable-angle-words.npz')==r['observable_words_sha256']
        bc=circuit_counts(folder/'rotation-bank.qpy')
        for key in bc:assert bc[key]==r['rotation_bank'][key]
        words=np.load(folder/'factor-angle-tables.npz');tt=0
        for i,t in enumerate(r['lookup_tables']):
            pop=int(np.bitwise_count(words['words_'+str(i)]).sum());assert pop==t['word_ones']
            m=t['N']//t['lambda_value'];lookupT=16*max(m-2,0)+7*r['precision_bits']*(t['lambda_value']-1)
            assert lookupT==t['literal_lookup_T'];tt+=2*lookupT+bc['T']
        assert 2*tt+16*n==r['controlled_walk_T_upper']
        assert r['full_matrix_shots_sufficient']==3*N*(N+1)//2*r['shots_per_entry']
        protein.append(dict(target=target,status='PASS_SCOPED_FORMULA_AND_RAW_DATA_AUDIT',full_circuit_compiled=False))
    candidate=receipts['select_swap'];dense=receipts['dense_norm_optimal'];ratio=candidate['cost']['expected_total_T']/dense['cost']['expected_total_T']
    result=dict(status='PASS_VERIFICATION',decision='FAIL_PREDECLARED_HALF_T_COST_GATE',candidate_to_strongest_evaluated_dense_expected_T_ratio=ratio,
        half_cost_gate_pass=ratio<=.5,scope='Verification passes independently of the scientific decision gate. Cost includes exact cached composition of actual small components; missing high-order/full-workspace simulation remains missing.',
        cases=cases,rotation_operator_checks_70_decimal_digits=rotation,protein=protein,wall_seconds=time.perf_counter()-start,source_sha256=sha(__file__))
    (ROOT/'independent-verification.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))

if __name__=='__main__':main()
