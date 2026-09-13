"""Bounded, local-only trace of unchanged scientific numerical statements."""
from pathlib import Path
import argparse, ast, contextlib, datetime, difflib, hashlib, io, json, os, platform, sys, time, types
import numpy as np
import scipy
from scipy.linalg import svd, eigh, qr
from scipy.sparse import csr_matrix
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'immutable'))
import operator_study as original
mm=original.mm

def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(p,x):
    Path(p).parent.mkdir(parents=True,exist_ok=True)
    Path(p).write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')
def save(p,**arrays):
    clean={}
    for k,v in arrays.items():
        a=np.asarray(v)
        if a.dtype.kind not in 'biufc':continue
        if not np.isfinite(a).all():raise ArithmeticError('Nonfinite '+k)
        clean[k]=a
    Path(p).parent.mkdir(parents=True,exist_ok=True);np.savez(p,**clean)
def maxdiff(a,b):return float(abs(a-b).max()) if a.size else 0.
def norm(a):return float(svd(a,compute_uv=False,check_finite=True)[0]) if a.size else 0.
def asarrays(x,names):return {k:x[k] for k in names.split() if k in x}
def gap(a,b):
    if not np.isfinite(a).all() or not np.isfinite(b).all():raise ValueError('Nonfinite subspace')
    if a.shape[0]!=b.shape[0]:raise ValueError('Different ambient dimensions')
    ua,sa,_=svd(a,full_matrices=False);ub,sb,_=svd(b,full_matrices=False)
    ta=np.finfo(float).eps*max(a.shape)*(float(sa[0]) if len(sa) else 0.)
    tb=np.finfo(float).eps*max(b.shape)*(float(sb[0]) if len(sb) else 0.)
    ra=int(np.sum(sa>ta));rb=int(np.sum(sb>tb))
    rank_info={'column_counts':[a.shape[1],b.shape[1]],'effective_ranks':[ra,rb],
      'singular_values_a':sa.tolist(),'singular_values_b':sb.tolist(),
      'diagnostic_rank_cutoffs':[ta,tb],'rank_deficient_input':ra<a.shape[1] or rb<b.shape[1],
      'rank_policy':'Diagnostic only: eps*max(shape)*largest singular value; scientific admission unchanged.'}
    if ra==0 or rb==0:
        z=0. if ra==rb else 1.
        return dict(rank_info,ranks=[ra,rb],same_rank=ra==rb,principal_cosines=[],angles_radians=[],directional_gaps=[float(rb>0),float(ra>0)],maximum_gap=z)
    qa=ua[:,:ra];qb=ub[:,:rb]
    overlap=mm(qa.T,qb)
    s=svd(overlap,compute_uv=False)
    ab=qb-mm(qa,overlap);ba=qa-mm(qb,overlap.T)
    sab=svd(ab,compute_uv=False);sba=svd(ba,compute_uv=False)
    # Residual singular values resolve tiny angles without arccos cancellation.
    small=np.sort(sab)[::-1][:min(ra,rb)]
    if ra!=rb:small=np.sqrt(np.maximum(0.,1-np.clip(s,0,1)**2))[::-1]
    return {**rank_info,'ranks':[ra,rb],'same_rank':ra==rb,
      'principal_cosines':s.tolist(),'angles_radians':np.arcsin(np.clip(small,0,1)).tolist(),
      'directional_gaps':[float(sab[0]) if len(sab) else 0.,float(sba[0]) if len(sba) else 0.],
      'maximum_gap':max(float(sab[0]) if len(sab) else 0.,float(sba[0]) if len(sba) else 0.)}

class Instrument(ast.NodeTransformer):
    def __init__(self):self.function=''
    def hook(self,tag):return ast.parse("_trace(%r, locals())"%tag).body[0]
    def visit_FunctionDef(self,node):
        prior=self.function;self.function=node.name
        node=self.generic_visit(node);self.function=prior;return node
    def visit_Return(self,node):
        if self.function in ['static_model','build','block_basis','projected','positive_quadratic_upper']:
            return [self.hook(self.function+'_return'),node]
        return node
    def visit_Assign(self,node):
        s=ast.unparse(node);tags=[]
        if self.function=='block_basis':
            if s=='V = left[:, keep]':tags=['seed']
            elif s=='block = H @ frontier':tags=['raw']
            elif s.startswith('original_scale = '):tags=['scale']
            elif isinstance(node.value,ast.Call) and isinstance(node.value.func,ast.Name) and node.value.func.id=='svd' and node.value.args and isinstance(node.value.args[0],ast.Name) and node.value.args[0].id=='block':tags=['svd']
            elif s=='left = left[:, keep]':tags=['eligible']
            elif s.startswith('left = left[:, :min('):tags=['prefix']
            elif s.startswith('V = np.column_stack'):tags=['admitted']
        return [node]+[self.hook(t) for t in tags] if tags else node
    def visit_For(self,node):
        node=self.generic_visit(node)
        if self.function=='block_basis' and ast.unparse(node.target)=='_' and ast.unparse(node.iter)=='range(2)':node.body.append(self.hook('orth'))
        if self.function=='projected' and ast.unparse(node.target)=='t':node.body.append(self.hook('time_bound'))
        return node

def strip_hooks(tree):
    class Strip(ast.NodeTransformer):
        def visit_Expr(self,node):
            if isinstance(node.value,ast.Call) and isinstance(node.value.func,ast.Name) and node.value.func.id=='_trace':return None
            return node
    return Strip().visit(tree)

class Trace:
    def __init__(self,out):self.out=Path(out);self.out.mkdir(parents=True,exist_ok=True);self.step=0;self.prev=None;self.project_count=0;self.quads=[];self.times=[];self.static=None;self.built=None;self.rows=[]
    def __call__(self,event,L):
        if event=='static_model_return':
            self.static={k:np.array(L[k],copy=True) for k in ['r0','edges','hess','vals','lam','B','sd']};return
        if event=='build_return':
            self.built={k:np.array(L[k],copy=True) for k in ['axis','delta','inds','energy','Ei','pi','means','F','times','sd']};return
        if event=='seed':
            save(self.out/'seed.npz',V=L['V'],s=L['s'],left=L['left'],right=L['_'],keep=L['keep']);return
        if event=='raw':
            self.step+=1;self.prev=L['block'].copy()
            save(self.out/f'step-{self.step:03d}-raw.npz',W0=L['block'],frontier=L['frontier'],incoming_rank=L['rank']);return
        if event=='scale':
            dump(self.out/f'step-{self.step:03d}.json',{'incoming_rank':L['rank'],'frontier_width':L['frontier'].shape[1],'original_scale':L['original_scale']});return
        if event=='orth':
            coefficient=mm(L['V'].T,self.prev)
            save(self.out/f'step-{self.step:03d}-orth{L["_"]+1}.npz',W=L['block'],projection_coefficients=coefficient)
            self.prev=L['block'].copy();return
        if event=='svd':
            rec=mm(L['left']*L['s'],L['_']);d=maxdiff(L['block'],rec)
            save(self.out/f'step-{self.step:03d}-svd.npz',left=L['left'],s=L['s'],right=L['_'],threshold=L['cutoff']*max(L['original_scale'],np.finfo(float).tiny),margin=L['s']-L['cutoff']*max(L['original_scale'],np.finfo(float).tiny),reconstruction_max_error=d);return
        if event in ['eligible','prefix']:
            save(self.out/f'step-{self.step:03d}-{event}.npz',left=L['left'],keep=L['keep']);return
        if event=='admitted':
            v=L['V'];m=mm(v.T,v)-np.eye(v.shape[1]);defect=L['F']-mm(v,mm(v.T,L['F']))
            piv=np.argmax(abs(L['left']),axis=0);ab=np.sort(abs(L['left']),axis=0)
            save(self.out/f'step-{self.step:03d}-admitted.npz',left=L['left'],pivots=piv,pivot_gap=ab[-1]-ab[-2])
            path=self.out/f'step-{self.step:03d}.json';row=json.loads(path.read_text());row.update(outgoing_rank=v.shape[1],admitted_width=L['left'].shape[1],eligible_width=int(L['keep'].sum()),prefix_cut=bool(int(L['keep'].sum())>L['left'].shape[1]),orthogonality_max=float(abs(m).max()),orthogonality_spectral=norm(m),seed_defect_frobenius=float(np.linalg.norm(defect)))
            dump(path,row);return
        if event=='positive_quadratic_upper_return':
            self.quads.append({k:L[k] for k in ['value','mass','allowance']});return
        if event=='time_bound':
            self.times.append({k:np.array(L[k],copy=True) for k in ['t','phi','psi','I','J','one','two','initial','static_error','response_bound','unif']});return
        if event=='projected_return':
            rank=L['V'].shape[1];self.project_count+=1
            p=self.out/f'checkpoint-{rank:03d}';p.mkdir(exist_ok=True)
            save(p/'matrices.npz',**asarrays(L,'V HV Hr coeff rot U lam lam_used Aused residual Rgram coupling mass_defect norms pnorms dnorms coefficient_norms initial static_error'),residual_parallel_matrix=mm(L['V'].T,L['residual']))
            save(p/'residual-svd.npz',singular_values=svd(L['residual'],compute_uv=False))
            for k,t in enumerate(self.times):save(p/f'time-{k}.npz',**t)
            dump(p/'quadratic-summations.json',self.quads);self.quads=[];self.times=[]
            save(p/'outputs.npz',**L['arrays']);dump(p/'row.json',L['row']);return
        if event=='block_basis_return':
            dump(self.out/'construction.json',{'selected':L['selected'],'records':L['records'],'steps':self.step,'final_rank':L['V'].shape[1]});return
        raise ValueError(event)

def instrument(trace):
    modules={}
    for name in ['reference_source','operator_study']:
        path=ROOT/'immutable'/(name+'.py');source=path.read_text();tree=Instrument().visit(ast.parse(source));ast.fix_missing_locations(tree)
        assert ast.dump(strip_hooks(ast.parse(ast.unparse(tree))),include_attributes=False)==ast.dump(ast.parse(source),include_attributes=False)
        counts={}
        for call in ast.walk(tree):
            if isinstance(call,ast.Call) and isinstance(call.func,ast.Name) and call.func.id=='_trace':
                tag=call.args[0].value;counts[tag]=counts.get(tag,0)+1
        expected={'static_model_return':1} if name=='reference_source' else {k:1 for k in ['build_return','positive_quadratic_upper_return','projected_return','time_bound','block_basis_return','seed','raw','scale','orth','svd','eligible','prefix','admitted']}
        assert counts==expected,(name,counts)
        module=types.ModuleType('traced_'+name);module.__file__=str(path);module.__dict__['_trace']=trace
        exec(compile(tree,str(path),'exec'),module.__dict__)
        modules[name]=module
    modules['operator_study'].reference=modules['reference_source'];return modules['operator_study']

def load_bundle(path):
    z=np.load(path,allow_pickle=False);h=csr_matrix((z['data'],z['indices'],z['indptr']),shape=tuple(z['shape']));return h,z['F'].copy(),z['times'].copy()
def case_info(index):return json.loads((ROOT/'protocol.json').read_text())['cases'][index]
def native_build(label,index):
    c=case_info(index);out=ROOT/'results'/label/c['name'];tr=Trace(out/'build-trace');study=instrument(tr);p=json.loads((ROOT/'immutable/preanalysis.json').read_text());g=next(x for x in p['geometries'] if x['id']==c['geometry'])
    h,f,t,sd,meta=study.build(p,g,c['kappa'],c['n'],4)
    save(out/'input.npz',data=h.data,indices=h.indices,indptr=h.indptr,shape=h.shape,**tr.built,**{'static_'+k:v for k,v in tr.static.items()})
    ho,fo,to,so,_=original.build(p,g,c['kappa'],c['n'],4)
    parity={'H_data':np.array_equal(h.data,ho.data),'H_indices':np.array_equal(h.indices,ho.indices),'H_indptr':np.array_equal(h.indptr,ho.indptr),'F':np.array_equal(f,fo),'times':np.array_equal(t,to),'sd':np.array_equal(sd,so)}
    assert all(parity.values());dump(out/'build.json',{'metadata':meta,'instrumentation_bitwise_parity':parity,'input_sha256':sha(out/'input.npz')})

def basis(label,index,shared=False,repair=False):
    c=case_info(index);inputlabel='pinned' if shared else label;kind='shared' if shared else 'native';out=ROOT/('trace-repair' if repair else 'results')/label/c['name']/kind;tr=Trace(out);study=instrument(tr)
    path=ROOT/'results'/inputlabel/c['name']/'input.npz';h,f,t=load_bundle(path);p=json.loads((ROOT/'immutable/preanalysis.json').read_text());snap,rows,selected,sec=study.block_basis(h,f,p,t)
    assert_trace(out)
    if repair:
        prior=ROOT/'results'/label/c['name']/kind
        equality={}
        for q in sorted(prior.rglob('*.npz')):
            newpath=out/q.relative_to(prior)
            if not newpath.exists() and q.name=='reference.npz':continue
            oldz=np.load(q);newz=np.load(newpath);assert set(oldz.files)==set(newz.files)
            equality[str(q.relative_to(prior))]=all(np.array_equal(oldz[k],newz[k]) for k in oldz.files)
        assert all(equality.values()),[k for k,v in equality.items() if not v]
        dump(out/'prior-trace-bitwise-parity.json',equality)
    if index==0 and not shared and not repair:
        old,oldrows,oldsel,_=original.block_basis(h,f,p,t);checks={str(rows[i]['rank'])+'_'+k:np.array_equal(a[k],old[i][k]) for i,a in enumerate(snap) for k in a}
        assert all(checks.values()) and oldsel==selected;dump(out/'instrumentation-parity.json',checks)
    # References are obtained only after the basis has stopped.
    ref,refsec=original.full_reference(h,f,t);save(out/'reference.npz',**ref)
    report={'input_sha256':sha(path),'selected':selected,'rank':rows[-1]['rank'],'checkpoints':[],'reference_seconds':refsec}
    for a,r in zip(snap,rows):
        err=abs(a['C']-ref['C']);report['checkpoints'].append({'rank':r['rank'],'max_C_error':float(err.max()),'max_bound':float(a['bound'].max()),'min_slack':float((a['bound']-err).min()),'selected_by_bound':r['bound_stop_pass']})
    dump(out/'evaluation.json',report)


def assert_trace(path):
    path=Path(path);j=json.loads((path/'construction.json').read_text());missing=[]
    for step in range(1,j['steps']+1):
        for kind in ['raw','orth1','orth2','svd','eligible','prefix','admitted']:
            p=path/('step-%03d-%s.npz'%(step,kind))
            if not p.exists():missing.append(str(p))
        if not missing:
            z=np.load(path/('step-%03d-svd.npz'%step));e=np.load(path/('step-%03d-eligible.npz'%step));assert np.array_equal(z['s']>z['threshold'],e['keep'])
            assert np.array_equal(z['left'][:,e['keep']],e['left'])
    for row in j['records']:
        p=path/('checkpoint-%03d'%row['rank'])
        for name in ['matrices.npz','residual-svd.npz','outputs.npz','row.json','quadratic-summations.json','time-0.npz','time-1.npz','time-2.npz']:
            if not (p/name).exists():missing.append(str(p/name))
        if not missing:assert len(json.loads((p/'quadratic-summations.json').read_text()))==18
    dump(path/'trace-completeness.json',{'complete':not missing,'steps':j['steps'],'checkpoints':len(j['records']),'missing':missing,'singular_values_and_eligible_vectors_match':not missing})
    assert not missing,missing

def eigensolvers(label,index):
    c=case_info(index);base=ROOT/'results'/'pinned'/c['name'];h,f,t=load_bundle(base/'input.npz')
    for cp in sorted((base/'native').glob('checkpoint-*')):
        v=np.load(cp/'matrices.npz')['V']
        for driver in ['evr','evd','ev']:
            out=ROOT/'results'/label/c['name']/'fixed-V-drivers'/cp.name/driver;tr=Trace(out);study=instrument(tr);study.eigh=lambda a,**kw:eigh(a,driver=driver,**kw)
            a,r=study.projected(h,f,v,t);save(out/'outputs.npz',**a);dump(out/'row.json',r)


def compare_runs(a,b,out):
    ad=Path(a);bd=Path(b);va=np.load(ad/'seed.npz')['V'];vb=np.load(bd/'seed.npz')['V'];rows=[{'step':0,'stage':'seed','byte_equal':np.array_equal(va,vb),'gap':gap(va,vb)}]
    indexes=sorted(int(p.name.split('-')[1]) for p in ad.glob('step-*-admitted.npz'))
    for idx in indexes:
        prefix=f'step-{idx:03d}'
        if not (bd/(prefix+'-admitted.npz')).exists():rows.append({'step':idx,'stage':'missing_admission'});break
        meta_a=json.loads((ad/(prefix+'.json')).read_text());meta_b=json.loads((bd/(prefix+'.json')).read_text())
        for name,key in [('raw','W0'),('orth1','W'),('orth2','W'),('eligible','left'),('prefix','left'),('admitted','left')]:
            aa=np.load(ad/(prefix+'-'+name+'.npz'))[key];bb=np.load(bd/(prefix+'-'+name+'.npz'))[key]
            item={'step':idx,'stage':name,'gap':gap(aa,bb),'shape_equal':aa.shape==bb.shape,'byte_equal':np.array_equal(aa,bb),'relative_frobenius_difference':float(np.linalg.norm(aa-bb)/max(np.linalg.norm(aa),np.finfo(float).tiny)) if aa.shape==bb.shape else None}
            rows.append(item)
        aa=np.load(ad/(prefix+'-admitted.npz'))['left'];bb=np.load(bd/(prefix+'-admitted.npz'))['left'];va=np.column_stack((va,aa));vb=np.column_stack((vb,bb))
        rows.append({'step':idx,'stage':'full_V','gap':gap(va,vb),'rank_a':va.shape[1],'rank_b':vb.shape[1],'metadata_a':meta_a,'metadata_b':meta_b})
    thresholds={}
    for threshold in [1e-12,1e-10,1e-8,1e-6]:
        thresholds[str(threshold)]=next(({'step':r['step'],'stage':r['stage'],'maximum_gap':r['gap']['maximum_gap']} for r in rows if 'gap' in r and r['gap']['maximum_gap']>threshold),None)
    dump(out,{'run_a':str(ad.relative_to(ROOT)),'run_b':str(bd.relative_to(ROOT)),'thresholds':thresholds,'rows':rows})

def comparisons():
    out=ROOT/'results/comparison';out.mkdir(exist_ok=True)
    for i in range(2):
        c=case_info(i);a=ROOT/'results/pinned'/c['name'];b=ROOT/'results/system'/c['name'];za=np.load(a/'input.npz');zb=np.load(b/'input.npz')
        checks={k:{'shape_equal':za[k].shape==zb[k].shape,'byte_equal':np.array_equal(za[k],zb[k]),'max_difference':maxdiff(za[k],zb[k]) if za[k].shape==zb[k].shape else None} for k in za.files}
        dump(out/(c['name']+'-input-comparison.json'),checks)
        compare_runs(a/'native',b/'native',out/(c['name']+'-native.json'))
        compare_runs(a/'native',b/'shared',out/(c['name']+'-shared.json'))
        compare_runs(b/'native',b/'shared',out/(c['name']+'-system-input-effect.json'))
        # Select the earliest shared-input residual span discrepancy above 1e-10.
        j=json.loads((out/(c['name']+'-shared.json')).read_text());r=next((r for r in j['rows'] if r['step']>0 and r['stage']=='eligible' and r['gap']['maximum_gap']>1e-10),None)
        if r is None:r={'step':1,'reason':'No eligible-span discrepancy above 1e-10; inspect first step only.'}
        dump(out/(c['name']+'-isolation-selection.json'),r)

def svd_isolation(label,index):
    c=case_info(index);sel=json.loads((ROOT/'results/comparison'/(c['name']+'-isolation-selection.json')).read_text());idx=sel['step'];p=ROOT/'results/pinned'/c['name']/'native';w=np.load(p/f'step-{idx:03d}-orth2.npz')['W'];meta=json.loads((p/f'step-{idx:03d}.json').read_text());threshold=1e-12*max(meta['original_scale'],np.finfo(float).tiny)
    for driver in ['gesdd','gesvd']:
        u,s,vh=svd(w,full_matrices=False,lapack_driver=driver);keep=s>threshold;eligible=u[:,keep];width=meta['admitted_width'];admitted=eligible[:,:width].copy()
        for col in admitted.T:
            if col[np.argmax(abs(col))]<0:col*=-1
        save(ROOT/'results'/label/c['name']/f'isolated-step-{idx:03d}-{driver}.npz',left=u,s=s,right=vh,keep=keep,eligible=eligible,admitted=admitted,reconstruction_max_error=maxdiff(w,mm(u*s,vh)))

def summarize():
    summaries=[]
    for i in range(2):
        c=case_info(i);cp=ROOT/'results/comparison';row={'case':c['name'],'native':json.loads((cp/(c['name']+'-native.json')).read_text())['thresholds'],'identical_input':json.loads((cp/(c['name']+'-shared.json')).read_text())['thresholds'],'evaluations':{},'fixed_V_driver_comparisons':[],'isolated_svd_comparisons':[]}
        for label,kind in [('pinned','native'),('system','native'),('system','shared')]:row['evaluations'][label+'_'+kind]=json.loads((ROOT/'results'/label/c['name']/kind/'evaluation.json').read_text())
        for p in sorted((ROOT/'results/pinned'/c['name']/'fixed-V-drivers').glob('checkpoint-*/evr/outputs.npz')):
            baseline=np.load(p);cpname=p.parents[1].name
            for label in ['pinned','system']:
                for driver in ['evr','evd','ev']:
                    q=ROOT/'results'/label/c['name']/'fixed-V-drivers'/cpname/driver/'outputs.npz';z=np.load(q)
                    row['fixed_V_driver_comparisons'].append({'label':label,'driver':driver,'rank':int(cpname.split('-')[1]),'C_max_difference':maxdiff(baseline['C'],z['C']),'bound_max_difference':maxdiff(baseline['bound'],z['bound'])})
        idx=json.loads((cp/(c['name']+'-isolation-selection.json')).read_text())['step'];base=np.load(ROOT/'results/pinned'/c['name']/f'isolated-step-{idx:03d}-gesdd.npz')
        for label in ['pinned','system']:
            for driver in ['gesdd','gesvd']:
                z=np.load(ROOT/'results'/label/c['name']/f'isolated-step-{idx:03d}-{driver}.npz');row['isolated_svd_comparisons'].append({'label':label,'driver':driver,'step':idx,'s_max_difference':maxdiff(base['s'],z['s']),'keep_equal':np.array_equal(base['keep'],z['keep']),'eligible_gap':gap(base['eligible'],z['eligible']),'admitted_gap':gap(base['admitted'],z['admitted'])})
        summaries.append(row)
    dump(ROOT/'results/summary.json',{'status':'DIAGNOSTIC COMPLETE; interpretation requires stage-specific comparison','cases':summaries})

def selftests(label):
    e=np.eye(7);v=e[:,:3];rot=np.array([[0.,-1.,0.],[1.,0.,0.],[0.,0.,-1.]])
    tests={'sign_gauge':gap(v,-v)['maximum_gap']<1e-12,'rotation_gauge':gap(v,mm(v,rot))['maximum_gap']<1e-12,'new_direction':gap(v,e[:,[0,1,4]])['maximum_gap']>.99,'unequal_nested':gap(v,v[:,:2])['maximum_gap']>.99}
    a=e[:,:2];b=mm(a,np.array([[2**-.5,-2**-.5],[2**-.5,2**-.5]]));tests['cluster_full_span']=gap(a,b)['maximum_gap']<1e-12;tests['cluster_prefix']=gap(a[:,:1],b[:,:1])['maximum_gap']>.7
    tests['rank_deficiency_reported']=gap(np.column_stack((e[:,0],e[:,0])),e[:,:1])['rank_deficient_input']
    tests['threshold_straddle']=bool(np.array_equal(np.array([2e-12,.5e-12])>1e-12,[True,False]))
    try:gap(np.full((7,2),np.nan),a);tests['nonfinite_rejection']=False
    except ValueError:tests['nonfinite_rejection']=True
    h=np.diag(np.arange(7,dtype=float));f=v[:,:2];before=original.projected(h,f,v,np.array([.1]))[0]['C'];after=original.projected(h+np.eye(7),f,v,np.array([.1]))[0]['C'];tests['fixed_basis_operator_change']=maxdiff(before,after)>1e-3
    assert all(tests.values());dump(ROOT/'results'/label/'selftests.json',tests)
    buf=io.StringIO()
    with contextlib.redirect_stdout(buf):
        print('NumPy configuration:', getattr(np.__config__, 'CONFIG', 'unavailable'))
        print('SciPy configuration:', getattr(scipy.__config__, 'CONFIG', 'unavailable'))
        if hasattr(np,'show_runtime'):np.show_runtime()
    pools=[]
    try:
        from threadpoolctl import threadpool_info
        pools=threadpool_info();assert all(x['num_threads']==1 for x in pools)
    except ImportError:pass
    dump(ROOT/'results'/label/'environment.json',{'python':platform.python_version(),'executable':sys.executable,'numpy':np.__version__,'scipy':scipy.__version__,'platform':platform.platform(),'machine':platform.machine(),'configuration':buf.getvalue(),'threadpools':pools,'float64':{'bits':np.finfo(np.float64).bits,'eps':float(np.finfo(np.float64).eps)},'longdouble':{'bits':np.finfo(np.longdouble).bits,'eps':float(np.finfo(np.longdouble).eps)}})

def generate_instrumentation_receipts():
    for name in ['operator_study','reference_source']:
        p=ROOT/'immutable'/(name+'.py');s=p.read_text();tree=Instrument().visit(ast.parse(s));ast.fix_missing_locations(tree);rendered=ast.unparse(tree)+'\n'
        assert ast.dump(strip_hooks(ast.parse(rendered)),include_attributes=False)==ast.dump(ast.parse(s),include_attributes=False)
        (ROOT/(name+'-instrumented.py.txt')).write_text(rendered)
        (ROOT/(name+'-instrumentation.diff')).write_text(''.join(difflib.unified_diff(s.splitlines(True),rendered.splitlines(True),fromfile=name+'.py',tofile=name+'-instrumented')))

if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('mode');a.add_argument('--label',default='pinned');a.add_argument('--case',type=int,default=0);a.add_argument('--shared',action='store_true');a.add_argument('--repair',action='store_true');args=a.parse_args();start=time.monotonic()
    if args.mode=='instrumentation':generate_instrumentation_receipts()
    elif args.mode=='inventory':selftests(args.label)
    elif args.mode=='build':native_build(args.label,args.case)
    elif args.mode=='basis':basis(args.label,args.case,args.shared,args.repair)
    elif args.mode=='eigh':eigensolvers(args.label,args.case)
    elif args.mode=='compare':comparisons()
    elif args.mode=='svd':svd_isolation(args.label,args.case)
    elif args.mode=='summary':summarize()
    else:raise ValueError(args.mode)
    print(json.dumps({'mode':args.mode,'label':args.label,'case':args.case,'seconds':time.monotonic()-start}),flush=True)
