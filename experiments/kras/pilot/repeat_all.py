"""Repeat preparation and every numeric pilot case in a separate local folder."""
from pathlib import Path
import json,shutil,subprocess,sys,datetime,time
import numpy as np
from prepare import ROOT,sha

def main():
    work=ROOT/'repeat-all';work.mkdir(exist_ok=True);(work/'model').mkdir(exist_ok=True);shutil.copytree(ROOT/'raw',work/'raw',dirs_exist_ok=True)
    for name in ['prepare.py','run_pilot.py','preanalysis-protocol.json','preanalysis-protocol.sha256']:shutil.copyfile(ROOT/name,work/name)
    begin=time.monotonic();runs=[]
    for script in ['prepare.py','run_pilot.py']:
        p=subprocess.run([sys.executable,'-W','error',script],cwd=work,capture_output=True,text=True,check=True);(work/(script+'.log')).write_text(p.stdout+p.stderr);runs.append({'script':script,'code_sha256':sha(ROOT/script)})
    results=[]
    for path in sorted((ROOT/'results').glob('*.npz')):
        a=np.load(path);b=np.load(work/'results'/path.name);assert set(a.files)==set(b.files)
        diffs={k:float(np.max(np.abs(a[k].astype(float)-b[k].astype(float)))) if a[k].size else 0. for k in a.files}
        exact=all(np.array_equal(a[k],b[k]) for k in a.files);assert max(diffs.values())<1e-12
        results.append({'file':path.name,'all_arrays_bitwise_identical':exact,'max_absolute_difference':max(diffs.values())})
    report={'repeated_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'seconds':time.monotonic()-begin,'scope':'fresh model preparation and all24 finite-grid cases, eight all-mode/retained-mode analytic cases and graph/geometric controls; same host and pinned Python environment','runs':runs,'result_files':results,'all_arrays_bitwise_identical':all(x['all_arrays_bitwise_identical'] for x in results),'cross_platform_guarantee':False}
    (ROOT/'full-repeat-verification.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({'files':len(results),'all_identical':report['all_arrays_bitwise_identical'],'seconds':report['seconds']}))

if __name__=='__main__':main()
