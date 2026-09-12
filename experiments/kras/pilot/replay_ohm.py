"""Rebuild the pinned GPL-3.0 comparator and reproduce its primary output.

Requires make, patch, zlib, and a GCC15+ compiler with C++23 support.
This does not install software or modify global environment settings.
"""
from pathlib import Path
import argparse,subprocess,tarfile,shutil,json,datetime,os
from prepare import ROOT,sha

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--compiler',default='g++-15');args=ap.parse_args();work=ROOT/'replay-ohm';source=work/'source';run=work/'run';work.mkdir(exist_ok=True);commands=[]
    def call(argv,cwd):
        p=subprocess.run(argv,cwd=cwd,check=True,capture_output=True,text=True);commands.append({'argv':argv,'cwd':str(cwd.relative_to(ROOT)),'stdout':p.stdout,'stderr':p.stderr})
    if not source.exists():
        source.mkdir()
        with tarfile.open(ROOT/'external/ohm-upstream-462a3b1.tar.gz') as tar:tar.extractall(source,filter='data')
        call(['patch','-p1','-i',str(ROOT/'external/ohm/gcc15-portability.patch')],source)
    # A relative prefix keeps upstream Makefile targets valid in folders with spaces.
    (source/'.prefix').write_text('install\n')
    call(['make','-j2','contacts','diffuse','CC='+args.compiler,'FLAGS=-std=c++23 -Wno-everything -lstdc++ -Wwrite-strings -lm -Isrc -MMD -DNDEBUG -O3'],source)
    run.mkdir(exist_ok=True);shutil.copyfile(ROOT/'external/ohm/kras.pdb',run/'kras.pdb')
    call([str(source/'install/bin/contacts'),'kras.pdb','-oi','kras.ind','-oc','kras.mat','-c','3.4'],run)
    sites='11+12+13+14+15+16+17+18+28+30+32+116+117+119+120+145+146+147'
    call([str(source/'install/bin/diffuse'),'aci','kras.mat',sites,'primary.nodes','primary.bonds','-n','10000','-c','0.05','-norm','-a','4.5','-seed','11','-mdist','primary.mdist'],run)
    comparisons={name:{'expected_sha256':sha(ROOT/'external/ohm'/name),'replayed_sha256':sha(run/name),'bitwise_equal':(ROOT/'external/ohm'/name).read_bytes()==(run/name).read_bytes()} for name in ['kras.mat','primary.nodes','primary.bonds','primary.mdist']}
    receipt={'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'upstream_commit':'462a3b1318e24ebe2061137fe88af719ef89c0ae','upstream_archive_sha256':sha(ROOT/'external/ohm-upstream-462a3b1.tar.gz'),'patch_sha256':sha(ROOT/'external/ohm/gcc15-portability.patch'),'comparisons':comparisons,'commands':commands,'caveat':'Exact numerical output may vary across compiler, standard library or architecture. The recorded successful repeat uses this host compiler; no cross-platform bitwise portability is claimed.'}
    (ROOT/'ohm-replay-verification.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(comparisons,indent=2))
    assert all(x['bitwise_equal'] for x in comparisons.values())

if __name__=='__main__':main()
