"""Build a self-contained viewer from frozen measured coordinates and scores."""
import json
from prepare import ROOT,sha
import numpy as np

def main():
    meta=json.loads((ROOT/'model/static-model.json').read_text());m=np.load(ROOT/'model/static-model.npz');p=np.load(ROOT/'results/biquadratic-d2-n33-e4.npz');top=p['order'][:5].tolist();rec=set(m['receiver'].tolist())
    nodes=[{'id':int(m['canonical'][i]),'name':x['resname'],'xyz':m['r0'][i].tolist(),'score':float(p['score'][i]),'receiver':i in rec,'eligible':bool(m['candidate'][i]),'topRank':top.index(i)+1 if i in top else None,'receiverDistance':float(m['receiver_distance'][i])} for i,x in enumerate(meta['node_map'])]
    data={'nodes':nodes,'input':'4OBE chain A','isoform':'KRAS4B / P01116-2','top5':[nodes[i]['id'] for i in top],'inputSha256':sha(ROOT/'raw/4OBE.pdb'),'predictionSha256':sha(ROOT/'results/biquadratic-d2-n33-e4.npz'),'protocolSha256':sha(ROOT/'preanalysis-protocol.json')}
    template=(ROOT/'viewer-template.html').read_text();html=template.replace('__VERIFIED_DATA__',json.dumps(data,separators=(',',':')));path=ROOT/'figures/kras-structure-viewer.html';path.write_text(html)
    (ROOT/'figures/viewer-provenance.json').write_text(json.dumps({'html_sha256':sha(path),'builder_sha256':sha(__file__),'template_sha256':sha(ROOT/'viewer-template.html'),'source':data,'geometry':'Only observed 4OBE C-alpha coordinates; consecutive residue positions are joined as a backbone trace. No reference ligand or calculated conformation is drawn. Rotation is a rigid orthographic projection.','network':'No CDN, fetch, external fonts or runtime dependencies.'},indent=2)+'\n');print(sha(path))

if __name__=='__main__':main()
