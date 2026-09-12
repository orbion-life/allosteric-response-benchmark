#!/usr/bin/env python3
"""Paper-native 493 × 216 pt comparison; self-contained verified data JSON."""
from pathlib import Path
import hashlib,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
P=Path(__file__).resolve().parent
data=json.loads((P/'protein-comparison-data.json').read_text())
ink='#16161A';violet='#6336A7';blue='#2B6DA6';paper='#FDFDFB';gray='#CCC8D1';muted='#56545A'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10.4,'svg.fonttype':'none','pdf.fonttype':42,'figure.facecolor':paper,'axes.facecolor':paper,'text.color':ink})
fig=plt.figure(figsize=(493/72,216/72));ax=fig.add_axes([0,0,1,1]);ax.set_xlim(0,493);ax.set_ylim(0,216);ax.set_axis_off()
def text(x,y,s,**kwargs):return ax.text(x,y,s,fontsize=10.4,**kwargs)
projected={};views={}
for name,t in data['targets'].items():
    xyz=np.array(t['coordinates_A']);center=xyz.mean(axis=0);_,_,v=np.linalg.svd(xyz-center,full_matrices=False)
    for j in range(3):
        if v[j,np.argmax(np.abs(v[j]))]<0:v[j]*=-1
    projected[name]=(xyz-center)@v[:2].T
    views[name]={'center_A':center.tolist(),'orthographic_view_rows':v[:2].tolist(),'method':'Principal-axis orthographic view, with deterministic largest-component signs; no conformation change.'}
# Both traces use exactly the same point-per-Å display scale.
maxwidth=max(np.ptp(x[:,0]) for x in projected.values());maxheight=max(np.ptp(x[:,1]) for x in projected.values());scale=min(88/maxwidth,119/maxheight)
for letter,name,cx in [('a','KRAS',57),('b','ABL',168)]:
    t=data['targets'][name];xy=projected[name];xy-=.5*(xy.max(axis=0)+xy.min(axis=0));xy=xy*scale+np.array([cx,112]);ids=np.array(t['canonical']);receiver=np.array(t['receiver_indices']);top=np.array([np.flatnonzero(ids==p)[0] for p in t['primary_top5']])
    text(cx-47,202,letter+'  '+name+' · '+t['input_pdb'],fontweight='bold',va='center')
    for i in range(len(ids)-1):
        if ids[i+1]==ids[i]+1:ax.plot(xy[i:i+2,0],xy[i:i+2,1],color=gray,lw=.65,zorder=1)
    ax.scatter(xy[receiver,0],xy[receiver,1],s=9,color=blue,linewidths=0,zorder=2)
    ax.scatter(xy[top,0],xy[top,1],s=18,color=violet,marker='D',linewidths=.3,edgecolors=paper,zorder=3)
text(239,202,'c  Reference labels in the top five',fontweight='bold',va='center')
text(365,180,'KRAS',ha='center');text(455,180,'ABL',ha='center')
row_ys=[160,141,122,103,84,65]
for y,row in zip(row_ys,data['rows']):
    text(239,y,row['method'],va='center',fontweight='bold' if row['method']=='Biquadratic' else 'normal')
    for target,start in [('KRAS',344),('ABL',434)]:
        for rank,state in enumerate(row[target]['states']):
            x=start+rank*10.5
            if state=='contact':ax.scatter(x,y,s=23,color=violet,edgecolors=violet,linewidths=.65,zorder=4)
            elif state=='observed_noncontact':ax.scatter(x,y,s=23,facecolors=paper,edgecolors=muted,linewidths=.7,zorder=4)
            elif state=='unknown':text(x,y,'?',ha='center',va='center',fontweight='bold')
            else:raise ValueError(state)
ax.scatter(13,35,s=13,color=blue);text(21,35,'Receiver',va='center')
ax.scatter(101,35,s=18,color=violet,marker='D');text(110,35,'Primary top five',va='center')
text(10,12,'Cα traces; common display scale.',va='center')
ax.scatter(244,35,s=23,color=violet);text(252,35,'≤5 Å',va='center')
ax.scatter(313,35,s=23,facecolors=paper,edgecolors=muted,linewidths=.7);text(321,35,'>5 Å',va='center')
text(386,35,'?',ha='center',va='center',fontweight='bold');text(395,35,'Unresolved',va='center')
text(239,12,'Left to right follows shortlist rank.',va='center')
fig.canvas.draw();renderer=fig.canvas.get_renderer();issues=[]
for t in ax.texts:
    box=t.get_window_extent(renderer)
    if box.x0<-.25 or box.y0<-.25 or box.x1>fig.bbox.width+.25 or box.y1>fig.bbox.height+.25:issues.append(t.get_text())
assert not issues,issues
for ext in ['pdf','svg','png']:fig.savefig(P/f'protein-comparison.{ext}',dpi=230)
views.update({'shared_points_per_A':scale,'figure_size_pt':[493,216],'minimum_text_font_pt':10.4,'no_text_clipped_at_figure_boundary':True,'structure_scope':'Static retained input nodes only; reference ligands are not drawn or overlaid.','label_scope':'Filled marker: a reference heavy-atom contact≤5Å. Open marker: observed beyond5Å. Question mark: unresolved reference position. The five markers follow prediction order.'})
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
(P/'protein-comparison-receipt.json').write_text(json.dumps({'views_and_checks':views,'source_data_sha256':sha(P/'protein-comparison-data.json'),'code_sha256':sha(Path(__file__)),'outputs':{f'protein-comparison.{ext}':sha(P/f'protein-comparison.{ext}') for ext in ['pdf','svg','png']}},indent=2)+'\n')
print(json.dumps(views,indent=2))
