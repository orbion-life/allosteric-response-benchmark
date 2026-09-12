#!/usr/bin/env python3
"""Data-only scientific figures; all geometry comes from the locked 1OPL model."""
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
from run import ROOT
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.labelsize':11,'axes.titlesize':12,'xtick.labelsize':10,'ytick.labelsize':10,'svg.fonttype':'none','pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False,'figure.facecolor':'#fffefb','axes.facecolor':'#fffefb','text.color':'#272330','axes.labelcolor':'#272330'})
violet='#6d4698';ink='#272330';gray='#99949d';teal='#3c7978'
model=np.load(ROOT/'model/static-model.npz');ev=json.loads((ROOT/'evaluation/evaluation.json').read_text());diag=json.loads((ROOT/'results/diagnostics.json').read_text());xyz=model['r0'];ids=model['canonical'];rec=model['receiver'];top=[r['canonical'] for r in ev['methods']['biquadratic_d2_n33']['top5']];ti=np.array([np.flatnonzero(ids==p)[0] for p in top])
fig=plt.figure(figsize=(12.8,4.8));gs=fig.add_gridspec(1,3,width_ratios=[1.1,1.38,1],left=.025,right=.985,bottom=.35,top=.9,wspace=.5)
ax=fig.add_subplot(gs[0],projection='3d');ax.plot(*xyz.T,color='#c4c0c7',lw=1);ax.scatter(*xyz[rec].T,color=teal,s=13,depthshade=False);ax.scatter(*xyz[ti].T,color=violet,s=32,depthshade=False)
for i in ti:ax.text(*xyz[i],str(ids[i]),fontsize=10,color=ink)
ax.set_axis_off();span=np.ptp(xyz,axis=0);ax.set_box_aspect(span);ax.view_init(18,-67);ax.set_title('a  Frozen input and shortlist',loc='left',pad=18)
ax.legend(handles=[Line2D([],[],marker='o',ls='',color=teal,label='Input P16 receiver'),Line2D([],[],marker='o',ls='',color=violet,label='Primary top five')],frameon=False,loc='lower center',bbox_to_anchor=(.48,-.17),fontsize=10)
names=['biquadratic_d2_n33','harmonic_d2_n33','Hookean_d2_n33','equilibrium_d2_n33','all_mode_harmonic','graph_heat_kernel','degree'];labels=['Biquadratic','Harmonic','Distance-Hookean','Equilibrium','All-mode harmonic','Graph diffusion','Degree']
bx=fig.add_subplot(gs[1]);ys=np.arange(len(names))[::-1]
for y,name in zip(ys,names):
    r=ev['methods'][name];q=r['percentile_enrichment_null']['null_q025_median_q975'];bx.plot([q[0],q[2]],[y,y],color='#bdb8c2',lw=2,zorder=1);bx.scatter(q[1],y,marker='|',s=90,color=gray);bx.scatter(r['mean_contact_percentile'],y,color=violet if name==names[0] else ink,s=35,zorder=3)
bx.set_yticks(ys,labels);bx.set_xlim(0,1);bx.set_xticks([0,.25,.5,.75,1]);bx.set_xlabel('Mean rank percentile of 20 AY7 contacts');bx.set_title('b  Known contacts rank poorly',loc='left',pad=18);bx.tick_params(axis='y',length=0);bx.spines['left'].set_visible(False);bx.set_ylim(-.7,6.7)
bx.text(.0,-.25,'Dots: observed. Gray: matched-label\nmedian and central 95% null interval.\nLarger ranks indicate higher priority.',transform=bx.transAxes,fontsize=10,va='top')
cx=fig.add_subplot(gs[2]);values=[diag['grid_receiver_block_max_C_difference'],diag['coordinate_receiver_block_max_C_difference']];limits=[.001,.002]
for y,value,limit in zip([1,0],values,limits):cx.plot([limit,value],[y,y],color=gray,lw=1.4);cx.scatter(value,y,color=violet,s=40);cx.scatter(limit,y,marker='|',color=ink,s=130);cx.annotate(f'{value:.2g}',(value,y),xytext=(0,12),textcoords='offset points',ha='center',fontsize=10)
cx.set_xscale('log');cx.set_xlim(4e-6,.6);cx.set_yticks([1,0],['Grid refinement','Mode truncation']);cx.set_ylim(-.65,1.65);cx.set_xlabel('Maximum receiver-block |ΔC|');cx.set_title('c  Grid passes; truncation fails',loc='left',pad=18);cx.tick_params(axis='y',length=0);cx.spines['left'].set_visible(False);cx.text(0,-.25,'Dots: error; marks: locked limits.\nGrid: 33² versus 65².\nModes: two versus 750.\nHarmonic dynamics for mode check.',transform=cx.transAxes,fontsize=10,va='top')
fig.text(.035,.025,'1OPL kinase-domain coordinates only; no superposed 5MO4 ligand is drawn. Residue 385 is absent from the reference.\nThe primary shortlist contains 0 known contacts among 4 observed residues and 1 unknown; contact recovery is not functional validation.',fontsize=10)
for ext in ['svg','pdf','png']:fig.savefig(ROOT/f'figures/abl-pilot-summary.{ext}',dpi=180)
plt.close(fig)
C=np.load(ROOT/'results/biquadratic-d2-n33-e4.npz')['C'];vmax=np.abs(C).max();fig,ax=plt.subplots(figsize=(7.5,6.5));cmap=LinearSegmentedColormap.from_list('pulsar',[violet,'#fffefb',teal]);im=ax.imshow(C,origin='lower',extent=[ids[0]-.5,ids[-1]+.5,ids[0]-.5,ids[-1]+.5],cmap=cmap,vmin=-vmax,vmax=vmax,interpolation='none');ax.set_xlabel('P00519 canonical residue');ax.set_ylabel('P00519 canonical residue');ax.set_title('ABL: complete two-mode response matrix, 33² grid',pad=14);fig.colorbar(im,ax=ax,label='C (all-mode harmonic SD normalization)',fraction=.046,pad=.04);fig.tight_layout()
for ext in ['svg','pdf','png']:fig.savefig(ROOT/f'figures/abl-response-matrix.{ext}',dpi=180)
plt.close(fig)
print('Wrote two SVG/PDF/PNG scientific figures.')
