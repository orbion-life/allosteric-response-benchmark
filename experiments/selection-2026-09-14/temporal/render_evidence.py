"""Standalone scientific figure from saved frozen-grid metrics only."""
from pathlib import Path
import json,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parent
p=json.loads((ROOT/'protocol.json').read_text());t=np.array(p['observation_times_over_tau']);keep=t>0
fig,axes=plt.subplots(1,2,figsize=(10,3.6),sharey=True)
colors=['#163E63','#C06529','#277466','#786682']
for target,color in zip(p['targets'],colors):
 s=json.loads((ROOT/'results'/target/'analysis.json').read_text());a=np.array(s['step']['error_by_time']);b=np.max([x['error_by_time'] for x in s['pulses']],axis=0)
 axes[0].plot(t[keep],a[keep],color=color,lw=1.5,label=target.upper());axes[1].plot(t[keep],b[keep],color=color,lw=1.5)
for ax,title in zip(axes,['Sustained perturbation','Largest error across three pulse durations']):
 ax.set_xscale('log');ax.set_yscale('log');ax.set_ylim(1e-9,.1);ax.axhline(.002,color='#333333',lw=.8,ls='--');ax.text(0.99,.002,' 0.002 threshold',transform=ax.get_yaxis_transform(),ha='right',va='bottom',fontsize=8);ax.set_title(title,fontsize=10,loc='left');ax.set_xlabel('Model time, t / τ',fontsize=9);ax.spines[['top','right']].set_visible(False);ax.tick_params(labelsize=8)
axes[0].set_ylabel('Maximum absolute normalized response error',fontsize=9);axes[0].legend(frameon=False,fontsize=8,loc='lower left',ncol=2)
fig.suptitle('Fixed Gaussian reductions fail the early-time response threshold',fontsize=12,x=.08,ha='left');fig.tight_layout();(ROOT/'figures').mkdir(exist_ok=True)
fig.savefig(ROOT/'figures/temporal-error.pdf',bbox_inches='tight');fig.savefig(ROOT/'figures/temporal-error.png',dpi=180,bbox_inches='tight')
