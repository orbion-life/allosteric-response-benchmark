from pathlib import Path
import json,numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=Path(__file__).resolve().parents[1];s=json.loads((P/'results/summary.json').read_text());out=P/'figures';out.mkdir(exist_ok=True)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.spines.top':False,'axes.spines.right':False,'axes.spines.left':False,'axes.edgecolor':'#999999','axes.labelcolor':'#303030','xtick.color':'#555555','ytick.color':'#303030','svg.fonttype':'none'})
labels={'gaussian':'Gaussian response','harmonic':'Harmonic response','ohm':'Ohm','degree':'Contact degree','negative_distance':'Receiver proximity'}
fig,axs=plt.subplots(1,2,figsize=(10.6,4.2),gridspec_kw={'width_ratios':[1,1.15]})
for i,r in enumerate(s['metrics']):
 lo,hi=r['rho_site_bootstrap95'];y=4-i;color='#126d78' if i==0 else '#777777';axs[0].plot([lo,hi],[y,y],color=color,lw=1.4);axs[0].scatter(r['rho'],y,color=color,s=26,zorder=3);axs[0].text(1.02,y,f"{r['rho']:.3f}",va='center',ha='left',color=color,transform=axs[0].get_yaxis_transform(),clip_on=False)
axs[0].set_yticks(range(5),[labels[r['method']] for r in s['metrics'][::-1]]);axs[0].set_xlim(-1,1);axs[0].set_xticks([-1,-.5,0,.5,1]);axs[0].set_ylim(-.65,4.65);axs[0].axvline(0,color='#aaaaaa',lw=.7,ls=':');axs[0].set_xlabel('Spearman correlation with binding sensitivity');axs[0].set_title('Observed association and 95% intervals',loc='left',pad=18)
for i,r in enumerate(s['gaussian_minus_controls']):
 lo,hi=r['paired_site_bootstrap98_75'];y=3-i;axs[1].plot([lo,hi],[y,y],color='#777777',lw=1.4);axs[1].scatter(r['delta_rho'],y,color='#126d78',s=26,zorder=3)
axs[1].set_yticks(range(4),[labels[r['comparator']] for r in s['gaussian_minus_controls'][::-1]]);axs[1].set_xlim(-.65,1.3);axs[1].set_xticks([-.5,0,.5,1]);axs[1].set_ylim(-.65,3.65);axs[1].axvline(0,color='#777777',lw=.8,ls=':');axs[1].set_xlabel('Gaussian − comparator correlation');axs[1].set_title('Paired differences and 98.75% intervals',loc='left',pad=18)
for ax in axs:ax.tick_params(axis='y',length=0,pad=8);ax.grid(axis='x',color='#eeeeee',lw=.5);ax.set_axisbelow(True)
fig.suptitle('GRB2: moderate association; incremental benefit remains unestablished',x=.02,ha='left',fontsize=12,fontweight='medium',y=.99)
fig.text(.02,.045,'21 distal residues · 20,000 paired residue-bootstrap draws · one inferred GAB2-binding endpoint',fontsize=8,color='#555555')
fig.text(.02,.01,'Intervals condition on observed site aggregates; spatial dependence and shared fit uncertainty limit interpretation.',fontsize=8,color='#555555')
fig.subplots_adjust(left=.17,right=.98,top=.76,bottom=.24,wspace=.65)
for ext in ['png','svg','pdf']:fig.savefig(out/f'grb2-control-comparison.{ext}',dpi=220,bbox_inches='tight')
plt.close(fig)
