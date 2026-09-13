"""Plot only archived CSV rows; never constructs a scientific model."""
from pathlib import Path
import csv,json,hashlib
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator,LogFormatterMathtext,NullLocator
ROOT=Path(__file__).resolve().parent
INK='#24242c';VIOLET='#6953a3';BLUE='#357584';GREY='#aaa4ad';PAPER='#fbfaf6';MUTED='#64636e'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11.2,'axes.titlesize':12,'axes.labelsize':11.2,'xtick.labelsize':10.5,'ytick.labelsize':11,'text.color':INK,'axes.labelcolor':INK,'xtick.color':MUTED,'ytick.color':MUTED,'axes.edgecolor':'#b5b1b8','axes.linewidth':.7,'pdf.fonttype':42,'svg.fonttype':'none'})
def rows(name):
 with (ROOT/name).open() as stream:return list(csv.DictReader(stream))
def clean(ax):
 for side in ['top','right']:ax.spines[side].set_visible(False)
 ax.set_facecolor(PAPER);ax.tick_params(which='both',length=3,pad=5)
fig=plt.figure(figsize=(8.8,10.0),facecolor=PAPER)
fig.text(.09,.964,'PROJECT PULSAR',fontsize=11.2,color=VIOLET,weight='bold')
fig.text(.09,.925,'Finite-model checks, numerical stability and complete costs',fontsize=14,weight='bold')
fig.text(.09,.892,'Computational evidence; physical convergence and biological accuracy remain open.',fontsize=11.2,color=MUTED)
a=fig.add_axes([.10,.605,.345,.220]);b=fig.add_axes([.61,.605,.345,.220]);c=fig.add_axes([.26,.19,.69,.22])
for ax in [a,b,c]:clean(ax)
# Harmonic error: all ten prescribed rows, with each fixed face kept separate.
h=rows('harmonic-all-cases.csv');styles={4:(BLUE,'o'),5:(GREY,'s'),6:(VIOLET,'o'),7:(INK,'x'),8:(INK,'+')}
for A in [4,5,6,7,8]:
 rr=sorted([r for r in h if int(r['A'])==A],key=lambda x:float(x['h']));color,marker=styles[A]
 a.plot([float(r['h']) for r in rr],[float(r['max_component_error']) for r in rr],color=color,marker=marker,lw=1.2,ms=5.5,label=f'A = {A}')
a.set_xscale('log',base=2);a.set_yscale('log');a.set_xlim(.012,.60);a.set_ylim(1e-5,6e-2);a.set_xticks([1/64,1/16,1/4],[r'$1/64$',r'$1/16$',r'$1/4$']);a.yaxis.set_major_locator(LogLocator(base=10,numticks=5));a.yaxis.set_major_formatter(LogFormatterMathtext());a.xaxis.set_minor_locator(NullLocator());a.yaxis.set_minor_locator(NullLocator());a.axhline(.001,color=GREY,lw=.8,ls=(0,(3,3)));a.text(.018,.00135,'0.001 screen',fontsize=10.5,color=MUTED)
a.set_xlabel('Grid spacing h');a.set_ylabel('Maximum normalized component error');a.set_title('A  Harmonic calibration',loc='left',pad=24,weight='bold');a.text(0,1.055,'Development, κ = 1; faces ±A',transform=a.transAxes,fontsize=10.5,color=MUTED)
a.legend(loc='upper left',bbox_to_anchor=(.095,.525),bbox_transform=fig.transFigure,ncol=3,frameon=False,fontsize=10.5,columnspacing=.8,handlelength=1.2,handletextpad=.4)
# Stability: native/shared final rank; fixed V is max over all checkpoints.
s=[r for r in rows('stability-all-cases.csv') if r['case']=='elongated-k100-n33-e4'];by={r['comparison']:r for r in s}
for i,key in enumerate(['native','shared','fixed_V']):
 rr=by[key];x=float(rr['C_difference']);z=float(rr['bound_difference']);b.plot([x,z],[2-i,2-i],color='#cbc7cd',lw=1.3);b.scatter([x],[2-i],color=VIOLET,s=34,zorder=3,label='Response ΔC' if i==0 else None);b.scatter([z],[2-i],color=INK,marker='D',s=30,zorder=3,label='Bound difference' if i==0 else None)
b.set_xscale('log');b.set_xlim(5e-16,1e-5);b.set_ylim(-.5,2.5);b.set_xticks([1e-15,1e-12,1e-9,1e-6]);b.xaxis.set_major_formatter(LogFormatterMathtext());b.xaxis.set_minor_locator(NullLocator());b.set_yticks([2,1,0],['Native H, F','Shared H, F','Fixed V']);b.axvline(1e-9,color=GREY,lw=.8,ls=(0,(3,3)));b.text(1.5e-9,2.38,'10⁻⁹',fontsize=10.5,color=MUTED);b.set_xlabel('Maximum absolute difference');b.set_title('B  Stability isolation',loc='left',pad=24,weight='bold');b.text(0,1.055,'Elongated, κ = 100; 33³',transform=b.transAxes,fontsize=10.5,color=MUTED);b.legend(loc='upper left',bbox_to_anchor=(.605,.525),bbox_transform=fig.transFigure,frameon=False,fontsize=10.5,handletextpad=.3)
# Full measured nonlinear invocation, decomposed without omitting overhead.
n=rows('nonlinear-all-cases.csv');keys=['preparation_and_array_serialization_seconds','Taylor_three_queries_seconds','Chebyshev_three_queries_seconds','other_validation_startup_shutdown_seconds','full_fresh_worker_seconds'];names=['Preparation + saving','Taylor: three queries','Chebyshev: three queries','Other checks + startup','Full worker invocation']
for j,row in enumerate(n):
 shift=.10 if j==0 else -.10;color=BLUE if j==0 else VIOLET;y=[4-i+shift for i in range(5)];values=[float(row[k]) for k in keys];c.scatter(values,y,s=35,color=color,marker='o' if j==0 else 's',label=f'h = {row["h"]}; {int(row["states"]):,} states',zorder=3)
 c.text(values[-1]*1.12,y[-1]-.025,f'{values[-1]:.2f} s',fontsize=10.5,color=color,va='center')
c.set_xscale('log');c.set_xlim(.009,90);c.set_ylim(-.55,4.5);c.set_yticks([4,3,2,1,0],names);c.set_xticks([.01,.1,1,10,100],['0.01','0.1','1','10','100']);c.xaxis.set_minor_locator(NullLocator());c.set_xlabel('Measured time (seconds; logarithmic scale)');c.set_title('C  Complete nonlinear profile',loc='left',pad=25,weight='bold');c.text(0,1.075,'Development, κ = 1; reflecting faces ±(4, 4, 17)',transform=c.transAxes,fontsize=10.5,color=MUTED);c.legend(loc='lower left',bbox_to_anchor=(-.02,-.38),frameon=False,fontsize=10.5,ncol=2,handletextpad=.3,columnspacing=1)
fig.text(.09,.067,'A: max error over G₀, K and C; harmonic tensor factorization only. A = 7, 8 overlap.',fontsize=10.5,color=MUTED)
fig.text(.09,.044,'B: native/shared use rank 48; fixed V spans every checkpoint and eigensolver driver.',fontsize=10.5,color=MUTED)
fig.text(.09,.021,'C: both methods share H, F, π and G₀. Only propagation is independently checked.',fontsize=10.5,color=MUTED)
for ext in ['pdf','png','svg']:fig.savefig(ROOT/('computational-evidence.'+ext),dpi=180,facecolor=PAPER)
plt.close(fig)
manifest={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(ROOT.glob('*.csv'))};manifest['plot_evidence.py']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest();(ROOT/'plot-input-hashes.json').write_text(json.dumps(manifest,indent=2)+'\n')
