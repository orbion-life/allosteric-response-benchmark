"""Render every primary finite instance from the immutable measured run.

This script only reads the protocol and results. It never runs or fits a model.
At the saved width of 174 mm, every text label is at least 10.4 pt.
"""
from pathlib import Path
import argparse, csv, hashlib, json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import FixedLocator, FixedFormatter, NullLocator

ROOT = Path(__file__).resolve().parent.parent
INK, VIOLET, BLUE, GRAY, MUTED = '#16161A', '#6336A7', '#2B6DA6', '#CCC8D1', '#56545A'
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10.4,
    'axes.labelsize': 10.4, 'xtick.labelsize': 10.4, 'ytick.labelsize': 10.4,
    'text.color': INK, 'axes.labelcolor': INK, 'xtick.color': INK,
    'pdf.fonttype': 42, 'svg.fonttype': 'none', 'axes.facecolor': 'white',
    'figure.facecolor': 'white', 'savefig.facecolor': 'white',
    'axes.unicode_minus': False})

def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()

def load_data(results):
    protocol = json.loads((ROOT/'preanalysis.json').read_text())
    run = json.loads((results/'run.json').read_text())
    assert run['status'] == 'COMPLETE'
    assert run['protocol_sha256'] == sha(ROOT/'preanalysis.json')
    source_files = [ROOT/'preanalysis.json', results/'run.json']
    rows, diagnostics = [], []
    for g in protocol['geometries']:
        for k in protocol['kappa']:
            for n in protocol['grid_n']:
                stem = f"{g['id']}-k{k}-n{n}-e4"
                jp, ap = results/(stem+'.json'), results/(stem+'.npz')
                item = json.loads(jp.read_text()); a = np.load(ap)
                source_files += [jp, ap]
                last = item['checkpoints'][-1]; r = last['rank']
                error = float(np.max(abs(a[f'rank{r}_C']-a['reference_C'])))
                bound = float(np.max(a[f'rank{r}_bound']))
                assert np.isclose(error, last['response_max_error'], atol=1e-15, rtol=0)
                assert np.isclose(bound, last['max_response_bound'], atol=1e-15, rtol=0)
                ratio = item['total_operator_seconds']/item['total_reference_seconds']
                assert np.isclose(ratio, item['total_time_ratio'], atol=1e-15, rtol=0)
                rows.append(dict(geometry=g['id'], kappa=k, grid_n=n, states=item['states'],
                    selected_or_last_rank=r, selected_by_bound=item['selected_by_bound'],
                    max_absolute_response_error=error, max_calculated_response_bound=bound,
                    total_operator_seconds=item['total_operator_seconds'],
                    total_reference_seconds=item['total_reference_seconds'], total_time_ratio=ratio,
                    finite_fidelity_pass=item['finite_fidelity_pass'], practical_pass=item['practical_pass']))
            arrays = []
            for n,e in [(17,4),(33,4),(41,5)]:
                ap = results/f"{g['id']}-k{k}-n{n}-e{e}.npz"
                source_files.append(ap); arrays.append(np.load(ap)['reference_C'])
            grid = float(np.max(abs(arrays[0]-arrays[1])))
            domain = float(np.max(abs(arrays[1]-arrays[2])))
            diagnostics.append(dict(geometry=g['id'], kappa=k, grid_response_difference=grid,
                domain_response_difference=domain, grid_pass=grid<=protocol['gates']['grid_response_difference'],
                domain_pass=domain<=protocol['gates']['domain_response_difference']))
    return protocol, run, rows, diagnostics, sorted(set(source_files))

def make_figure(results, output, replay_status):
    output.mkdir(parents=True, exist_ok=True)
    p, run, rows, diagnostics, sources = load_data(results)
    assert len(rows)==18 and len(diagnostics)==9
    if replay_status=='review complete':
        verification=ROOT/'audit/parent-direct-verification.json'
        review=ROOT/'audit/results-independent-review.md'
        assert verification.exists() and review.exists()
        assert json.loads(verification.read_text())['status']=='PASS'
        sources += [verification,review]
    finite = sum(x['finite_fidelity_pass'] for x in rows)
    practical = sum(x['practical_pass'] for x in rows)
    cost_only = sum(x['total_time_ratio']<=1 for x in rows)
    failed_grid = sum(not x['grid_pass'] for x in diagnostics)
    failed_domain = sum(not x['domain_pass'] for x in diagnostics)
    fig = plt.figure(figsize=(174/25.4, 118/25.4))
    ax = fig.add_axes([.255,.265,.398,.565])
    bx = fig.add_axes([.745,.265,.222,.565], sharey=ax)
    ys = [0,1,2,3.7,4.7,5.7,7.4,8.4,9.4]
    labels = [('development','Development'),('compact','Compact'),('elongated','Elongated')]
    for xaxis in [ax,bx]:
        xaxis.set_ylim(10.2,-.75)
        xaxis.set_yticks([])
        for spine in ['left','right','top']: xaxis.spines[spine].set_visible(False)
        xaxis.spines['bottom'].set_color(GRAY); xaxis.spines['bottom'].set_linewidth(.65)
        xaxis.tick_params(axis='x',length=3,width=.6,pad=5)
        xaxis.set_xscale('log'); xaxis.xaxis.set_minor_locator(NullLocator())
        for y in [2.85,6.55]: xaxis.axhline(y,color=GRAY,lw=.45,zorder=0)
    ax.set_xlim(3e-10,2e-2)
    ax.xaxis.set_major_locator(FixedLocator([1e-9,1e-6,1e-3]))
    ax.xaxis.set_major_formatter(FixedFormatter(['10⁻⁹','10⁻⁶','10⁻³']))
    bx.set_xlim(.07,1.6)
    bx.xaxis.set_major_locator(FixedLocator([.1,.3,1]))
    bx.xaxis.set_major_formatter(FixedFormatter(['0.1','0.3','1']))
    ax.axvline(.002,color=INK,lw=.8,ls=(0,(3,3)),zorder=0)
    bx.axvline(1,color=INK,lw=.8,ls=(0,(3,3)),zorder=0)
    for ig,(gid,gname) in enumerate(labels):
        center=ys[ig*3+1]
        fig_y=fig.transFigure.inverted().transform(ax.transData.transform((1e-9,center)))[1]
        fig.text(.012,fig_y,gname,fontsize=10.4,va='center',weight='medium')
        for ik,k in enumerate([1,10,100]):
            y=ys[ig*3+ik]
            fy=fig.transFigure.inverted().transform(ax.transData.transform((1e-9,y)))[1]
            fig.text(.229,fy,str(k),ha='right',va='center',fontsize=10.4)
            for n,offset,color,marker in [(17,-.19,BLUE,'o'),(33,.19,VIOLET,'s')]:
                row=next(r for r in rows if (r['geometry'],r['kappa'],r['grid_n'])==(gid,k,n))
                z=y+offset; e=row['max_absolute_response_error']; b=row['max_calculated_response_bound']
                ax.plot([e,b],[z,z],color=color,lw=.65,alpha=.75,zorder=2)
                ax.scatter(e,z,s=20,marker=marker,c=color,edgecolors=color,linewidths=.7,zorder=4)
                ax.scatter(b,z,s=28,marker=marker,facecolors='white',edgecolors=color,linewidths=1,zorder=5)
                bx.scatter(row['total_time_ratio'],z,s=24,marker=marker,c=color,edgecolors=color,linewidths=.7,zorder=4)
    fig.text(.255,.957,'a  Finite-response accuracy',weight='bold',fontsize=11.2,va='top')
    fig.text(.745,.957,'b  Total time ratio',weight='bold',fontsize=11.2,va='top')
    fig.text(.012,.860,'Geometry',fontsize=10.4,weight='bold')
    fig.text(.229,.860,'κ',ha='right',fontsize=10.4,weight='bold')
    fig.text(.255,.860,'Filled: error. Open: bound.',fontsize=10.4)
    fig.text(.745,.860,'Projection / full',fontsize=10.4)
    handles=[Line2D([],[],color=BLUE,marker='o',linestyle='none',markersize=4.8,label='17³ grid'),
             Line2D([],[],color=VIOLET,marker='s',linestyle='none',markersize=4.8,label='33³ grid')]
    fig.legend(handles=handles,loc='upper left',bbox_to_anchor=(.0,.927),ncol=2,
               frameon=False,borderaxespad=0,handletextpad=.55,columnspacing=1.25,fontsize=10.4)
    ax.text(.002,1.015,'0.002',transform=ax.get_xaxis_transform(),ha='center',va='bottom',fontsize=10.4)
    bx.text(1,1.015,'1',transform=bx.get_xaxis_transform(),ha='center',va='bottom',fontsize=10.4)
    ax.set_xlabel('Maximum across entries and 3 times',labelpad=8)
    bx.set_xlabel('Cold elapsed time',labelpad=8)
    fig.text(.012,.120,f'{finite}/18 meet finite-fidelity gates; {practical}/18 also meet the cost gate.',fontsize=10.4)
    fig.text(.012,.074,f'All {failed_grid} grid and {failed_domain} domain checks fail. Continuum accuracy is unresolved.',fontsize=10.4)
    fig.text(.012,.028,'Primary run. Independent replay '+replay_status+'.',fontsize=10.4,color=MUTED)
    fig.canvas.draw()
    renderer=fig.canvas.get_renderer(); overflow=[]
    for text in fig.findobj(matplotlib.text.Text):
        if not text.get_text() or not text.get_visible():continue
        assert text.get_fontsize()>=10.2, (text.get_text(),text.get_fontsize())
        box=text.get_window_extent(renderer)
        if box.x0 < -1 or box.y0 < -1 or box.x1 > fig.bbox.width+1 or box.y1 > fig.bbox.height+1:
            overflow.append(text.get_text())
    assert not overflow, overflow
    metadata={'Creator':'Project Pulsar; reproducible Matplotlib plot','Title':'Finite response-operator accuracy and total cost',
              'CreationDate':None,'ModDate':None}
    fig.savefig(output/'response-operator.pdf',metadata=metadata)
    fig.savefig(output/'response-operator.png',dpi=240)
    fig.savefig(output/'response-operator.svg',metadata={'Date':None})
    plt.close(fig)
    for name,items in [('instance-data.csv',rows),('convergence-data.csv',diagnostics)]:
        with (output/name).open('w',newline='') as f:
            writer=csv.DictWriter(f,fieldnames=list(items[0]));writer.writeheader();writer.writerows(items)
    receipt={'protocol_sha256':sha(ROOT/'preanalysis.json'),'result_run_sha256':sha(results/'run.json'),
      'plot_script_sha256':sha(Path(__file__)),'width_mm':174,'height_mm':118,'minimum_label_pt':10.4,
      'finite_instances':len(rows),'times_per_instance':3,'finite_fidelity_pass':finite,
      'cost_ratio_only_pass':cost_only,'combined_practical_pass':practical,
      'grid_diagnostic_failures':failed_grid,'domain_diagnostic_failures':failed_domain,
      'report_variant':{'file':'operator-preservation.pdf','width_mm':174,'height_mm':58,'minimum_label_pt':10.4,'scope_in_caption':True},
      'independent_replay_status':replay_status,'timing_scope':'Primary run; no timing confidence interval or repeated performance benchmark.',
      'run_wall_seconds':run['wall_seconds'],'matplotlib_version':matplotlib.__version__,'numpy_version':np.__version__,
      'source_sha256':{str(x.relative_to(ROOT)):sha(x) for x in sources}}
    (output/'figure-provenance.json').write_text(json.dumps(receipt,indent=2)+'\n')
    compact_figure(rows,output)
    print(json.dumps({k:receipt[k] for k in ['finite_instances','finite_fidelity_pass','cost_ratio_only_pass','combined_practical_pass','grid_diagnostic_failures','domain_diagnostic_failures']}))


def compact_figure(rows, output):
    """All 18 instances at the report width; detailed scope belongs in caption."""
    fig=plt.figure(figsize=(174/25.4,58/25.4))
    ax=fig.add_axes([.255,.192,.398,.615])
    bx=fig.add_axes([.745,.192,.222,.615],sharey=ax)
    for a in [ax,bx]:
        a.set_ylim(8.6,-.6);a.set_yticks([]);a.set_xscale('log')
        for spine in ['left','right','top']:a.spines[spine].set_visible(False)
        a.spines['bottom'].set_color(GRAY);a.spines['bottom'].set_linewidth(.5)
        a.tick_params(axis='x',length=2.5,width=.6,pad=3)
        a.xaxis.set_minor_locator(NullLocator())
        for y in [2.5,5.5]:a.axhline(y,color=GRAY,lw=.4,zorder=0)
    ax.set_xlim(3e-10,2e-2)
    ax.xaxis.set_major_locator(FixedLocator([1e-9,1e-6,.002]))
    ax.xaxis.set_major_formatter(FixedFormatter(['10⁻⁹','10⁻⁶','0.002']))
    bx.set_xlim(.07,1.6)
    bx.xaxis.set_major_locator(FixedLocator([.1,.3,1]))
    bx.xaxis.set_major_formatter(FixedFormatter(['0.1','0.3','1']))
    ax.axvline(.002,color=INK,lw=.7,ls=(0,(2,3)),zorder=0)
    bx.axvline(1,color=INK,lw=.7,ls=(0,(2,3)),zorder=0)
    for ig,(gid,gname) in enumerate([('development','Development'),('compact','Compact'),('elongated','Elongated')]):
        center=3*ig+1
        fy=fig.transFigure.inverted().transform(ax.transData.transform((1e-9,center)))[1]
        fig.text(.01,fy,gname,fontsize=10.4,va='center')
        for ik,k in enumerate([1,10,100]):
            y=3*ig+ik
            fy=fig.transFigure.inverted().transform(ax.transData.transform((1e-9,y)))[1]
            fig.text(.23,fy,str(k),fontsize=10.4,ha='right',va='center')
            for n,offset,color,marker in [(17,-.21,BLUE,'o'),(33,.21,VIOLET,'s')]:
                r=next(x for x in rows if (x['geometry'],x['kappa'],x['grid_n'])==(gid,k,n))
                z=y+offset;e=r['max_absolute_response_error'];b=r['max_calculated_response_bound']
                ax.plot([e,b],[z,z],color=color,lw=.55,alpha=.75,zorder=2)
                ax.scatter(e,z,s=8,marker=marker,c=color,edgecolors=color,linewidths=.5,zorder=3)
                ax.scatter(b,z,s=11,marker=marker,facecolors='white',edgecolors=color,linewidths=.7,zorder=4)
                bx.scatter(r['total_time_ratio'],z,s=10,marker=marker,c=color,edgecolors=color,linewidths=.5,zorder=3)
    fig.text(.255,.985,'a  Finite-response accuracy',fontsize=10.4,weight='bold',va='top')
    fig.text(.745,.985,'b  Total time ratio',fontsize=10.4,weight='bold',va='top')
    fig.text(.23,.835,'κ',fontsize=10.4,ha='right',weight='bold')
    fig.text(.255,.835,'Filled: error; open: bound',fontsize=10.4)
    fig.text(.745,.835,'Projection / full',fontsize=10.4)
    handles=[Line2D([],[],color=BLUE,marker='o',linestyle='none',markersize=3.8,label='17³'),
             Line2D([],[],color=VIOLET,marker='s',linestyle='none',markersize=3.8,label='33³')]
    fig.legend(handles=handles,loc='upper left',bbox_to_anchor=(.002,.986),ncol=1,frameon=False,
               handletextpad=.4,borderaxespad=0,labelspacing=.15,fontsize=10.4)
    fig.text(.255,.013,'Maximum over entries and 3 times',fontsize=10.4,va='bottom')
    fig.text(.745,.013,'Cold elapsed time',fontsize=10.4,va='bottom')
    fig.canvas.draw();renderer=fig.canvas.get_renderer()
    overflow=[]
    for t in fig.findobj(matplotlib.text.Text):
        if not t.get_text() or not t.get_visible():continue
        assert t.get_fontsize()>=10.2
        box=t.get_window_extent(renderer)
        if box.x0 < -1 or box.y0 < -1 or box.x1 > fig.bbox.width+1 or box.y1 > fig.bbox.height+1:overflow.append(t.get_text())
    assert not overflow,overflow
    meta={'Creator':'Project Pulsar; reproducible Matplotlib plot','Title':'Finite response-operator preservation and cost','CreationDate':None,'ModDate':None}
    fig.savefig(output/'operator-preservation.pdf',metadata=meta)
    fig.savefig(output/'operator-preservation.png',dpi=300)
    fig.savefig(output/'operator-preservation.svg',metadata={'Date':None})
    plt.close(fig)

if __name__=='__main__':
    cli=argparse.ArgumentParser();cli.add_argument('--results',type=Path,default=ROOT/'results')
    cli.add_argument('--output',type=Path,default=Path(__file__).resolve().parent)
    cli.add_argument('--replay-status',default='review pending',choices=['review pending','review complete'])
    args=cli.parse_args();make_figure(args.results.resolve(),args.output.resolve(),args.replay_status)
