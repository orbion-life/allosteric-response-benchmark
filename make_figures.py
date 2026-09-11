"""Publication-style figures generated only from recorded numerical data."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PURPLE, BLUE, INK, PAPER = "#6336A7", "#2B6DA6", "#16161A", "#FDFDFB"

def make_figures(output, destination):
    destination.mkdir(parents=True, exist_ok=True)
    raw=json.loads((output/"results.json").read_text())
    data=json.loads((output/"sensitivity-results.json").read_text())
    plt.rcParams.update({"font.family":"DejaVu Serif", "font.size":10, "axes.spines.top":False,
        "axes.spines.right":False, "axes.facecolor":PAPER, "figure.facecolor":PAPER,
        "text.color":INK, "axes.labelcolor":INK, "svg.hashsalt":"pulsar-benchmark-v1"})
    def save(fig, name):
        for extension in ("svg", "png"):
            fig.savefig(destination/f"{name}.{extension}", bbox_inches="tight", dpi=180,
                        metadata={"Date":None} if extension=="svg" else None)
        plt.close(fig)
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.4), gridspec_kw={"width_ratios":[1,1.35]})
    r=np.array(raw["parameters"]["positions"])
    # A labelled linear projection; edges/coordinates come from the fixture.
    xy=np.column_stack((r[:,0]+.38*r[:,1],r[:,2]+.45*r[:,1]))
    ax=axes[0]
    for i,j in raw["parameters"]["edges_zero_based"]:
        ax.plot(xy[[i,j],0],xy[[i,j],1], color="#A8A5AF", lw=1.2, zorder=1)
    for i,(x,y) in enumerate(xy):
        ax.scatter(x,y,s=340,color=PURPLE if i==3 else BLUE if i==4 else INK,zorder=2)
        ax.text(x,y,str(i),ha="center",va="center",color="white",fontsize=11,zorder=3)
    ax.set_title("a  Five nodes · nine contacts",loc="left",fontweight="bold")
    ax.text(.5,-.10,"3 → 4: no direct contact; three shared neighbours",transform=ax.transAxes,ha="center",fontsize=8.5)
    ax.margins(.24);ax.set_aspect("equal");ax.axis("off")
    ax=axes[1]
    for model,color,label in [("harmonic",BLUE,"Harmonic"),("biquadratic",PURPLE,"Biquadratic")]:
        case=next(x for x in data["strength"] if x["kappa"]==100 and x["model"]==model)
        curve=case["time_curve"]
        ax.plot([x["t_over_tau"] for x in curve],[x["C43"] for x in curve],"o-",color=color,ms=3,label=label)
    ax.axvline(1,color="#A8A5AF",ls=":",lw=.8)
    ax.set(xlabel="Time / harmonic relaxation time",ylabel="Normalized response C₄₃",xlim=(0,10))
    ax.set_title("b  Same intervention, different response",loc="left",fontweight="bold")
    ax.legend(frameon=False,fontsize=9)
    fig.suptitle("PROJECT PULSAR   |   Synthetic finite-time response",x=.08,ha="left",fontsize=13)
    fig.tight_layout();save(fig,"overview")
    fig, axes = plt.subplots(1,3,figsize=(12,3.7))
    for model,color,label in [("harmonic",BLUE,"Harmonic"),("biquadratic",PURPLE,"Biquadratic")]:
        cases=[x for x in data["grid_kappa100"] if x["model"]==model]
        axes[0].plot([x["n"] for x in cases],[x["C43"] for x in cases],"o-",color=color,ms=4,label=label)
        cases=[x for x in data["fixed_spacing_domain_kappa100"] if x["model"]==model]
        axes[1].plot([x["extent"] for x in cases],[x["C43"] for x in cases],"o-",color=color,ms=4)
        for d,marker in [(2,"o"),(3,"s"),(4,"^")]:
            cases=[x for x in data["dimension"] if x["model"]==model and x["dimensions"]==d]
            axes[2].plot([x["n"] for x in cases],[x["C43_common_d4"] for x in cases],marker+"-",color=color,alpha=1 if d==2 else .7,ms=4,label=f"{label}, d={d}")
    for ax,title in zip(axes,["a  Grid spacing","b  Domain extent","c  Retained coordinates"]):
        ax.set_title(title,loc="left",fontweight="bold");ax.tick_params(labelsize=8)
    axes[0].set(xlabel="Nodes per coordinate",ylabel="Normalized response C₄₃",xscale="log",xticks=[4,8,16,32,64],xticklabels=[4,8,16,32,64]);axes[0].legend(frameon=False,fontsize=8)
    axes[1].set(xlabel="Half-width / slow harmonic σ",ylabel="Normalized response C₄₃",xticks=[3,4,5,6])
    axes[2].set(xlabel="Nodes per coordinate",ylabel="C₄₃ with common d=4 normalization");axes[2].legend(frameon=False,fontsize=6.5,ncol=2)
    fig.suptitle("Accuracy has three separate questions",x=.06,ha="left",fontsize=13)
    fig.tight_layout();save(fig,"convergence")

if __name__=="__main__":
    import argparse
    parser=argparse.ArgumentParser();parser.add_argument("--input",type=Path,default=Path("results"));parser.add_argument("--output",type=Path,default=Path("figures"));args=parser.parse_args()
    make_figures(args.input,args.output)
