#!/usr/bin/env python3
"""Picard gene-body coverage plot (5'->3'): poly(dT) vs random priming, human (dark) & mouse (light),
mean + 95% CI band. Reads the per-sample Picard metric files under data_derived/genebody_metrics/.

Layout of data_derived/genebody_metrics/:
  random_priming/<sample>.<human|mouse>.rnaseq_metrics.txt   (included in this deposit; 9 samples)
  polydt/<sample>.<human|mouse>.rnaseq_metrics.txt           (add from the main species-mixing deposit
                                                              to reproduce the comparison)
  <group>_counts.tsv  (optional, columns run/sample + human_reads + mouse_reads) enforces the
                       paired >=MIN_READS/species inclusion used in the manuscript.
The metric files were produced by picard_genebody.sh (inserts <=500, RP removed, MINIMUM_LENGTH=1000).

Usage: python3 plot_genebody_profile.py [--min-reads 1000] [--out-dir ../figures_reference]
Requires numpy, matplotlib; scipy optional (t-based CI, else normal approx).
"""
import os, argparse, glob, numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
matplotlib.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','Helvetica','DejaVu Sans'],
                            'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none'})
try:
    from scipy import stats as _st
    def tcrit(n): return float(_st.t.ppf(0.975,n-1)) if n>1 else 0.0
except Exception:
    def tcrit(n): return 1.96 if n>1 else 0.0
HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.dirname(HERE)
MDIR=os.path.join(ROOT,"data_derived","genebody_metrics")
GROUPS=[("random priming","#AD1457","#F06292","random_priming"),   # (label, human colour, mouse colour, dir)
        ("poly(dT) 3'",   "#1B5E20","#66BB6A","polydt")]
def hist(path):
    if not os.path.exists(path): return None
    cov=[];inh=False
    for l in open(path).read().splitlines():
        if l.startswith('## HISTOGRAM'): inh=True; continue
        if inh:
            if l.startswith('normalized_position') or not l.strip(): continue
            q=l.split('\t')
            if q[0].isdigit(): cov.append(float(q[1]))
    return cov if len(cov)==101 else None
def load_counts(group_dir):
    p=os.path.join(MDIR,f"{group_dir}_counts.tsv"); c={}
    if os.path.exists(p):
        for i,l in enumerate(open(p)):
            if i==0: continue
            f=l.rstrip("\n").split("\t")
            # accept either run\tbc\th\tm  or  run\tbarcode\th\tm
            c[(f[0],f[1])]=(int(f[-2]),int(f[-1]))
    return c
def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--min-reads",type=int,default=1000)
    ap.add_argument("--out-dir",default=os.path.join(ROOT,"figures_reference")); a=ap.parse_args()
    X=list(range(101)); fig,ax=plt.subplots(figsize=(9.5,7.2)); MIN=a.min_reads
    for glabel,hcol,mcol,gdir in GROUPS:
        d=os.path.join(MDIR,gdir)
        if not os.path.isdir(d):
            print(f"[skip] {glabel}: {d} not present"); continue
        counts=load_counts(gdir)
        samples=sorted({os.path.basename(f).rsplit('.',3)[0] for f in glob.glob(os.path.join(d,"*.rnaseq_metrics.txt"))})
        keep=[]
        for s in samples:
            cH=hist(os.path.join(d,f"{s}.human.rnaseq_metrics.txt")); cM=hist(os.path.join(d,f"{s}.mouse.rnaseq_metrics.txt"))
            ok=cH is not None and cM is not None
            if ok and counts:
                key=next((k for k in counts if k[0] in s or k[1] in s or f"{k[0]}__{k[1]}"==s), None)
                if key and (counts[key][0]<MIN or counts[key][1]<MIN): ok=False
            if ok: keep.append((cH,cM))
        if not keep: continue
        for sidx,species,col in ((0,"human",hcol),(1,"mouse",mcol)):
            arr=np.array([k[sidx] for k in keep]); n=len(arr); mean=arr.mean(0)
            sem=arr.std(0,ddof=1)/np.sqrt(n) if n>1 else np.zeros_like(mean); half=tcrit(n)*sem
            ax.fill_between(X,mean-half,mean+half,color=col,alpha=0.16,lw=0,zorder=2)
            ax.plot(X,mean,color=col,lw=2.8,zorder=4,label=f"{glabel} — {species} (n={n})")
    ax.axhline(1,color='black',lw=1.0,zorder=3); ax.set_box_aspect(1)
    ax.set_xlabel("gene-body percentile (5'→3')"); ax.set_ylabel("normalized coverage")
    for s in ('top','right'): ax.spines[s].set_visible(False)
    ax.legend(fontsize=9,frameon=False,loc='upper left',bbox_to_anchor=(1.02,1.0))
    fig.suptitle("Gene-body coverage — poly(dT) vs random priming, human (dark) vs mouse (light)   (shaded = 95% CI)\n"
                 f"inserts ≤500 bp · ribosomal proteins removed · transcripts ≥1000 bp · ≥{MIN} reads/species",
                 fontsize=11,fontweight='bold',x=0.02,ha='left')
    os.makedirs(a.out_dir,exist_ok=True)
    for ext in ("png","pdf","svg"): fig.savefig(f"{a.out_dir}/genebody_profile_polydt_vs_random.{ext}",dpi=175,bbox_inches='tight')
    print("wrote genebody_profile_polydt_vs_random.png/.pdf/.svg")
if __name__=="__main__": main()
