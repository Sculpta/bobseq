#!/usr/bin/env python3
"""Per-sample 5'->3' coverage balance boxplot: random priming (magenta) vs poly(dT) (green),
human + mouse, single panel. Balance = 2 x coverage centroid of the gene-body metagene:
1 = evenly balanced, <1 = 5'-leaning, >1 = 3'-leaning (bounded [0,2]; robust for 3'-piled poly(dT)).

Reads the same Picard metric files as plot_genebody_profile.py
(data_derived/genebody_metrics/{random_priming,polydt}/<sample>.<species>.rnaseq_metrics.txt;
inserts <=500, RP removed, MINIMUM_LENGTH=1000). Optional <group>_counts.tsv enforces the paired
>=MIN reads/species inclusion. Requires numpy + matplotlib.
"""
import os, argparse, glob, numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
matplotlib.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','Helvetica','DejaVu Sans'],
                            'pdf.fonttype':42,'svg.fonttype':'none'})
HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.dirname(HERE)
MDIR=os.path.join(ROOT,"data_derived","genebody_metrics")
# (label, dir, human colour, mouse colour) — random priming = magenta, poly(dT) = green
GROUPS=[("random priming","random_priming","#AD1457","#F06292"),
        ("poly(dT) 3'","polydt","#1B5E20","#66BB6A")]
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
def balance(cov):
    c=np.array(cov); pos=np.linspace(0,1,101); return 2.0*(pos*c).sum()/c.sum()
def load_counts(gdir):
    p=os.path.join(MDIR,f"{gdir}_counts.tsv"); c={}
    if os.path.exists(p):
        for i,l in enumerate(open(p)):
            if i==0: continue
            f=l.rstrip("\n").split("\t"); c[(f[0],f[1])]=(int(f[-2]),int(f[-1]))
    return c
def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--min-reads",type=int,default=1000)
    ap.add_argument("--out-dir",default=os.path.join(ROOT,"figures_reference")); a=ap.parse_args(); MIN=a.min_reads
    series=[]  # (label, colour, values)  in order: RP-human, RP-mouse, polydT-human, polydT-mouse
    for glabel,gdir,hcol,mcol in GROUPS:
        d=os.path.join(MDIR,gdir)
        if not os.path.isdir(d): print(f"[skip] {glabel}: {d} not present"); continue
        counts=load_counts(gdir)
        samples=sorted({os.path.basename(f).rsplit('.',3)[0] for f in glob.glob(os.path.join(d,"*.rnaseq_metrics.txt"))})
        bH=[];bM=[]
        for s in samples:
            cH=hist(os.path.join(d,f"{s}.human.rnaseq_metrics.txt")); cM=hist(os.path.join(d,f"{s}.mouse.rnaseq_metrics.txt"))
            if cH is None or cM is None: continue
            if counts:
                key=next((k for k in counts if k[0] in s or k[1] in s or f"{k[0]}__{k[1]}"==s),None)
                if key and (counts[key][0]<MIN or counts[key][1]<MIN): continue
            bH.append(balance(cH)); bM.append(balance(cM))
        series.append((f"{glabel}",hcol,mcol,np.array(bH),np.array(bM)))
    # flatten to 4 boxes: [RP-h, RP-m, polydT-h, polydT-m]
    boxes=[]
    for glabel,hcol,mcol,bH,bM in series:
        boxes.append((glabel,"human",hcol,bH)); boxes.append((glabel,"mouse",mcol,bM))
    positions=[1,1.5,2.2,2.7][:len(boxes)]
    fig,ax=plt.subplots(figsize=(3.8,6.8)); rng=np.random.default_rng(3)
    for (glab,sp,col,vals),pos in zip(boxes,positions):
        if len(vals)==0: continue
        ax.boxplot([vals],positions=[pos],widths=0.3875,patch_artist=True,showfliers=False,
                   medianprops=dict(color=col,lw=2.2),boxprops=dict(facecolor=col,alpha=0.28,edgecolor='none'),
                   whiskerprops=dict(color=col,lw=1.2),capprops=dict(color=col,lw=1.2))
        x=rng.normal(pos,0.03,size=len(vals)); ax.scatter(x,vals,s=23,color=col,alpha=0.7,edgecolor='white',linewidth=0.25,zorder=3)
    ax.axhline(1,color='black',lw=1.2,zorder=2); ax.set_ylim(0,2); ax.set_yticks([0,0.5,1.0,1.5,2.0])
    ax.set_xlim(0.55,3.85); ax.set_xticks(positions); ax.set_xticklabels([b[1] for b in boxes],fontsize=9.5)
    if len(positions)>=2:
        tr=ax.get_xaxis_transform()
        ax.text((positions[0]+positions[1])/2,-0.085,"random priming",transform=tr,ha='center',va='top',fontsize=10,fontweight='bold')
        if len(positions)>=4:
            ax.text((positions[2]+positions[3])/2,-0.085,"poly(dT) 3'",transform=tr,ha='center',va='top',fontsize=10,fontweight='bold')
    ax.set_ylabel("5′→3′ balance  (2×coverage centroid;  1 = even, >1 = 3′-leaning)",fontsize=10)
    for s in ('top','right'): ax.spines[s].set_visible(False)
    fig.suptitle("Per-sample 5′→3′ coverage balance — random priming vs poly(dT)\n"
                 f"inserts ≤500 bp · RP removed · transcripts ≥1000 bp · ≥{MIN} reads/species",
                 fontsize=11,fontweight='bold',x=0.02,ha='left')
    os.makedirs(a.out_dir,exist_ok=True)
    for ext in ("png","pdf","svg"): fig.savefig(f"{a.out_dir}/balance_box_random_vs_polydt.{ext}",dpi=175,bbox_inches='tight')
    for glab,sp,col,vals in boxes:
        if len(vals): print(f"  {glab} {sp}: n={len(vals)} median={np.median(vals):.3f}")
    print("wrote balance_box_random_vs_polydt.png/.pdf/.svg")
if __name__=="__main__": main()
