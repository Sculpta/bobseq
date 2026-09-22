#!/usr/bin/env python3
"""Genome-browser view of one gene for the two example samples, with the gene model expanded to
all annotated isoforms and per-base coverage of inserts <=500 bp.

The two samples (as used in the manuscript RACK1 panel):
  * 260730 bc12  — poly(dT)-primed 3'TSO bobcode  (run 260730_RT_conditions_BC11-13, barcode12)   [comparison]
  * 260827 BC15  — RANDOM-priming 3'TSO bobcode   (run 260827_BC14_16, barcode15; 1x poly(A))       [this study]
Both reads its coordinate-sorted, region-sliced BAM in data_derived/browser_bams/.

Usage:
  python3 genome_browser_gene.py --gene RACK1 --gtf /path/to/combined_genome.gtf
  # any single-copy human gene present on a HUMAN_-prefixed contig works: RACK1, ACTB, GAPDH, ...
Requires: samtools on PATH; Python 3 with numpy + matplotlib.
"""
import os, sys, argparse, subprocess, re, numpy as np, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Patch
matplotlib.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','Helvetica','DejaVu Sans'],
                            'pdf.fonttype':42,'svg.fonttype':'none'})
HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.dirname(HERE)
BAMDIR=os.path.join(ROOT,"data_derived","browser_bams")
# (label, colour, sliced BAM) — poly(dT) comparison first, random-priming study sample second
TRACKS=[("260730 bc12  ·  poly(dT)","#1B5E20", os.path.join(BAMDIR,"260730_bc12_polydT.rack1_actb_gapdh.bam")),
        ("260827 BC15  ·  random priming (1x polyA)","#AD1457", os.path.join(BAMDIR,"260827_bc15_randomPriming.rack1_actb_gapdh.bam"))]
MAXINS=500
INS_AWK=('BEGIN{OFS="\\t"} /^@/{print;next}{c=$6;ql=0;n="";for(i=1;i<=length(c);i++){ch=substr(c,i,1);'
         'if(ch>="0"&&ch<="9")n=n ch;else{if(ch=="M"||ch=="I"||ch=="="||ch=="X")ql+=n;n=""}} if(ql<=%d)print}' % MAXINS)
BT_ORDER={'protein_coding':0,'protein_coding_CDS_not_defined':1,'retained_intron':2,'nonsense_mediated_decay':3}
BT_COL={'protein_coding':'#1a1a1a','protein_coding_CDS_not_defined':'#666','retained_intron':'#9a9a9a','nonsense_mediated_decay':'#C1666B'}
tid_re=re.compile(r'transcript_id "([^"]+)"'); bt_re=re.compile(r'transcript_biotype "([^"]+)"')

def gene_lines(gene,gtf):
    out=subprocess.run(["grep",f'gene_name "{gene}"',gtf],capture_output=True,text=True).stdout.splitlines()
    return [l for l in out if l.startswith("HUMAN_")]

def depth(bam,gs,ge,region):
    v=subprocess.Popen(["samtools","view","-h","-F","0x904",bam,region],stdout=subprocess.PIPE,text=True)
    aw=subprocess.Popen(["awk",INS_AWK],stdin=v.stdout,stdout=subprocess.PIPE,text=True); v.stdout.close()
    b=subprocess.run(["samtools","view","-"],stdin=aw.stdout,capture_output=True,text=True); aw.stdout.close(); aw.wait()
    d=np.zeros(ge-gs+1)
    for line in b.stdout.splitlines():           # per-base depth of the filtered reads
        f=line.split('\t'); pos=int(f[3]); cig=f[5]
        for n,op in re.findall(r'(\d+)([MIDNSHP=X])',cig):
            n=int(n)
            if op in 'M=X':
                a=max(pos-gs,0); z=min(pos-gs+n,ge-gs+1)
                if z>a: d[a:z]+=1
                pos+=n
            elif op in 'DN': pos+=n
    return d

def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--gene",required=True); ap.add_argument("--gtf",required=True)
    ap.add_argument("--out-dir",default=os.path.join(ROOT,"figures_reference")); a=ap.parse_args()
    L=gene_lines(a.gene,a.gtf); gl=[l for l in L if l.split('\t')[2]=='gene']
    if not gl: sys.exit(f"{a.gene}: not found on a HUMAN_ contig in {a.gtf}")
    p=gl[0].split('\t'); chrom,gs,ge,strand=p[0],int(p[3]),int(p[4]),p[6]
    pad=int((ge-gs)*0.03)+100; gs-=pad; ge+=pad; region=f"{chrom}:{gs}-{ge}"
    tx={}
    for l in L:
        p=l.split('\t')
        if p[0]!=chrom: continue
        if p[2]=='transcript':
            t=tid_re.search(p[8]).group(1); b=bt_re.search(p[8]); tx.setdefault(t,{})['biotype']=b.group(1) if b else '?'; tx[t].setdefault('exons',[])
        elif p[2]=='exon':
            t=tid_re.search(p[8]).group(1); tx.setdefault(t,{}).setdefault('exons',[]).append((int(p[3]),int(p[4])))
    for t in tx: tx[t]['exons'].sort(); tx[t].setdefault('biotype','?')
    order=sorted(tx,key=lambda t:(BT_ORDER.get(tx[t]['biotype'],9),-(tx[t]['exons'][-1][1]-tx[t]['exons'][0][0])))
    X=np.arange(gs,ge+1); covs=[(lab,col,depth(bam,gs,ge,region)) for lab,col,bam in TRACKS]; nt=len(order)
    nz=np.zeros(ge-gs+1,bool)
    for _,_,d in covs: nz|=d>0
    if nz.any():
        lo=int(np.argmax(nz)); hi=len(nz)-1-int(np.argmax(nz[::-1])); w=hi-lo; pd=max(int(w*0.06),250)
        wlo=gs+max(0,lo-pd); whi=gs+min(len(nz)-1,hi+pd)
    else: wlo,whi=gs,ge
    zoomed=(whi-wlo)<0.9*(ge-gs)
    fig=plt.figure(figsize=(9.5,3.4+0.16*nt)); gspec=GridSpec(3,1,height_ratios=[1.3,1.3,0.16*nt+0.4],hspace=0.22)
    for i,(lab,col,d) in enumerate(covs):
        ax=fig.add_subplot(gspec[i]); ax.fill_between(X,d,color=col,lw=0,alpha=0.9); mx=d.max()
        ax.set_ylim(0,max(mx*1.15,1)); ax.set_yticks([int(mx)] if mx>0 else [0]); ax.set_xlim(wlo,whi)
        ax.set_ylabel(lab,rotation=0,ha='right',va='center',fontsize=8.5); ax.tick_params(axis='y',labelsize=7,length=2)
        ax.tick_params(axis='x',length=0,labelbottom=False)
        for s in ('top','right','bottom'): ax.spines[s].set_visible(False)
    gm=fig.add_subplot(gspec[2]); gm.set_xlim(wlo,whi); gm.set_ylim(-1,nt)
    for yi,t in enumerate(order):
        y=nt-1-yi; ex=tx[t]['exons']; col=BT_COL.get(tx[t]['biotype'],'#bbb')
        gm.plot([ex[0][0],ex[-1][1]],[y,y],color=col,lw=0.5,zorder=1)
        for s,e in ex: gm.add_patch(plt.Rectangle((s,y-0.34),max(e-s,(whi-wlo)*0.0015),0.68,color=col,lw=0,zorder=2))
        gm.text(wlo-(whi-wlo)*0.012,y,t.replace("ENST00000",""),ha='right',va='center',fontsize=5.0,color=col)
    gm.set_yticks([])
    for s in ('top','right','left'): gm.spines[s].set_visible(False)
    five="minus strand · 3′ end at left" if strand=='-' else "plus strand · 3′ end at right"
    zn="   ·   zoomed to expressed body" if zoomed else ""
    gm.set_xlabel(f"{chrom} position (bp)   ·   {a.gene} {five}   ·   {nt} isoforms{zn}",fontsize=9)
    gm.ticklabel_format(axis='x',style='plain',useOffset=False); gm.tick_params(axis='x',labelsize=8)
    leg=[Patch(color=BT_COL[b],label=b.replace('_',' ')) for b in ('protein_coding','protein_coding_CDS_not_defined','retained_intron','nonsense_mediated_decay')]
    gm.legend(handles=leg,fontsize=6.5,frameon=False,loc='lower left',bbox_to_anchor=(0,1.0),ncol=4,handlelength=1.0,columnspacing=1.2)
    fig.suptitle(f"{a.gene} — two samples over all isoforms  ·  coverage inserts ≤500 bp",fontsize=11,fontweight='bold',x=0.02,ha='left')
    fig.subplots_adjust(left=0.20,right=0.985,top=0.955,bottom=0.05)
    os.makedirs(a.out_dir,exist_ok=True)
    for ext in ("png","pdf","svg"): fig.savefig(f"{a.out_dir}/{a.gene.lower()}_two_samples_isoforms.{ext}",dpi=175,bbox_inches='tight')
    print(f"wrote {a.gene.lower()}_two_samples_isoforms.png/.pdf/.svg  ({nt} isoforms; "+", ".join(f'{lab.split(chr(183))[0].strip()} max={int(d.max())}' for lab,col,d in covs)+")")

if __name__=="__main__": main()
