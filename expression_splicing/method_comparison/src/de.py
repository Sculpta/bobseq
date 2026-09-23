"""Stage de: differential gene expression between the methods with DESeq2 (src/deseq2_contrast.R; needs R with DESeq2 and apeglm)
on the thinned counts of stage `mrna` or `tables`. Six pairwise contrasts (bobseq_vs_truseq, primeseq_vs_truseq, brbseq_vs_truseq,
bobseq_vs_primeseq, bobseq_vs_brbseq, primeseq_vs_brbseq; log2FC > 0 = higher in the first-named method), design ~ condition,
Wald test, apeglm-shrunk log2FC, genes with >= 10 counts over the samples, DE = padj < 0.001, no fold-change cutoff.
Length normalisation as DESeq2 offsets: for a contrast that involves BOBseq or TruSeq the R script receives one length per gene
per METHOD (the mean effective length over that method's samples; 1 for prime-seq / BRB-seq), stored as the avgTxLength assay,
from which DESeq2 builds per-gene normalisation factors as it does for tximport input, so a BOBseq-vs-prime-seq log2FC compares
TPM with CPM. A per-sample length would cancel within a method only up to Salmon's per-sample isoform noise, which doubles the
dispersion. Run twice: on the full gene set (results/de, figures/de) and on the CAPPED common set (results/de_capped,
figures/de_capped): the CAP genes with the highest mean log2(CPM + 1) in their least-expressing method
(results/de_capped_gene_set.tsv), DESeq2 re-fitted on the subset, so that no method shows more DE merely because it detects more
genes. Per contrast: de_<contrast>.tsv (baseMean, shrunk and raw log2FC, lfcSE, stat, p, padj, log2 normalised counts per sample),
a volcano (plain and with the top 5 up + 5 down named), a bar chart and a heatmap of those genes; per run de_summary.tsv,
de_top10.tsv, de_counts.svg and, when mRNA lengths are available (common.mrna_length), de_length_dependence.tsv."""
import os, re, subprocess, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import HERE, RESULTS, FIGURES, MRNA_THIN, EFFLEN, COL, TPM_METHODS, LENGTH_EDGES, LENGTH_LABELS, mrna_length, style
from scipy.stats import spearmanr
plt = style()
from scipy.cluster.hierarchy import linkage, leaves_list
SC = 0.625                              # publication size: 0.625 x the working size, fonts unchanged
MAX_PADJ, MIN_LFC = 0.001, 0.0          # DE = padj < 0.001, no log2FC cutoff
RD, FD = os.path.join(RESULTS, "de"), os.path.join(FIGURES, "de")
CONTRASTS = {"bobseq_vs_truseq": ("BOBseq", "TruSeq"), "primeseq_vs_truseq": ("prime-seq", "TruSeq"), "brbseq_vs_truseq": ("BRB-seq", "TruSeq"),
             "bobseq_vs_primeseq": ("BOBseq", "prime-seq"), "bobseq_vs_brbseq": ("BOBseq", "BRB-seq"), "primeseq_vs_brbseq": ("prime-seq", "BRB-seq")}
RNAME = {"BOBseq": "BOBseq", "prime-seq": "primeseq", "BRB-seq": "BRBseq", "TruSeq": "TruSeq"}       # syntactically valid R factor levels

def dodge(desired, gap, lo, hi):
    """1-D label placement: keep labels near `desired` (sorted descending) with >= gap between
    neighbours inside [lo, hi]; top-down pass, then a bottom-up push if the chain ran past lo."""
    y = []
    for d in desired: y.append(min(d, hi) if not y else min(d, y[-1] - gap))
    if y and y[-1] < lo:
        y[-1] = lo
        for i in range(len(y) - 2, -1, -1): y[i] = max(y[i], y[i + 1] + gap)
    return y

def deseq(tag, test, ref, M, Lm, method, out_dir):
    """Run one pairwise DESeq2 contrast (reference first in coldata, dose 0/1 orders the levels); with a length matrix (the method's
    mean effective length per gene for the TPM methods' samples, 1 for the others) when either method is length-normalised.
    Returns the coldata, the size factors (or the normalisation-factor geometric means) and whether lengths were used."""
    cols = [c for c in M.columns if method[c] == ref] + [c for c in M.columns if method[c] == test]
    cts = M[cols].round().astype(int); cts.index.name = "gene_id"
    cd = pd.DataFrame({"well": cols, "condition": [RNAME[method[c]] for c in cols], "dose": [int(method[c] == test) for c in cols]})
    cp, mp, lp, out = (os.path.join(out_dir, f"_{tag}_{x}.tsv") for x in ("counts", "coldata", "lengths", "")); out = os.path.join(out_dir, f"de_{tag}.tsv")
    cts.to_csv(cp, sep="\t"); cd.to_csv(mp, sep="\t", index=False); cmd = ["Rscript", os.path.join(HERE, "deseq2_contrast.R"), cp, mp, "pairwise", out]
    with_len = test in TPM_METHODS or ref in TPM_METHODS
    if with_len:
        lens = pd.DataFrame({c: (Lm[method[c]] if method[c] in TPM_METHODS else 1.0) for c in cols}, index=M.index); lens.index.name = "gene_id"; lens.to_csv(lp, sep="\t"); cmd.append(lp)
    r = subprocess.run(cmd, capture_output=True, text=True)
    for f in (cp, mp, lp):
        if os.path.exists(f): os.remove(f)
    if r.returncode: raise SystemExit(r.stderr[-3000:])
    print(f"{tag} {r.stdout.strip()}"); sf = re.search(r"size factors (.*)$", r.stdout.strip()).group(1)
    return cd, dict(zip(cols, sf.split())), with_len

def run(M, sym, Lm, method, glen, RD, FD, label):
    """The six contrasts, tables and figures on one count matrix (rows = the gene set of the run) into RD / FD."""
    summary, top_all, ldep = [], [], []
    for tag, (test, ref) in CONTRASTS.items():
        rd, fd = os.path.join(RD, tag), os.path.join(FD, tag); os.makedirs(rd, exist_ok=True); os.makedirs(fd, exist_ok=True)
        cd, sf, with_len = deseq(tag, test, ref, M, Lm, method, rd); order = [ref, test]; title = f"{test} vs {ref}"
        offsets = ("; ".join(m for m in (test, ref) if m in TPM_METHODS) + " (method mean per gene)") if with_len else "none"
        d = pd.read_csv(os.path.join(rd, f"de_{tag}.tsv"), sep="\t"); d.insert(1, "symbol", d.gene_id.map(sym))
        d["name"] = d.symbol.str.replace(r" \[ENSG\d+\]$", "", regex=True).where(d.symbol != "", d.gene_id)     # plot label: tx2gene's "SYMBOL [ENSG..]" disambiguation dropped
        d["hit"] = (d.padj < MAX_PADJ).astype(int); d = d.sort_values(["hit", "padj"], ascending=[False, True])
        d.to_csv(os.path.join(rd, f"de_{tag}.tsv"), sep="\t", index=False)
        t = d[d.padj.notna()].copy(); t["len"] = glen.reindex(t.gene_id).to_numpy(); t = t[t.len.notna()]
        if len(t):                                                     # the length columns need results/mrna_length.tsv (common.mrna_length)
            b = pd.cut(10 ** t.len, LENGTH_EDGES, labels=LENGTH_LABELS, right=False)
            ldep.append({"contrast": tag, "length_offsets": offsets, "genes": len(t), "spearman_rho_log2fc_vs_log10_length": round(float(spearmanr(t.log2fc, t.len).correlation), 3),
                         **{f"median_log2fc_{lab}": round(float(t.log2fc[(b == lab).to_numpy()].median()), 3) for lab in LENGTH_LABELS}})
        up, dn = d[(d.hit == 1) & (d.log2fc > 0)].sort_values("padj"), d[(d.hit == 1) & (d.log2fc < 0)].sort_values("padj")
        top = pd.concat([up.head(5), dn.head(5)])
        summary.append({"contrast": tag, "title": title, "test": test, "reference": ref, "model": "pairwise", "length_offsets": offsets, "samples": len(cd), "genes_tested": int(d.padj.notna().sum()),
                        "padj_lt_0.05": int((d.padj < 0.05).sum()), "de_genes": int(d.hit.sum()), "up": len(up), "down": len(dn), "size_factors": " ".join(f"{w}={s}" for w, s in sf.items())})
        # ---- volcano, two versions: plain, and labelled with the top 5 up + 5 down on a rail in the right
        # margin, dodged vertically and joined to their point by a leader line. padj underflowing to 0 is drawn at 1e-300.
        x, y, h = d.log2fc.to_numpy(), -np.log10(d.padj.fillna(1).clip(lower=1e-300).to_numpy()), d.hit.to_numpy() == 1
        XL = max(3.0, np.abs(x[np.isfinite(x)]).max() * 1.05); YL = max(5.0, y.max() * 1.1)
        for labelled in (False, True):
            fig, ax = plt.subplots(figsize=((5.6 if not labelled else 7.0) * SC, 5.4 * SC))
            ax.scatter(x[~h], y[~h], s=4, color="#bdbdbd", lw=0, alpha=0.7, rasterized=True); ax.scatter(x[h & (x > 0)], y[h & (x > 0)], s=9, color="#d62728", lw=0, alpha=0.7, zorder=3, rasterized=True, label=f"up ({len(up)})")
            ax.scatter(x[h & (x < 0)], y[h & (x < 0)], s=9, color="#1f77b4", lw=0, alpha=0.7, zorder=3, rasterized=True, label=f"down ({len(dn)})")
            ax.axhline(-np.log10(MAX_PADJ), ls=":", color="#999", lw=0.8)
            ax.set_xlim(-XL, XL); ax.set_ylim(0, YL); ax.set_box_aspect(1)
            ax.set_xlabel(f"log2FC {test} vs {ref}"); ax.set_ylabel("-log10 padj")
            ax.set_title(f"{title}\n{len(up) + len(dn)} DE genes (padj<{MAX_PADJ})\nlength offsets: {', '.join(m for m in (test, ref) if m in TPM_METHODS) or 'none'}", fontsize=8, loc="left")
            fig.tight_layout(rect=[0, 0, 0.82 if labelled else 1, 1])        # lay out BEFORE the rail labels, which would otherwise push the y-label off-canvas
            if labelled and len(top):
                pts = sorted([(r.log2fc, -np.log10(max(r.padj, 1e-300)), r["name"], r.log2fc > 0) for _, r in top.iterrows()], key=lambda t: -t[1])
                ax.scatter([p[0] for p in pts], [p[1] for p in pts], s=22, facecolor="none", edgecolor="#111", lw=0.7, zorder=4)
                ys = dodge([p[1] / YL for p in pts], gap=0.075, lo=0.03, hi=0.97)
                for (x0, y0, g, isup), yl in zip(pts, ys):
                    ax.annotate(g, xy=(x0, y0), xycoords="data", xytext=(1.04, yl), textcoords="axes fraction", fontsize=7, ha="left", va="center",
                                color="#d62728" if isup else "#1f77b4", annotation_clip=False, arrowprops=dict(arrowstyle="-", color="#888", lw=0.5, shrinkA=0, shrinkB=2))
                ax.legend(frameon=False, fontsize=7, loc="upper left", bbox_to_anchor=(1.04, 0.0), handletextpad=0.2, borderpad=0.2)
            else:
                ax.legend(frameon=False, fontsize=7, loc="upper left", bbox_to_anchor=(1.01, 1.0), handletextpad=0.2, borderpad=0.2)
            name = f"volcano_{tag}" + ("_labelled" if labelled else "")
            fig.savefig(os.path.join(fd, f"{name}.svg")); plt.close(fig)
        if len(top) == 0: continue
        wells_ = [c[9:] for c in d.columns if c.startswith("log2norm_")]; cond_of = {w: method[w] for w in wells_}; worder = [w for c in order for w in wells_ if cond_of[w] == c]
        for _, r in top.iterrows(): top_all.append({"contrast": tag, "gene_id": r.gene_id, "symbol": r.symbol, "direction": "up" if r.log2fc > 0 else "down", "log2fc": r.log2fc, "padj": r.padj, **{f"log2norm_{w}": r[f"log2norm_{w}"] for w in worder}})
        # ---- bar charts: mean +- SD, dots = samples
        fig, axes = plt.subplots(2, 5, figsize=(11.5 * SC, 7.5 * SC))
        for ax, (_, r) in zip(axes.flat, top.iterrows()):
            for i, c in enumerate(order):
                v = np.array([r[f"log2norm_{w}"] for w in wells_ if cond_of[w] == c])
                ax.bar(i, v.mean(), yerr=v.std(ddof=1) if len(v) > 1 else 0, color=COL[c], alpha=0.6, width=0.7, capsize=3, error_kw={"lw": 1})
                ax.scatter([i + (j - (len(v) - 1) / 2) * 0.08 for j in range(len(v))], v, s=10, color="#111", zorder=3)
            ax.set_xticks(range(len(order))); ax.set_xticklabels(order, fontsize=6.5, rotation=45, ha="right"); ax.set_ylabel("log2(norm+1)", fontsize=6.5); ax.tick_params(axis="y", labelsize=6.5)
            ax.set_title(f"{r['name']}\n{r.log2fc:+.1f}, q {'<1e-300' if r.padj == 0 else f'{r.padj:.0e}'}", fontsize=7, loc="left")
        for ax in list(axes.flat)[len(top):]: ax.axis("off")
        fig.suptitle(f"{title}: top 5 up / 5 down; bar = mean, error = SD, dots = samples", fontsize=8)
        fig.tight_layout(); fig.savefig(os.path.join(fd, f"bars_top10_{tag}.svg")); plt.close(fig)
        # ---- heatmap of the same genes x samples
        Mz = top[[f"log2norm_{w}" for w in worder]].to_numpy(); Z = (Mz - Mz.mean(1, keepdims=True)) / (Mz.std(1, keepdims=True) + 1e-9)
        ordr = leaves_list(linkage(Z, "average", metric="correlation")) if len(top) > 2 else np.arange(len(top))
        fig, ax = plt.subplots(figsize=(5.5 * SC, 4.4 * SC)); im = ax.imshow(Z[ordr], aspect="auto", cmap="RdBu_r", vmin=-2, vmax=2)
        ax.set_yticks(range(len(top))); ax.set_yticklabels(top["name"].to_numpy()[ordr], fontsize=6.5); ax.set_xticks(range(len(worder))); ax.set_xticklabels(worder, rotation=90, fontsize=6)
        for k, w in enumerate(worder): ax.get_xticklabels()[k].set_color(COL[cond_of[w]])
        for b in np.cumsum([sum(cond_of[w] == c for w in worder) for c in order])[:-1]: ax.axvline(b - 0.5, color="#111", lw=1)
        ax.set_title(f"{title}: top DE genes, z-score", fontsize=7.5, loc="left"); cb = fig.colorbar(im, ax=ax, shrink=0.6); cb.set_label("z-score", fontsize=7); cb.ax.tick_params(labelsize=6.5)
        fig.tight_layout(); fig.savefig(os.path.join(fd, f"heatmap_top10_{tag}.svg")); plt.close(fig)
    S = pd.DataFrame(summary); S.to_csv(os.path.join(RD, "de_summary.tsv"), sep="\t", index=False); pd.DataFrame(top_all).to_csv(os.path.join(RD, "de_top10.tsv"), sep="\t", index=False)
    if ldep: pd.DataFrame(ldep).to_csv(os.path.join(RD, "de_length_dependence.tsv"), sep="\t", index=False)
    # ---- DE gene counts per contrast (up / down)
    fig, ax = plt.subplots(figsize=(7.0 * SC, 4.6 * SC)); xs = np.arange(len(S)); pad = 0.015 * max(S["up"].max(), S["down"].max())
    ax.bar(xs, S["up"], color="#d62728", alpha=0.85, label="up"); ax.bar(xs, -S["down"], color="#1f77b4", alpha=0.85, label="down")
    for i, r in S.iterrows(): ax.text(i, r.up + pad, str(r.up), ha="center", fontsize=6.5); ax.text(i, -r.down - pad, str(r.down), ha="center", va="top", fontsize=6.5)
    ax.axhline(0, color="#333", lw=0.8); ax.set_xticks(xs); ax.set_xticklabels(S["title"], rotation=35, ha="right", fontsize=6.5); ax.set_ylim(-1.18 * S["down"].max(), 1.15 * S["up"].max())
    ax.set_ylabel(f"DE genes (padj < {MAX_PADJ})"); ax.legend(frameon=False, fontsize=7, loc="upper left", bbox_to_anchor=(1.01, 1.0)); ax.set_title(f"DE genes per method contrast (DESeq2; {' / '.join(TPM_METHODS)} with gene-length offsets){label}", fontsize=8.5, loc="left")
    fig.tight_layout(); fig.savefig(os.path.join(FD, "de_counts.svg")); plt.close(fig)
    print(S.drop(columns=["title", "size_factors"]).to_string(index=False))
    if ldep: print("\nlog2FC vs mRNA length (tested genes):"); print(pd.DataFrame(ldep).to_string(index=False))

CAP = 7500                              # the capped run: the CAP genes with the highest expression in their least-expressing method

def capped_gene_set(M, method):
    """Genes 'expressed in all studies' for the capped DE run: rank every gene by the MINIMUM over the four methods of its mean
    log2(CPM + 1) on the thinned counts (CPM over the mRNA gene set) and keep the top CAP -- a gene qualifies only through its
    least-expressing method, so the set cannot be carried by the method that detects most genes: a common gene set of
    roughly 7,500 expressed genes so that a method does not show more DE just because it detects more genes."""
    cpm = M / M.sum() * 1e6; score = pd.DataFrame({m: np.log2(cpm[[c for c in M.columns if method[c] == m]] + 1).mean(1) for m in ORDER}).min(1).sort_values(ascending=False)
    keep = score.index[:CAP]; cut = float(score.iloc[CAP - 1])
    return keep, cut, score

ORDER = ["BOBseq", "prime-seq", "BRB-seq", "TruSeq"]

def main():
    M = pd.read_csv(MRNA_THIN, index_col=0); sym = M["SYMBOL"].fillna(""); M = M.drop(columns=["SYMBOL"]); L = pd.read_csv(EFFLEN, index_col=0).drop(columns=["SYMBOL"]).loc[M.index]
    samples = pd.read_csv(os.path.join(RESULTS, "sample_table.tsv"), sep="\t"); method = dict(zip(samples["sample"], samples["method"]))
    Lm = {m: L[[c for c in M.columns if method[c] == m]].mean(1) for m in TPM_METHODS}        # one length per gene per method: the mean over its samples
    ml = mrna_length(); glen = np.log10(ml) if len(ml) else ml
    print(f"== full gene set: {len(M):,} genes"); run(M, sym, Lm, method, glen, RD, FD, "")
    keep, cut, score = capped_gene_set(M, method); Mc = M.loc[keep]
    pd.DataFrame({"gene_id": keep, "symbol": sym.loc[keep].to_numpy(), "min_over_methods_mean_log2cpm": score.loc[keep].round(4).to_numpy(), "rank": np.arange(1, len(keep) + 1)}).to_csv(os.path.join(RESULTS, "de_capped_gene_set.tsv"), sep="\t", index=False)
    print(f"\n== capped gene set: top {CAP:,} genes by the minimum over methods of the mean log2(CPM+1) (cutoff {cut:.2f} = CPM {2 ** cut - 1:.1f}); "
          f"{int((Mc >= 1).all(1).sum()):,} of them have >= 1 thinned count in every sample; they carry " + ", ".join(f"{m} {100 * Mc[[c for c in M.columns if method[c] == m]].sum().sum() / M[[c for c in M.columns if method[c] == m]].sum().sum():.1f} %" for m in ORDER) + " of the thinned counts")
    run(Mc, sym, Lm, method, glen, RD + "_capped", FD + "_capped", f"; capped to the {CAP:,}-gene common set")

if __name__ == "__main__":
    main()
