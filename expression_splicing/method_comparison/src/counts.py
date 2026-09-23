"""Stage counts (the quant.sf route; stage `tables` builds the same results from the supplementary data tables): gene counts,
Salmon TPM and gene effective lengths of the 29 samples from Salmon output, <SALMON_QUANT>/<source>/{quant.sf, aux_info/meta_info.json},
with the tximport summarizeToGene rule (the sample-arms: the reads of the UMI-deduplicated BAM; the TruSeq runs: the trimmed raw
reads, no UMIs): count = sum of the transcripts' NumReads, TPM = sum of their TPM, effective length =
TPM-weighted mean of their EffectiveLength (NaN where the gene has no reads; the unweighted mean kept as the fallback); tx2gene =
TX2GENE (Ensembl 113, version-less transcript ids). count / effective length equals gene TPM up to a per-sample constant, checked
per sample. Writes results/gene_counts.csv.gz, gene_tpm.csv.gz, gene_efflen.csv.gz, gene_efflen_unweighted.csv.gz (GENEID, SYMBOL,
29 samples; every gene of the Salmon index) and results/sample_table.tsv (sample, method, group, unit, source, Salmon library type
and fragment-length distribution, reads, genes detected)."""
import json, os, sys
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import RESULTS, ANALYSIS, UNIT, SALMON_QUANT, TX2GENE, write_tsv

# The 29 samples: BOBseq = its 3 untreated Hek controls only; every other method = all its control samples. prime-seq: all
# 8 untreated wells (one extraction condition, "Incubation + ProtK", 10k cells, per the E-MTAB-10142 SDRF). BRB-seq: all 8
# batch-5 no-Sendai control wells (2 "ctrl" + 6 empty-vector) of GSE334309. TruSeq: all 10 WT HEK293T runs of the three
# studies (group = laboratory / GEO series). `source` is the column name of the sample in the supplementary data tables
# (S7 / S8 / S9; the TruSeq runs are their truseq_public__<run> columns) and the quant.sf directory under SALMON_QUANT.
SAMPLES = [
    ("BOBseq",    "BOBseq_ctl_1",    "local", "bobseq_pe_native__Hek_control_1", "control"),
    ("BOBseq",    "BOBseq_ctl_2",    "local", "bobseq_pe_native__Hek_control_2", "control"),
    ("BOBseq",    "BOBseq_ctl_3",    "local", "bobseq_pe_native__Hek_control_3", "control"),
    ("prime-seq", "prime-seq_1",     "local", "primeseq_native__HEK2", "E-MTAB-10142"),
    ("prime-seq", "prime-seq_2",     "local", "primeseq_native__HEK3", "E-MTAB-10142"),
    ("prime-seq", "prime-seq_3",     "local", "primeseq_native__HEK9", "E-MTAB-10142"),
    ("prime-seq", "prime-seq_4",     "local", "primeseq_native__HEK20", "E-MTAB-10142"),
    ("prime-seq", "prime-seq_5",     "local", "primeseq_native__HEK31", "E-MTAB-10142"),
    ("prime-seq", "prime-seq_6",     "local", "primeseq_native__HEK42", "E-MTAB-10142"),
    ("prime-seq", "prime-seq_7",     "local", "primeseq_native__HEK53", "E-MTAB-10142"),
    ("prime-seq", "prime-seq_8",     "local", "primeseq_native__HEK64", "E-MTAB-10142"),
    ("BRB-seq",   "BRB-seq_1",       "local", "brbseq_native__ctl_TACGTTATTCCGAA", "ctrl well"),       # well A01, ctrl
    ("BRB-seq",   "BRB-seq_2",       "local", "brbseq_native__ctl_ACGAGCAGATGCAG", "ctrl well"),       # well A04, ctrl
    ("BRB-seq",   "BRB-seq_3",       "local", "brbseq_native__ctl_CCGCAAGAGGTTGG", "vector well"),       # well B01, vector
    ("BRB-seq",   "BRB-seq_4",       "local", "brbseq_native__ctl_TAACACCTCCAACG", "vector well"),       # well B02, vector
    ("BRB-seq",   "BRB-seq_5",       "local", "brbseq_native__ctl_CAGACATGATAGAG", "vector well"),       # well B03, vector
    ("BRB-seq",   "BRB-seq_6",       "local", "brbseq_native__ctl_AACAACCGAAGTAA", "vector well"),       # well B04, vector
    ("BRB-seq",   "BRB-seq_7",       "local", "brbseq_native__ctl_ATGCAATTACACCG", "vector well"),       # well F09, vector
    ("BRB-seq",   "BRB-seq_8",       "local", "brbseq_native__ctl_CCTCCATCCATGGA", "vector well"),       # well F10, vector
    ("TruSeq",    "TruSeq_1",        "gcs",   "SRR33482324", "GSE296724"),                          # GSE296724 HEK293T_WT1
    ("TruSeq",    "TruSeq_2",        "gcs",   "SRR33482323", "GSE296724"),                          # GSE296724 HEK293T_WT2
    ("TruSeq",    "TruSeq_3",        "gcs",   "SRR33482322", "GSE296724"),                          # GSE296724 HEK293T_WT3
    ("TruSeq",    "TruSeq_4",        "gcs",   "SRR33776933", "GSE298560"),                          # GSE298560 WT 293T rep1
    ("TruSeq",    "TruSeq_5",        "gcs",   "SRR33776932", "GSE298560"),                          # GSE298560 WT 293T rep2
    ("TruSeq",    "TruSeq_6",        "gcs",   "SRR33776931", "GSE298560"),                          # GSE298560 WT 293T rep3
    ("TruSeq",    "TruSeq_7",        "gcs",   "SRR15425654", "GSE181978"),                          # GSE181978 HEK293T WT rep 1
    ("TruSeq",    "TruSeq_8",        "gcs",   "SRR15425655", "GSE181978"),                          # GSE181978 HEK293T WT rep 2
    ("TruSeq",    "TruSeq_9",        "gcs",   "SRR15425656", "GSE181978"),                          # GSE181978 HEK293T WT rep 3
    ("TruSeq",    "TruSeq_10",       "gcs",   "SRR15425657", "GSE181978"),                          # GSE181978 HEK293T WT rep 4
]

def quant_path(kind, key):
    """Salmon output directory of one sample: <SALMON_QUANT>/<key>/ with quant.sf and aux_info/meta_info.json."""
    d = os.path.join(SALMON_QUANT, key)
    if not os.path.exists(os.path.join(d, "quant.sf")):
        raise SystemExit(f"{d}/quant.sf not found (set BOBSEQ_SALMON_QUANT, or use stage `tables`)")
    return d

def tx2gene():
    t = pd.read_csv(TX2GENE, sep="\t")                                            # TXNAME GENEID SYMBOL DESCRIPTION
    return dict(zip(t.TXNAME, t.GENEID)), dict(zip(t.GENEID, t.SYMBOL.fillna("")))

def gene_tables(path, t2g):
    """tximport summarizeToGene equivalent (ignoreTxVersion): per gene the sum of NumReads, the sum of TPM, the TPM-weighted mean
    EffectiveLength (NaN without reads) and the unweighted mean EffectiveLength; plus reads without a gene and total reads."""
    q = pd.read_csv(path, sep="\t", usecols=["Name", "EffectiveLength", "TPM", "NumReads"])
    q["GENEID"] = q.Name.str.split(".").str[0].map(t2g)
    miss = q.GENEID.isna(); g = q.dropna(subset=["GENEID"]); by = g.groupby("GENEID")
    counts, tpm, ulen = by.NumReads.sum(), by.TPM.sum(), by.EffectiveLength.mean()
    wlen = ((g.TPM * g.EffectiveLength).groupby(g.GENEID).sum() / tpm).where(tpm > 0)
    # identity check: count / weighted length == TPM x (sum over all transcripts of NumReads / EffectiveLength) / 1e6
    k = float((q.NumReads / q.EffectiveLength).sum()) / 1e6; ok = counts >= 10
    dev = float(((counts[ok] / wlen[ok]) / (tpm[ok] * k) - 1).abs().max())
    return counts, tpm, wlen, ulen, float(q.loc[miss, "NumReads"].sum()), float(q.NumReads.sum()), dev

def salmon_meta(d):
    """Library type and fragment-length distribution Salmon used (aux_info/meta_info.json)."""
    j = json.load(open(os.path.join(d, "aux_info", "meta_info.json")))
    return {"salmon_version": j["salmon_version"], "library_type": "/".join(j["library_types"]), "frag_length_mean": round(j["frag_length_mean"], 1), "frag_length_sd": round(j["frag_length_sd"], 1),
            "salmon_mapped_pct": round(j["percent_mapped"], 2)}

def main():
    os.makedirs(RESULTS, exist_ok=True); t2g, sym = tx2gene(); tabs = {"counts": {}, "tpm": {}, "efflen": {}, "efflen_unweighted": {}}; rows = []
    for method, name, kind, key, group in SAMPLES:
        d = quant_path(kind, key); c, t, wl, ul, unassigned, total, dev = gene_tables(os.path.join(d, "quant.sf"), t2g)
        tabs["counts"][name], tabs["tpm"][name], tabs["efflen"][name], tabs["efflen_unweighted"][name] = c, t, wl, ul
        rows.append({"sample": name, "method": method, "group": group, "expression_unit": UNIT[method], "source": key,
                     "unit": "dedup reads" if kind == "local" else "reads (no UMI)", **salmon_meta(d), "salmon_reads": round(total), "reads_no_gene": round(unassigned),
                     "gene_counts": round(float(c.sum())), "genes_ge1": int((c >= 1).sum()), "genes_ge10": int((c >= 10).sum()),
                     "efflen_count_weighted_mean": round(float((c * wl).sum() / c.sum())), "tpm_identity_max_rel_dev": f"{dev:.1e}"})
        print(f"  {name:<15} {key:<34} {rows[-1]['library_type']:<4} frag {rows[-1]['frag_length_mean']:>6.1f} +- {rows[-1]['frag_length_sd']:>5.1f}  {total:>12,.0f} reads -> {c.sum():>12,.0f} gene counts, "
              f"{int((c>=10).sum()):>6,} genes >= 10; count/length vs TPM max rel. dev. {dev:.1e}")
    for k, series in tabs.items():
        M = pd.DataFrame(series); M = M.fillna(0.0) if k in ("counts", "tpm") else M                    # lengths keep NaN (= no reads in that sample)
        M.insert(0, "SYMBOL", [sym.get(g, "") for g in M.index]); M.index.name = "GENEID"; M.to_csv(os.path.join(RESULTS, f"gene_{k}.csv.gz"))
    write_tsv(os.path.join(RESULTS, "sample_table.tsv"), rows, list(rows[0].keys()))
    print(f"\ngene_counts / gene_tpm / gene_efflen / gene_efflen_unweighted .csv.gz: {M.shape[0]:,} genes x {len(SAMPLES)} samples")

if __name__ == "__main__":
    main()
