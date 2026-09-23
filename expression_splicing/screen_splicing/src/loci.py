"""Stage loci: genome coordinates for every top-10 junction (the input of the splicing_tracks figures).

Per row of results/ds_top10.tsv: the junction (intron, 1-based inclusive), the full LSV it belongs to (min start .. max end
over all member junctions of the table), and a browser window = the LSV extent padded by max(PAD_MIN, PAD_FRAC x
span). Contigs get UCSC names (1 -> chr1, MT -> chrM). Writes results/top10_loci.tsv and results/top10_junctions.bed
(BED6, one line per top-10 junction, name = <contrast>|<rank>|<gene>, score = 1000 x dPSI)."""
import os, sys
from collections import defaultdict
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import *
PAD_MIN, PAD_FRAC = 500, 0.2

def main():
    rn = junction_ids(COUNTS_FRAG); ext = {}
    for r in rn:
        d = parse_row(r); lo, hi = ext.get(d["lsv"], (d["start"], d["end"])); ext[d["lsv"]] = (min(lo, d["start"]), max(hi, d["end"]))
    members = defaultdict(int)
    for r in rn: members[r.split("|")[0]] += 1
    out, bed = [], []
    for r in read_tsv("ds_top10.tsv"):
        ch = ucsc_chrom(r["chrom"]); js, je = int(r["start"]), int(r["end"]); ls, le = ext[r["lsv"]]; span = le - ls; pad = int(max(PAD_MIN, PAD_FRAC * span))
        out.append({"contrast": r["contrast"], "rank": r["rank"], "gene": r["gene"], "hit": r["hit"], "dpsi": r["dpsi"], "padj": f"{float(r['padj']):.3g}",
                    "junction_id": r["junction_id"], "lsv": r["lsv"], "lsv_members": members[r["lsv"]], "ucsc_chrom": ch, "strand": r["strand"],
                    "junction_start_1based": js + 1, "junction_end": je, "junction_position": f"{ch}:{js + 1}-{je}",
                    "lsv_start_1based": ls + 1, "lsv_end": le, "lsv_position": f"{ch}:{ls + 1}-{le}",
                    "window_position": f"{ch}:{max(1, ls + 1 - pad)}-{le + pad}", "pad_bp": pad})
        bed.append(f"{ch}\t{js}\t{je}\t{r['contrast']}|{r['rank']}|{r['gene']}\t{int(round(1000 * float(r['dpsi'])))}\t{r['strand']}")
    write_tsv(os.path.join(RESULTS, "top10_loci.tsv"), out, list(out[0].keys()))
    with open(os.path.join(RESULTS, "top10_junctions.bed"), "w") as fh:
        fh.write('track name="screen_splicing_top10" description="top-10 DS junctions per contrast (name = contrast|rank|gene; score = 1000 x dPSI)" useScore=0\n' + "\n".join(bed) + "\n")
    print(f"loci: {len(out)} rows -> top10_loci.tsv, top10_junctions.bed; window = LSV extent + max({PAD_MIN} bp, {PAD_FRAC:.0%} of span)")
    for o in out[:3] + out[-2:]: print(f"  {o['contrast']:<11}{o['rank']:>3} {o['gene']:<10} {o['junction_position']:<32} LSV {o['lsv_members']} members {o['lsv_position']:<32} window {o['window_position']}")

if __name__ == "__main__":
    main()
