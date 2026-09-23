"""The positive-control splicing events drawn by splicing_tracks: five exons the SMA literature reports as risdiplam
off-targets (STRN3, FOXM1, APLP2, MADD, SLC25A17) and three NMD-coupled poison exons that respond to cycloheximide
(CASP2, SRSF3, TRA2B), one event per gene.

Risdiplam is an SMN2 exon-7 5'ss modifier; the same chemistry stabilises the U1:5'ss duplex at other exons whose 5'ss carries
the same GA|GUAAGU context. Coordinates are GRCh38 / Ensembl 113, 1-based inclusive for exons; junctions are 0-based half-open
intron intervals (BED convention), as in the junction tables.
  incl  junctions that can only come from a transcript CONTAINING the target exon / extension
  skip  the junction that bypasses it (shares an anchor with an incl junction, so the two compete directly)
  ss5   the target exon's own 5' splice site (3 exonic | 8 intronic, transcript orientation)
"""
# MANE Select transcript per gene (Ensembl 113) and the locus to slice (gene span + 2 kb margins)
GENES = {
    "STRN3":    dict(tx="ENST00000357479", gene="ENSG00000196792", chrom="14", strand="-", span=(30893799, 31026401)),
    "FOXM1":    dict(tx="ENST00000359843", gene="ENSG00000111206", chrom="12", strand="-", span=(2857680, 2877174)),
    "APLP2":    dict(tx="ENST00000338167", gene="ENSG00000084234", chrom="11", strand="+", span=(130068147, 130144811)),
    "MADD":     dict(tx="ENST00000706887", gene="ENSG00000110514", chrom="11", strand="+", span=(47269161, 47330031)),
    "SLC25A17": dict(tx="ENST00000435456", gene="ENSG00000100372", chrom="22", strand="-", span=(40769630, 40819399)),
    # three NMD-coupled poison exons, readouts of the CHX arm
    "CASP2":    dict(tx="ENST00000310447", gene="ENSG00000106144", chrom="7", strand="+", span=(143288215, 143307696)),
    "SRSF3":    dict(tx="ENST00000373715", gene="ENSG00000112081", chrom="6", strand="+", span=(36594353, 36605911)),
    "TRA2B":    dict(tx="ENST00000453386", gene="ENSG00000136527", chrom="3", strand="-", span=(185914558, 185938103)),
}
# ---------------------------------------------------------------------------------------------------------------
# The target event per gene: the junction usage that changes between control and Ris 25 mM (or CHX) in the data,
# annotated against Ensembl 113 and frozen here. Coordinates are 1-based inclusive for exons, 0-based half-open (BED) for the
# introns of each junction, matching the CIGAR-derived junctions of stage `fragments`.
#   incl  junctions that can only come from a transcript CONTAINING the target exon / extension
#   skip  the junction that bypasses it (shares an anchor with an incl junction, so the two compete directly)
#   ss5   the target exon's own 5' splice site (3 exonic | 8 intronic, transcript orientation): risdiplam acts by
#         stabilising U1 snRNP on a weak 5'ss whose last exonic base is A and whose intron starts GUAAG (as in SMN2
#         exon 7, GGA|GUAAGUCU) -- so the motif is the mechanistic check that an event can be a risdiplam target.
EVENTS = {
 "FOXM1": dict(exon=(2861299, 2861412), exon_nt=114, label="exon 9", ss5="TCA|GTAAGTCT", motif=True,
               note="cassette exon 9 of FOXM1-201 (protein-coding), absent from the MANE transcript. It leaves through "
                    "either of two donors: its own 5'ss at 2,861,299 (TCA|GTAAGTCT, the full 114 nt exon) or an internal "
                    "one at 2,861,379 (9S, TGA|GTAAGTTC, the 34 nt form). Both are canonical, both carry the A|GUAAG "
                    "risdiplam motif and both pair with the same acceptor, so they are KEPT APART and the schematic draws "
                    "the 80 nt between them as its own box, as it does for MADD's exon-13 extension. Numbering follows "
                    "FOXM1-201: ex8 - ex9 - ex10 (last).",
               incl=[("ex8>ex9", (2861412, 2864319)), ("ex9S>ex10", (2859663, 2861378)), ("ex9>ex10", (2859663, 2861298))],
               skip=[("ex8>ex10", (2859663, 2864319))]),
 "SLC25A17": dict(exon=(40797284, 40797336), exon_nt=53, label="poison exon (intron 2)", ss5="TGA|GTAAGATT", motif=True,
               note="NMD-transcript exon 3 of SLC25A17-213/-214 (53 nt; a 83-nt form starting 40,797,254 exists in -209/-211 "
                    "but its 5'ss, AAG|GTTAAACC, does not carry the motif). Inclusion puts a premature stop in the mRNA.",
               incl=[("ex2>PE", (40797336, 40799022)), ("PE>ex3", (40794580, 40797283)), ("PE83>ex3", (40794580, 40797253))],
               skip=[("ex2>ex3", (40794580, 40799022))]),
 "STRN3": dict(exon=(30929201, 30929311), exon_nt=111, label="exon 8", ss5="AGA|GTAAGTGC", motif=True,
               note="cassette exon 8 (in MANE); the exon the SMA literature reports as a risdiplam / branaplam target.",
               incl=[("ex7>ex8", (30929311, 30935162)), ("ex8>ex9", (30919106, 30929200))],
               skip=[("ex7>ex9", (30919106, 30935162))]),
 "APLP2": dict(exon=(130123612, 130123779), exon_nt=168, label="exon 7 (KPI domain)", ss5="TGA|GTAAGTCC", motif=True,
               note="cassette exon 7 (in MANE), the Kunitz protease-inhibitor domain of APLP2.",
               incl=[("ex6>ex7", (130122513, 130123611)), ("ex7>ex8", (130123779, 130126699))],
               skip=[("ex6>ex8", (130122513, 130126699))]),
 "MADD": dict(exon=(47285066, 47285194), exon_nt=129, label="exon 13 extension (distal 5'ss)", ss5="TCG|GTGAGAGC", motif=False,
               note="alternative 5'ss of exon 13, not a cassette: the distal donor (47,285,194) extends exon 13 by 129 nt, the "
                    "proximal one (47,285,065) ends it early. NEITHER donor carries the A|GUAAG motif -- see README.",
               incl=[("ex13long>ex14", (47285194, 47285450))],
               skip=[("ex13short>ex14", (47285065, 47285450))]),
 # --- three NMD-coupled exons; exon numbering follows the isoform that carries the exon (the MANE transcript lacks it,
 # so MANE's exon n becomes n+1 here), as for FOXM1 above.
 "SRSF3": dict(exon=(36599821, 36600276), exon_nt=456, label="exon 4 (poison exon)", ss5="CTA|GTAAGTTT", motif=True,
               note="exon 4 of SRSF3-203 (biotype nonsense_mediated_decay), the classic SR-protein poison exon: inclusion "
                    "brings a premature stop and the transcript is degraded by NMD. Numbering follows SRSF3-203: "
                    "ex3 - ex4 - ex5 (= MANE exon 4).",
               incl=[("ex3>ex4", (36598983, 36599820)), ("ex4>ex5", (36600276, 36601151))],
               skip=[("ex3>ex5", (36598983, 36601151))]),
 "CASP2": dict(exon=(143300380, 143300440), exon_nt=61, label="exon 9", ss5="TCT|GTAAGTGT", motif=False,
               note="the 61 nt exon 9 of CASP2-206; including it gives Casp2S, which carries a premature stop and is "
                    "degraded by NMD. Numbering follows CASP2-206: ex8 - ex9 - ex10 (= MANE exon 9).",
               incl=[("ex8>ex9", (143300294, 143300379)), ("ex9>ex10", (143300440, 143303783))],
               skip=[("ex8>ex10", (143300294, 143303783))]),
 "TRA2B": dict(exon=(185931577, 185931852), exon_nt=276, label="exon 2 (poison exon)", ss5="TAA|GTAATTAC", motif=False,
               note="exon 2 of TRA2B-205 (biotype nonsense_mediated_decay), the TRA2B poison exon. Numbering follows "
                    "TRA2B-205: ex1 - ex2 - ex3 (= MANE exon 2) - ex4 (= MANE exon 3); both junctions that leave exon 1 "
                    "without the poison exon count as skips.",
               incl=[("ex1>ex2", (185931852, 185937824)), ("ex2>ex3", (185926734, 185931576))],
               skip=[("ex1>ex3", (185926734, 185937824)), ("ex1>ex4", (185925626, 185937824))]),
}
for _g, _d in EVENTS.items(): _d["chrom"], _d["strand"] = GENES[_g]["chrom"], GENES[_g]["strand"]
