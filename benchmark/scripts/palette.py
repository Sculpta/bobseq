"""One palette for every preprint panel: a green for BOBseq, a crimson for DRUG-seq, and two hues built at the same OKLCH
lightness and chroma for the other two methods.
"""

COL = {'DRUG-seq': '#A8325A', 'BRB-seq': '#2260A9', 'prime-seq': '#8E4E00', 'BOBseq': '#2E6B32'}
METH = [
    'DRUG-seq',
    'prime-seq',
    'BOBseq',
]  # BRB-seq is not a benchmark method (per-sample RNA purification, not a cells-in-well multiplexed method); its colour is kept for extra/ material
# read-composition classes of the preprint's classification rules (composition_rules.py -> composition_panels.py), stacked-bar order
CLS = [
    ('mRNA', 'mRNA (protein-coding)', '#4C9A4C'),
    ('ribosomal-protein', 'ribosomal-protein mRNA', '#6E6E6E'),
    ('exonic other biotype', 'exonic, other biotype', '#9E9E9E'),
    ('intronic', 'intronic', '#C6C6C6'),
    ('intergenic', 'intergenic', '#B5379B'),
    ('rRNA', 'rRNA', '#E6B422'),
    ('mitochondrial', 'mitochondrial', '#2260A9'),
]  # intergenic magenta, mRNA green, rRNA yellow, MT blue, others greys
