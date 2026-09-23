"""Track colours = the 24-plex condition palette of the screen analyses (screen_splicing, screen_pca): control dark grey,
CHX 1 µg/mL light orange, CHX 50 µg/mL dark orange, Ris 25 mM blue. Replicate tracks take their condition's colour.
"""
COL = {'control': '#444444', 'chx_low': '#fdae6b', 'chx_high': '#d94801', 'ris_low': '#1f77b4'}
ARC = {'control': '#1a1a1a', 'chx_low': '#c2661a', 'chx_high': '#7f2a00', 'ris_low': '#0b4a7a'}   # junction arcs and their count labels, a darker shade
GENECOL = '#2F6DB5'                                                       # transcript model (as in the benchmark script)
EVENTCOL = '#b2182b'                                                      # the event exon when it is not in the canonical transcript
