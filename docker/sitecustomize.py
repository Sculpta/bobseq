"""Font alias for the container: the figure scripts set rcParams['font.family'] = 'Arial'. Debian has no Arial;
Liberation Sans (fonts-liberation) is metric-compatible with it, so its faces are registered under the additional
name 'Arial' in matplotlib's font list. matplotlib matches on a font's internal family name and, with font.family
set explicitly, never consults the sans-serif fallback list, so installing the font alone is not enough.
Fails open: any problem here leaves matplotlib exactly as it was."""

try:
    import matplotlib
    matplotlib.use('Agg')
    from matplotlib import font_manager as _fm

    if 'Arial' not in {f.name for f in _fm.fontManager.ttflist}:
        for _f in [f for f in _fm.fontManager.ttflist if f.name == 'Liberation Sans']:
            _fm.fontManager.ttflist.append(_fm.FontEntry(
                fname=_f.fname, name='Arial', style=_f.style, variant=_f.variant,
                weight=_f.weight, stretch=_f.stretch, size=_f.size))
except Exception:
    pass
