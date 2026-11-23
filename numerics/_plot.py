"""Minimal dependency-free SVG plotter + text table writer.

Only stdlib. Produces self-contained SVG files (no external fonts required).
Used by the Yang-Mills RG-handoff numerics scripts.
"""
import math

PALETTE = ["#1f77b4", "#d62728", "#2ca02c", "#9467bd", "#ff7f0e",
           "#17becf", "#8c564b", "#e377c2", "#bcbd22", "#7f7f7f"]


class Curve:
    def __init__(self, label, pts, color=None, width=2.0, dash=None, marker=False):
        self.label = label
        self.pts = [(float(x), float(y)) for (x, y) in pts if _finite(x) and _finite(y)]
        self.color = color
        self.width = width
        self.dash = dash          # e.g. "6,4"
        self.marker = marker      # draw small circles at points


def _finite(v):
    try:
        return math.isfinite(float(v))
    except (TypeError, ValueError):
        return False


def _nice_lin(vmin, vmax, n=6):
    if vmax <= vmin:
        vmax = vmin + 1.0
    span = vmax - vmin
    raw = span / max(n, 1)
    mag = 10.0 ** math.floor(math.log10(raw)) if raw > 0 else 1.0
    norm = raw / mag
    if norm < 1.5:
        step = mag
    elif norm < 3.0:
        step = 2.0 * mag
    elif norm < 7.0:
        step = 5.0 * mag
    else:
        step = 10.0 * mag
    start = math.ceil(vmin / step - 1e-9) * step
    ticks = []
    v = start
    while v <= vmax + 1e-9 * abs(step):
        ticks.append(v)
        v += step
    return ticks, step


def _log_ticks(vmin, vmax):
    lo = math.floor(math.log10(vmin))
    hi = math.ceil(math.log10(vmax))
    return [10.0 ** k for k in range(lo, hi + 1)]


def _fmt_tick(v, logmode):
    if logmode:
        if v != 0 and (abs(math.log10(abs(v)) - round(math.log10(abs(v)))) < 1e-6):
            e = int(round(math.log10(abs(v))))
            return "1e%d" % e
        return "%.2g" % v
    if v == 0:
        return "0"
    av = abs(v)
    if av >= 1e4 or av < 1e-3:
        return "%.2g" % v
    if av >= 100:
        return "%.0f" % v
    if av >= 10:
        return "%.1f" % v
    return "%.2f" % v


def svg_plot(path, curves, xlabel, ylabel, title,
             xlog=False, ylog=False,
             width=780, height=500,
             vlines=None, hlines=None, annotations=None,
             legend_loc="upper right", invert_x=False):
    """Render curves to an SVG file.

    curves       : list of Curve
    vlines/hlines: list of (value, label, color) or list of values
    annotations  : list of (x, y, text, color)
    """
    vlines = vlines or []
    hlines = hlines or []
    annotations = annotations or []

    # normalise vlines/hlines entries
    def _norm(items):
        out = []
        for it in items:
            if isinstance(it, (list, tuple)):
                v = it[0]
                lab = it[1] if len(it) > 1 else ""
                col = it[2] if len(it) > 2 else "#888888"
            else:
                v, lab, col = it, "", "#888888"
            out.append((float(v), lab, col))
        return out
    vlines = _norm(vlines)
    hlines = _norm(hlines)

    allx, ally = [], []
    for c in curves:
        for (x, y) in c.pts:
            allx.append(math.log10(x) if (xlog and x > 0) else x)
            ally.append(math.log10(y) if (ylog and y > 0) else y)
    if not allx:
        allx, ally = [0, 1], [0, 1]

    for v, _, _ in vlines:
        allx.append(math.log10(v) if (xlog and v > 0) else v)
    for v, _, _ in hlines:
        ally.append(math.log10(v) if (ylog and v > 0) else v)
    for (x, y, *_rest) in annotations:
        allx.append(math.log10(x) if (xlog and x > 0) else x)
        ally.append(math.log10(y) if (ylog and y > 0) else y)

    xmin, xmax = min(allx), max(allx)
    ymin, ymax = min(ally), max(ally)
    xpad = (xmax - xmin) * 0.06 or 1.0
    ypad = (ymax - ymin) * 0.08 or 1.0
    xmin -= xpad; xmax += xpad
    ymin -= ypad; ymax += ypad

    ml, mr, mt, mb = 78, 30, 52, 62
    pw = width - ml - mr
    ph = height - mt - mb

    def X(x):
        tx = math.log10(x) if (xlog and x > 0) else x
        f = (tx - xmin) / (xmax - xmin) if xmax > xmin else 0.5
        if invert_x:
            f = 1.0 - f
        return ml + f * pw

    def Y(y):
        ty = math.log10(y) if (ylog and y > 0) else y
        f = (ty - ymin) / (ymax - ymin) if ymax > ymin else 0.5
        return height - mb - f * ph

    out = []
    out.append('<?xml version="1.0" encoding="UTF-8"?>')
    out.append('<svg xmlns="http://www.w3.org/2000/svg" width="%d" height="%d" '
               'viewBox="0 0 %d %d" font-family="DejaVu Sans, Arial, sans-serif">'
               % (width, height, width, height))
    out.append('<rect width="%d" height="%d" fill="#ffffff"/>' % (width, height))
    # title
    out.append('<text x="%d" y="30" font-size="17" font-weight="bold" '
               'text-anchor="middle" fill="#111">%s</text>' % (width // 2, _esc(title)))

    # plot frame
    out.append('<rect x="%d" y="%d" width="%d" height="%d" fill="none" stroke="#333" '
               'stroke-width="1"/>' % (ml, mt, pw, ph))

    # ticks  (X/Y already apply log transforms when needed)
    if xlog:
        xtval = [t for t in _log_ticks(10 ** xmin, 10 ** xmax)
                 if xmin <= math.log10(t) <= xmax]
    else:
        xtval, _ = _nice_lin(xmin, xmax, 7)
    for tv in xtval:
        px = X(tv)
        out.append('<line x1="%.1f" y1="%d" x2="%.1f" y2="%d" stroke="#e2e2e2"/>'
                   % (px, mt, px, mt + ph))
        out.append('<line x1="%.1f" y1="%d" x2="%.1f" y2="%d" stroke="#333"/>'
                   % (px, mt + ph, px, mt + ph + 5))
        out.append('<text x="%.1f" y="%d" font-size="11" text-anchor="middle" fill="#222">%s</text>'
                   % (px, mt + ph + 18, _fmt_tick(tv, xlog)))

    if ylog:
        ytval = [t for t in _log_ticks(10 ** ymin, 10 ** ymax)
                 if ymin <= math.log10(t) <= ymax]
    else:
        ytval, _ = _nice_lin(ymin, ymax, 7)
    for tv in ytval:
        py = Y(tv)
        out.append('<line x1="%d" y1="%.1f" x2="%d" y2="%.1f" stroke="#e2e2e2"/>'
                   % (ml, py, ml + pw, py))
        out.append('<line x1="%d" y1="%.1f" x2="%d" y2="%.1f" stroke="#333"/>'
                   % (ml - 5, py, ml, py))
        out.append('<text x="%d" y="%.1f" font-size="11" text-anchor="end" fill="#222">%s</text>'
                   % (ml - 8, py + 4, _fmt_tick(tv, ylog)))

    # axis labels
    out.append('<text x="%d" y="%d" font-size="13" text-anchor="middle" fill="#111">%s</text>'
               % (ml + pw // 2, height - 16, _esc(xlabel)))
    out.append('<text x="18" y="%d" font-size="13" text-anchor="middle" fill="#111" '
               'transform="rotate(-90 18 %d)">%s</text>'
               % (mt + ph // 2, mt + ph // 2, _esc(ylabel)))

    # vlines / hlines
    for v, lab, col in vlines:
        if xlog and v <= 0:
            continue
        px = X(v)
        if not (ml - 1 <= px <= ml + pw + 1):
            continue
        out.append('<line x1="%.1f" y1="%d" x2="%.1f" y2="%d" stroke="%s" '
                   'stroke-width="1.2" stroke-dasharray="6,4"/>'
                   % (px, mt, px, mt + ph, col))
        if lab:
            out.append('<text x="%.1f" y="%d" font-size="10" fill="%s" '
                       'text-anchor="end">%s</text>' % (px - 3, mt + 13, col, _esc(lab)))
    for v, lab, col in hlines:
        if ylog and v <= 0:
            continue
        py = Y(v)
        if not (mt - 1 <= py <= mt + ph + 1):
            continue
        out.append('<line x1="%d" y1="%.1f" x2="%d" y2="%.1f" stroke="%s" '
                   'stroke-width="1.2" stroke-dasharray="6,4"/>'
                   % (ml, py, ml + pw, py, col))
        if lab:
            out.append('<text x="%d" y="%.1f" font-size="10" fill="%s">%s</text>'
                       % (ml + 6, py - 4, col, _esc(lab)))

    # curves
    for i, c in enumerate(curves):
        col = c.color or PALETTE[i % len(PALETTE)]
        if not c.pts:
            continue
        pts = sorted(c.pts)
        # split into contiguous runs for NaNs handled by _finite filter above
        d = "M " + " L ".join("%.2f,%.2f" % (X(x), Y(y)) for (x, y) in pts)
        dash = (' stroke-dasharray="%s"' % c.dash) if c.dash else ""
        out.append('<path d="%s" fill="none" stroke="%s" stroke-width="%.1f"%s/>'
                   % (d, col, c.width, dash))
        if c.marker:
            for (x, y) in pts:
                out.append('<circle cx="%.2f" cy="%.2f" r="2.2" fill="%s"/>'
                           % (X(x), Y(y), col))

    # annotations (point labels)
    for ann in annotations:
        x, y, text, color = ann[0], ann[1], ann[2], (ann[3] if len(ann) > 3 else "#222")
        px, py = X(x), Y(y)
        out.append('<circle cx="%.2f" cy="%.2f" r="3.4" fill="%s"/>' % (px, py, color))
        out.append('<text x="%.2f" y="%.2f" font-size="10.5" fill="%s">%s</text>'
                   % (px + 6, py - 6, color, _esc(text)))

    # legend
    if len(curves) > 1:
        ly = mt + 14
        if legend_loc == "upper left":
            lx = ml + 14
            anchor = "start"
            sw_x1, sw_x2 = lx, lx + 22
            tx = lx + 28
        else:
            lx = ml + pw - 14
            anchor = "end"
            sw_x1, sw_x2 = lx - 26, lx - 4
            tx = lx
        for i, c in enumerate(curves):
            col = c.color or PALETTE[i % len(PALETTE)]
            yy = ly + i * 18
            dash = (' stroke-dasharray="%s"' % c.dash) if c.dash else ""
            out.append('<line x1="%d" y1="%d" x2="%d" y2="%d" stroke="%s" '
                       'stroke-width="2.4"%s/>' % (sw_x1, yy, sw_x2, yy, col, dash))
            out.append('<text x="%d" y="%d" font-size="11" text-anchor="%s" fill="#111">%s</text>'
                       % (tx, yy + 4, anchor, _esc(c.label)))

    out.append('</svg>')
    with open(path, "w", encoding="utf-8") as f:
        f.write("\n".join(out))


def _esc(s):
    return (str(s).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;"))


def write_table(path, title, headers, rows, note=None):
    """Write a plain-text aligned table."""
    cols = list(zip(*([headers] + rows))) if rows else [[h] for h in headers]
    widths = []
    for col in cols:
        widths.append(max(len(str(x)) for x in col))
    lines = []
    lines.append(title)
    lines.append("=" * min(100, max(40, sum(widths) + 2 * len(widths))))
    fmt = "  ".join("{:>%d}" % w for w in widths)
    lines.append(fmt.format(*[str(h) for h in headers]))
    lines.append("-" * min(100, max(40, sum(widths) + 2 * len(widths))))
    for r in rows:
        lines.append(fmt.format(*[str(x) for x in r]))
    if note:
        lines.append("")
        lines.append(note)
    with open(path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
