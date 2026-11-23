#!/usr/bin/env python3
"""Independent re-analysis of results_*.csv (written by o3_low/high.lomt).

Recomputes C_conn, a thresholded log-linear fit, a nearest-neighbour ratio
estimate, and the second-moment xi from the raw C_raw columns.
"""
import sys, math
from collections import defaultdict

path = sys.argv[1] if len(sys.argv) > 1 else "results_low_beta.csv"

rows = defaultdict(list)
summ = {}
for line in open(path):
    line = line.strip()
    if not line:
        continue
    if line.startswith("# summary"):
        parts = line[len("# summary"):].strip().split(",")
        if len(parts) >= 6:
            b = float(parts[0])
            summ[b] = dict(accept=float(parts[1]), absm=float(parts[2]),
                           chi=float(parts[3]), xi2=float(parts[4]), xiexp=float(parts[5]))
        continue
    if line.startswith("beta,"):
        continue
    parts = line.split(",")
    if len(parts) != 4:
        continue
    try:
        b = float(parts[0]); r = int(parts[1]); craw = float(parts[2])
    except ValueError:
        continue
    rows[b].append((r, craw))

print(f"{'beta':>6} {'plateau':>9} {'xi_ratio':>8} {'xi_fit':>7} {'xi2':>7} "
      f"{'|m|':>7} {'chi':>8} {'acc':>6} {'nfit':>5}")
for b in sorted(rows):
    rs = sorted(rows[b])
    if not rs:
        continue
    plateau = rs[-1][1]
    conn = [(r, c - plateau) for r, c in rs]
    c0 = conn[0][1]
    xir = float("nan")
    if len(conn) > 1 and conn[1][1] > 0 and c0 > 0:
        ratio = conn[1][1] / c0
        if 0 < ratio < 1:
            xir = -1.0 / math.log(ratio)
    thr = 0.01 * c0 if c0 > 0 else 0.0
    xs = [(r, math.log(c)) for r, c in conn
          if 1 <= r <= min(8, rs[-1][0] - 1) and c > max(thr, 1e-9)]
    xif = float("nan")
    if len(xs) >= 2:
        n = len(xs)
        sx = sum(x for x, _ in xs); sy = sum(y for _, y in xs)
        sxx = sum(x * x for x, _ in xs); sxy = sum(x * y for x, y in xs)
        den = n * sxx - sx * sx
        if den != 0:
            slope = (n * sxy - sx * sy) / den
            if slope < 0:
                xif = -1.0 / slope
    num = sum(r * r * c for r, c in conn)
    den2 = sum(c for _, c in conn)
    xi2 = math.sqrt(num / den2) if den2 > 0 and num > 0 else float("nan")
    s = summ.get(b, {})
    print(f"{b:6.2f} {plateau:9.4f} {xir:8.3f} {xif:7.3f} {xi2:7.3f} "
          f"{s.get('absm', float('nan')):7.3f} {s.get('chi', float('nan')):8.3f} "
          f"{s.get('accept', float('nan')):6.3f} {len(xs):5d}")
