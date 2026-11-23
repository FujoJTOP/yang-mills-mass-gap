import math

rows = []
summary = []
with open('results.csv') as f:
    for line in f:
        line = line.strip()
        if not line:
            continue
        if line.startswith('# summary'):
            parts = [p for p in line[len('# summary'):].strip().split(',') if p != '']
            if parts:
                summary.append([float(x) for x in parts])
        else:
            parts = line.split(',')
            if len(parts) == 4:
                try:
                    rows.append((float(parts[0]), int(parts[1]),
                                 float(parts[2]), float(parts[3])))
                except ValueError:
                    continue

betas = sorted({r[0] for r in rows})
out = []
out.append('# O(3) L=32 results\n')
out.append('## Summary\n')
out.append('| beta | acceptance | <|m|> | chi | xi_2nd | xi_exp |')
out.append('|---|---|---|---|---|---|')
for s in summary:
    if len(s) >= 6:
        out.append('| %.2f | %.3f | %.4f | %.3f | %.3f | %.3f |' %
                   (s[0], s[1], s[2], s[3], s[4], s[5]))

out.append('\n## Correlation C_conn(r)\n')
for b in betas:
    sub = sorted([r for r in rows if r[0] == b], key=lambda t: t[1])
    out.append('### beta = %.2f\n' % b)
    out.append('| r | C_raw | C_conn |')
    out.append('|---|---|---|')
    for _, r, craw, cconn in sub:
        out.append('| %d | %.5f | %.5f |' % (r, craw, cconn))
    # log-slope over positive C_conn for r>=1
    pts = [(r, cconn) for _, r, _, cconn in sub if r >= 1 and cconn > 1e-6]
    if len(pts) >= 2:
        xs = [p[0] for p in pts]
        ys = [math.log(p[1]) for p in pts]
        n = len(xs)
        mx = sum(xs) / n
        my = sum(ys) / n
        den = sum((x - mx) ** 2 for x in xs)
        slope = sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / den if den else 0.0
        out.append('\nlog-slope = %.3f, xi_log = %.3f (n=%d)\n' %
                   (slope, -1.0 / slope if slope < 0 else float('inf'), n))
    out.append('')

open('SUMMARY.md', 'w', encoding='utf-8').write('\n'.join(out))
print('\n'.join(out))
