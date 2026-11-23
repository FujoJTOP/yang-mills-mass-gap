"""Markdown/CSV summary tables from results.json."""
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
res = json.load(open(os.path.join(HERE, 'results.json')))
res.sort(key=lambda r: (r['L'], r['beta']))

hdr = ("| beta | L | b | k3(0,0) alpha=1 | c3_corr | c3(>2sig) | alpha* | "
       "k3(alpha*) | field k3 | max|V| | top1% k3 |")
print(hdr)
print("|" + "---|" * 11)
for r in res:
    for b in (1, 2):
        if f'b{b}_k3_000' not in r:
            continue
        print(f"| {r['beta']:.2f} | {r['L']} | {b} | "
              f"{r[f'b{b}_a1.0_k3']:+.4f}({r[f'b{b}_a1.0_ek3']:.4f}) | "
              f"{r[f'b{b}_c3_corr']:+.3f} | {r[f'b{b}_c3_2sig']:.3f} | "
              f"{r[f'b{b}_alpha_min']:.2f} | {r[f'b{b}_k3_min']:+.4f} | "
              f"{r[f'b{b}_field_k3']:+.2e} | {r[f'b{b}_maxabsV']:.1f} | "
              f"{r[f'b{b}_k3_top1_share']:.3f} |")

print("\n### decay / correlation length")
print("| beta | L | b | C2(0) | C2(1) | C2(1)/C2(0) | white ratio b2 |")
print("|---|---|---|---|---|---|---|")
for r in res:
    for b in (1, 2):
        c2 = r.get(f'b{b}_C2')
        if not c2:
            continue
        print(f"| {r['beta']:.2f} | {r['L']} | {b} | {c2[0]:.4f} | "
              f"{c2[1]:+.4f} | {c2[1]/c2[0]:+.5f} | "
              f"{r.get('white_ratio','') if b==2 else ''} |")

print("\n### alpha scan of k3 (b=1)   [alpha: k3 +- err]")
for r in res:
    row = " ".join(f"{a}:{r[f'b1_a{a}_k3']:+.3f}" for a in
                   (0.0, 0.25, 0.5, 0.6, 0.75, 1.0, 1.25, 1.5)
                   if f'b1_a{a}_k3' in r)
    print(f"| {r['beta']:.2f} L={r['L']} | {row} |")
