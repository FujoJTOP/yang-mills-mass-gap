#!/usr/bin/env python3
"""Analyze the 4D Z2 parallel-tempering raw samples.

raw.csv columns: sweep,r,beta,plaq,W11,W12,W13,W21,W22,W23,W31,W32,W33
Rows with a leading '#' are comments/header.
"""
import csv
import math
import sys
from collections import defaultdict

OBS = ["plaq", "W11", "W12", "W13", "W21", "W22", "W23", "W31", "W32", "W33"]


def mean(xs):
    return sum(xs) / len(xs)


def naive_err(xs):
    n = len(xs)
    if n < 2:
        return float("nan")
    m = mean(xs)
    var = sum((x - m) ** 2 for x in xs) / (n - 1)
    return math.sqrt(var / n)


def blocking_err(xs, block=16):
    n = len(xs)
    nb = n // block
    if nb < 4:
        return naive_err(xs)
    bm = [mean(xs[i * block:(i + 1) * block]) for i in range(nb)]
    return naive_err(bm)


def main(path="raw.csv"):
    by_beta = defaultdict(lambda: defaultdict(list))
    with open(path, newline="") as f:
        for line in f:
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            parts = line.strip().split(",")
            if parts[0] != "S":
                continue
            beta = float(parts[2])
            vals = [float(x) for x in parts[3:]]
            if len(vals) != len(OBS):
                continue
            for name, v in zip(OBS, vals):
                by_beta[beta][name].append(v)

    betas = sorted(by_beta)
    print(f"# betas: {', '.join(f'{b:.4f}' for b in betas)}")
    print(f"# samples per beta: {len(by_beta[betas[0]]['plaq'])}")
    print()
    header = ["beta"] + [f"{o}_mean" for o in OBS[:4]] + [f"{o}_err" for o in OBS[:4]]
    print(",".join(header))

    stats = {}
    for b in betas:
        stats[b] = {}
        for o in OBS:
            xs = by_beta[b][o]
            stats[b][o] = (mean(xs), blocking_err(xs), naive_err(xs))
        row = [f"{b:.4f}"] + [f"{stats[b][o][0]:.6f}" for o in OBS[:4]] + \
              [f"{stats[b][o][1]:.6f}" for o in OBS[:4]]
        print(",".join(row))

    print()
    print("# Wilson loops W(R,T)")
    print("beta,R,T,W,err_block,err_naive")
    for b in betas:
        for R in (1, 2, 3):
            for T in (1, 2, 3):
                m, eb, en = stats[b][f"W{R}{T}"]
                print(f"{b:.4f},{R},{T},{m:.6f},{eb:.6f},{en:.6f}")

    print()
    print("# Creutz ratios chi(R,T) = -ln[W(R,T)W(R-1,T-1)/(W(R-1,T)W(R,T-1))]")
    print("beta,R,T,chi,note")

    def w(b, R, T):
        if R == 0 or T == 0:
            return 1.0
        return stats[b][f"W{R}{T}"][0]

    for b in betas:
        for R in (2, 3):
            for T in (2, 3):
                num = w(b, R, T) * w(b, R - 1, T - 1)
                den = w(b, R - 1, T) * w(b, R, T - 1)
                if num <= 0 or den <= 0:
                    print(f"{b:.4f},{R},{T},nan,nonpositive_W")
                else:
                    print(f"{b:.4f},{R},{T},{-math.log(num / den):.6f},ok")

    # Machine-readable results.csv
    with open("results.csv", "w", newline="") as f:
        wr = csv.writer(f)
        wr.writerow(["beta", "observable", "mean", "err_block", "err_naive", "n"])
        for b in betas:
            for o in OBS:
                m, eb, en = stats[b][o]
                wr.writerow([f"{b:.4f}", o, f"{m:.8f}", f"{eb:.8f}", f"{en:.8f}",
                             len(by_beta[b][o])])
        wr.writerow([])
        wr.writerow(["beta", "R", "T", "creutz", "note"])
        for b in betas:
            for R in (2, 3):
                for T in (2, 3):
                    num = w(b, R, T) * w(b, R - 1, T - 1)
                    den = w(b, R - 1, T) * w(b, R, T - 1)
                    if num <= 0 or den <= 0:
                        wr.writerow([f"{b:.4f}", R, T, "nan", "nonpositive_W"])
                    else:
                        wr.writerow([f"{b:.4f}", R, T, f"{-math.log(num / den):.8f}", "ok"])
    print("\n# wrote results.csv")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "raw.csv")
