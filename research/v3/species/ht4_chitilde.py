#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ht4_chitilde.py -- HT4: chi~(k_min;L) = chi(k_min;L) - C(0) at beta=2.5.

ZERO new MC.  Reads ONLY stored npz:
  partA_su2_N4000.npz   SU(2) beta=2.5 L=4  48 fixed transverse plane waves, a=<v,log U>
  partA_su2_N800.npz    SU(2) beta=2.5 L=6  48 modes
  k0_su2_b2.4_L8_n2.npz SU(2) beta=2.4 L=8  12 modes (n=1: k=pi/4, n=2: k=pi/2), N=400
  rerun_su2_L4.npz / rerun_su2_L6.npz  independent 6-mode cross-checks

Modes are UNIT-normalised in the (4,L^4,3) real algebra space (gpp_core.plane_mode:79),
so  Var(a) = (1/3)<|A(k)|^2> = chi(k)  up to O(1/N)  (cos-mode identity, k != 0,pi).

C(0) = (1/3)E|log U|^2 = (4pi^2-6)/9 = 3.719824178 (HR1-HAAR theorem, ALL beta, ALL L).
Errors: delete-block jackknife over configurations (20 blocks), applied to the
mode-averaged variance.  Deterministic; no MC.
"""
import numpy as np, json, os

HERE = os.path.dirname(os.path.abspath(__file__))
W7 = os.path.abspath(os.path.join(HERE, "..", "..", "unconditional", "w7"))
OUT = []
def P(*a):
    s = " ".join(str(x) for x in a)
    OUT.append(s); print(s, flush=True)

C0 = (4 * np.pi ** 2 - 6) / 9.0
P("C(0) = (4pi^2-6)/9 = %.9f   [HR1-HAAR: exact, all beta, all L]" % C0)

def load_modes(path):
    """partA-style: one <mode>_a key per mode;  k0-style: single 'a' (N,Nmode)."""
    z = np.load(path, allow_pickle=True)
    keys = sorted([k for k in z.files if k.endswith("_a")])
    if keys:
        A = np.stack([z[k] for k in keys], axis=1)      # (N, Nmode)
        names = [k[:-2] for k in keys]
    else:
        A = np.asarray(z["a"], float)
        names = [str(x) for x in z["names"]] if "names" in z.files else [f"m{i}" for i in range(A.shape[1])]
    meta = {k: z[k] for k in z.files if not k.endswith(("_a", "_H", "_G"))}
    return A, names, meta

def jack_var(A, nblk=20):
    """mean-over-modes variance + delete-block jackknife SE (blocks of configs)."""
    N, M = A.shape
    B = A - A.mean(0)                                    # mean-subtracted per mode
    per = (B ** 2).mean(0)                               # per-mode Var (population)
    est = per.mean()
    edges = np.linspace(0, N, nblk + 1).astype(int)
    jk = []
    for i in range(nblk):
        m = np.ones(N, bool); m[edges[i]:edges[i + 1]] = False
        jk.append(((B[m] ** 2).mean(0)).mean())
    jk = np.array(jk)
    se = np.sqrt((nblk - 1) / nblk * ((jk - jk.mean()) ** 2).sum())
    return est, se, per

def report(tag, path, kmode_filter=None):
    A, keys, meta = load_modes(path)
    N, M = A.shape
    if kmode_filter is not None:
        sel = [i for i, k in enumerate(keys) if kmode_filter(k)]
        A = A[:, sel]; keys = [keys[i] for i in sel]
    est, se, per = jack_var(A)
    q = np.percentile(per, [0, 25, 50, 75, 100])
    P(f"{tag:34s} N={N:5d} M={M:3d}  <a>={A.mean():+.4f}  "
      f"chi=<Var(a)>={est:.4f}+-{se:.4f}  (noise={se:.4f})")
    P(f"{'':34s} per-mode Var quartiles [{q[0]:.3f},{q[1]:.3f},{q[2]:.3f},{q[3]:.3f},{q[4]:.3f}]")
    P(f"{'':34s} chi~ = chi - C(0) = {est - C0:+.4f} +- {se:.4f}   ({abs(est-C0)/se:.2f} sigma)")
    return dict(N=N, M=M, chi=est, se=se, chit=est - C0, keys=keys,
                mean=float(A.mean()), q=[float(x) for x in q])

res = {}
P("")
P("=" * 78)
P("[1] beta = 2.5, L = 4,6  (48-mode PartA sets, k_min = 2pi/L)")
P("=" * 78)
res["L4"] = report("SU2 b2.5 L4 k=pi/2", os.path.join(W7, "partA_su2_N4000.npz"))
res["L6"] = report("SU2 b2.5 L6 k=pi/3", os.path.join(W7, "partA_su2_N800.npz"))

P("")
P("=" * 78)
P("[2] beta = 2.4, L = 8  (k0 run, n=1 -> k=pi/4 = k_min; n=2 -> k=pi/2)")
P("=" * 78)
z = np.load(os.path.join(W7, "k0_su2_b2.4_L8_n2.npz"), allow_pickle=True)
nm = [str(x) for x in z["names"]]
A8 = z["a"]
sel1 = [i for i, n in enumerate(nm) if n.endswith("_n1")]
sel2 = [i for i, n in enumerate(nm) if n.endswith("_n2")]
P(f"  stored names: {nm}")
P(f"  kvals: {np.round(z['kvals'],5)}  beta={float(z['beta'])}  L={int(z['L'])}  N={A8.shape[0]}")
res["L8"] = report("SU2 b2.4 L8 k=pi/4 (k_min)", os.path.join(W7, "k0_su2_b2.4_L8_n2.npz"),
                   kmode_filter=lambda k: k.endswith("_n1"))
res["L8_pi2"] = report("SU2 b2.4 L8 k=pi/2", os.path.join(W7, "k0_su2_b2.4_L8_n2.npz"),
                       kmode_filter=lambda k: k.endswith("_n2"))

P("")
P("=" * 78)
P("[3] independent 6-mode cross-checks (rerun sets, beta=2.5)")
P("=" * 78)
res["r4"] = report("rerun SU2 L4 (6 modes)", os.path.join(W7, "rerun_su2_L4.npz"))
res["r6"] = report("rerun SU2 L6 (6 modes)", os.path.join(W7, "rerun_su2_L6.npz"))

P("")
P("=" * 78)
P("[4] scale judgement  chi~ = A * L^p      (log-log LSQ on the amplitude)")
P("=" * 78)

def fit(pts, label):
    if len(pts) < 2:
        P(f"  {label}: <2 points, no fit"); return None
    L = np.array([p[0] for p in pts], float)
    y = np.array([p[1] for p in pts], float)
    e = np.array([p[2] for p in pts], float)
    if np.any(y <= 0):
        P(f"  {label}: non-positive chi~, log-fit undefined (values: "
          f"{['%.4f'%v for v in y]})")
        # linear-in-log fit with sign kept: y = a + b x  (absolute scale)
        x = np.log(L); b, a = np.polyfit(x, y, 1, w=1 / np.maximum(e, 1e-9))
        P(f"     LSQ y = a + b log L : b = {b:+.4f} +- ? (unweighted-sign fit, "
          f"y(4,6,8) = {['%.3f'%v for v in y]})")
        return None
    x = np.log(L); w = 1 / np.maximum(e, 1e-9)
    b, a = np.polyfit(x, y, 1, w=w)
    yhat = a + b * x
    chi2 = float(np.sum(((y - yhat) / np.maximum(e, 1e-9)) ** 2))
    P(f"  {label}: p = {b:+.3f}   A = {np.exp(a):.4f}   chi2 = {chi2:.2f}/{len(pts)-2}")
    for l, v, ee in pts:
        P(f"      L={l:2d}  chi~ = {v:+.4f} +- {ee:.4f}")
    return b

p25 = [(4, res["L4"]["chit"], res["L4"]["se"]), (6, res["L6"]["chit"], res["L6"]["se"])]
fit(p25, "beta=2.5 only (L=4,6)")
p_all = p25 + [(8, res["L8"]["chit"], res["L8"]["se"])]
fit(p_all, "L=4,6 (b2.5) + L=8 (b2.4)")

# each fixed exponent, amplitude from the L=4 point: predicted vs measured
P("")
P("  fixed-exponent test, amplitude anchored at L=4 (chi~ = A L^p, A = chi~_4/4^p):")
for p in (-2.0, 0.0, 2.0):
    A = res["L4"]["chit"] / 4.0 ** p
    line = f"    p={p:+.0f}: A={A:+.5f}"
    for L, key in ((4, "L4"), (6, "L6"), (8, "L8")):
        pred = A * L ** p
        obs, se = res[key]["chit"], res[key]["se"]
        line += f" | L={L}: pred {pred:+.4f} vs obs {obs:+.4f}+-{se:.4f} ({abs(pred-obs)/se:.1f}s)"
    P(line)

P("")
P("[5] weighted mean & consistency with 0")
jk = [(res[k]["chit"], res[k]["se"]) for k in ("L4", "L6", "L8")]
w = np.array([1 / e ** 2 for _, e in jk]); y = np.array([v for v, _ in jk])
mu, smu = float((w * y).sum() / w.sum()), float(1 / np.sqrt(w.sum()))
chi2 = float((w * (y - mu) ** 2).sum())
P(f"  chi~(k_min) weighted mean = {mu:+.4f} +- {smu:.4f}   "
  f"chi2/dof(vs 0) = {float((w*y**2).sum()):.2f}/3   chi2/dof(vs const) = {chi2:.2f}/2")

P("")
P("[6] KL-atom bound: if chi~(k_min) = f^2/khat^2 (light atom, residue f^2),")
P("    f^2 = chi~ * khat^2,  khat^2 = 4 sin^2(pi/L) = 2 / 1 / 0.5858 at L=4/6/8")
for L, key in ((4, "L4"), (6, "L6"), (8, "L8")):
    kh2 = 4 * np.sin(np.pi / L) ** 2
    f2, sf2 = res[key]["chit"] * kh2, res[key]["se"] * kh2
    P(f"    L={L}: f^2 = {f2:+.4f} +- {sf2:.4f}    (2sigma upper: {f2 + 2*sf2:.4f})")

P("")
P("[7] free-massless log-carrier row: chi = 1/(4sin^2(pi/L)) => chi~ = chi - C(0)")
for L, key in ((4, "L4"), (6, "L6"), (8, "L8")):
    ch = 1.0 / (4 * np.sin(np.pi / L) ** 2)
    pred = ch - C0
    obs, se = res[key]["chit"], res[key]["se"]
    P(f"    L={L}: free chi={ch:.4f} chi~_free={pred:+.4f} | obs {obs:+.4f}+-{se:.4f}"
      f" -> excluded at {abs(pred-obs)/se:.0f} sigma")
P("")
P("[8] noise floor (se at L=4,6,8): "
  f"{res['L4']['se']:.4f}, {res['L6']['se']:.4f}, {res['L8']['se']:.4f}")

P("")
P("[9] k=0 projector at beta=2.5 (M1/GY3 stored JSON), chi~(0) = chi_flat - C(0)")
for tag, path in (("M1", r"research/v3/w7/m1_4d_run_results.json"),
                  ("GY3", r"research/v3/w7/gy3_4d_results.json")):
    d = json.load(open(os.path.join("D:/Dev/math", path)))
    rows = []
    for key, v in d.items():
        if not isinstance(v, dict) or "L" not in v or abs(v.get("beta", -1) - 2.5) > 1e-9:
            continue
        rows.append((int(v["L"]), int(v["N"]), float(v["chi_flat"]), float(v["chi_flat_se"])))
    for L, N, ch, se in sorted(rows):
        P(f"    {tag} L={L:2d} N={N:5d}: chi(k=0)={ch:.4f}+-{se:.4f}  "
          f"chi~(0)={ch - C0:+.4f}+-{se:.4f}  ({abs(ch-C0)/se:.2f}s from 0)")

with open(os.path.join(HERE, "ht4_chitilde.log"), "w", encoding="utf-8") as fh:
    fh.write("\n".join(OUT) + "\n")
print("\n# log -> ht4_chitilde.log")
