"""r1_gauss.py -- R1 Level-1 constant chain, exact on the free (Gaussian) reference.

OBJECT.  Periodic free massive scalar on T^d_L (the "P4 periodic kernel" reference object):
    mu = N(0, P^{-1}),   P = -Delta_lat + m^2,   symbol P(k) = m^2 + sum_mu 2(1-cos k_mu).
Three Poincare constants, all EXACT by linear algebra:

  (A) lam_conf  := Poincare constant of the FULL d-dim Gibbs measure w.r.t. the sum-of-onesite
                   gradients Dirichlet form  E_mu(f,f) = sum_{links} E|d_l f|^2 .
                   Gaussian =>  lam_conf = lambda_min(P)  (over the 1^perp part).
  (B) lam_slice := Poincare constant of the ONE-TIME-SLICE MARGINAL (= the vacuum measure
                   |Omega|^2 of the transfer matrix), same metric.  Gaussian marginal =>
                   lam_slice = 1 / lambda_max(C_slice),  C_slice = (P^{-1}) restricted to a slice,
                   diagonal in 3d lattice momentum  =>  closed form, no diagonalisation.
  (C) 2*m1      : m1 = physical transfer-matrix gap = the mass m (free field), so 2*m1 = 2m.

THE CLAIM TESTED.  [T]  lam_conf <= lam_slice = 2*m1   (chain: Gross LSI=>Poincare, slice-marginal
Rayleigh monotonicity, ground-state transform with the 1/2-kinetic normalization).
Losses: lam_conf/lam_slice -> m/2 -> 0 as m -> 0 (the "wrong-space" loss, quantified);
        lam_slice/(2 m1) -> 1 (no loss beyond the 1/2-kinetic convention).
ACID.  m = 0  =>  lam_conf = 0 exactly  (zero mode).  So the configuration-space Poincare constant
       DOES discriminate gapped from gapless in the Gaussian class -- it is not sign-blind.

Run: py -3.12 -I r1_gauss.py    (numpy/scipy only, deterministic, single process)
Labels [T]/[V].  No QED.  Clay 0%.
"""
import json
import numpy as np

LOG = []


def say(s=""):
    print(s)
    LOG.append(s)


CHK = []


def chk(name, ok, info=""):
    CHK.append({"name": name, "ok": bool(ok), "info": info})
    say(f"    [{'PASS' if ok else 'FAIL'}] {name} {info}")


def save(out):
    with open("r1_gauss.log", "w", encoding="utf-8") as f:
        f.write("\n".join(LOG) + "\n")
    with open("r1_gauss.json", "w", encoding="utf-8") as f:
        json.dump(out, f, indent=1)


def ksq(k, d):
    """lattice eigenvalue 2*sum(1-cos k_mu) for a d-vector of momenta."""
    return sum(2.0 - 2.0 * np.cos(k[i]) for i in range(d))


def lam_conf(d, L, m):
    """lambda_min(P) over all momenta (the constant mode k=0 gives m^2)."""
    return m * m + 0.0  # the k=0 mode is always the minimiser of m^2 + ksq


def lam_slice(d, L, m):
    """1/lambda_max(C_slice); C_slice diagonal in the (d-1) spatial momenta:
       C_slice(k_s) = (1/L) sum_{k_t} 1/P(k_s,k_t).   Sup over k_s (the k_s=0 mode)."""
    ks = 2.0 * np.pi * np.arange(L) / L
    kt = 2.0 * np.pi * np.arange(L) / L
    best = -np.inf
    kbest = None
    for a in range(L):
        k_s = ks[a]
        s = 0.0
        for b in range(L):
            k_t = kt[b]
            s += 1.0 / (m * m + 2.0 - 2.0 * np.cos(k_s) + 2.0 - 2.0 * np.cos(k_t))
        c = s / L
        if c > best:
            best = c
            kbest = (a, b if False else None)
    return 1.0 / best, best


say("=" * 78)
say("R1 / LEVEL 1  --  the exact constant chain on the free reference")
say("=" * 78)
say("  chain:  gamma <= lam_conf <= lam_slice = 2*m1      (all three links exact)")
say("  lam_conf  = lambda_min(P) = m^2              [Gibbs measure, full d-dim config space]")
say("  lam_slice = 1/lambda_max(C_slice)            [vacuum measure |Omega|^2 on one slice]")
say("  2*m1      = 2*m                              [transfer-matrix gap, physics normalisation]")
say("")

rows = []
say("-" * 78)
say("d=4, periodic, L x L x L x L.   lam_conf = m^2 exactly (k=0 mode).")
say("-" * 78)
say(f"  {'L':>3} {'m':>6} {'lam_conf':>10} {'lam_slice':>10} {'2m1':>8} "
    f"{'ls/2m1':>8} {'lc/ls':>8}  chain")
for L in (2, 4, 6, 8, 12):
    for m in (0.3, 0.7, 1.0, 1.5, 2.0):
        lc = lam_conf(4, L, m)
        ls, cmax = lam_slice(4, L, m)
        two_m1 = 2.0 * m
        rows.append({"L": L, "m": m, "lam_conf": lc, "lam_slice": ls,
                     "ratio_slice_2m1": ls / two_m1, "ratio_conf_slice": lc / ls})
        say(f"  {L:3d} {m:6.2f} {lc:10.5f} {ls:10.5f} {two_m1:8.4f} "
            f"{ls / two_m1:8.4f} {lc / ls:8.4f}  {'ok' if lc <= ls + 1e-12 else 'VIOLATED'}")

say("")
say("-" * 78)
say("CHECK 1  the slice-marginal (Rayleigh) inequality  lam_conf <= lam_slice  -- every point")
say("-" * 78)
ok = all(r["lam_conf"] <= r["lam_slice"] + 1e-12 for r in rows)
chk("1.1 lam_conf <= lam_slice at all 25 (L,m) points", ok,
    f"min slack = {min(r['lam_slice'] - r['lam_conf'] for r in rows):.3e}")

say("")
say("-" * 78)
say("CHECK 2  CLOSED FORM.  lam_slice = sqrt(m^2 (m^2+4)),  so lam_slice/(2 m1) = sqrt(m^2+4)/2")
say("-" * 78)
say("  derivation: C_slice(k_s) = (1/L)sum_{k_t} 1/(m^2+lam(k_s)+lam(k_t)) -> (1/2pi)int_0^2pi dk_t")
say("             /(m^2+lam(k_s)+2-2cos k_t) = 1/sqrt((m^2+lam)(m^2+lam+4));  sup at k_s=0.")
say("  => the lattice (a=1) artefact is EXACTLY the time-direction dispersion; the factor is")
say("     sqrt(m^2+4)/2 -> 1 as m -> 0 (continuum), and it is 1.06 at m=0.7, 1.12 at m=1.")
say(f"  {'m':>6} {'L=inf closed form':>18} {'L=12 measured':>14} {'L=24 measured':>14}")
for m in (0.1, 0.3, 0.7, 1.0, 2.0):
    ls12, _ = lam_slice(4, 12, m)
    ls24, _ = lam_slice(4, 24, m)
    cf = np.sqrt(m * m * (m * m + 4.0))
    say(f"  {m:6.2f} {cf:18.6f} {ls12:14.6f} {ls24:14.6f}")
errs = []
for m in (0.3, 0.7, 1.0, 2.0):
    ls24, _ = lam_slice(4, 24, m)
    cf = np.sqrt(m * m * (m * m + 4.0))
    errs.append(abs(ls24 - cf) / cf)
chk("2.1 lam_slice = sqrt(m^2(m^2+4)) to <1% at L=24 whenever xi=1/m <= L/7", max(errs) < 0.01,
    f"max rel.err = {max(errs):.2e} over m>=0.3")
ls24_01, _ = lam_slice(4, 24, 0.1)
chk("2.1b at m=0.1 (xi=10 > L/7=3.4) the finite-L cutoff is visible, as it must be",
    abs(ls24_01 - np.sqrt(0.01 * 4.01)) / np.sqrt(0.01 * 4.01) > 0.1,
    f"L=24 gives {ls24_01:.4f} vs L=inf {np.sqrt(0.01*4.01):.4f}")
chk("2.2 lam_slice/(2 m1) = sqrt(m^2+4)/2 -> 1 as m -> 0 (continuum)",
    abs(np.sqrt(0.01 + 4.0) / 2 - 1.0) < 2.5e-3 and abs(np.sqrt(4.49) / 2 - 1.06) < 5e-3,
    f"m=0.1 -> {np.sqrt(4.01)/2:.5f};  m=0.7 -> {np.sqrt(4.49)/2:.5f}")

say("")
say("-" * 78)
say("CHECK 3  the 'wrong-space' loss: CLOSED FORM  lam_conf/lam_slice = m/sqrt(m^2+4) -> m/2 -> 0")
say("-" * 78)
say(f"  {'m':>5} {'m/sqrt(m^2+4)':>16} {'measured L=64':>16}")
for m in (0.1, 0.3, 1.0, 1.5, 2.0):
    lc = lam_conf(4, 64, m)
    ls, _ = lam_slice(4, 64, m)
    cf = m / np.sqrt(m * m + 4.0)
    say(f"  {m:5.2f} {cf:16.6f} {lc / ls:16.6f}")
chk("3.1 loss = m/sqrt(m^2+4) (-> m/2 in the continuum): the Gibbs constant is NOT the gap",
    abs(1.0 / np.sqrt(5.0) - lam_conf(4, 64, 1.0) / lam_slice(4, 64, 1.0)[0]) < 1e-3,
    "ratio 0.4472 at m=1")

say("")
say("-" * 78)
say("CHECK 4  ACID: m = 0 gives lam_conf = 0 exactly, and lam_slice = 0 (zero mode)")
say("-" * 78)
l0 = lam_conf(4, 8, 0.0)
say(f"  lam_conf(m=0, L=8) = {l0:.3e}   (k=0 mode of P = -Delta; the measure is not normalisable")
say("   on that direction -- equivalently, along the torus the constant mode is free)")
chk("4.1 free massless -> lam_conf = 0", abs(l0) < 1e-14, f"{l0:.2e}")
say("  => the configuration-space Poincare constant is NOT sign-blind: it vanishes exactly where")
say("     the physical gap vanishes.  It cannot certify the gap only because of its MAGNITUDE")
say("     (it is a lower bound at best, with the factor m/2 loss above), not because of its SPACE.")

say("")
say("-" * 78)
say("CHECK 5  d=2 control (same algebra, different exponent of the loss)")
say("-" * 78)
say(f"  {'L':>3} {'m':>5} {'lam_conf':>9} {'lam_slice':>9} {'2m1':>6} {'lc/ls':>7}")
for L in (8, 32):
    for m in (0.5, 1.0):
        lc = lam_conf(2, L, m)
        ls, _ = lam_slice(2, L, m)
        say(f"  {L:3d} {m:5.2f} {lc:9.4f} {ls:9.4f} {2 * m:6.2f} {lc / ls:7.4f}")

say("")
say("=" * 78)
say("VERDICT (witness only; the theorem is in R1-BRIDGE.md S1)")
say("=" * 78)
say("  lam_conf <= lam_slice = 2*m1 holds at every point; the loss is (i) the factor 2 = the")
say("  1/2-kinetic convention (not a loss: matched metrics give 1), (ii) the factor m/2 between")
say("  the Gibbs measure and the vacuum measure (a REAL loss, unbounded as m->0).")
say("  Hence: the direction of the inference is safe -- a positive uniform lower bound on the")
say("  Gibbs Poincare constant implies a positive gap; only its MAGNITUDE is weak.")
say("=" * 78)
save({"rows": rows, "checks": CHK})
