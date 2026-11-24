"""
wi_modified_functional.py -- WI/IR flood: a MODIFIED analytic functional robust to a
background continuum (witness for research/unconditional/WI-MODIFIED-FUNCTIONAL.md).

THE FAILURE BEING FIXED
-----------------------
The correlator is Z(tau) = w0 + sum_k c_k e^{-E_k tau}, tau complex, w0 = |<O>|^2 > 0.
The RAW functional  A = pi/|Im tau|  of the NEAREST zero fails once a continuum is
present: the continuum creates spurious zeros at SMALL |Im|, so the nearest zero is
NOT the atom's (WI-ADMISSIBLE, T2).

THE CANDIDATE ROBUST FUNCTIONAL: the REWEIGHTING-PERSISTENT RUNG
-----------------------------------------------------------------
An isolated ATOM (one level E1 + VEV) puts zeros on ONE vertical line
    tau_n = (1/E1)[ ln(c1/w0) - i(2n+1)pi ],   |Im tau_n| = (2n+1)pi/E1.
Under operator reweighting only the ATOM WEIGHT c1 changes (levels E_k are fixed):
  * |Im tau_n| -- the RUNG POSITION -- is INDEPENDENT of c1;
  * Re tau_n = ln(c1/w0)/E1 -- MOVES LOG-LINEARLY in c1.
A CONTINUUM has no such c1: reweighting it (changing the density, support fixed) moves
the zeros in BOTH Re and Im, and generically leaves no c1-stable |Im| rung.

So the modified functional is: RUNG(atom) := the |Im| value that is INVARIANT under a
family of operator reweightings, together with the check that Re tau along it scales as
ln(c1).  It is non-local (a property of the whole zero ladder + a reweighting family),
unlike the single nearest zero.

MODELS (w0=1, atom E1=1):
  ATOM      : c1 e^{-tau}                              -- reference (isolated)
  ATOM+GAP  : atom + small continuum with edge mc>E1    -- PHYSICAL atom (isolated by a gap)
  ATOM+THR  : atom + small continuum with edge mc=E1    -- atom AT threshold (not isolated)
  CUT       : continuum only (no atom)                  -- pure cut (no c1 to reweight)

TEST: for each model, reweight (c1 in {0.4, 1.0, 2.5} for atoms; density scale for the cut),
find the |Im| rungs present in ALL reweightings, and report:
  * whether a rung at ~pi/E1 persists (ATOM signal),
  * the slope of Re tau vs ln(c1) along it (ATOM -> slope = 1/E1; CUT -> no such rung).

Run:  py -3.12 -I wi_modified_functional.py
"""
import numpy as np

np.set_printoptions(precision=5, suppress=True, linewidth=120)


def make_Z(E, c, w0):
    E = np.asarray(E, float); c = np.asarray(c, float)
    def Z(t):   return w0 + np.sum(c * np.exp(-E * t))
    def dZ(t):  return np.sum(-E * c * np.exp(-E * t))
    return Z, dZ


def zeros_newton(Z, dZ, seeds, tol=1e-11, maxit=300):
    roots = []
    for s in seeds:
        t = complex(s); ok = True
        for _ in range(maxit):
            f = Z(t); fp = dZ(t)
            if not np.isfinite(f) or not np.isfinite(fp) or abs(fp) < 1e-300:
                ok = False; break
            step = f / fp
            if not np.isfinite(step): ok = False; break
            t = t - step
            if abs(step) < tol: break
        else:
            ok = False
        if ok and abs(t.real) < 30 and 0.05 < abs(t.imag) < 40 and abs(Z(t)) < 1e-6:
            roots.append(t)
    # dedupe by (Re, |Im|) proximity
    uniq = []
    for r in sorted(roots, key=lambda z: abs(z.imag)):
        if all(abs(r.real - u.real) > 1e-3 or abs(abs(r.imag) - abs(u.imag)) > 1e-3
               for u in uniq):
            uniq.append(r)
    return uniq


def seeds_for(E1=1.0, w0r=1.0, kmax=16):
    out = []
    for k in range(kmax):
        out.append((np.log(w0r) - 1j * (2 * k + 1) * np.pi) / E1)
    for re in np.linspace(-3.0, 4.0, 15):
        for im in np.linspace(-0.2, 30, 70):
            out.append(re - 1j * im)
    return out


def build(kind, reweight, w0=1.0, E1=1.0, mc=2.0, eps=0.01, Nc=48, dE=0.05, decay=1.0):
    """reweight = atom weight c1 (ATOM/ATOM+GAP/ATOM+THR) or density scale (CUT)."""
    if kind == "CUT":
        Ec = mc + (np.arange(Nc) + 0.5) * dE
        wc = reweight * dE * np.exp(-(Ec - mc) / decay)
        return make_Z(list(Ec), list(wc), w0)
    E1c = E1
    if kind == "ATOM":
        E = [E1c]; c = [reweight]
    elif kind == "ATOM+GAP":
        Ec = mc + (np.arange(Nc) + 0.5) * dE
        E = [E1c] + list(Ec); c = [reweight] + list(eps * dE * np.exp(-(Ec - mc) / decay))
    elif kind == "ATOM+THR":
        Ec = E1 + (np.arange(Nc) + 0.5) * dE
        E = [E1c] + list(Ec); c = [reweight] + list(eps * dE * np.exp(-(Ec - E1) / decay))
    else:
        raise ValueError(kind)
    return make_Z(E, c, w0)


def zeros_at(kind, rw, **kw):
    Z, dZ = build(kind, rw, **kw)
    return zeros_newton(Z, dZ, seeds_for())


def persistent_rungs(models, tol_im=0.06):
    """|Im| rungs present (within tol_im) in EVERY reweighting."""
    rungs = []
    for z0 in models[0]:
        im0 = abs(z0.imag)
        ok = True
        for zp in models[1:]:
            if not any(abs(abs(z.imag) - im0) < tol_im for z in zp):
                ok = False; break
        if ok:
            rungs.append(im0)
    rungs.sort()
    # merge nearby
    out = []
    for r in rungs:
        if not out or r - out[-1] > tol_im:
            out.append(r)
    return out


def re_at_rung(model, rung, tol_im=0.06):
    """mean Re tau of roots whose |Im| is within tol of the rung (None if absent)."""
    rs = [z.real for z in model if abs(abs(z.imag) - rung) < tol_im]
    return np.mean(rs) if rs else np.nan


def banner(s):
    print("\n" + "=" * 80); print(s); print("=" * 80)


if __name__ == "__main__":
    w0 = 1.0; E1 = 1.0
    ATOMWS = [0.4, 1.0, 2.5]
    CUTWS  = [0.4, 1.0, 2.5]      # density scale for the cut

    banner("RAW nearest-zero functional (the failing one) -- pi/|Im| of smallest-|Im| root")
    print(f"{'model':10s} " + "  ".join(f"rw={w:<4}" for w in ATOMWS))
    for kind, ws in [("ATOM", ATOMWS), ("ATOM+GAP", ATOMWS),
                     ("ATOM+THR", ATOMWS), ("CUT", CUTWS)]:
        vals = []
        for rw in ws:
            zs = zeros_at(kind, rw, mc=(1.0 if kind == "ATOM+THR" else 2.0))
            vals.append(np.pi / abs(zs[0].imag) if zs else np.nan)
        print(f"{kind:10s} " + "  ".join(f"{v:7.4f}" for v in vals))

    banner("MODIFIED functional: reweighting-PERSISTENT rungs (|Im| stable across c1)")
    for kind, ws in [("ATOM", ATOMWS), ("ATOM+GAP", ATOMWS),
                     ("ATOM+THR", ATOMWS), ("CUT", CUTWS)]:
        models = [zeros_at(kind, rw, mc=(1.0 if kind == "ATOM+THR" else 2.0)) for rw in ws]
        rungs = persistent_rungs(models)
        print(f"\n{kind}:  persistent rungs |Im| = {[f'{r:.3f}' for r in rungs]}")
        for rung in rungs:
            rezs = [re_at_rung(m, rung) for m in models]
            if any(np.isnan(x) for x in rezs):
                continue
            # slope of Re vs ln(c1) fit
            xs = np.log(np.array(ws)); ys = np.array(rezs)
            slope = np.polyfit(xs, ys, 1)[0] if len(set(xs)) > 1 else np.nan
            print(f"   rung |Im|={rung:6.3f}  pi/|Im|={np.pi/rung:6.4f}  "
                  f"Re(c1)={['%.3f' % r for r in rezs]}  dRe/dln(c1)={slope:6.3f} "
                  f"(atom: 1/E1=1.000)")

    banner("VERDICT table: does a rung at pi/E1 (=3.1416) persist, and is its Re slope ~1/E1?")
    print(f"{'model':10s} {'rung@~pi':>10s} {'pi/rung':>8s} {'dRe/dln(c1)':>12s} {'verdict':>18s}")
    for kind, ws in [("ATOM", ATOMWS), ("ATOM+GAP", ATOMWS),
                     ("ATOM+THR", ATOMWS), ("CUT", CUTWS)]:
        models = [zeros_at(kind, rw, mc=(1.0 if kind == "ATOM+THR" else 2.0)) for rw in ws]
        rungs = persistent_rungs(models)
        near_pi = min(rungs, key=lambda r: abs(r - np.pi)) if rungs else np.nan
        if not np.isnan(near_pi) and abs(near_pi - np.pi) < 0.15:
            rezs = [re_at_rung(m, near_pi) for m in models]
            if not any(np.isnan(x) for x in rezs):
                xs = np.log(np.array(ws)); slope = np.polyfit(xs, np.array(rezs), 1)[0]
                pi_over = np.pi / near_pi
                verdict = ("ATOM (robust)" if abs(slope - 1.0) < 0.15 and abs(pi_over - 1) < 0.02
                           else "atom-ish")
                print(f"{kind:10s} {near_pi:10.3f} {pi_over:8.4f} {slope:12.3f} {verdict:>18s}")
            else:
                print(f"{kind:10s} rung {near_pi:.3f} but Re missing")
        else:
            print(f"{kind:10s} {'none':>10s} {'--':>8s} {'--':>12s} {'NO atom rung':>18s}")

    print("""
READING:
  * ATOM: rung at |Im|=pi, pi/|Im|=1.000, dRe/dln(c1)=1/E1=1 -> ATOM.
  * ATOM+GAP: if the same rung at |Im|~pi persists (pi/|Im|~1) with the SAME log-slope,
    the modified functional is ROBUST to a GAPPED continuum -- it separates the
    physical atom (isolated by a gap) from a pure cut, unlike the raw nearest zero.
  * ATOM+THR: atom at threshold -> rung smeared, expected failure (not isolated).
  * CUT: no rung at pi/E1 with atom log-slope -> NO atom (correctly).
""")
