"""
wi_reslope_lattice.py -- witness for WI-RESLOPE-LATTICE.md (WI line 2, real-tau lattice protocol).

Compares the two finite-L atom tests on the SAME synthetic correlator family
Z_r(tau) = w0 + sum_k c_k^r e^{-E_k tau},  reweighting = change c_1 only (levels E_k fixed):

  (WI-SU3) RAW:  A_r = pi/|Im of NEAREST zero|;  decision = spread_r A_r  (atom~0, cut large).
  (Re-slope)  :  persistent |Im tau| RUNG with Re tau log-linear in c_1 at slope 1/E1;
                 decision = R = #rungs with slope ~ 1/E1  (atom >=3, cut 0).

Claim tested: for an ATOM riding on a continuum (the physical case) the RAW nearest-zero
spread FAILS (continuum creates spurious small-|Im| zeros -> WI-ADMISSIBLE T2 / WI-EMBEDDED),
while the Re-slope R stays >=3.  For a pure CUT, R=0.

Run: py -3.12 -I wi_reslope_lattice.py
"""
import numpy as np
np.set_printoptions(precision=5, suppress=True, linewidth=120)
np.seterr(over='ignore', invalid='ignore')   # divergent Newton seeds overflow exp(); harmless

def make_Z(E, c, w0):
    E = np.asarray(E, float); c = np.asarray(c, float)
    Z  = lambda t: w0 + np.sum(c * np.exp(-E * t))
    dZ = lambda t: np.sum(-E * c * np.exp(-E * t))
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
    uniq = []
    for r in sorted(roots, key=lambda z: abs(z.imag)):
        if all(abs(r.real - u.real) > 1e-3 or abs(abs(r.imag)-abs(u.imag)) > 1e-3 for u in uniq):
            uniq.append(r)
    return uniq

def seeds_ladder(E1, w0, c1, kmax=14, extra=True):
    out = [(np.log(c1/w0) - 1j*(2*k+1)*np.pi)/E1 for k in range(kmax)]
    if extra:  # catch spurious continuum zeros at small |Im|
        for re in np.linspace(-2, 2, 5):
            for im in (0.3, 0.6, 1.0, 1.6, 2.4):
                out.append(re - 1j*im)
    return out

def spectral(E1, contE, contC, w0, c1):  # levels & weights for a given atom weight c1
    E = np.concatenate([[E1], contE]); c = np.concatenate([[c1], contC])
    return E, c

def raw_test(E1, contE, contC, w0, c1s):
    A = []
    for c1 in c1s:
        E, c = spectral(E1, contE, contC, w0, c1)
        Z, dZ = make_Z(E, c, w0)
        zs = zeros_newton(Z, dZ, seeds_ladder(E1, w0, c1))
        if not zs: A.append(np.inf); continue
        A.append(np.pi/min(abs(z.imag) for z in zs))     # nearest-zero functional
    return np.array(A)

def re_slope_test(E1, contE, contC, w0, c1s, tol_rung=0.35, tol_slope=0.12):
    """R = number of rungs (2n+1)pi/E1 with Re slope ~ 1/E1 across the reweighting family."""
    rung_Re = {n: [] for n in range(1, 8)}    # rung index n -> list of Re per c1
    for c1 in c1s:
        E, c = spectral(E1, contE, contC, w0, c1)
        Z, dZ = make_Z(E, c, w0)
        zs = zeros_newton(Z, dZ, seeds_ladder(E1, w0, c1))
        ims = np.array([abs(z.imag) for z in zs])
        for n in range(1, 8):
            target = (2*n+1)*np.pi/E1
            if ims.size == 0: rung_Re[n].append(np.nan); continue
            j = int(np.argmin(np.abs(ims - target)))
            if abs(ims[j]-target) < tol_rung*(np.pi/E1):
                rung_Re[n].append(zs[j].real)
            else:
                rung_Re[n].append(np.nan)
    lnc = np.log(np.array(c1s))
    R = 0; slopes = {}
    for n in range(1, 8):
        re = np.array(rung_Re[n])
        if np.any(np.isnan(re)): continue
        s = np.polyfit(lnc, re, 1)[0]
        slopes[n] = s
        if abs(s - 1.0/E1) < tol_slope:   # slope ~ 1/E1
            R += 1
    return R, slopes

def dense_continuum(edge, span, total, p=40):
    xs = np.linspace(edge, edge+span, p); w = np.full(p, total/p); return xs, w

def main():
    w0 = 1.0; E1 = 1.0; c1s = [0.4, 1.0, 2.5]
    print("="*74); print("WI-RESLOPE-LATTICE witness  (w0=1, atom E1=1, reweight c1 in", c1s, ")"); print("="*74)
    cases = {}
    # pure atom
    cases["ATOM      "] = (np.array([]), np.array([]))
    # atom + gapped continuum (edge 2 > E1), 1% total weight
    eE, eC = dense_continuum(2.0, 2.0, 0.01); cases["ATOM+GAP  "] = (eE, eC)
    # atom + AT-THRESHOLD continuum (edge = E1)
    tE, tC = dense_continuum(1.0, 2.0, 0.01); cases["ATOM+THR  "] = (tE, tC)
    # pure cut (no atom): c1 = 0
    cutE, cutC = dense_continuum(1.0, 3.0, 1.0)
    print("\ncase        |  RAW(WI-SU3): A_r per c1 -> spread   |  RE-SLOPE: R atom-slope rungs")
    print("-"*74)
    for name,(eE,cC) in cases.items():
        if name.startswith("ATOM") and "GAP" not in name and "THR" not in name:
            pass
        A = raw_test(E1, eE, cC, w0, c1s)
        spread = np.ptp(A) if np.all(np.isfinite(A)) else np.inf
        R, sl = re_slope_test(E1, eE, cC, w0, c1s)
        print(f"{name}  |  A={A}  spread={spread:8.4f}  |  R={R}   slopes={ {k:round(v,3) for k,v in sl.items()} }")
    # pure cut: no atom weight c1; reweight the continuum density scale instead
    print("-"*74)
    A = []
    for s in (0.5, 1.0, 2.0):
        Z, dZ = make_Z(cutE, s*cutC, 1.0)
        zs = zeros_newton(Z, dZ, seeds_ladder(E1, 1.0, 1.0))
        A.append(np.pi/min(abs(z.imag) for z in zs) if zs else np.inf)
    A = np.array(A)
    # Re-slope for cut: reweight density scale; expect NO persistent atom-slope rung
    rung_Re = {n: [] for n in range(1,8)}; lnc = np.log([0.5,1.0,2.0])
    for s in (0.5,1.0,2.0):
        Z,dZ = make_Z(cutE, s*cutC, 1.0)
        zs = zeros_newton(Z, dZ, seeds_ladder(E1,1.0,1.0)); ims=np.array([abs(z.imag) for z in zs])
        for n in range(1,8):
            tgt=(2*n+1)*np.pi; j=int(np.argmin(np.abs(ims-tgt)))
            rung_Re[n].append(zs[j].real if abs(ims[j]-tgt)<0.35*np.pi else np.nan)
    Rcut=0
    for n in range(1,8):
        re=np.array(rung_Re[n])
        if np.any(np.isnan(re)): continue
        if abs(np.polyfit(lnc,re,1)[0]-1.0)<0.12: Rcut+=1
    print(f"CUT        |  A={A}  spread={np.ptp(A):8.4f}  |  R={Rcut}")
    print("="*74)
    print("READ: ATOM -> R>=3, spread~0.  ATOM+GAP/THR -> RAW spread LARGE (continuum spoils nearest zero)")
    print("      but RE-SLOPE R>=3 (deep rungs survive).  CUT -> R=0. => Re-slope robust to continuum.")

if __name__ == "__main__":
    main()
