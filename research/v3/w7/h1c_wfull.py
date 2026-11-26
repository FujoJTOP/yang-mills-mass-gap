#!/usr/bin/env python3
r"""h1c_wfull.py -- H1C: does W_full have a quantitative lower bound?  (py -3.12, numpy only)

DEFINITIONS (copied from docs/paper-v2, not from memory):
  ch05-theorems.tex eq:wgeo / eq:wgeo-quant:
      W_geo(v;a) := E[Phi''|a] - Var(Phi'|a)  >=  lam(v) >= (N/12)||v||^2      (flat chart,
      Haar dU = J dX, Phi = -log J; SU(2): Phi = phi0(rho) = 2log(rho/2) - 2log sin(rho/2))
  app.tex compendium (= ch05 eq:split, ch08 line 18, ch16 eq:nacounter):
      W_full := W_geo + X,   X = -2 Cov(G, Phi'),   G = d_v S,  Phi' = d_v Phi
      G''(0) = E[H+Phi''] - Var(G+Phi') = Delta_dyn + W_geo + X = Delta_dyn + W_full.
  Hence  W_full >= W_geo  <=>  X >= 0  <=>  Cov(G,Phi'|a) <= 0  =  hypothesis [NA] (app.tex l.800),
  which the corpus lists as an *unproved* hypothesis, not a theorem.

MODELS
 (A) single plaquette SU(2): one free link, B = exp(i theta sigma_3/2), mode v = e3 (||v||=1);
     exact 3D flat-chart quadrature (sm1_exact.moment_decomp_1x1).
 (B) P4 kernel: periodic 4-plaquette constant-mode transfer chain (p4_periodic_transfer.py model),
     v_l = e3/sqrt(L), ||v||=1; exact z-quadrature at a0 = 0 with z_1 = -(z_2+z_3+z_4), and the
     transverse (two flat directions per link) integrated exactly -> per-link conditional moments
     M1(z)=E[p|z], M2(z)=E[p^2|z], Q(z)=E[q|z] with p,q the per-link pieces of Phi', Phi''.
     Full 3D-chart Jacobian is kept, so the N/12 floor of Lemma B is the one being tested.

Output: W_geo, X, W_full, Var(G), and the ratio  R = W_full / ((N/12)||v||^2)  = 6*W_full (N=2, ||v||=1).
"""
import numpy as np
from scipy.special import sici
from sm1_exact import pr, prr, phi0, moment_decomp_1x1, rho_v_1x1, fd_curve

np.seterr(invalid='ignore', divide='ignore')
N, FLOOR = 2, 2 / 12.0          # N/12 * ||v||^2 ,  ||v|| = 1
TWOPI = 2 * np.pi
_GX, _GW = np.polynomial.legendre.leggauss(3000)     # one grid for all calls


def ratio(w):
    return w / FLOOR


# --------------------------------------------------------------- exact 1-link moment audit
def moments_1x1(a0, b, t):
    """Same quadrature as sm1_exact.moment_decomp_1x1 (v=e3), returning every piece needed for
    the exact identity  W_full = E[Phi''] + Var(G) - Var(G+Phi')."""
    xg, wg = _GX, _GW
    v = np.array([0.0, 0.0, 1.0])
    Rm = np.sqrt(max(4 * np.pi ** 2 - a0 * a0, 1e-14))
    rr = 0.5 * Rm * (xg + 1.0); wrr = 0.5 * Rm * wg
    ea = np.array([1.0, 0, 0])
    w = 2 * np.pi * rr * wrr
    xi = a0 * v[None, :] + rr[:, None] * ea[None, :]
    rho = np.linalg.norm(xi, axis=1)
    Sv = b * (1 - np.cos(rho / 2) * np.cos(t / 2)) + phi0(rho) + b * np.sin(rho / 2) * np.sin(t / 2) * \
         (np.dot(xi, np.array([0, 0, 1.0])) / rho)
    wt = np.exp(-Sv) * w
    Z = wt.sum()
    xiv = xi @ v
    Phig = pr(rho / 2) * xiv / rho
    Phigg = prr(rho / 2) * (xiv / rho) ** 2 + pr(rho / 2) * (1.0 - (xiv / rho) ** 2) / rho

    def S_of(a):
        xi2 = a * v[None, :] + rr[:, None] * ea[None, :]
        rho2 = np.linalg.norm(xi2, axis=1)
        return b * (1 - np.cos(rho2 / 2) * np.cos(t / 2)) + b * np.sin(rho2 / 2) * np.sin(t / 2) * \
               (np.dot(xi2, np.array([0, 0, 1.0])) / rho2)
    da = 1e-5
    Gs = (S_of(a0 + da) - S_of(a0 - da)) / (2 * da)
    Hs = (S_of(a0 + da) - 2 * S_of(a0) + S_of(a0 - da)) / (da * da)

    def E(x): return (wt * x).sum() / Z
    def V(x): 
        m = E(x); return E(x * x) - m * m
    Gt = Gs + Phig
    return dict(EG=E(Gs), VarG=V(Gs), VarPhi=V(Phig), EGpp=E(Phigg), VarPhiG=V(Gt),
                CovGP=E(Gs * Phig) - E(Gs) * E(Phig), EH=E(Hs))


# --------------------------------------------------------------- (A) single plaquette, exact 3D
def single_plaquette(betas=(0.0, 1.0, 2.5, 5.7), thetas=(0.0, np.pi / 2, np.pi),
                     a0s=np.linspace(-3.4, 2.6, 19)):
    print("\n" + "=" * 104)
    print("[A] SINGLE PLAQUETTE  SU(2), one free link, v=e3, ||v||=1.  floor (N/12)||v||^2 = %.6f" % FLOOR)
    print("    W_geo=E[Phi'']-Var(Phi'),  X=-2Cov(G,Phi'),  W_full=W_geo+X;  R_Wf = W_full/floor,")
    print("    R_Wg = W_geo/floor.  [NA] <=> X>=0 <=> W_full>=W_geo.")
    print("=" * 104)
    rows = []
    hdr = "  %6s %7s | %9s %9s %9s | %7s %7s | %8s" % (
        "beta", "theta", "W_geo", "X", "W_full", "R_Wg", "R_Wf", "min(R_Wf)")
    print(hdr)
    for b in betas:
        for t in thetas:
            Wg = np.empty(len(a0s)); Xc = np.empty(len(a0s)); Dd = np.empty(len(a0s))
            Id = np.empty(len(a0s)); Vp = np.empty(len(a0s)); Df = np.empty(len(a0s))
            EGpp = np.empty(len(a0s))
            for i, a0 in enumerate(a0s):
                Dd[i], Wg[i], Xc[i] = moment_decomp_1x1(a0, b, t, 'e3')
                m = moments_1x1(a0, b, t)
                Vp[i] = m['VarPhi']; EGpp[i] = m['EGpp']
                Df[i] = m['VarPhiG'] - m['VarG']                 # Var(G+Phi') - Var(G)
                Id[i] = m['EGpp'] + m['VarG'] - m['VarPhiG']     # identity value of W_full
            Wf = Wg + Xc
            aS = a0s[int(np.argmax([rho_v_1x1(a, b, t, 'e3') for a in a0s]))]
            j = int(np.argmin(np.abs(a0s - aS)))
            rows.append(dict(beta=b, theta=t, a0s=a0s, D=Dd, Wg=Wg, X=Xc, Wf=Wf, a0star=aS,
                             Id=Id, Egp=EGpp, Vp=Vp, Df=Df))
            print("  %6.2f %7.3f | %+9.5f %+9.5f %+9.5f | %7.3f %7.3f | %8.3f   (a0*=%.2f: W_geo=%+.4f"
                  " X=%+.4f W_full=%+.4f R=%.3f)  min_a0 X=%+.4f  max_id_err=%.2e" %
                  (b, t, Wg[j], Xc[j], Wf[j], ratio(Wg[j]), ratio(Wf[j]), Wf.min() / FLOOR,
                   aS, Wg[j], Xc[j], Wf[j], ratio(Wf[j]), Xc.min(),
                   np.max(np.abs(Id - Wf))))
    # explicit negativity witnesses: the a0 that minimises X at each (beta, theta)
    print("\n  [A-neg] points with X<0 (=> W_full < W_geo, so the (a) route FAILS there):")
    anyneg = False
    for r in rows:
        if r['X'].min() < 0:
            anyneg = True
            i = int(np.argmin(r['X']))
            print("    beta=%.2f theta=%.3f  a0=%+7.3f  W_geo=%+8.5f  X=%+8.5f  W_full=%+8.5f"
                  "  deficit W_geo-W_full=%+8.5f  R_Wf=%.3f  ratio (W_full/floor) %s 1" %
                  (r['beta'], r['theta'], r['a0s'][i], r['Wg'][i], r['X'][i], r['Wf'][i],
                   -r['X'][i], ratio(r['Wf'][i]), ">=" if r['Wf'][i] >= FLOOR else "<"))
    if not anyneg:
        print("    none on the scanned box.")
    # does W_full ever drop below the floor, and what does the identity route say?
    wmin = min((r['Wf'].min(), r['beta'], r['theta']) for r in rows)
    print("  [A-floor] min over box of W_full = %+.5f (beta=%.2f, theta=%.3f) -> min R = %.3f  %s 1"
          % (wmin[0], wmin[1], wmin[2], ratio(wmin[0]), ">=" if wmin[0] >= FLOOR else "<"))
    print("  [A-ident] W_full = E[Phi''] + Var(G) - Var(G+Phi')   (exact; E[Phi''] >= floor = Lemma B)")
    print("      %6s %7s %8s | %9s %9s %9s | %9s" %
          ("beta", "theta", "a0", "E[Phi'']", "Var(Phi')", "Def", "E[Phi'']-floor"))
    best = None
    for r in rows:
        i = int(np.argmin(r['Wf']))
        d = dict(beta=r['beta'], theta=r['theta'], a0=r['a0s'][i], Egp=r['Egp'][i],
                 Vp=r['Vp'][i], Df=r['Df'][i], Wf=r['Wf'][i])
        if best is None or d['Wf'] < best['Wf']:
            best = d
    for r in rows:
        i = int(np.argmin(np.abs(r['a0s'] - r['a0star'])))
        d = dict(beta=r['beta'], theta=r['theta'], a0=r['a0s'][i], Egp=r['Egp'][i],
                 Vp=r['Vp'][i], Df=r['Df'][i], Wf=r['Wf'][i])
        print("      %6.2f %7.3f %8.3f | %+9.5f %+9.5f %+9.5f | %+9.5f   [a0*]" %
              (d['beta'], d['theta'], d['a0'], d['Egp'], d['Vp'], d['Df'], d['Egp'] - FLOOR))
    print("      tightest point overall (min W_full): beta=%.2f theta=%.3f a0=%+.3f  E[Phi'']=%+.5f"
          " Var(Phi')=%+.5f Def=%+.5f  E[Phi'']-floor=%+.5f  W_full=%+.5f R=%.3f" %
          (best['beta'], best['theta'], best['a0'], best['Egp'], best['Vp'], best['Df'],
           best['Egp'] - FLOOR, best['Wf'], ratio(best['Wf'])))
    return rows


# --------------------------------------------------------------- (B) P4 kernel (periodic 4-chain)
def rho1(z):
    """exact one-link Cartan marginal, pi(Cin(2pi)-Cin|z|)  (un-normalised; ratios are invariant)."""
    return np.pi * (sici(TWOPI)[0] - sici(np.abs(z))[0])


def link_moments(L, nz=900):
    """M1(z)=E[p|z], M2(z)=E[p^2|z], Q(z)=E[q|z] for v_l=e3/sqrt(L), full 3D chart transverse avg.
    transverse disk: 2pi r dr * J(rho),  rho=sqrt(z^2+r^2),  J=exp(-phi0(rho))."""
    z = np.linspace(-TWOPI, TWOPI, nz)
    xg, wg = np.polynomial.legendre.leggauss(600)
    M1 = np.empty(nz); M2 = np.empty(nz); Q = np.empty(nz)
    for i, zz in enumerate(z):
        R = np.sqrt(max(4 * np.pi ** 2 - zz * zz, 1e-14))
        r = 0.5 * R * (xg + 1.0); wr = 0.5 * R * wg
        rho = np.sqrt(zz * zz + r * r)
        J = np.exp(-phi0(rho))
        w = 2 * np.pi * r * wr * J                       # transverse measure
        wv = zz / np.sqrt(L)                             # xi . v_l
        p = pr(rho / 2) * wv / rho
        q = prr(rho / 2) * wv ** 2 / rho ** 2 + pr(rho / 2) * (1.0 / L - wv ** 2 / rho ** 2) / rho
        Z = w.sum()
        if Z <= 0 or not np.isfinite(Z):          # endpoints |z|=2pi: R=0, no measure
            M1[i] = M2[i] = Q[i] = 0.0
            continue
        M1[i] = (w * p).sum() / Z
        M2[i] = (w * p * p).sum() / Z
        Q[i] = (w * q).sum() / Z
    return z, M1, M2, Q


def p4_core(L, b, n, zt=None, tab=None):
    """all constrained moments on the periodic chain at a0=0; returns (W_geo, X, Var Phi', Var G)."""
    if tab is None:
        tab = link_moments(L)
    zt_, M1t, M2t, Qt = tab
    g = np.linspace(-TWOPI, TWOPI, n)
    h = g[1] - g[0]
    Z2, Z3, Z4 = np.meshgrid(g, g, g, indexing='ij')
    zs = [-Z2 - Z3 - Z4, Z2, Z3, Z4]                      # z_1 = -sum_{2..4} for a0=0
    inside = np.abs(zs[0]) < TWOPI
    S = np.zeros_like(Z2)
    for l in range(L):
        S += b * (1.0 - np.cos(0.5 * (zs[(l + 1) % L] - zs[l])))
    W = np.exp(-S)
    for zz in zs:
        W = W * rho1(zz)
    W = np.where(inside, W, 0.0)
    Zt = W.sum() * h ** 3
    def E(x): return (W * x).sum() * h ** 3 / Zt
    M1 = [np.interp(zz, zt_, M1t) for zz in zs]
    M2 = [np.interp(zz, zt_, M2t) for zz in zs]
    Q = [np.interp(zz, zt_, Qt) for zz in zs]
    G = np.zeros_like(Z2)
    for l in range(L):
        dl = 0.5 * (zs[l] - zs[(l - 1) % L])
        dl1 = 0.5 * (zs[(l + 1) % L] - zs[l])
        G += 0.5 * b * (np.sin(dl) - np.sin(dl1))        # telescoping: G == 0 identically
    G /= np.sqrt(L)
    EPg = sum(E(m) for m in M1)
    EPg2 = sum(E(m2) for m2 in M2)
    for l in range(L):
        for k in range(L):
            if l != k:
                EPg2 += E(M1[l] * M1[k])
    EPgg = sum(E(q) for q in Q)
    EG = E(G); VarG = E(G * G) - EG * EG
    CovGP = sum(E(G * m) - EG * E(m) for m in M1)
    Vg = EPg2 - EPg * EPg
    Wg = EPgg - Vg
    return Wg, -2 * CovGP, Vg, VarG, EPgg


def p4_kernel(L=4, betas=(0.0, 1.0, 2.5, 5.7), n=120, nchk=60):
    print("\n" + "=" * 104)
    print("[B] P4 KERNEL: periodic %d-plaquette constant-mode chain (p4 period), a0=0, ||v||=1." % L)
    print("    S=beta*sum(1-cos((z_{l+1}-z_l)/2)), measure prod rho1(z_l), z_1=-(z_2+..+z_L).")
    print("    H=0 identically (translation flat) AND G=d_a S == 0 identically (telescoping sines)")
    print("    => W_full = W_geo exactly on this kernel; G'' = W_geo - Var(G) = W_geo.")
    print("=" * 104)
    tab = link_moments(L)
    out = []
    for b in betas:
        Wg, Xc, Vg, VarG, EPgg = p4_core(L, b, n, tab=tab)
        Wf = Wg + Xc
        Gpp = -VarG + Wg + Xc
        out.append(dict(L=L, beta=b, Wgeo=Wg, X=Xc, Wfull=Wf, VarPhi=Vg, VarG=VarG,
                        EH=0.0, Gpp=Gpp, R=ratio(Wf), Rg=ratio(Wg), Egp=EPgg))
        print("  beta=%5.2f | W_geo=%+9.5f  X=%+9.5f  W_full=%+9.5f | Var(Phi')=%8.5f  Var(G)=%8.2e"
              " | G''=%+9.5f | R_Wg=%.3f R_Wf=%.3f | E[Phi'']=%+.5f" %
              (b, Wg, Xc, Wf, Vg, VarG, Gpp, ratio(Wg), ratio(Wf), EPgg))
    conv = [abs(p4_core(L, b, nchk, tab=tab)[0] / p4_core(L, b, n, tab=tab)[0] - 1.0)
            for b in (0.0, 2.5)]
    print("  [B-conv] W_geo(n=%d) vs W_geo(n=%d): rel.diff(beta=0, 2.5) = %.2e %.2e"
          % (n, nchk, conv[0], conv[1]))
    return out


def main():
    print("=" * 104)
    print("H1C  W_full vs W_geo : is W_full >= (N/12)||v||^2 free?")
    print("W_full = W_geo + X,  X = -2Cov(G,Phi')  [app.tex l.796; ch05 eq:split]")
    print("W_full >= W_geo  <=>  X >= 0  <=>  Cov(G,Phi'|a) <= 0  <=>  hypothesis [NA] (app.tex l.800)")
    print("=" * 104)
    ra = single_plaquette()
    rb = p4_kernel()
    tot = np.concatenate([np.asarray(r['Wf']).ravel() for r in ra] +
                         [np.array([r['Wfull']]) for r in rb])
    xs = np.concatenate([np.asarray(r['X']).ravel() for r in ra])
    print("\n" + "=" * 104)
    print("VERDICT DATA")
    print("  min W_full over all measured points : %+.5f   (floor %.5f, R=%.3f)" %
          (tot.min(), FLOOR, ratio(tot.min())))
    print("  min (W_full - W_geo) = X over all     : %+.5f   -> negative somewhere: %s" %
          (xs.min(), bool((xs < 0).any())))
    print("  fraction of (beta,theta,a0) points with W_full < W_geo: %.3f" % np.mean(xs < 0))
    print("  fraction of all points with W_full < (N/12)||v||^2    : %.3f" % np.mean(tot < FLOOR))
    print("  min R = W_full/floor over all points  : %.3f" % ratio(tot.min()))
    print("=" * 104)


if __name__ == "__main__":
    main()
