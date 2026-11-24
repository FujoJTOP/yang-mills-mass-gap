"""
h2_defect_bound.py  --  H2 witness (deterministic, numpy only, py -3.12 -I)

Bounds the measure-carried defect D at the head of the UV tower.

Exact objects (INVENT-VARIANCE-INEQ / SLOT-SD-FIXPOINT / SAVVIDY-DTERM):
    G''(0) = S - D,     S = E[h_n - <K^{-1}c,c>|0],   D = Var(f|0) - E[<K^{-1}c,c>|0].
Rank-one Savvidy split:  K = K_+ - |lam_-| phi phi^T  (phi rank one, |lam_-|=g|B|,
    phi in ker K_+ by construction),  => <K^{-1}c,c> = <K_+^+ c,c> - |<phi,c>|^2/(g|B|)
    => D = D_+ + D_Sav,   D_+ = Var(f|0)-E<K_+^+ c,c>  (<=0 by BL on log-concave fiber),
                            D_Sav = E[|<phi,c>|^2]/(g|B|)  (>=0).

Checks:
 (1) rank-one split D = D_+ + D_Sav with D_+ <= 0, D_Sav >= 0   [exact algebra]
 (2) Fisher-information identity  E_{a~rho}[G''(a)] = I(rho) = E[(d_a log rho)^2] >= 0
     (the measure-weighted statement); free massless (rho flat) => 0   [acid test]
 (3) rank-one bound  D_Sav <= 2 g^2 / lambda  (exact Var(phi^2) <= 2/lambda), compared to
     S at beta=2.5 -> the bound overshoots (marginal, fails to certify D<S)
 (4) measured per-mode D/S on k0_su2_b2.4_L8_n2.npz: fraction D<S, margin E[G'']
 (5) overlap suppression |<phi,c>|^2 ~ (k/sqrt(gB))^2 deep IR vs O(1) at head
"""
import numpy as np, os

np.set_printoptions(precision=5, suppress=True)
rng = np.random.default_rng(20261124)

print("="*74)
print("H2-DEFECT-BOUND witness : measure-carried defect D at the UV head")
print("="*74)

# ============================================================================
# (1) EXACT RANK-ONE SPLIT
# ============================================================================
print("\n(1) EXACT RANK-ONE SPLIT  D = D_+ + D_Sav   [T]")
n = 6
Q, _ = np.linalg.qr(rng.standard_normal((n, n)))
phi = Q[:, 0]                                       # rank-one negative mode
pos = np.abs(rng.standard_normal(n)) + 0.3
pos[0] = 0.0                                        # phi in ker(K_+): the positive part
Kplus = (Q * pos) @ Q.T; Kplus = 0.5*(Kplus + Kplus.T)
lam_minus = 1.7                                     # |lam_-| = g|B|
K = Kplus - lam_minus * np.outer(phi, phi)          # phi eigenvalue = -lam_minus < 0

c = rng.standard_normal(n)
Kplus_pinv = np.linalg.pinv(Kplus, rcond=1e-12)     # pseudo-inverse (phi -> 0)
Kinv = np.linalg.inv(K)
lhs = c @ Kinv @ c
rhs = c @ Kplus_pinv @ c - (phi @ c)**2 / lam_minus
print(f"    <K^-1 c,c>                     = {lhs:+.6f}")
print(f"    <K_+^+ c,c> - <phi,c>^2/lam_-  = {rhs:+.6f}")
print(f"    split error                    = {abs(lhs-rhs):.2e}   (exact algebra)")

def fiber_Dplus(B, k4, h=0.3, npts=400001):
    y = np.linspace(-20, 20, npts)
    w = np.exp(-(0.5*B*y**2 + 0.25*k4*y**4)); w /= np.trapezoid(w, y)
    f = h*y**2; K_ = B + 3.0*k4*y**2; cc = 2.0*h*y
    Var_f = np.trapezoid(w*f**2, y) - np.trapezoid(w*f, y)**2
    EcK = np.trapezoid(w*cc**2/K_, y)
    return Var_f - EcK
print("    D_+ on log-concave positive fiber (BL => <=0):")
for B in (0.5, 1.0, 2.0):
    print(f"      B={B:>4}:  D_+ = {fiber_Dplus(B, 1.0):+.5f}")
print("    => D_Sav = E|<phi,c>|^2/(g|B|) >= 0 is the only positive part of D   [T]")

# ============================================================================
# (2) FISHER-INFORMATION IDENTITY + ACID TEST  (use analytic log-density)
# ============================================================================
print("\n(2) MEASURE-WEIGHTED:  E_{a~rho}[G''(a)] = I(rho) = E[(d_a log rho)^2] >= 0   [T]")
def fisher(logrho, dlog, a, da):
    e = logrho(a); e = e - e.max()
    rho = np.exp(e); Z = np.trapezoid(rho, a); rho = rho/Z
    d1 = dlog(a)                                  # analytic score
    d2 = np.gradient(d1, da)                      # d_a of the score
    Gpp = -d2
    E_Gpp = np.trapezoid(rho*Gpp, a)
    I = np.trapezoid(rho*d1**2, a)
    return E_Gpp, I
a = np.linspace(-12, 12, 800001); da = a[1]-a[0]
cases = {
 "Gaussian (m=1)       ": (lambda x: -0.5*x**2,              lambda x: -x),
 "Gaussian (m=2)       ": (lambda x: -0.5*(2*x)**2,         lambda x: -4*x),
 "non-Gaussian (well)  ": (lambda x: -0.5*x**2-0.05*x**4,   lambda x: -x-0.2*x**3),
 "FREE massless (flat) ": (lambda x: 0*x,                   lambda x: 0*x),   # acid test
}
for name,(lr,dl) in cases.items():
    E_Gpp, I = fisher(lr, dl, a, da)
    print(f"    {name}:  E[G''] = {E_Gpp:+.6f}   I(rho) = {I:+.6f}")
print("    => typical-field AVERAGE of G'' is the Fisher information of the coarse")
print("       marginal; >=0 automatically, =0 for the free massless field (acid test).")

# ============================================================================
# (3) RANK-ONE BOUND  D_Sav <= 2 g^2 / lambda  vs S  at beta=2.5
# ============================================================================
print("\n(3) RANK-ONE BOUND  D_Sav <= 2 g^2/lambda  (exact Var(phi^2)<=2/lambda)   [T]")
beta = 2.5; g2 = 4.0/beta; S_meas = 1.40
print(f"    beta={beta}, g^2=4/beta={g2:.3f},  S_bare ~ {S_meas}")
print(f"    {'lambda/g^2':>10} {'D_Sav bound=2g^2/lam':>22} {'< S ?':>16}")
for r in (0.5, 1.0, 1.5, 2.0):
    bound = 2*g2/(r*g2)
    print(f"    {r:>10.2f} {bound:>22.3f} {('YES' if bound<S_meas else 'NO (overshoot)'):>16}")
print("    => natural lambda~g^2 gives bound ~2 > S~1.4: too lossy to certify D<S at")
print("       the head, though MEASURED D/S<1.  The bound-vs-measured gap is the defect.")
def var_phi2(c, lam, npts=600001):
    x = np.linspace(-40, 40, npts)
    e = 0.5*c*x**2 - 0.25*lam*x**4; e -= e.max()
    w = np.exp(e); w /= np.trapezoid(w, x)
    m2 = np.trapezoid(w*x**2, x); m4 = np.trapezoid(w*x**4, x)
    return m4 - m2**2
vp = [var_phi2(c,1.0) for c in (0.1,1.0,5.0,100.0)]
print("    exact Var(phi^2), lambda=1, c=g|B| =",
      [f"{v:.4f}" for v in vp], " -> 2/lambda=2 (bounded, c-independent)")

# ============================================================================
# (4) MEASURED per-mode D/S  (ensemble identity per mode)
# ============================================================================
print("\n(4) MEASURED per-mode D/S at beta=2.4, L=8   [V-num]")
path = os.path.normpath(r"D:\Dev\math\research\unconditional\w7\k0_su2_b2.4_L8_n2.npz")
if os.path.exists(path):
    d = np.load(path)
    H,G,P1,P2 = d['H'],d['G'],d['P1'],d['P2']           # (Ncfg, Nmode)
    S    = H.mean(0)                                    # E[H] per mode
    VarG = G.var(0)                                     # Var(G) per mode
    EP2  = P2.mean(0); VarP1 = P1.var(0)
    Cov  = ((G-G.mean(0))*(P1-P1.mean(0))).mean(0)
    W    = EP2 - VarP1 - 2*Cov
    D    = VarG - W
    Gpp  = S - D
    DS   = D/S
    print(f"    Ncfg={H.shape[0]}  Nmode={H.shape[1]}")
    print(f"    D/S : mean={DS.mean():.4f}  min={DS.min():.4f}  max={DS.max():.4f}")
    print(f"    G'' : mean={Gpp.mean():+.4f}  min={Gpp.min():+.4f}  max={Gpp.max():+.4f}")
    print(f"    modes with D<S (G''>0): {int((Gpp>0).sum())}/{len(Gpp)}"
          f"   with D>S: {int((Gpp<0).sum())}/{len(Gpp)}")
    print(f"    typical-field margin  E[G''] = mean(S)-mean(D) = {Gpp.mean():+.4f}"
          f"  = {1-DS.mean():.3f}*mean(S)")
else:
    print("    [npz not found; lattice-cell doc: D/S mean 0.969, 5/12 modes G''<0]")

# ============================================================================
# (5) OVERLAP SUPPRESSION channel
# ============================================================================
print("\n(5) OVERLAP SUPPRESSION  |<phi,c>|^2/(gB) across k/sqrt(gB)  [T]/[O]")
for r in (0.05, 0.1, 0.3, 0.5, 1.0, 2.0):
    ov = r**2/(1.0+r**2)
    tag = "deep IR -> 0" if r < 0.2 else ("HEAD: O(1)" if r >= 0.5 else "")
    print(f"    k/sqrt(gB)={r:>5}:  overlap^2 = {ov:>.4f}   {tag}")
print("    => deep IR: D_Sav/S -> 0 (uniform, floor survives).")
print("       head k~sqrt(gB): D_Sav/S = O(1) -> marginal O(1) competition = H1.")

print("\n" + "="*74)
print("VERDICT: deep IR D<S is uniform [T]; head D<S is an O(1) marginal competition [O].")
print("Typical-field average E[G''] = Fisher info of the coarse marginal [T], =0 free field.")
print("Clay 0%; no QED.")
print("="*74)
