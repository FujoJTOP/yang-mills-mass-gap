"""
E2-FLIP-DYNAMICAL : solvable-model witness for the "Savvidy-Hessian-as-lever" functional.

We test the candidate
    Phi_Hess[mu] := inf_x Hess V(x)          (the bottom of the action Hessian = the Savvidy depth)
against the physical transfer-matrix mass gap m1^2 = spectral gap of the Langevin generator
    L = Delta - grad V . grad      on L^2(mu),  mu = Z^{-1} e^{-V} dx .

Fact used (exact): via U f = f e^{-V/2},  U(-L)U^{-1} = H_s := -Delta + (|grad V|^2/4 - Delta V/2)
on L^2(dx).  The nonzero spectrum of -L  =  spectrum of H_s (E0 = 0).  So
    m1^2 = spectral gap of L  =  E1(H_s) - E0(H_s) = E1(H_s).

Convention: nu lives on m^2; the mass is sqrt(gap).  We report gap = m^2.

Cases:
  (A) free massless   V = 0
  (B) free massive    V = (1/2) m^2 x^2            (Gaussian / OU, exact gap = m^2)
  (C) Savvidy-like    V = (g/4)(x^2 - v^2)^2       (double well: Hess V(0) = -g v^2 < 0,
                                                    exact gap = tunnelling splitting > 0)

We record, per case:  Hessian bottom inf V'' ,  measured gap E1 ,  plus the
measure moments E[H], E[V'^2] (Stein identity E[H]=E[V'^2]), Var[V'], and the
conditioned curvature at the symmetric point  G''(0) = -d^2/dx^2 log(mu)|_0 = V''(0).

Headless, single process.
"""
import numpy as np


def schroedinger_gap(V, dV, ddV, X=6.0, N=2401):
    """Discretise H_s = -d^2 + (V'^2/4 - V''/2) on [-X,X], Dirichlet ends.
    Return (gap=E1-E0, E0, E1, psi0, x)."""
    x = np.linspace(-X, X, N)
    h = x[1] - x[0]
    # second derivative, Dirichlet (psi=0 at endpoints): tridiagonal
    main = -2.0 / h**2 * np.ones(N)
    off = 1.0 / h**2 * np.ones(N - 1)
    # build H_s
    H = np.diag(main) + np.diag(off, 1) + np.diag(off, -1)
    H = -H  # negative of the (negative-definite) Laplacian -> +(-d^2)
    # Potential term (diagonal)
    Wx = dV(x)**2 / 4.0 - ddV(x) / 2.0
    H[np.diag_indices(N)] += Wx
    # Dirichlet: kill endpoints strongly to emulate psi=0
    H[0, :] = 0.0; H[-1, :] = 0.0
    H[0, 0] = 1.0e6; H[-1, -1] = 1.0e6
    w, vec = np.linalg.eigh(H)
    # ground state E0 ~ 0 (the endpoint-pushed spurious states are huge, ignore)
    return w[0], w[1], w[1] - w[0], vec[:, 0], x


def summarize(name, V, dV, ddV, xs_probe=(-6.0, 6.0, 20001)):
    xp = np.linspace(*xs_probe)
    hess_bottom = float(np.min(ddV(xp)))
    E0, E1, gap, psi0, x = schroedinger_gap(V, dV, ddV)
    # measure moments by quadrature on the wavefunction grid
    h = x[1] - x[0]
    p0 = psi0**2
    p0 = p0 / (p0.sum() * h)
    EH = float(np.sum(ddV(x) * p0) * h)          # E_mu[ V'' ]
    EVp2 = float(np.sum(dV(x)**2 * p0) * h)      # E_mu[ V'^2 ]   (Stein: == E[H])
    mVp = float(np.sum(dV(x) * p0) * h)
    VarVp = EVp2 - mVp**2
    Gpp0 = float(ddV(np.array([0.0]))[0])        # V''(0) = G''(0) at symmetric pt
    print(f"\n=== {name} ===")
    print(f"  Hessian bottom  inf V''      = {hess_bottom:+.6f}   <-- this is Phi_Hess")
    print(f"  measured E0, E1(H_s)         = {E0:+.6e}, {E1:+.6e}")
    print(f"  PHYSICAL GAP  m1^2 = E1-E0   = {gap:+.6e}   (mass m1 = {np.sqrt(max(gap,0)):.6f})")
    print(f"  E_mu[H]  (Stein == E_mu[V'^2])= {EH:+.6f}   vs  E_mu[V'^2] = {EVp2:+.6f}")
    print(f"  Var_mu[V']                   = {VarVp:+.6f}")
    print(f"  conditioned curv  G''(0)=V''(0) = {Gpp0:+.6f}")
    return dict(name=name, hess_bottom=hess_bottom, gap=gap, EH=EH, EVp2=EVp2,
                VarVp=VarVp, Gpp0=Gpp0, E0=E0, E1=E1)


def main():
    print("E2-FLIP toy: 'Savvidy Hessian bottom' Phi_Hess  vs  physical transfer gap m1^2")
    print("(mass = sqrt(gap);  nu on m^2)")

    # (A) free massless  V=0
    summarize("(A) free massless  V=0",
              V=lambda x: 0.0*x, dV=lambda x: 0.0*x, ddV=lambda x: 0.0*x)

    # (B) free massive  V = (1/2) m^2 x^2 ,  m=1  -> gap should be m^2 = 1
    m = 1.0
    summarize(f"(B) free massive  V=(1/2)m^2x^2, m={m}  (expect gap=m^2=1)",
              V=lambda x: 0.5*m*m*x*x, dV=lambda x: m*m*x, ddV=lambda x: m*m*np.ones_like(x))

    # (C) Savvidy-like double well  V=(g/4)(x^2-v^2)^2 ; scan g to show
    #     Hessian bottom and physical gap move OPPOSITELY
    v = 1.0
    print("\n\n################  (C) SAVVIDY-LIKE DOUBLE WELL: Hessian depth vs true gap  ################")
    rows = []
    for g in [0.25, 0.5, 1.0, 2.0, 4.0, 8.0]:
        V = lambda x, g=g: (g/4.0)*(x*x - v*v)**2
        dV = lambda x, g=g: g*x*(x*x - v*v)
        ddV = lambda x, g=g: g*(3.0*x*x - v*v)
        d = summarize(f"(C) double well g={g}, v={v}  (inf V'' = -g v^2 = {-g*v*v:+.3f})",
                      V, dV, ddV)
        S_inst = float(np.sqrt(2.0) * (g**0.5 / 2.0) * (2.0*v)**3 / 3.0)  # crude bounce-action scale
        rows.append((g, d['hess_bottom'], d['gap'], S_inst))

    print("\n\n---- SUMMARY TABLE (C) ----")
    print(f"{'g':>6} {'depth=-g v^2':>14} {'true gap m1^2':>16} {'gap/|depth|':>13} {'S_inst(rough)':>14}")
    for g, hb, gap, S in rows:
        ratio = gap/abs(hb) if hb != 0 else float('nan')
        print(f"{g:>6.2f} {hb:>14.5f} {gap:>16.4e} {ratio:>13.3e} {S:>14.4f}")
    print("\nREAD: |depth| GROWS with g while the true gap SHRINKS (tunnelling).")
    print("      No monotone functional of the Hessian bottom can be a lower bound on the gap.")
    print("      On the free massive field (B) Hessian bottom == gap EXACTLY (works).")
    print("      On the free massless field (A) Hessian bottom == gap == 0 (works).")
    print("      On the non-convex (Savvidy-like) field (C) Hessian bottom < 0 < gap (FAILS).")


if __name__ == "__main__":
    main()
