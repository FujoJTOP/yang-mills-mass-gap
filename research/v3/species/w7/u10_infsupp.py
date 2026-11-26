"""U10 - inf supp mu on the lattice: is E_eff(t) := -(1/t)log|C(t)| an upper or lower bound?

C(t) = int e^{-E t} dmu(E),  mu >= 0 normalised (C(0)=1),  E_eff(t) := -(1/t) log|C(t)|.

Part A  exact direction test on model measures (thin atom, gapless continuum, gapped atomless)
Part B  SU(2) beta=2.5 L=2 lattice: 0++ plaquette-plaquette correlator, E_eff(1) = E_eff(L/2)

py -3.12 u10_infsupp.py > u10_infsupp.log
"""
import numpy as np


# ------------------------------------------------------------------ Part A
def E_eff(pts, ws, t):
    return -np.log(abs(np.sum(np.asarray(ws) * np.exp(-np.asarray(pts) * t)))) / t


def part_A():
    print("== Part A: exact direction test on model measures ==")
    ts = [0.25, 0.5, 1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0]

    w, eps, S = 0.3, 1e-2, 1.5
    print(f"\n[A1] thin atom {w}*d_{eps}+{1-w}*d_{S}   inf supp = {eps}")
    print("        t      E_eff(t)     E_eff - inf supp")
    for t in ts:
        e = E_eff([eps, S], [w, 1 - w], t)
        print(f"  {t:7.2f} {e:12.6f}   {e - eps:+.3e}")

    s = np.linspace(0, 1, 200001)[1:]
    nu = np.sqrt(s) / np.trapezoid(np.sqrt(s), s)          # CFT-like, edge at 0
    print(f"\n[A2] gapless continuum dnu ~ s^(1/2) on [0,1]   inf supp = 0")
    print("        t      E_eff(t)")
    for t in ts:
        C = np.trapezoid(nu * np.exp(-s * t), s)
        print(f"  {t:7.2f} {-np.log(abs(C)) / t:12.6f}")

    m2 = 1.0
    ss = m2 + np.linspace(0, 3, 300001)[1:]
    rho = np.sqrt(ss - m2) / np.trapezoid(np.sqrt(ss - m2), ss)   # 2-particle sqrt edge
    print(f"\n[A3] gapped atomless sqrt edge at {m2}   inf supp = {m2}")
    print("        t      E_eff(t)")
    for t in ts:
        C = np.trapezoid(rho * np.exp(-ss * t), ss)
        print(f"  {t:7.2f} {-np.log(abs(C)) / t:12.6f}")

    e_gap = -np.log(abs(np.trapezoid(nu * np.exp(-s), s)))
    e_gapped = -np.log(abs(np.trapezoid(rho * np.exp(-ss), ss)))
    print(f"\n[A4] ACID TEST: E_eff(1) gapless = {e_gap:.6f} ; gapped(atomless,m^2=1) = {e_gapped:.6f}")
    print("     -> one finite-t number does not separate inf supp=0 from inf supp>0;")


# ------------------------------------------------------------------ Part B
def plaquette_spatial(U, L):
    """P[t,x1,x2,x3] = (1/3) sum_{spatial planes} [1 - 1/2 Re tr U_p]."""
    acc = np.zeros((L, L, L, L))
    for i, j in ((1, 2), (1, 3), (2, 3)):
        Ui = U[i]
        Ujf = np.roll(U[j], -1, axis=i)
        Uj = U[j]
        Uif = np.roll(U[i], -1, axis=j)
        M = np.matmul(np.matmul(Ui, Ujf),
                      np.matmul(Uj.conj().swapaxes(-1, -2), Uif.conj().swapaxes(-1, -2)))
        acc += 1.0 - 0.5 * M[..., 0, 0].real
    return acc / 3.0


def run_lattice(L=2, beta=2.5, nmeas=4000, gap=2, seed=20261125, nblk=40):
    """SU(2) sampler with the repo's beta=2.5 convention (m1_4d_run.make_stream, kappa=beta/2)."""
    import sys, os
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'w7'))
    import m1_4d_run as M
    U, sweep = M.make_stream(L, beta, seed, cold=True)
    for _ in range(800):
        sweep()

    A = np.empty((nmeas, L))
    for n in range(nmeas):
        for _ in range(gap):
            sweep()
        A[n] = M.plaq_operator_t(U)

    nt = L // 2 + 1

    def C_of(idx):
        a = A[idx]
        abar = a.mean()
        return np.array([np.mean(np.roll(a, -tau, axis=1) * a) - abar ** 2 for tau in range(nt)])

    C = C_of(np.arange(nmeas))
    blocks = np.array_split(np.arange(nmeas), nblk)
    jk = np.empty((nt - 1, nblk))
    for kb in range(nblk):
        idx = np.concatenate([bl for k, bl in enumerate(blocks) if k != kb])
        Cb = C_of(idx)
        jk[:, kb] = -np.log(abs(Cb[1:] / Cb[0]))

    print(f"\n== SU(2) beta={beta}, L={L} ({L}^4 periodic), nmeas={nmeas}, "
          f"sweeps/meas={gap}, seed={seed} ==")
    print("  t   C(t)=<O(t)O(0)>_c (zeromom)    C(t)/C(0)      E_eff(t)=-log|C(t)/C(0)|")
    for tau in range(nt):
        if tau == 0:
            print(f"  {tau}   {C[0]:+.6e}          1.000000       (E_eff undefined at t=0)")
        else:
            r = C[tau] / C[0]
            e = -np.log(abs(r))
            err = np.sqrt((nblk - 1) / nblk * np.sum((jk[tau - 1] - jk[tau - 1].mean()) ** 2))
            tag = "   <-- t=L/2 (wrap: fwd+bwd both counted)" if tau == L // 2 else ""
            print(f"  {tau}   {C[tau]:+.6e}          {r:+.6f}       {e:.6f} +- {err:.6f}{tag}")
            if r < 0:
                print("        WARNING: ratio<0 -> |.| in play; mu>=0 violated or noise dominated")
    return C


def part_B():
    run_lattice(L=2, nmeas=8000, gap=2)
    run_lattice(L=4, nmeas=900, gap=5)


if __name__ == "__main__":
    part_A()
    part_B()
