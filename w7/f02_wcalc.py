"""
F02 W-CALC -- Savvidy negative mode vs. positive multiplicative weights.
Pointwise cancellation: possible or not?

Exact identities tested here (f = W^{1/2}, M = -Delta + V self-adjoint, A = f M f^{-1}):

  T1  f(-Delta)f^{-1} = -Delta + 2(grad log f).grad + (Delta log f - |grad log f|^2)
  T2  1/2(A + A*) = M - (1/4)|grad log W|^2                (<= M everywhere: correction is <= 0)
  T3  spec(A) = spec(M)   (similarity)  -> lam_min unchanged for EVERY W
  T4  <psi|A|psi> = <psi|V|psi> - (1/4)|| grad log(W psi^2) ||^2_{psi}   <= <psi|V|psi>
      => for the Savvidy mode (V_LLL = 2 s gB = -2gB, <D psi>^2 = gB) this reads
         <psi0|A psi0> = - (1/4)|| grad log(W psi0^2) ||^2_{psi0} <= 0 ,
         equality iff W = kappa / psi0^2  (the "oracle": spectral, singular, not a certificate).

§1 = exact Savvidy location (Landau ladder, oscillator basis).
§2-§4 = conjugation audit on a lattice model carrying a localised negative mode of the same
        shape (LLL-like localisation, eigenvalue -gB), pure numpy/scipy.
"""
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as sla

OUT = []


def say(*a):
    s = " ".join(str(x) for x in a)
    OUT.append(s)
    print(s, flush=True)


# ----------------------------------------------------------------------------------------
def sec1():
    say("=" * 84)
    say("§1  WHERE the Savvidy negative mode is (exact): colour-3 transverse magnetic B")
    say("=" * 84)
    gB, N = 1.0, 64
    a = np.diag(np.sqrt(np.arange(1, N)), 1)
    num = a.T @ a                                   # a^dag a
    lap = gB * (2.0 * num + np.eye(N))              # -D_perp^2 = gB(2 a^dag a + 1)
    for s in (+1, -1):
        ev = np.linalg.eigvalsh(lap + 2.0 * s * gB * np.eye(N))
        say(f"   helicity s={s:+d}:  lam_min = {ev.min():+.12f}")
    say("   => lam_{n,s} = k_par^2 + (2n+1) g|B| + 2 s g|B| ; bottom at n=0, s=-1, k_3=k_4=0:")
    say("      lam_min = -g|B|.  Direction: transverse magnetic (colour 3 background,")
    say("      eps_- = (e_1 - i e_2)/sqrt2), LLL profile in the 1-2 plane, constant in 3,4.")
    say("      (matches PROOF-SAVVIDY-EXACT-SIGN; lattice: O06 'pointwise min -1.893', W4 7.1% neg)")
    say("")
    return gB


def laplace2d(n, h):
    N = n * n
    T = sp.lil_matrix((N, N))
    for i in range(n):
        for j in range(n):
            p = i * n + j
            for di, dj in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                ii, jj = i + di, j + dj
                if 0 <= ii < n and 0 <= jj < n:
                    T[p, p] += 1.0 / h ** 2
                    T[p, ii * n + jj] += -1.0 / h ** 2
    return T.tocsr()


def grad_h(f, n, h):
    X = f.reshape(n, n)
    g1 = np.zeros_like(X); g2 = np.zeros_like(X)
    g1[1:-1, :] = (X[2:, :] - X[:-2, :]) / (2 * h)
    g2[:, 1:-1] = (X[:, 2:] - X[:, :-2]) / (2 * h)
    g1[0, :] = (X[1, :] - X[0, :]) / h; g1[-1, :] = (X[-1, :] - X[-2, :]) / h
    g2[:, 0] = (X[:, 1] - X[:, 0]) / h; g2[:, -1] = (X[:, -1] - X[:, -2]) / h
    return g1.ravel(), g2.ravel()


def sec2(gB):
    say("=" * 84)
    say("§2  Lattice model with a localised negative mode of the Savvidy shape")
    say("    M = -Delta_h + V,  V = -c inside a disk of radius R (depth tuned to lam_min ~ -gB)")
    say("=" * 84)
    n, L, R = 64, 6.0, 2.4
    h = L / (n - 1)
    x = np.linspace(0, L, n)
    X, Y = np.meshgrid(x, x, indexing="ij")
    r = np.hypot(X - L / 2, Y - L / 2)
    disk = (r <= R).astype(float).ravel()
    T = laplace2d(n, h)
    best = None
    for c in (1.6, 1.8, 2.0, 2.2):
        M = T - sp.diags(c * disk)
        w = sla.eigsh(M, k=1, which="SA", return_eigenvectors=False)[0]
        if best is None or abs(w + gB) < abs(best[1] + gB):
            best = (c, w)
    c = best[0]
    M = (T - sp.diags(c * disk)).tocsr()
    w, v = sla.eigsh(M, k=2, which="SA")
    lam_min = w[0]
    psi = v[:, 0]
    k = np.argmax(np.abs(psi)); psi = psi * np.exp(-1j * np.angle(psi[k]))
    psi = psi / np.linalg.norm(psi)
    V = -c * disk
    p = np.abs(psi) ** 2
    say(f"    grid {n}x{n}, h={h:.5f}, L={L}, depth c={c}, target lam_min=-gB={-gB}")
    say(f"    lambda_min(M) = {lam_min:+.8f}   (next {w[1]:+.6f})")
    say(f"    <psi|V|psi>   = {np.sum(p*V):+.8f}  (= <psi|M|psi> - ||grad psi||^2)")
    say(f"    localisation: max site weight {p.max():.3e}, sites with w>1e-3: {(p>1e-3).sum()}")
    say("    (transverse LLL is localised exactly this way: width 1/sqrt(gB) << box)")
    say("")
    return dict(n=n, h=h, M=M, psi=psi, lam_min=lam_min, V=V, gB=gB, c=c)


def sec3(mod):
    n, h, M, psi, lam_min, V, gB = (mod[k] for k in ("n", "h", "M", "psi", "lam_min", "V", "gB"))
    say("=" * 84)
    say("§3  Conjugation audit   A = W^{1/2} M W^{-1/2}   (candidate designs)")
    say("=" * 84)
    x = np.linspace(0, 6.0, n)
    X, Y = np.meshgrid(x, x, indexing="ij")
    X = X.ravel(); Y = Y.ravel()
    dens = np.abs(psi) ** 2
    Phi = 0.5 * gB * gB * ((X - 3) ** 2 + (Y - 3) ** 2)
    Fp = gB ** 2 * (1.0 + 0.1 * ((X - 3) ** 2 + (Y - 3) ** 2))
    rng = np.random.default_rng(7)
    Ws = {
        "e^{-c Phi}, c=0.02": np.exp(-0.02 * Phi),
        "e^{-c Phi}, c=0.10": np.exp(-0.10 * Phi),
        "(1+|F|^2)^k, k=1/2": Fp ** 0.5,
        "(1+|F|^2)^k, k=1": Fp ** 1.0,
        "det(Hess Phi)^a (const)": np.ones(n * n),
        "random smooth": np.abs(2 + np.cos(0.7 * X) * np.cos(0.9 * Y)
                                + 0.1 * rng.standard_normal(n * n)) + 0.5,
        "oracle kappa/psi^2": 1.0 / (dens + 1e-13),
    }
    say(f"    {'W':26s} {'||A s.psi - lam min s.psi||':>26s} {'<psi|A|psi>':>13s} "
        f"{'pred=lam_min-¼<|grad logW|²>':>28s} {'lam_min(sym)':>13s}")
    rows = []
    for name, W in Ws.items():
        s = np.sqrt(W)
        A = (sp.diags(s) @ M @ sp.diags(1.0 / s)).tocsr()
        ph = s * psi; ph = ph / np.linalg.norm(ph)
        resid = np.linalg.norm(A @ ph - lam_min * ph)
        wit = float(np.real(np.vdot(psi, A @ psi)))
        g1, g2 = grad_h(np.log(W), n, h)
        pen = 0.25 * float(np.sum(dens * (g1 ** 2 + g2 ** 2)))
        pred = lam_min - pen
        sym = (0.5 * (A + A.getH())).tocsr()
        lsym = float(sla.eigsh(sym, k=1, which="SA", return_eigenvectors=False)[0])
        say(f"    {name:26s} {resid:26.3e} {wit:+13.6f} {pred:+28.6f} {lsym:+13.6f}")
        rows.append(dict(W=name, resid=resid, witness=wit, pred=pred, lsym=lsym))
    say("")
    say("    col2 ~ 0  : W^{1/2}psi is an EXACT eigenvector of A with eigenvalue lam_min(M) -- for every W")
    say("    col3 < 0  : the negative mode witnesses against A >= Delta^2 > 0 -- for every designed W")
    say("    col4 = col3: <psi|A|psi> = lam_min - ¼<|grad log W|^2>_psi  <= lam_min  (identity T4)")
    say("                 => the best weight in the whole class is the TRIVIAL one (grad log W = 0 on supp psi);")
    say("                    every designed gradient makes the negative mode strictly worse.")
    say("    col5 <= lam_min(M): the self-adjoint part is pushed DOWN by -¼|grad log W|²")
    say("")
    return rows


def sec4(mod):
    say("=" * 84)
    say("§4  The conjugation term EXACTLY: bonds get multiplied by cosh(δ log f) >= 1, diagonal untouched")
    say("=" * 84)
    n, h, M = mod["n"], mod["h"], mod["M"]
    x1 = np.linspace(0, 6.0, n)
    X1, Y1 = np.meshgrid(x1, x1, indexing="ij")
    lw = (0.3 * np.sin(X1) * np.cos(0.8 * Y1)).ravel()
    s = np.sqrt(np.exp(lw))
    A = (sp.diags(s) @ M @ sp.diags(1.0 / s)).tocsr()
    sym = (0.5 * (A + A.getH())).tocsr()
    # exact discrete penalty  P = (1/h^2) sum_mu (cosh(delta_mu log f) - 1)(S_mu + S_-mu)
    P = sp.lil_matrix((n * n, n * n))
    cmin, cmax = 1e9, -1e9
    for i in range(n):
        for j in range(n):
            p = i * n + j
            f0 = lw[p]
            for di, dj in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                ii, jj = i + di, j + dj
                if 0 <= ii < n and 0 <= jj < n:
                    d = lw[ii * n + jj] - f0
                    c = np.cosh(0.5 * d) - 1.0     # delta log f = (1/2) delta log W
                    cmin = min(cmin, c); cmax = max(cmax, c)
                    P[p, ii * n + jj] += c / h ** 2
    P = P.tocsr()
    D = (sym - (M - P)).tolil()
    D[0, :] = 0; D[-1, :] = 0; D[:, 0] = 0; D[:, -1] = 0
    say(f"    max|½(A+A*) - (M - P_cosh)| (interior) = {np.max(np.abs(D.tocsr().data)):.3e}   [exact identity]")
    pen = 0.25 * (grad_h(lw, n, h)[0] ** 2 + grad_h(lw, n, h)[1] ** 2)
    say(f"    continuum form: P -> ¼|grad logW|² ;  max(P_cosh)={cmax/h**2:.4e}, "
        f"max(¼|grad logW|²)={pen.max():.4e}")
    dg = np.abs(np.real(sym.diagonal()) - np.real(M.diagonal()))
    say(f"    diagonal correction |diag(A)-diag(M)| : max = {dg.max():.3e}  (exactly zero: A_ii = M_ii)")
    say(f"    every bond weight cosh(delta log f) >= 1 : min = {1+cmin:.6f}, max = {1+cmax:.6f}"
        f"  -> hopping is strictly ENLARGED, never reduced")
    say("    => no positive zero-order term is created anywhere; the only effect is to spread the")
    say("       operator over larger hoppings, which lowers the bottom (spectral flow downward).")
    say("")
    return float(np.max(np.abs(D.tocsr().data))), float(1 + cmin)


def sec5():
    say("=" * 84)
    say("§5  General theorem: spectrum is similarity-invariant for ANY positive invertible W")
    say("=" * 84)
    rng = np.random.default_rng(11)
    d = 80
    Q = rng.standard_normal((d, d)); Q = Q + Q.T
    M = Q @ Q.T / d + np.diag(rng.standard_normal(d)) * 2.0
    M = 0.5 * (M + M.T)
    ev = np.linalg.eigvalsh(M)
    say(f"    random self-adjoint M ({d}x{d}): lam_min={ev.min():+.8f}  lam_max={ev.max():+.8f}")
    err = 0.0
    for _ in range(6):
        S = np.eye(d) + 0.4 * rng.standard_normal((d, d))
        S = S @ S.T + 0.05 * np.eye(d)
        ea = np.sort(np.linalg.eigvals(S @ M @ np.linalg.inv(S)))
        err = max(err, float(np.max(np.abs(ea.real - ev)) + np.max(np.abs(ea.imag))))
    say(f"    max |spec(S M S^-1) - spec(M)| over 6 random SPD S = {err:.3e}")
    say("    => no W of any kind (multiplicative or not) can move the bottom of the spectrum.")
    say("")
    return err


def main():
    gB = sec1()
    mod = sec2(gB)
    rows = sec3(mod)
    errs, minbond = sec4(mod)
    e5 = sec5()
    say("=" * 84)
    say("§6  VERDICT")
    say("=" * 84)
    say(f"    (a) negative-mode location      : lam_min = -g|B| (n=0, s=-1, k_3=k_4=0)  [exact]")
    say(f"    (b) similarity residual, worst W: {max(r['resid'] for r in rows):.2e} -> lam_min(A)=lam_min(M)")
    say(f"    (c) witness <psi0|A|psi0>       : {max(r['witness'] for r in rows):+.6f}"
        f" (max over designs; = lam_min - ¼<|grad logW|²> <= lam_min)")
    say(f"    (d) conjugation term            : exact lattice identity max err {errs:.2e};"
        f" every bond weight >= {minbond:.4f} (>=1, i.e. hopping only grows)")
    say(f"    (e) arbitrary SPD W             : spec shift {e5:.2e}")
    say("    => W-CALC (pointwise cancellation of the Savvidy negative mode by a positive weight)")
    say("       DOES NOT EXIST.  Verdict (c): structural theorem. D1/W-SHIFT is a similarity.")
    with open("f02_wcalc.log", "w", encoding="utf-8") as fh:
        fh.write("\n".join(OUT) + "\n")


if __name__ == "__main__":
    main()
