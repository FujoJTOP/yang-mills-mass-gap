"""Correctness checks for su3.py (no physics production here)."""
import numpy as np
import su3 as S


def test_staple():
    rng = np.random.default_rng(1234)
    L, beta = 4, 2.2
    U = S.random_su3(rng, (4, L, L, L, L))
    mu = 2
    x = (1, 3, 0, 2)
    Ssum = S.staple_sum(U, mu)
    Slink = Ssum[x]
    # random small left rotation
    xi = 1e-3 * rng.standard_normal(8)
    V = S.unitaries_from_xi(xi, np.eye(3))
    Unew = U.copy()
    Unew[mu][x] = V @ U[mu][x]
    dS = S.wilson_action(Unew, beta) - S.wilson_action(U, beta)
    pred = -(beta / 3.0) * np.real(np.trace((V - np.eye(3)) @ U[mu][x] @ Slink))
    print(f"[staple] dS={dS:.10e}  prediction={pred:.10e}  rel.err={abs(dS-pred)/abs(pred):.2e}")
    assert abs(dS - pred) < 1e-6 * abs(pred)


def test_loop_eq_plaquette():
    rng = np.random.default_rng(7)
    L = 4
    U = S.random_su3(rng, (4, L, L, L, L))
    w = S.wilson_loop(U, 1, 1)
    p = S.mean_plaquette(U)
    print(f"[W(1,1)] {w:.12f}   plaquette {p:.12f}")
    assert abs(w - p) < 1e-12


def test_gauge_invariance():
    rng = np.random.default_rng(11)
    L = 4
    U = S.random_su3(rng, (4, L, L, L, L))
    g = S.random_su3(rng, (L, L, L, L))
    Ug = S.gauge_transform(U, g)
    print(f"[gauge] dP={abs(S.mean_plaquette(U)-S.mean_plaquette(Ug)):.2e}  "
          f"dW22={abs(S.wilson_loop(U,2,2)-S.wilson_loop(Ug,2,2)):.2e}")
    assert abs(S.mean_plaquette(U) - S.mean_plaquette(Ug)) < 1e-12
    assert abs(S.wilson_loop(U, 2, 2) - S.wilson_loop(Ug, 2, 2)) < 1e-12


def test_su2_sampler():
    """Compare sampled x-marginal with direct numerical integration."""
    rng = np.random.default_rng(3)
    for lam in (0.3, 2.0, 5.0, 12.0):
        x = S._sample_w0(np.full(400000, lam), rng)
        # theory normalisation and first moments by quadrature
        xs = np.linspace(-1, 1, 20001)
        w = np.sqrt(1 - xs ** 2) * np.exp(lam * xs)
        Z = np.trapezoid(w, xs)
        m1 = np.trapezoid(xs * w, xs) / Z
        m2 = np.trapezoid(xs ** 2 * w, xs) / Z
        e1, e2 = x.mean(), (x ** 2).mean()
        s1 = x.std() / np.sqrt(len(x))
        print(f"[su2 x] lam={lam:5.1f}  <x>={e1:+.5f} (th {m1:+.5f})"
              f"  <x^2>={e2:.5f} (th {m2:.5f})  err={abs(e1-m1)/s1:.2f} sigma")
        assert abs(e1 - m1) < 5 * s1 + 1e-3
        assert abs(e2 - m2) < 0.01


def test_heatbath_vs_strong_coupling():
    """At beta=1.0 the leading strong-coupling result is <P> = beta/18."""
    rng = np.random.default_rng(99)
    L, beta = 4, 1.0
    lat = S.Lattice(L, beta, rng, cold=True)
    for _ in range(300):
        lat.sweep()
    vals = [lat.sweep() or S.mean_plaquette(lat.U) for _ in range(400)]
    m = np.mean(vals)
    se = np.std(vals) / np.sqrt(len(vals))
    print(f"[beta=1.0 L=4] <P>={m:.5f} +- {se:.5f}   beta/18={beta/18:.5f}   "
          f"sampler eff={lat.sampler_efficiency():.3f}")
    assert abs(m - beta / 18) < 0.01




# ----------------------------------------------------------------------
def metropolis_run(L, beta, n_sweeps, rng, U=None, step=0.25):
    """Independent SU(2)-subgroup Metropolis (cross-check of the heat bath)."""
    if U is None:
        U = S.random_su3(rng, (4, L, L, L, L))
    kappa = beta / 3.0
    pairs = ((0, 1), (0, 2), (1, 2))
    acc = 0
    prop = 0
    for _ in range(n_sweeps):
        Uf = U.reshape(4, -1, 3, 3)
        for c in (0, 1):
            idx = None
            # same parity checkerboard as the heat bath
            x0, x1, x2, x3 = np.meshgrid(*[np.arange(L)] * 4, indexing='ij')
            par = ((x0 + x1 + x2 + x3) % 2).reshape(-1)
            idx = np.nonzero(par == c)[0]
            for mu in range(4):
                Sst = S.staple_sum(U, mu).reshape(-1, 3, 3)[idx]
                A = Uf[mu][idx]
                M = A @ Sst
                for (i, j) in pairs:
                    # proposal: SU(2) rotation with random axis, |angle| ~ step
                    ax = rng.standard_normal((len(idx), 3))
                    ax /= np.linalg.norm(ax, axis=1, keepdims=True)
                    th = step * rng.standard_normal(len(idx))
                    c0 = np.cos(th / 2)
                    s0 = np.sin(th / 2)
                    v = np.empty((len(idx), 2, 2), complex)
                    v[:, 0, 0] = c0 - 1j * s0 * ax[:, 2]
                    v[:, 0, 1] = -1j * s0 * (ax[:, 0] - 1j * ax[:, 1])
                    v[:, 1, 0] = -1j * s0 * (ax[:, 0] + 1j * ax[:, 1])
                    v[:, 1, 1] = c0 + 1j * s0 * ax[:, 2]
                    B00 = M[:, i, i]; B01 = M[:, i, j]
                    B10 = M[:, j, i]; B11 = M[:, j, j]
                    # Delta Re tr(v B) = Re tr((v-1) B)
                    trvB = v[:, 0, 0] * B00 + v[:, 0, 1] * B10 + v[:, 1, 0] * B01 + v[:, 1, 1] * B11
                    trB = B00 + B11
                    dS = -kappa * np.real(trvB - trB)
                    r = rng.random(len(idx))
                    ok = r < np.exp(-np.clip(dS, -700, 700))
                    acc += int(ok.sum()); prop += len(idx)
                    if ok.any():
                        vv = v[ok]
                        sel = idx[ok]
                        A2 = A[ok]
                        M2 = M[ok]
                        Mi = M2[:, i, :].copy(); Mj = M2[:, j, :].copy()
                        M2[:, i, :] = vv[:, 0, 0, None] * Mi + vv[:, 0, 1, None] * Mj
                        M2[:, j, :] = vv[:, 1, 0, None] * Mi + vv[:, 1, 1, None] * Mj
                        Ai = A2[:, i, :].copy(); Aj = A2[:, j, :].copy()
                        A2[:, i, :] = vv[:, 0, 0, None] * Ai + vv[:, 0, 1, None] * Aj
                        A2[:, j, :] = vv[:, 1, 0, None] * Ai + vv[:, 1, 1, None] * Aj
                        A[ok] = A2
                        M[ok] = M2
                Uf[mu][idx] = A
    return U, acc / max(prop, 1)


def test_heatbath_vs_metropolis():
    """Metropolis started from an equilibrated heat-bath configuration must stay
    at the heat-bath plaquette (stationarity / absence of algorithmic bias)."""
    for beta in (1.0, 2.2):
        rng1 = np.random.default_rng(2024)
        rng2 = np.random.default_rng(2025)
        lat = S.Lattice(4, beta, rng1, cold=False)
        for _ in range(1500):
            lat.sweep()
        U = lat.U.copy()
        hb = []
        for _ in range(600):
            lat.sweep()
            hb.append(S.mean_plaquette(lat.U))
        hb = np.array(hb)
        U, ar = metropolis_run(4, beta, 200, rng2, U=U, step=0.8)
        mh = []
        for k in range(400):
            U, _ = metropolis_run(4, beta, 1, rng2, U=U, step=0.8)
            mh.append(S.mean_plaquette(U))
        mh = np.array(mh)

        def err(x):
            return x.std(ddof=1) / np.sqrt(len(x))

        d = abs(hb.mean() - mh.mean())
        print(f"[HB vs Metropolis beta={beta}]  HB={hb.mean():.5f}+-{err(hb):.5f}  "
              f"Met={mh.mean():.5f}+-{err(mh):.5f}  met_acc={ar:.3f}  diff={d:.5f}")
        assert d < 3 * np.hypot(err(hb), err(mh)) + 2e-3


def _main():
    test_staple()
    test_loop_eq_plaquette()
    test_gauge_invariance()
    test_su2_sampler()
    test_heatbath_vs_strong_coupling()
    test_heatbath_vs_metropolis()
    print('ALL CHECKS PASSED')


if __name__ == '__main__':
    _main()
