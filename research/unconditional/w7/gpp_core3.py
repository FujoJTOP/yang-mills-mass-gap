"""SU(3) counterpart of gpp_core for the exact marginal curvature G''(0).

Force/curvature are taken directly from the Wilson action along the ray
U -> exp(i s v) U  (first/second differences of su3.wilson_action), so no
gradient routine is needed and the SU(2)/SU(3) estimators are identical.
"""
import sys
import numpy as np

sys.path.insert(0, 'D:/Dev/math/research/lattice-cpb')
sys.path.insert(0, 'D:/Dev/math/research/unconditional/w7')
import su3 as S3          # noqa: E402
import gpp_core as C      # reuse veff_curvature, ev_estimator  # noqa: E402

TG = S3.T_GEN             # (8,3,3), tr(t^a t^b)=1/2 delta


def su3_log_x(U):
    """Principal log: x_a = 2 Im tr(t^a log U), U=exp(i x_a t^a)."""
    w, V = np.linalg.eig(U)
    th = np.angle(w)
    Vin = np.linalg.inv(V)
    logU = np.einsum('...ij,...j,...jk->...ik', V, 1j * th, Vin)
    trt = np.einsum('aij,...ji->...a', TG, logU)   # complex
    return 2.0 * trt.imag


def apply_exp(U, comp):
    """U -> exp(i comp.t) U, comp: (4,L,L,L,L,8)."""
    A = np.einsum('...a,aij->...ij', comp, TG)     # Hermitian
    w, V = np.linalg.eigh(A)
    E = np.einsum('...ik,...k,...jk->...ij', V, np.exp(1j * w), np.conj(V))
    return E @ U


def plane_mode(L, axis, pol, color, n=1, kind='cos'):
    xs = np.arange(L)
    xa = xs.reshape([L if i == axis else 1 for i in range(4)])
    ph = (np.cos if kind == 'cos' else np.sin)(2.0 * np.pi * n * xa / L)
    ph = np.broadcast_to(ph, (L, L, L, L)).astype(np.float64)
    v = np.zeros((4, L, L, L, L, 8))
    v[pol, :, :, :, :, color] = ph
    v /= np.linalg.norm(v)
    return v


def ray_action(U, comp, beta, s):
    return S3.wilson_action(apply_exp(U, s * comp), beta)


def mode_obs(U, beta, modes, h=1e-3):
    x = su3_log_x(U)
    S0 = S3.wilson_action(U, beta)
    a, Hvv, Gv = [], [], []
    for v in modes:
        a.append(float(np.sum(v * x)))
        Sp = ray_action(U, v, beta, h)
        Sm = ray_action(U, v, beta, -h)
        Gv.append((Sp - Sm) / (2.0 * h))
        Hvv.append((Sp - 2.0 * S0 + Sm) / (h * h))
    return np.array(a), np.array(Hvv), np.array(Gv)


def sample_cfgs(L, beta, ntherm, ncfg, gapevery, seed):
    rng = np.random.default_rng(seed)
    lat = S3.Lattice(L, beta, rng, cold=True)
    for _ in range(ntherm):
        lat.sweep()
    out = []
    for _ in range(ncfg):
        for _ in range(gapevery):
            lat.sweep()
        out.append(lat.U.copy())
    return out


def panel(L, colors=(0, 3, 7), nmax=1):
    modes = {}
    for axis in range(4):
        for pol in range(4):
            if pol == axis:
                continue
            for n in range(1, nmax + 1):
                for kind in ('cos',):
                    for color in colors:
                        key = f"ax{axis}_p{pol}_n{n}_{kind}_c{color}"
                        modes[key] = plane_mode(L, axis, pol, color, n=n, kind=kind)
    return modes
