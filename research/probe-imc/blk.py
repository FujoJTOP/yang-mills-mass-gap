"""b=2 link blocking of SU(3) Wilson configurations, and blocked-lattice
observables (plaquette / 6-link / 8-link densities, character coefficients,
2-point correlators, 4-point ratio).

Blocking (Teper/Swendsen type, scale-2):
  V_mu(y) = P_SU3 [ U_mu(2y) U_mu(2y+mu)
                    + sum_{nu != mu} ( B1 + B2 ) ]
  B1 = U_nu(2y) U_mu(2y+nu) U_mu(2y+nu+mu) U_nu^dag(2y+2mu)
  B2 = U_nu^dag(2y-nu) U_mu(2y-nu) U_mu(2y-nu+mu) U_nu(2y-nu+2mu)
all seven paths weighted equally, then projected onto SU(3).

Conventions match su3.py: U shape (4,L,L,L,L,3,3), U[mu][x] links.
"""
import numpy as np

import os
import sys
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                              '..', 'lattice-cpb'))

import su3 as S


# ---------------------------------------------------------------- blocking
def _ev(U, mu, off):
    """U[mu] at fine sites 2y+off, off a dict axis->{0,1}. Shape (n^4,3,3)."""
    sl = []
    for d in range(4):
        sl.append(slice(off.get(d, 0), None, 2))
    return U[mu][tuple(sl)]


def block_links(U, wstaple=1.0):
    """b=2 blocked links V, shape (4,n,n,n,n,3,3), n = L/2.

    wstaple: weight of the six 4-link staple paths relative to the straight
    double link (default 1.0 = the equal-weight blocking used in the report).
    """
    L = U.shape[1]
    n = L // 2
    V = np.empty((4, n, n, n, n, 3, 3), dtype=np.complex128)
    for mu in range(4):
        acc = _ev(U, mu, {}) @ _ev(U, mu, {mu: 1})
        for nu in range(4):
            if nu == mu:
                continue
            B1 = (_ev(U, nu, {}) @ _ev(U, mu, {nu: 1})
                  @ _ev(U, mu, {nu: 1, mu: 1})
                  @ S.adj(np.roll(_ev(U, nu, {}), -1, axis=mu)))
            B2 = (S.adj(np.roll(_ev(U, nu, {nu: 1}), 1, axis=nu))
                  @ np.roll(_ev(U, mu, {nu: 1}), 1, axis=nu)
                  @ np.roll(_ev(U, mu, {nu: 1, mu: 1}), 1, axis=nu)
                  @ np.roll(np.roll(_ev(U, nu, {nu: 1}), -1, axis=mu), 1, axis=nu))
            acc = acc + wstaple * (B1 + B2)
        V[mu] = S.project_su3_cuberoot(acc)
    return V


# ---------------------------------------------------------------- loops
def _off(*pairs):
    o = [0, 0, 0, 0]
    for d, k in pairs:
        o[d] += k
    return tuple(o)


def mat_loop(V, steps):
    """Product of links along steps=[(dir, offset4, conj), ...];  (N,3,3)."""
    A = None
    for d, off, conj in steps:
        M = V[d]
        for ax in range(4):
            if off[ax]:
                M = np.roll(M, -off[ax], axis=ax)
        if conj:
            M = S.adj(M)
        A = M if A is None else A @ M
    return A


def _plaquette_steps(mu, nu):
    return [(mu, _off(), False), (nu, _off((mu, 1)), False),
            (mu, _off((nu, 1)), True), (nu, _off(), True)]


def _rect_steps(mu, nu, nmu, nnu):
    """Rectangle with nmu links along mu and nnu along nu."""
    st = []
    for k in range(nmu):
        st.append((mu, _off((mu, k)), False))
    for k in range(nnu):
        st.append((nu, _off((mu, nmu), (nu, k)), False))
    for k in range(nmu - 1, -1, -1):                    # return along -mu
        st.append((mu, _off((mu, k), (nu, nnu)), True))
    for k in range(nnu - 1, -1, -1):                    # return along -nu
        st.append((nu, _off((nu, k)), True))
    return st


PAIRS = [(mu, nu) for mu in range(4) for nu in range(mu + 1, 4)]


def loop_matrices(V, kind):
    """Stack of loop matrices (nloop, N,3,3) for a loop family."""
    if kind == 'P':
        steps = [_plaquette_steps(mu, nu) for mu, nu in PAIRS]
    elif kind == 'R':          # 1x2 rectangles, both orientations
        steps = []
        for mu, nu in PAIRS:
            steps.append(_rect_steps(mu, nu, 2, 1))
            steps.append(_rect_steps(mu, nu, 1, 2))
    elif kind == 'S':          # 2x2 squares
        steps = [_rect_steps(mu, nu, 2, 2) for mu, nu in PAIRS]
    elif kind == 'R3':         # 1x3 rectangles, both orientations
        steps = []
        for mu, nu in PAIRS:
            steps.append(_rect_steps(mu, nu, 3, 1))
            steps.append(_rect_steps(mu, nu, 1, 3))
    else:
        raise ValueError(kind)
    return np.stack([mat_loop(V, st) for st in steps])


def loop_density(V, kind):
    """Mean over loops of (1/3) Re tr, shape-agnostic scalar."""
    M = loop_matrices(V, kind)
    return float(np.real(np.trace(M, axis1=-2, axis2=-1)).mean() / 3.0)


# ---------------------------------------------------------------- observables
def char_coeffs(Mp):
    """Mp (nloop,N,3,3) plaquette matrices -> c3,c6,c8,c10 (means)."""
    t = np.trace(Mp, axis1=-2, axis2=-1)
    t2 = np.trace(Mp @ Mp, axis1=-2, axis2=-1)
    t3 = np.trace(Mp @ Mp @ Mp, axis1=-2, axis2=-1)
    c3 = float(np.real(t).mean())
    c6 = float(np.real(0.5 * (t * t + t2)).mean())
    c8 = float(np.real(np.abs(t) ** 2 - 1.0).mean())
    c10 = float(np.real((t ** 3 + 3.0 * t * t2 + 2.0 * t3) / 6.0).mean())
    return c3, c6, c8, c10


def field_P(V):
    """Site field P(x) = (1/6) sum_{mu<nu} (1/3) Re tr U_p."""
    Mp = loop_matrices(V, 'P')            # (6,N,3,3)
    tr = np.real(np.trace(Mp, axis1=-2, axis2=-1)) / 3.0
    return tr.mean(axis=0)                # (N,)


def corr_P(P, n, orbits):
    """Connected plaquette correlator for displacement orbits (list of 4-tuples)."""
    P2 = P.reshape(n, n, n, n)
    mu = P2.mean()
    out = {}
    for r in orbits:
        acc = 0.0
        cnt = 0
        for sign in (1, -1):
            sh = tuple(sign * x for x in r)
            A = np.roll(P2, [-x for x in sh], axis=(0, 1, 2, 3))
            acc += float(np.mean(P2 * A))
            cnt += 1
        out[r] = acc / cnt - mu * mu
    return out


def corr_link(V, n, rmax):
    """W2(r) = <(1/3) Re tr V_mu(x) V_mu(x+r mu_hat)>, straight 2-link line."""
    out = {}
    for r in range(1, rmax + 1):
        acc = 0.0
        cnt = 0
        for d in range(4):
            M = V[d]
            B = np.roll(M, -r, axis=d)     # V_d(x + r d)
            W = np.real(np.trace(M @ S.adj(B), axis1=-2, axis2=-1)) / 3.0
            acc += float(W.mean())
            cnt += 1
        out[r] = acc / cnt
    return out


def r4(P, n):
    """R4 = kappa4/kappa2^2 of the blocked plaquette field."""
    v = P - P.mean()
    k2 = float(np.mean(v * v))
    k4 = float(np.mean(v ** 4) - 3.0 * k2 * k2)
    return k4 / (k2 * k2) if k2 > 0 else 0.0


HIST_EDGES = np.linspace(-1.0, 1.0, 81)


def measure_all(V, orbits):
    """All blocked observables of one configuration -> dict of scalars."""
    n = V.shape[1]
    rec = {}
    for kind in ('P', 'R', 'S', 'R3'):
        rec['op_' + kind] = loop_density(V, kind)
    Mp = loop_matrices(V, 'P')
    Pv = np.real(np.trace(Mp, axis1=-2, axis2=-1)) / 3.0
    rec['hist'] = np.histogram(Pv.ravel(), bins=HIST_EDGES)[0].astype(float)
    rec['P2'] = float(np.mean(Pv ** 2))
    rec['P3'] = float(np.mean(Pv ** 3))
    pv = Pv.ravel()
    mv = float(pv.mean())
    for tag, th in (('0', 0.0), ('01', 0.01), ('02', 0.02), ('03', 0.03),
                    ('05', 0.05)):
        rec['n_' + tag] = float(np.mean(pv < th))
    rec['n_hh'] = float(np.mean(pv < 0.5 * mv))
    rec['n_h75'] = float(np.mean(pv < 0.75 * mv))
    c3, c6, c8, c10 = char_coeffs(Mp)
    rec['c3'] = c3
    rec['c6'] = c6
    rec['c8'] = c8
    rec['c10'] = c10
    P = field_P(V)
    rec['R4'] = r4(P, n)
    for r, val in corr_P(P, n, orbits).items():
        rec['CP_' + ''.join(map(str, r))] = val
    for r, val in corr_link(V, n, min(2, n // 2 + 1)).items():
        rec['CL_%d' % r] = val
    return rec


def orbits_for(n):
    """Displacement orbits available on an n^4 torus (n=3 or 4)."""
    if n == 3:
        return [(1, 0, 0, 0), (1, 1, 0, 0)]
    return [(1, 0, 0, 0), (2, 0, 0, 0), (1, 1, 0, 0), (1, 1, 1, 0), (2, 1, 0, 0)]
