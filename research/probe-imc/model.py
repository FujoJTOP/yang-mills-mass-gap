"""Metropolis sampler for a blocked-lattice action
    S_eff = sum_kind w_kind * sum_loops (1 - (1/3) Re tr loop)
with loop families  P (4-link), R (1x2, 6-link), S (2x2, 8-link), R3 (1x3, 8-link).

The action is not linear in a single link, so heat bath does not apply; we use
Metropolis with Haar-random SU(2)-subgroup proposals (Cabibbo-Marinari style)
and the full local Delta S.

Loop table (built once in numpy, consumed by the numba kernel):
  every loop is an ordered list of (flat link, conj) pairs; for every link it
  contains we register a *term* (loop, position).
"""
import numpy as np
from numba import njit

import blk


# ---------------------------------------------------------------- loop table
def _loop_coords(kind, mu, nu):
    """Ordered (offset4, dir, conj) list for a loop at base site 0."""
    if kind == 'P':
        return [(blk._off(), mu, 0), (blk._off((mu, 1)), nu, 0),
                (blk._off((nu, 1)), mu, 1), (blk._off(), nu, 1)]
    nmu, nnu = {'R': (2, 1), 'S': (2, 2), 'R3': (3, 1)}[kind]
    st = []
    for k in range(nmu):
        st.append((blk._off((mu, k)), mu, 0))
    for k in range(nnu):
        st.append((blk._off((mu, nmu), (nu, k)), nu, 0))
    for k in range(nmu - 1, -1, -1):
        st.append((blk._off((mu, k), (nu, nnu)), mu, 1))
    for k in range(nnu - 1, -1, -1):
        st.append((blk._off((nu, k)), nu, 1))
    return st


def build_table(L, kinds):
    """kinds = list of family names; returns numpy arrays for the kernel."""
    N = L ** 4
    loop_defs = []                     # (kind_index, [(link,conj),...])
    CHAR_OF = {'C6': 6, 'C8': 8, 'C10': 10}
    for ki, kind in enumerate(kinds):
        if kind in ('P', 'C6', 'C8', 'C10'):
            specs = [('P', mu, nu) for mu in range(4) for nu in range(mu + 1, 4)]
        elif kind in ('R', 'R3'):
            specs = []
            for mu in range(4):
                for nu in range(mu + 1, 4):
                    specs.append((kind, mu, nu))
                    specs.append((kind, nu, mu))
        else:                          # 'S'
            specs = [('S', mu, nu) for mu in range(4) for nu in range(mu + 1, 4)]
        for kindname, mu, nu in specs:
            steps_local = _loop_coords(kindname, mu, nu)
            # when mu>nu (second orientation) _loop_coords still builds the
            # correct rectangle with extents along mu/nu
            for s in range(N):
                x0 = s // L ** 3
                x1 = (s // L ** 2) % L
                x2 = (s // L) % L
                x3 = s % L
                base = np.array([x0, x1, x2, x3])
                links = []
                for off, d, conj in steps_local:
                    pos = tuple((base + np.array(off)) % L)
                    sf = ((pos[0] * L + pos[1]) * L + pos[2]) * L + pos[3]
                    links.append((sf * 4 + d, conj))
                loop_defs.append((ki, links))

    nloops = len(loop_defs)
    maxnt = max(len(l) for _, l in loop_defs)
    nterms = sum(len(l) for _, l in loop_defs)
    t_first = np.empty(nterms, np.int64)
    t_nt = np.empty(nterms, np.int64)
    t_pos = np.empty(nterms, np.int64)
    t_kind = np.empty(nterms, np.int64)
    t_code = np.empty(nterms, np.int64)
    slot_link = np.empty(nterms * maxnt, np.int64)
    slot_conj = np.empty(nterms * maxnt, np.int64)
    ti = 0
    si = 0
    per_link = [[] for _ in range(4 * N)]
    for ki, links in loop_defs:
        nt = len(links)
        for pos in range(nt):
            t_first[ti] = si
            t_nt[ti] = nt
            t_pos[ti] = pos
            t_kind[ti] = ki
            t_code[ti] = CHAR_OF.get(kinds[ki], 0)
            for k, (lk, cj) in enumerate(links):
                slot_link[si] = lk
                slot_conj[si] = cj
                si += 1
            per_link[links[pos][0]].append(ti)
            ti += 1
    slot_link = slot_link[:si].copy()
    slot_conj = slot_conj[:si].copy()
    link_off = np.zeros(4 * N + 1, np.int64)
    for l in range(4 * N):
        link_off[l + 1] = link_off[l] + len(per_link[l])
    link_term = np.empty(link_off[-1], np.int64)
    for l in range(4 * N):
        for k, t in enumerate(per_link[l]):
            link_term[link_off[l] + k] = t

    # checkerboard order: (parity, direction) -> links
    order_idx = np.full((8, 4 * N), -1, np.int64)
    order_cnt = np.zeros(8, np.int64)
    for l in range(4 * N):
        s = l // 4
        d = l % 4
        x0 = s // L ** 3
        x1 = (s // L ** 2) % L
        x2 = (s // L) % L
        x3 = s % L
        c = (x0 + x1 + x2 + x3) % 2
        row = c * 4 + d
        order_idx[row, order_cnt[row]] = l
        order_cnt[row] += 1

    return dict(t_first=t_first, t_nt=t_nt, t_pos=t_pos, t_kind=t_kind,
                t_code=t_code,
                slot_link=slot_link, slot_conj=slot_conj, maxnt=maxnt,
                link_off=link_off, link_term=link_term,
                order_idx=order_idx, order_cnt=order_cnt, nlinks=4 * N)


# ---------------------------------------------------------------- 3x3 helpers
@njit(cache=True, inline='always', fastmath=True)
def _mul(A, B, C):
    for i in range(3):
        for j in range(3):
            s = 0.0 + 0.0j
            for k in range(3):
                s += A[i, k] * B[k, j]
            C[i, j] = s


@njit(cache=True, inline='always', fastmath=True)
def _adj_mul(A, B, C):      # C = A^dag @ B
    for i in range(3):
        for j in range(3):
            s = 0.0 + 0.0j
            for k in range(3):
                s += np.conj(A[k, i]) * B[k, j]
            C[i, j] = s


@njit(cache=True, inline='always', fastmath=True)
def _mul_adj(A, B, C):      # C = A @ B^dag
    for i in range(3):
        for j in range(3):
            s = 0.0 + 0.0j
            for k in range(3):
                s += A[i, k] * np.conj(B[j, k])
            C[i, j] = s


@njit(cache=True, inline='always', fastmath=True)
def _tr3(A):
    return (A[0, 0] + A[1, 1] + A[2, 2]).real / 3.0


@njit(cache=True, inline='always', fastmath=True)
def _conj3(A, C):
    for i in range(3):
        for j in range(3):
            C[i, j] = np.conj(A[j, i])


@njit(cache=True, inline='always', fastmath=True)
def _left_apply(v, i, j, M):
    """M <- g M with g the 2x2 SU(2) v embedded at rows/cols (i,j)."""
    for k in range(3):
        a = M[i, k]
        b = M[j, k]
        M[i, k] = v[0, 0] * a + v[0, 1] * b
        M[j, k] = v[1, 0] * a + v[1, 1] * b


@njit(cache=True, inline='always', fastmath=True)
def _emb_left(v, i, j, M, C):
    """C = g M, g embedded."""
    for k in range(3):
        a = M[i, k]
        b = M[j, k]
        C[i, k] = v[0, 0] * a + v[0, 1] * b
        C[j, k] = v[1, 0] * a + v[1, 1] * b
    for k in range(3):
        for r in range(3):
            if r != i and r != j:
                C[r, k] = M[r, k]


@njit(cache=True, inline='always', fastmath=True)
def _emb_right_adj(v, i, j, M, C):
    """C = M g^dag, g the 2x2 SU(2) v embedded at rows/cols (i,j)."""
    for r in range(3):
        a = M[r, i]
        b = M[r, j]
        C[r, i] = a * np.conj(v[0, 0]) + b * np.conj(v[1, 0])
        C[r, j] = a * np.conj(v[0, 1]) + b * np.conj(v[1, 1])
    for r in range(3):
        for k in range(3):
            if k != i and k != j:
                C[r, k] = M[r, k]


@njit(cache=True, inline='always')
def _rand_su2(width):
    """SU(2) proposal: Haar for width=1.0, narrow around identity for width<1."""
    b = np.random.normal(0.0, 1.0) * width
    c = np.random.normal(0.0, 1.0) * width
    d = np.random.normal(0.0, 1.0) * width
    r2 = b * b + c * c + d * d
    if r2 < 1.0:
        a = np.sqrt(1.0 - r2)
    else:
        n = np.sqrt(r2)
        a = 0.0
        b /= n; c /= n; d /= n
    v = np.empty((2, 2), np.complex128)
    v[0, 0] = a + 1j * d
    v[0, 1] = c + 1j * b
    v[1, 0] = -c + 1j * b
    v[1, 1] = a - 1j * d
    return v



@njit(cache=True, inline='always', fastmath=True)
def _term_value(A, code, s1, s2):
    """Term value: (1/3)Re tr for a loop family (code 0), character chi_R otherwise."""
    if code == 0:
        return _tr3(A)
    t = A[0, 0] + A[1, 1] + A[2, 2]
    _mul(A, A, s1)                       # s1 = A^2
    t2 = s1[0, 0] + s1[1, 1] + s1[2, 2]
    if code == 6:
        return (0.5 * (t * t + t2)).real
    if code == 8:
        return (t * np.conj(t)).real - 1.0
    _mul(A, s1, s2)                      # s2 = A^3 (no aliasing)
    t3 = s2[0, 0] + s2[1, 1] + s2[2, 2]
    return ((t * t * t + 3.0 * t * t2 + 2.0 * t3) / 6.0).real


@njit(cache=True, inline='always', fastmath=True)
def _trace_term(Uf, t, base, nt, pos, X, slot_link, slot_conj, tmpA, tmpB, tmpC,
                code):
    """(1/3) Re tr of the loop product with the target link replaced by X.

    The loop product is M_0 M_1 ... M_{nt-1} (left to right along the path);
    at the target slot the matrix U is replaced by X (X^dag if the slot is a
    conjugated link).
    """
    for i in range(3):
        for j in range(3):
            tmpA[i, j] = 1.0 if i == j else 0.0
    for k in range(nt):
        if k == pos:
            if slot_conj[base + k] == 1:
                _mul_adj(tmpA, X, tmpB)      # tmpA @ X^dag
            else:
                _mul(tmpA, X, tmpB)          # tmpA @ X
            for i in range(3):
                for j in range(3):
                    tmpA[i, j] = tmpB[i, j]
        else:
            Lk = slot_link[base + k]
            if slot_conj[base + k] == 1:
                _mul_adj(tmpA, Uf[Lk], tmpB) # tmpA @ U^dag
            else:
                _mul(tmpA, Uf[Lk], tmpB)
            for i in range(3):
                for j in range(3):
                    tmpA[i, j] = tmpB[i, j]
    return _term_value(tmpA, code, tmpB, tmpC)


@njit(cache=True, fastmath=True)
def sweep(Uf, w, t_first, t_nt, t_pos, t_kind, t_code, slot_link, slot_conj,
          link_off, link_term, order_idx, order_cnt, nacc, width):
    """One full Metropolis sweep: 2 parities x 4 directions, sequential."""
    maxt = 64
    Tcur = np.empty(maxt)
    Tnew = np.empty(maxt)
    tloc = np.empty(maxt, np.int64)
    tmpA = np.empty((3, 3), np.complex128)
    tmpB = np.empty((3, 3), np.complex128)
    tmpC = np.empty((3, 3), np.complex128)
    X = np.empty((3, 3), np.complex128)
    GX = np.empty((3, 3), np.complex128)
    acc = 0
    tot = 0
    for c in range(2):
        for d in range(4):
            row = c * 4 + d
            nk = order_cnt[row]
            for kk in range(nk):
                l = order_idx[row, kk]
                m = 0
                for ti in range(link_off[l], link_off[l + 1]):
                    t = link_term[ti]
                    tloc[m] = t
                    Tcur[m] = _trace_term(Uf, t, t_first[t], t_nt[t], t_pos[t],
                                          Uf[l], slot_link, slot_conj,
                                          tmpA, tmpB, tmpC, t_code[t])
                    m += 1
                for sg in range(3):
                    if sg == 0:
                        i = 0; j = 1
                    elif sg == 1:
                        i = 0; j = 2
                    else:
                        i = 1; j = 2
                    v = _rand_su2(width)
                    dS = 0.0
                    for a in range(m):
                        t = tloc[a]
                        base = t_first[t]
                        pos = t_pos[t]
                        _emb_left(v, i, j, Uf[l], X)         # X = g U (new link)
                        Tnew[a] = _trace_term(Uf, t, base, t_nt[t], pos, X,
                                              slot_link, slot_conj,
                                              tmpA, tmpB, tmpC, t_code[t])
                        dS += -w[t_kind[t]] * (Tnew[a] - Tcur[a])
                    tot += 1
                    if dS <= 0.0 or np.random.random() < np.exp(-dS):
                        _left_apply(v, i, j, Uf[l])
                        for a in range(m):
                            Tcur[a] = Tnew[a]
                        acc += 1
    nacc[0] = acc
    nacc[1] = tot


@njit(cache=True)
def _ml(x, d):
    return x % d


@njit(cache=True, fastmath=True)
def loop_densities(Uf, w, t_first, t_nt, t_pos, t_kind, t_code, slot_link,
                   slot_conj, link_off, link_term, kinds):
    """Mean (1/3)Re tr per loop, for each loop family (in order of kinds)."""
    nk = len(kinds)
    acc = np.zeros(nk)
    cnt = np.zeros(nk, np.int64)
    tmpA = np.empty((3, 3), np.complex128)
    tmpB = np.empty((3, 3), np.complex128)
    tmpC = np.empty((3, 3), np.complex128)
    # unique loops: take only the terms with position 0
    for t in range(len(t_first)):
        if t_pos[t] != 0:
            continue
        k = t_kind[t]
        X = Uf[slot_link[t_first[t]]]
        T = _trace_term(Uf, t, t_first[t], t_nt[t], 0, X,
                        slot_link, slot_conj, tmpA, tmpB, tmpC, t_code[t])
        acc[k] += T
        cnt[k] += 1
    out = np.empty(nk)
    for k in range(nk):
        out[k] = acc[k] / cnt[k]
    return out


@njit(cache=True, fastmath=True)
def debug_delta(Uf, l, v, i, j, w, t_first, t_nt, t_pos, t_kind, t_code,
                slot_link, slot_conj, link_off, link_term):
    """Delta S of the proposal Uf[l] -> g Uf[l] for one SU(2) subgroup."""
    tmpA = np.empty((3, 3), np.complex128)
    tmpB = np.empty((3, 3), np.complex128)
    tmpC = np.empty((3, 3), np.complex128)
    X = np.empty((3, 3), np.complex128)
    _emb_left(v, i, j, Uf[l], X)
    dS = 0.0
    for ti in range(link_off[l], link_off[l + 1]):
        t = link_term[ti]
        base = t_first[t]
        pos = t_pos[t]
        T0 = _trace_term(Uf, t, base, t_nt[t], pos, Uf[l], slot_link,
                         slot_conj, tmpA, tmpB, tmpC, t_code[t])
        T1 = _trace_term(Uf, t, base, t_nt[t], pos, X, slot_link,
                         slot_conj, tmpA, tmpB, tmpC, t_code[t])
        dS += -w[t_kind[t]] * (T1 - T0)
    return dS
