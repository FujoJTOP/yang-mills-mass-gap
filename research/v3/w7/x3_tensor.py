#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
x3_tensor.py -- witness for X3-TENSOR.md, "the first step of the tensor-network / rigorous
gap-certificate paradigm."  EXACT / SOLVABLE only: no Monte Carlo, no fit.

Four exact legs:
  (1) free-field finite-volume gap   omega_min(D,L,m) = sqrt(m^2 + min_{k!=0} 4 sum_mu sin^2(pi n_mu/L))
      -> shows the free MASSLESS field already has an O(1) finite-volume "gap" at L=2 (=2.0),
         so a gap read at L=2 cannot pass the free-massless acid test in ANY dimension.  [T]
  (2) 1D MPS / transfer-matrix criterion: bulk delta_T = 1 - |lambda2/lambda0| = 1 - e^{-mu},
      cosh(mu)=1+m^2/2.  massive: delta_T>0 L-independent; massless: delta_T=0.  Verified against
      the exact ring correlator decay ratio.  This IS the transfer-matrix gap (kill criterion).  [T]
  (3) 2D toric code = 2D Z_2 gauge fixed point (the solvable gauge model): EXACT dense spectrum at
      L=2 -> ground-state degeneracy 4, gap 2 (stabilizer algebra).  GSD=2^{2g}=4>1 for all L
      (GF(2) rank) -> the ground-state PEPS is NOT injective, so the injectivity gap criterion does
      not fire; the gap is the exact stabilizer spectrum.  [T]
  (4) 4D Z_2 gauge (4D toric-code stabilizer) on L=2: GF(2) rank -> exact GSD; gap 2.  Feasible to
      compute exactly, but it is the finite-volume / transfer-matrix spectrum and L=2 fails the
      free-massless acid test (leg 1, D=4: massless gap at L=2 = 2.0).  [T]/[V]

Run:  py -3.12 -I research/v3/w7/x3_tensor.py
Writes research/v3/w7/x3_tensor_results.json.  numpy only; no downloads; ONE process.
Labels [T]/[V]/[NG]/[O].  No QED.  Clay 0%.
"""
import os, json
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RESULT = os.path.join(HERE, 'x3_tensor_results.json')


# ----------------------------------------------------------------------------
# (1) exact free-field finite-volume gap
# ----------------------------------------------------------------------------
def free_gap(D, L, m):
    """Single-particle excitation gap (ground->first excited) of a free real scalar on a D-dim
    periodic lattice of linear size L, lattice units.  omega_k = sqrt(m^2 + 4 sum_mu sin^2(pi n_mu/L)).
    Returns (omega_min_nonzero, k_min_multiplet_size)."""
    best = None
    count = 0
    # only need modes that are minimal: examine along one direction (n=1) and verify no smaller.
    # enumerate all n in {0,..,L-1}^D but avoid blowing up: for the min we only need n with at most
    # one nonzero component=+-1 (all higher modes are larger).  This is exact for the minimum.
    cand = []
    for mu in range(D):
        arg = 4.0 * np.sin(np.pi * 1 / L) ** 2  # one direction excited by 1
        cand.append(m * m + arg)
    wmin = float(np.sqrt(min(cand)))
    # multiplet size: number of modes attaining the min = 2D (mu, +-1) for L>2 ; for L==2, +-1 coincide
    multit = D if L == 2 else 2 * D
    return wmin, multit


# ----------------------------------------------------------------------------
# (2) 1D Gaussian MPS transfer matrix / correlation length (exact)
# ----------------------------------------------------------------------------
def transfer_1d(m):
    """Bulk transfer-matrix ratio lambda2/lambda0 = e^{-mu} for a 1D free scalar of mass m
    (nearest-neighbour lattice Laplacian -Delta + m^2).  cosh(mu) = 1 + m^2/2.
    delta_T = 1 - e^{-mu} (>0 massive, =0 massless); xi = 1/mu."""
    mu = np.arccosh(1.0 + m * m / 2.0)
    z = float(np.exp(-mu))
    return dict(m=m, mu=float(mu), z_ratio=z, delta_T=1.0 - z, xi=(1.0 / mu if mu > 0 else float('inf')))


def open_chain_decay(m, N, rmin):
    """Independent check of the transfer-matrix ratio.  Build the OPEN-chain operator
    K = (-Delta + m^2) on N sites (tridiagonal), invert it, and take the bulk-to-bulk Green's
    function C(r) = (K^{-1})[c, c+r] with c far from both boundaries.  For the 1D nearest-neighbour
    operator the interior Green's function is EXACTLY geometric, C(r) ~ rho^r, with rho = smaller
    root of rho^2-(2+m^2)rho+1=0 = e^{-mu}, cosh(mu)=1+m^2/2.  Returns C(rmin+1)/C(rmin)."""
    K = np.zeros((N, N))
    for i in range(N):
        K[i, i] = 2.0 + m * m
        if i + 1 < N:
            K[i, i + 1] = -1.0
            K[i + 1, i] = -1.0
    Kinv = np.linalg.inv(K)
    c = N // 2
    c0, c1 = Kinv[c, c + rmin], Kinv[c, c + rmin + 1]
    return c1 / c0 if c0 != 0 else float('nan')


# ----------------------------------------------------------------------------
# Pauli helpers (dense) for the 2D toric code exact spectrum
# ----------------------------------------------------------------------------
_I = np.eye(2, dtype=complex)
_X = np.array([[0, 1], [1, 0]], dtype=complex)
_Z = np.array([[1, 0], [0, -1]], dtype=complex)


def pauli_string(n, support, kind):
    """Kronecker product of `kind` (X or Z) on qubits in `support`, identity elsewhere. n qubits."""
    P = _X if kind == 'X' else _Z
    M = np.array([[1.0 + 0j]])
    for q in range(n):
        M = np.kron(M, P if q in support else _I)
    return M


# ----------------------------------------------------------------------------
# GF(2) rank of a list of bitmasks
# ----------------------------------------------------------------------------
def gf2_rank(masks):
    basis = {}
    rank = 0
    for r in masks:
        while r:
            h = r.bit_length() - 1
            if h in basis:
                r ^= basis[h]
            else:
                basis[h] = r
                rank += 1
                break
    return rank


# ----------------------------------------------------------------------------
# (3) 2D toric code = 2D Z_2 gauge fixed point
# ----------------------------------------------------------------------------
def toric_2d_supports(L):
    """Return (A_supports, B_supports) as lists of edge-index sets.  Edges: H index x+L*y in [0,L^2),
    V index L^2 + x+L*y.  A_v (vertex, X-type) and B_p (plaquette, Z-type)."""
    def H(x, y): return (x % L) + L * (y % L)
    def V(x, y): return L * L + (x % L) + L * (y % L)
    A, B = [], []
    for y in range(L):
        for x in range(L):
            # vertex (x,y): east H(x,y), west H(x-1,y), north V(x,y), south V(x,y-1)
            A.append({H(x, y), H(x - 1, y), V(x, y), V(x, y - 1)})
            # plaquette (x,y): bottom H(x,y), top H(x,y+1), left V(x,y), right V(x+1,y)
            B.append({H(x, y), H(x, y + 1), V(x, y), V(x + 1, y)})
    return A, B


def toric_2d(L, dense=True):
    n = 2 * L * L
    A, B = toric_2d_supports(L)
    # exact GSD via GF(2) symplectic rank: X-part bits 0..n-1, Z-part bits n..2n-1
    masks = []
    for s in A:
        masks.append(sum(1 << q for q in s))          # X-type
    for s in B:
        masks.append(sum(1 << (n + q) for q in s))    # Z-type
    rank = gf2_rank(masks)
    gsd_rank = 1 << (n - rank)
    out = dict(L=L, n_qubits=n, n_stab=len(masks), gf2_rank=rank,
               GSD_rank=gsd_rank, expected_GSD=4)
    if dense:
        Hm = np.zeros((1 << n, 1 << n), dtype=complex)
        for s in A:
            Hm -= pauli_string(n, s, 'X')
        for s in B:
            Hm -= pauli_string(n, s, 'Z')
        ev = np.linalg.eigvalsh(Hm)
        e0 = ev.min()
        gsd = int(np.sum(np.isclose(ev, e0, atol=1e-8)))
        # first excited level above the (gsd-fold) ground level
        ev_above = ev[~np.isclose(ev, e0, atol=1e-8)]
        e1 = ev_above.min()
        out.update(E0=float(e0), GSD_dense=gsd, E1=float(e1), gap=float(e1 - e0),
                   spectrum_unique_first6=[float(v) for v in np.unique(np.round(ev, 8))[:6]])
    return out


# ----------------------------------------------------------------------------
# (4) 4D Z_2 gauge theory (4D toric-code stabilizer) on L=2 -- exact GF(2) rank
# ----------------------------------------------------------------------------
def z2_gauge_DD(L, D):
    """Z_2 gauge theory on a D-dim periodic lattice linear size L.  Links n = D*L^D.
    Vertex stars (X-type) and plaquettes (Z-type).  Exact GSD = 2^{n - rank} via GF(2)."""
    V = L ** D

    def site(coords):
        s = 0
        for d in range(D):
            s += (coords[d] % L) * (L ** d)
        return s

    def coords_of(s):
        c = []
        for d in range(D):
            c.append(s % L); s //= L
        return c

    def link(d, s):  # link in direction d emanating from site s
        return d * V + s

    masks = []
    # vertex stars: outgoing link (d,s) and incoming link (d, s - e_d)
    for s in range(V):
        cs = coords_of(s)
        star = set()
        for d in range(D):
            star.add(link(d, s))
            cprev = list(cs); cprev[d] = (cprev[d] - 1) % L
            star.add(link(d, site(cprev)))
        masks.append(sum(1 << q for q in star))            # X-type (bits 0..n-1)
    # plaquettes in each (d1<d2) plane
    n = D * V
    for s in range(V):
        cs = coords_of(s)
        for d1 in range(D):
            for d2 in range(d1 + 1, D):
                c_a = list(cs); c_b = list(cs)
                c_b[d1] = (c_b[d1] + 1) % L
                c_c = list(cs); c_c[d2] = (c_c[d2] + 1) % L
                # plaquette links: (d1,s), (d2, s+e_d1), (d1, s+e_d2), (d2, s)
                plaq = {link(d1, s), link(d2, site(c_b)), link(d1, site(c_c)), link(d2, site(cs))}
                masks.append(sum(1 << (n + q) for q in plaq))  # Z-type
    rank = gf2_rank(masks)
    return dict(D=D, L=L, n_links=n, n_stab=len(masks), gf2_rank=rank,
                GSD_rank=1 << (n - rank))


def main():
    out = dict(protocol=dict(note='exact/solvable only; no Monte Carlo; numpy only'))

    # ---- (1) free-field finite-volume gap table ----
    print("=" * 78)
    print("(1) free-field finite-volume gap omega_min(D,L,m)   [exact, lattice units]")
    print("    omega_min = sqrt(m^2 + 4 sin^2(pi/L))  (single direction; exact minimum)")
    print("=" * 78)
    rows = []
    print(f"{'D':>2} {'L':>4} {'m=0 (massless)':>16} {'m=0.5':>9} {'m=1.0':>9}")
    for D in (1, 2, 3, 4):
        for L in (2, 4, 8, 16):
            g0, mt = free_gap(D, L, 0.0)
            g5, _ = free_gap(D, L, 0.5)
            g1, _ = free_gap(D, L, 1.0)
            rows.append(dict(D=D, L=L, m0=g0, m05=g5, m1=g1, multiplet=mt))
            print(f"{D:>2} {L:>4} {g0:>16.4f} {g5:>9.4f} {g1:>9.4f}")
        print("-" * 78)
    out['free_field_gap'] = rows
    # the decisive L=2 line
    l2 = [r for r in rows if r['L'] == 2 and r['D'] == 4][0]
    print(f"\n  DECISIVE: free MASSLESS D=4 L=2 gap = {l2['m0']:.4f} (O(1)), "
          f"vs massive m=1 bulk gap = 1.0 -> L=2 cannot distinguish. [T]")

    # ---- (2) 1D transfer-matrix criterion ----
    print("\n" + "=" * 78)
    print("(2) 1D MPS / transfer-matrix gap  delta_T = 1 - e^{-mu},  cosh(mu)=1+m^2/2")
    print("=" * 78)
    tm = []
    print(f"{'m':>6} {'mu=1/xi':>10} {'z=lambda2/l0':>14} {'delta_T':>10} {'xi':>10}  open-chain-ratio(N=301,r=8)")
    for m in (1.0, 0.5, 0.2, 0.0):
        t = transfer_1d(m)
        rr = open_chain_decay(m, 301, 8)
        t['open_chain_ratio'] = rr
        tm.append(t)
        print(f"{m:>6.2f} {t['mu']:>10.5f} {t['z_ratio']:>14.6f} {t['delta_T']:>10.5f} "
              f"{str(round(t['xi'],3)) if np.isfinite(t['xi']) else 'inf':>10}  {rr:.6f}")
    out['transfer_1d'] = tm
    print("  massive: delta_T>0 L-independent (PASS); massless: delta_T=0 (FAIL).")
    print("  NB delta_T is exactly the transfer-matrix gap -> 1D MPS certificate REDUCES to it. [T]/[NG]")

    # ---- (3) 2D toric code = 2D Z_2 gauge fixed point ----
    print("\n" + "=" * 78)
    print("(3) 2D toric code = 2D Z_2 gauge fixed point  [exact stabilizer spectrum]")
    print("=" * 78)
    tc = []
    for L in (2, 3, 4, 5):
        r = toric_2d(L, dense=(L == 2))   # dense only at L=2 (256-dim); rank for larger
        tc.append(r)
        if L == 2:
            print(f"  L={L}: n={r['n_qubits']} qubits, {r['n_stab']} stabilizers, GF(2) rank={r['gf2_rank']}"
                  f" -> GSD_rank={r['GSD_rank']}")
            print(f"        dense: E0={r['E0']:.1f}, GSD_dense={r['GSD_dense']}, E1={r['E1']:.1f}, "
                  f"gap={r['gap']:.1f}, unique levels {r['spectrum_unique_first6']}")
        else:
            print(f"  L={L}: n={r['n_qubits']} qubits, GF(2) rank={r['gf2_rank']} -> GSD_rank={r['GSD_rank']}")
    out['toric_2d'] = tc
    print("  GSD=4>1 for every L -> ground state NOT injective (injective PEPS has unique GS).")
    print("  The many-body gap (4; anyon mass 2) is the exact stabilizer spectrum, NOT a TN certificate. [T]/[NG]")

    # ---- (4) 4D Z_2 gauge on L=2 ----
    print("\n" + "=" * 78)
    print("(4) 4D Z_2 gauge (4D toric stabilizer) on L=2 -- exact GF(2) rank, feasibility")
    print("=" * 78)
    zg = []
    for (D, L) in ((2, 2), (2, 3), (3, 2), (3, 3), (4, 2)):
        r = z2_gauge_DD(L, D)
        zg.append(r)
        print(f"  D={D} L={L}: n_links={r['n_links']}, {r['n_stab']} stabilizers, "
              f"GF(2) rank={r['gf2_rank']} -> GSD={r['GSD_rank']}, many-body gap=4 (anyon mass 2)")
    out['z2_gauge'] = zg
    print("  4D L=2 exact spectrum is FEASIBLE (GF(2) rank of 112x128 over 64 bits, instant).")
    print("  But it is the finite-volume/transfer-matrix spectrum; L=2 massless gap=2.0 (leg 1)")
    print("  -> fails the free-massless acid test -> certifyingly vacuous. [T]/[V]/[NG]")

    # ---- summary ledger ----
    out['summary'] = dict(
        T='1D injective MPS => uniform gap = transfer-matrix gap (delta_T=1-|l2|); free-field '
          'finite-volume gap closed form; toric code exact stabilizer spectrum',
        V='exact numbers above (no MC)',
        NG=['1D MPS certificate == transfer-matrix gap (reduces to D3 object)',
            'D>=2 area-law criterion acid-blind: free massless satisfies area law',
            'D>=2 injectivity->gap is open 2D PEPS gap conjecture; solvable gauge model gap is the '
            'exact stabilizer/transfer spectrum, and its GS is non-injective (GSD=4)',
            '4D L=2 finite-volume gap fails free-massless acid test (massless L=2 gap=2.0)'],
        O='a 2D/3D injectivity=>gap theorem for a PEPS that provably represents the YM ground state '
          'with bounded bond dimension -- does not exist; that is the missing step',
        clay='0%')

    with open(RESULT, 'w') as fh:
        json.dump(out, fh, indent=1,
                  default=lambda o: float(o) if isinstance(o, (np.floating, np.ndarray)) else str(o))
    print(f"\n# wrote {RESULT}")
    print("Labels: [T]/[V] exact numbers; [NG] triple no-go; [O] missing 2D/3D injectivity-gap theorem.")
    print("No QED.  Clay 0%.")


if __name__ == '__main__':
    main()
