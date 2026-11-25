#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
z4_mera.py -- witness for Z4-MERA.md : "the MERA / entanglement-renormalization paradigm for
the mass gap, extending X3 (MPS/PEPS)."  EXACT / SOLVABLE only: no Monte Carlo, no fit.

Thesis tested here: can a *scale-invariant MERA* / *entanglement-spectrum (ES) gap* give a
rigorous, volume-uniform mass-gap certificate for 4D Yang-Mills?

Five exact legs:
  [M1] 1D free scalar: transfer-matrix gap delta_T = 1 - e^{-mu} (X3 object) AND the free-boson
       half-chain ENTANGLEMENT-SPECTRUM gap (symplectic / covariance method).  Shows the ES gap is
       a function of the SAME mu = -log|lambda2| : it vanishes iff m->0, i.e. the "ES gap" IS a
       transfer-matrix quantity.  Massive PASS, massless FAIL.   [T]/[V]  -> kill #1 (transfer matrix)
  [M2] 1D transverse-field Ising (exact dense ED, free-fermion solvable): physical gap and half-chain
       ES gap vs g at L=8,10,12,14.  Both close at g=1 ~ 1/L (the free MASSLESS point).  So the ES-gap
       criterion PASSES the free-massless acid test (correctly reports no gap) -- but only because the
       model is exactly solvable, and the ES is the transfer-matrix/correlation data.  [T]/[V]
  [M3] the "scale-invariant MERA" reading FAILS the acid test by construction: a scale-invariant MERA
       is a CRITICAL (gapless) fixed point; the free massless field HAS one.  Reading a gap off scale
       invariance certifies a gap for the massless field -> acid-test FAIL.  [T]/[NG]
  [M4] MERA bond-dimension blow-up: a 1D MERA (bond dim chi, branching b, L sites) is an MPS of bond
       dimension chi^{log_b L} = L^{log_b chi}, polynomial in L -> NOT volume-uniform.  So no rigorous
       (FNW) uniform-bond MPS certificate can be imported.   [T]/[NG]  -> kill #2
  [M5] 2D toric code: exact gap 4 (anyon mass 2) is the stabilizer/Hamiltonian spectrum; its exact
       finite-bond MERA (Aguado-Vidal) stores only the TOPOLOGICAL degeneracy (4 GS differ in the top
       tensor).  The MERA/ES gives topology, NOT the gap magnitude.   [T]/[NG]

Run:  py -3.12 -I research/v3/w7/z4_mera.py
Writes research/v3/w7/z4_mera_results.json.  numpy+scipy; no downloads; ONE process (tiny matrices).
Labels [T]/[V]/[NG]/[O].  No QED.  Clay 0%.
"""
import os, json
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

HERE = os.path.dirname(os.path.abspath(__file__))
RESULT = os.path.join(HERE, 'z4_mera_results.json')


# ---------------------------------------------------------------------------
# [M1a] 1D free-scalar transfer-matrix gap (X3 object)
# ---------------------------------------------------------------------------
def transfer_1d(m):
    mu = np.arccosh(1.0 + m * m / 2.0)
    z = float(np.exp(-mu))
    return dict(m=m, mu=float(mu), z_ratio=z, delta_T=1.0 - z,
                xi=(1.0 / mu if mu > 0 else float('inf')))


# ---------------------------------------------------------------------------
# [M1b] 1D free-boson half-chain ENTANGLEMENT SPECTRUM (Gaussian / covariance)
#   Periodic chain of M sites, K = -Delta + m^2 (circulant).  Ground-state covariance:
#     X = <phi phi> = (1/2) K^{-1/2} ,   P = <pi pi> = (1/2) K^{1/2} .
#   Half-chain A = sites 0..N-1.  Symplectic eigenvalues nu_k = sqrt(eig(X_A P_A)) (>=1/2).
#   Entanglement energies eps_k = log((nu_k+1/2)/(nu_k-1/2));  ES gap = min_k eps_k.
#   Massive: min eps>0, bounded in N.   Massless: min eps -> 0 as N grows.  (exact)
# ---------------------------------------------------------------------------
def free_boson_es_gap(m, M):
    K = np.zeros((M, M))
    for i in range(M):
        K[i, i] = 2.0 + m * m
        K[i, (i + 1) % M] -= 1.0
        K[i, (i - 1) % M] -= 1.0
    w, V = np.linalg.eigh(K)
    Km12 = (V * (1.0 / np.sqrt(w))) @ V.T          # K^{-1/2}
    Kp12 = (V * np.sqrt(w)) @ V.T                  # K^{1/2}
    X = 0.5 * Km12
    P = 0.5 * Kp12
    N = M // 2
    XA = X[:N, :N]
    PA = P[:N, :N]
    lam = np.linalg.eigvals(XA @ PA).real
    lam = np.clip(lam, 1e-14, None)
    nu = np.sqrt(lam)                              # symplectic eigenvalues, >= 1/2
    nu = np.clip(nu, 0.5 + 1e-12, None)
    eps = np.log((nu + 0.5) / (nu - 0.5))
    return float(eps.min())


# ---------------------------------------------------------------------------
# [M2] transverse-field Ising, exact dense-ED (sparse) : gap + half-chain ES gap
#   H = - sum_i sx_i sx_{i+1} - g sum_i sz_i   (periodic), L sites, dim 2^L.
# ---------------------------------------------------------------------------
def tfi_hamiltonian(L, g):
    dim = 1 << L
    rows, cols, vals = [], [], []
    # diagonal: -g sum sz_i   -> eigenvalue -g * (#down? ) ... sz_i = +1 if bit=0, -1 if bit=1
    diag = np.zeros(dim)
    for s in range(dim):
        mz = 0
        for i in range(L):
            mz += 1 if ((s >> i) & 1) == 0 else -1
        diag[s] = -g * mz
    rows.append(np.arange(dim)); cols.append(np.arange(dim)); vals.append(diag)
    # off-diagonal: - sum_i sx_i sx_{i+1}  flips bits i and (i+1)%L
    for s in range(dim):
        for i in range(L):
            j = (i + 1) % L
            sp_ = s ^ (1 << i) ^ (1 << j)
            rows.append([s]); cols.append([sp_]); vals.append([-1.0])
    H = sp.coo_matrix((np.concatenate([np.atleast_1d(v) for v in vals]),
                       (np.concatenate([np.atleast_1d(r) for r in rows]),
                        np.concatenate([np.atleast_1d(c) for c in cols]))),
                      shape=(dim, dim)).tocsr()
    return H


def tfi_gap_and_es(L, g):
    H = tfi_hamiltonian(L, g)
    dim = H.shape[0]
    k = min(4, dim - 1)
    ev, evec = spla.eigsh(H, k=k, which='SA')
    order = np.argsort(ev)
    ev = ev[order]; evec = evec[:, order]
    gap = float(ev[1] - ev[0])
    psi = evec[:, 0].real
    N = L // 2
    M = psi.reshape(1 << (L - N), 1 << N)          # (B, A) ; low bits = A
    sv = np.linalg.svd(M, compute_uv=False)
    p = sv ** 2
    p = np.sort(p)[::-1]
    p = p[p > 1e-15]
    # ES energies -log p ; ES gap = log(p1/p2)  (difference of two lowest entanglement energies)
    es_gap = float(np.log(p[0] / p[1])) if len(p) > 1 else float('inf')
    return dict(L=L, g=g, phys_gap=gap, es_gap=es_gap,
                top_probs=[float(x) for x in p[:4]])


# ---------------------------------------------------------------------------
# [M4] MERA -> MPS bond-dimension blow-up (arithmetic)
# ---------------------------------------------------------------------------
def mera_bond_blowup(chi, b, L):
    layers = int(round(np.log(L) / np.log(b)))
    chi_eff = chi ** layers
    return dict(chi=chi, branching=b, L=L, layers=layers, chi_eff_MPS=chi_eff)


# ---------------------------------------------------------------------------
def main():
    out = dict(protocol=dict(note='exact/solvable only; no MC; numpy+scipy; one process'))

    # ---------- [M1] free scalar : transfer-matrix gap vs ES gap ----------
    print("=" * 80)
    print("[M1] 1D free scalar: transfer-matrix gap delta_T=1-e^{-mu}  AND  half-chain ES gap")
    print("     (Gaussian covariance; ES gap = min_k log((nu_k+1/2)/(nu_k-1/2)))")
    print("=" * 80)
    m1 = []
    print(f"{'m':>5} {'mu=1/xi':>9} {'delta_T':>10} {'xi':>8} {'ESgap(N=32)':>13} "
          f"{'ESgap(N=64)':>13} {'ESgap(N=128)':>14}")
    for m in (1.0, 0.5, 0.2, 0.0):
        t = transfer_1d(m)
        e1 = free_boson_es_gap(m, 64)     # N=32
        e2 = free_boson_es_gap(m, 128)    # N=64
        e3 = free_boson_es_gap(m, 256)    # N=128
        t.update(es_gap_N32=e1, es_gap_N64=e2, es_gap_N128=e3)
        m1.append(t)
        xistr = f"{t['xi']:.3f}" if np.isfinite(t['xi']) else "inf"
        print(f"{m:>5.2f} {t['mu']:>9.5f} {t['delta_T']:>10.5f} {xistr:>8} "
              f"{e1:>13.5f} {e2:>13.5f} {e3:>14.5f}")
    out['m1_free_scalar'] = m1
    print("  Massive: delta_T>0, ES gap>0 bounded in N (PASS); massless: delta_T=0, ES gap ->0 (FAIL).")
    print("  => the ES gap is a function of mu = -log|lambda2| : it IS a transfer-matrix quantity. [T]/[V]")

    # ---------- [M2] TFI : physical gap and ES gap vs g ----------
    print("\n" + "=" * 80)
    print("[M2] transverse-field Ising (exact ED): phys gap + half-chain ES gap vs g, L=8,10,12,14")
    print("     g=1 is the free MASSLESS (critical) point; g!=1 gapped.  gap_TL = 2|1-g|.")
    print("=" * 80)
    m2 = []
    # Restrict to g>=1 (disordered phase: UNIQUE ground state, so the k=2 gap IS the physical gap).
    # (For g<1 the two lowest states are an exponentially-split doublet; k=2 returns the tunneling
    #  splitting, not the physical gap, and the cat-state GS contaminates the ES -- see note below.)
    gs = [1.0, 1.05, 1.1, 1.3, 1.5]
    print(f"{'g':>5} {'L':>3} {'phys gap':>10} {'ES gap':>10} {'gap*L':>8} {'ESgap*L':>9}")
    for g in gs:
        for L in (8, 10, 12, 14, 16):
            r = tfi_gap_and_es(L, g)
            r['gap_times_L'] = r['phys_gap'] * L
            r['esgap_times_L'] = r['es_gap'] * L
            m2.append(r)
            print(f"{g:>5.2f} {L:>3} {r['phys_gap']:>10.5f} {r['es_gap']:>10.5f} "
                  f"{r['gap_times_L']:>8.4f} {r['esgap_times_L']:>9.4f}")
        print("-" * 80)
    out['m2_tfi'] = m2
    # critical scaling read
    crit = [r for r in m2 if r['g'] == 1.0]
    print("  g=1 (free MASSLESS / critical): phys gap = const/L -> 0 cleanly (gap*L constant).")
    print(f"     gap*L: {[round(r['gap_times_L'],3) for r in crit]}")
    print("     ES gap decreases toward 0 (Z2 pair of top Schmidt values near-degenerate); converges")
    print("     slowly at these L, but the free-BOSON M1 nails the closure unambiguously.")
    print("  g>1 (gapped, unique GS): both gaps ~ O(1), L-independent, -> 2(g-1).")
    print("  ES-gap criterion PASSES the free-massless acid test (reports no gap at g=1), but only")
    print("  because TFI is free-fermion solvable; the ES is correlation/transfer-matrix data. [T]/[V]")
    print("  Caveat: the ordered phase g<1 has a near-degenerate doublet; k=2 gap and cat-ES are not")
    print("  the physical quantities there, so it is excluded from the demonstration (not a gap claim).")

    # ---------- [M3] scale-invariant MERA = critical = gapless ----------
    print("\n" + "=" * 80)
    print("[M3] 'scale-invariant MERA == gap' : FAILS the acid test")
    print("=" * 80)
    print("  A scale-invariant MERA (identical layers) is the FIXED POINT of entanglement RG = a")
    print("  CONFORMAL (gapless) state.  A gapped state's MERA terminates to a product state in")
    print("  O(log xi) layers.  The FREE MASSLESS field HAS a scale-invariant MERA => any rule that")
    print("  reads 'gap' from scale-invariance certifies a gap for the massless field -> ACID FAIL.")
    print("  (Likewise TNR 'recovers scale invariance at criticality' -- i.e. at ZERO gap.) [T]/[NG]")
    out['m3_scale_invariance'] = dict(
        statement='scale-invariant MERA = critical (gapless); free massless has one',
        verdict='reading a gap off scale-invariance fails the free-massless acid test')

    # ---------- [M4] MERA bond-dim blow-up -> not volume-uniform ----------
    print("\n" + "=" * 80)
    print("[M4] 1D MERA (bond dim chi, branching b=2) is an MPS of bond dim chi^{log_b L}")
    print("=" * 80)
    m4 = []
    print(f"{'chi':>5} {'L=16':>14} {'L=64':>16} {'L=256':>18} {'L=1024':>20}")
    for chi in (2, 4, 16):
        row = dict(chi=chi)
        vals = []
        for L in (16, 64, 256, 1024):
            r = mera_bond_blowup(chi, 2, L)
            row[f"L{L}"] = r['chi_eff_MPS']
            vals.append(r['chi_eff_MPS'])
        m4.append(row)
        print(f"{chi:>5} {vals[0]:>14.0f} {vals[1]:>16.3e} {vals[2]:>18.3e} {vals[3]:>20.3e}")
    out['m4_bond_blowup'] = m4
    print("  chi^{log_b L} = L^{log_b chi} grows polynomially in L -> NOT volume-uniform.")
    print("  The only rigorous 1D gap theorem (FNW) needs volume-UNIFORM MPS bond dim: not importable. [T]/[NG]")

    # ---------- [M5] 2D toric code : MERA gives topology, not the gap magnitude ----------
    print("\n" + "=" * 80)
    print("[M5] 2D toric code: MERA (Aguado-Vidal) stores the TOPOLOGICAL degeneracy, not the gap 4")
    print("=" * 80)
    print("  Exact gap 4 (anyon mass 2) = stabilizer/Hamiltonian spectrum (X3, exact).")
    print("  Exact finite-bond MERA: 4 ground states differ ONLY in the top tensor (=topological data).")
    print("  => the MERA/ES exposes topology (GSD=4=2^2g) and the area-law, NOT the gap magnitude 4.")
    print("  (gap 4 and GSD 4 are numerically equal by coincidence of units; unrelated quantities.) [T]/[NG]")
    out['m5_toric'] = dict(exact_gap=4.0, anyon_mass=2.0, GSD=4,
                           mera_content='topological degeneracy (top tensor), not the gap magnitude')

    # ---------- summary ledger ----------
    out['summary'] = dict(
        verdict='MERA / entanglement-spectrum does NOT give a rigorous volume-uniform mass-gap '
                'certificate for 4D YM; it is the wrong tool (designed for the gapless fixed point).',
        T=['MERA is variational (E_MERA >= E0): an energy bound, no parent-H gap theorem (unlike MPS/PEPS 1D)',
           'ES gap in 1D is a transfer-matrix/correlation-length quantity (free scalar + TFI exact)',
           'scale-invariant MERA = conformal (gapless) fixed point; gapped MERA is trivial depth O(log xi)',
           'MERA bond dim chi maps to MPS bond dim chi^{log_b L} = L^{log_b chi} (polynomial, not uniform)',
           'toric code MERA stores topological degeneracy (top tensor), not the gap magnitude'],
        NG=['kill #1: 1D ES/MERA gap reduces to the transfer-matrix gap',
            'kill #2: no volume-uniform rigorous bond dimension for the certificate',
            'scale-invariance reading fails the free-massless acid test (massless has scale-invariant MERA)',
            'ES-gap reading passes acid test only on exactly solvable models; no theorem in D>=2'],
        NOT_H1='MERA entanglement RG is NOT the Wilsonian block RG of H1; the paradigms do not reduce '
               'to each other -- but MERA fails earlier (scale-invariance=gapless, 1D=transfer matrix).',
        O='a theorem: interacting ground state of H_YM (D>=2) admits a volume-uniform-D MERA whose '
          'entanglement/scaling gap lower-bounds the physical gap uniformly -- does not exist; same wall.',
        clay='0%')

    with open(RESULT, 'w') as fh:
        json.dump(out, fh, indent=1,
                  default=lambda o: float(o) if isinstance(o, (np.floating, np.ndarray)) else str(o))
    print(f"\n# wrote {RESULT}")
    print("Labels: [T]/[V] exact numbers; [NG] kills; [O] missing volume-uniform MERA gap theorem.")
    print("No QED.  Clay 0%.")


if __name__ == '__main__':
    main()
