#!/usr/bin/env python3
"""na2_selftest.py -- unit tests of the analytic machinery in na2_su3_na.py.
   (1) d_a U, d_a^2 U  (Loewner divided differences)  vs explicit expm finite differences
   (2) Phi', Phi''  (eigenvalue perturbation theory)   vs finite differences of Phi
   (3) eigenvalue-flip test: rebuild A from (Lam,V) and check tr(A lam^a)=ksi^a, a=ksi_3
Run: py -3.12 -I na2_selftest.py
"""
import numpy as np
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))   # -I implies -P
import na2_su3_na as M

def expm_antiherm(A, h=0.0):
    """exp(iA) for Hermitian traceless A (explicit, via eigh)"""
    w, V = np.linalg.eigh(A)
    return (V * np.exp(1j * (w + h))[None, :]) @ V.conj().T

def test(N, seed=1, ns=6):
    rng = np.random.default_rng(seed)
    lam = M.lambdas(N); lam3 = lam[2]
    print(f"--- N={N} ---", flush=True)
    worst = {}
    for trial in range(ns):
        # random flat-chart point: random Hermitian traceless A, random small-ish
        w = rng.standard_normal(N); w -= w.mean()
        X = rng.standard_normal((N, N)) + 1j * rng.standard_normal((N, N))
        U0, _ = np.linalg.qr(X)
        A = (U0 * w[None, :]) @ U0.conj().T
        A = 0.5 * (A + A.conj().T)
        w, V = np.linalg.eigh(A)                       # ascending
        o = np.argsort(-w); w = w[o]; V = V[:, o]
        # ---- (1) dU, ddU
        Mt = V.conj().T @ (lam3 / 2.0) @ V
        d = w[:, None] - w[None, :]
        same = np.abs(d) < 1e-12
        dd = np.where(same, 1.0, d)
        e = np.exp(1j * w)
        D1 = 2j * np.exp(0.5j * (w[:, None] + w[None, :])) * (np.sin(dd / 2) / dd)
        D1 = np.where(same, 1j * e[:, None], D1)
        fy = 1j * e[None, :]
        D2 = np.where(same, -0.5 * np.exp(0.5j * (w[:, None] + w[None, :])), (D1 - fy) / dd)
        dU = V @ (Mt * D1) @ V.conj().T
        ddU = 2.0 * V @ ((Mt @ Mt) * D2) @ V.conj().T
        for h in (1e-3, 1e-4, 1e-5):
            fp = expm_antiherm(A + h * lam3 / 2)
            fm = expm_antiherm(A - h * lam3 / 2)
            f0 = expm_antiherm(A)
            dU_fd = (fp - fm) / (2 * h)
            ddU_fd = (fp - 2 * f0 + fm) / h ** 2
            e1 = np.max(np.abs(dU - dU_fd)) / max(np.max(np.abs(dU)), 1e-30)
            e2 = np.max(np.abs(ddU - ddU_fd)) / max(np.max(np.abs(ddU)), 1e-30)
            worst.setdefault('dU', []).append((h, e1))
            worst.setdefault('ddU', []).append((h, e2))
        # ---- (2) Phi
        iu = np.triu_indices(N, 1)
        def Phi_of(Amat):
            w2 = np.linalg.eigvalsh(Amat)
            D2v = w2[:, None] - w2[None, :]
            iu2 = np.triu_indices(N, 1)
            Dv = D2v[iu2]
            return np.sum(2 * np.log(np.abs(Dv) / 2) - 2 * np.log(np.abs(np.sin(Dv / 2))))
        Dl = w[:, None] - w[None, :]
        Dl = Dl[iu]
        g = 2 * np.log(np.abs(Dl) / 2) - 2 * np.log(np.abs(np.sin(Dl / 2)))
        gp = 2.0 / Dl - 1.0 / np.tan(Dl / 2)
        gpp = -2.0 / Dl ** 2 + 0.5 / np.sin(Dl / 2) ** 2
        dPhi = np.zeros(N); d2Phi = np.zeros((N, N))
        for idx, (j, k) in enumerate(zip(iu[0], iu[1])):
            dPhi[j] += gp[idx]; dPhi[k] -= gp[idx]
            d2Phi[j, j] += gpp[idx]; d2Phi[k, k] += gpp[idx]
            d2Phi[j, k] -= gpp[idx]; d2Phi[k, j] -= gpp[idx]
        m1 = np.real(np.diag(Mt))
        dinf = w[:, None] - w[None, :]; dinf[np.arange(N), np.arange(N)] = np.inf
        s2 = 2.0 * np.real(np.sum(np.abs(Mt) ** 2 / dinf, axis=1))
        P1 = dPhi @ m1
        P2 = dPhi @ s2 + np.einsum('j,k,jk->', m1, m1, d2Phi)
        for h in (1e-3, 1e-4, 1e-5):
            P1f = (Phi_of(A + h * lam3 / 2) - Phi_of(A - h * lam3 / 2)) / (2 * h)
            P2f = (Phi_of(A + h * lam3 / 2) - 2 * Phi_of(A) + Phi_of(A - h * lam3 / 2)) / h ** 2
            worst.setdefault('Phi1', []).append((h, abs(P1 - P1f) / max(abs(P1), 1e-30)))
            worst.setdefault('Phi2', []).append((h, abs(P2 - P2f) / max(abs(P2), 1e-30)))
        # ---- (3) flat-chart coordinate bookkeeping
        ksi = np.array([np.real(np.trace(A @ l)) for l in lam])
        assert abs(ksi[2] - np.real(np.trace(A @ lam3))) < 1e-12
    for k, v in worst.items():
        for h in (1e-3, 1e-4, 1e-5):
            errs = [e for (hh, e) in v if hh == h]
            print(f"     {k:6s} h={h:.0e}  max rel err = {max(errs):.3e}", flush=True)

test(2); test(3)
print("selftest done", flush=True)
