#!/usr/bin/env python3
r"""audit_a1_sun_hessian.py -- INDEPENDENT recomputation of Hess Phi|_0 for su(N).

Charts.  (P) physicist: U = exp(i xi^a T^a), T^a hermitian, tr T^a T^b = delta^{ab}/2.
         (K) Killing-orthonormal: b_a = i sqrt2 T^a, antihermitian, -tr b_a b_b = delta_ab.
Phi(X) = -log J(X),  J = det((1-e^{-ad X})/ad X) = prod_{alpha>0} (sin(alpha(X)/2)/(alpha(X)/2))^2.

Checks:
  (1) finite-difference Hessian of Phi at X=0 in the K basis  -> eigenvalues (isotropy = Ad-invariance)
  (2) Rayleigh quotient 2 Phi(t v)/t^2 on random unit v and on a pure Cartan v
  (3) exact Cartan closed form:  X = i diag(x),  Phi(X) = sum_{k<l} g(x_k - x_l),
      g(t) = -2 log(sin(t/2)/(t/2));  g''(0) = 1/6;  sum_{k<l}(x_k-x_l)^2 = N|x|^2
  (4) coroot/root coordinate:  X = t*H/alpha(H) with alpha(X)=t  ->  Phi'' = N/12  (all N)
  (5) reproduce the ORIGINAL hc_spherical.py Hessian and localise its error
"""
import numpy as np


def basis_antiherm(n):
    B = []
    for j in range(n):
        for k in range(j + 1, n):
            E = np.zeros((n, n), complex); E[j, k] = 1.0; E[k, j] = -1.0
            B.append(E / np.sqrt(2.0))
    for j in range(n):
        for k in range(j + 1, n):
            E = np.zeros((n, n), complex); E[j, k] = 1j; E[k, j] = 1j
            B.append(E / np.sqrt(2.0))
    for l in range(1, n):
        E = np.zeros((n, n), complex)
        for j in range(l):
            E[j, j] = 1.0
        E[l, l] = -float(l)
        B.append(1j * E / np.sqrt(l * (l + 1.0)))        # NOTE the factor i: Cartan generators are antihermitian
    return B


def ad(n, X, B):
    d = len(B)
    A = np.zeros((d, d))
    for i, bi in enumerate(B):
        Y = X @ bi - bi @ X
        for j, bj in enumerate(B):
            A[j, i] = -np.trace(Y @ bj).real
    return A


def phi(coords, B):
    X = sum(c * b for c, b in zip(coords, B))
    z = np.linalg.eigvals(ad(len(B[0]), X, B))
    small = np.abs(z) < 1e-13
    f = np.where(small, 1.0 + 0j, (1.0 - np.exp(-z)) / np.where(small, 1.0, z))
    return -np.log(np.abs(np.prod(f)))


def hess_fd(B, h=1e-3):
    d = len(B)
    H = np.zeros((d, d))
    for i in range(d):
        for j in range(i, d):
            ei = np.zeros(d); ei[i] = h
            ej = np.zeros(d); ej[j] = h
            v = (phi(ei + ej, B) - phi(ei - ej, B) - phi(-ei + ej, B) + phi(-ei - ej, B)) / (4 * h * h)
            H[i, j] = H[j, i] = v
    return H


def g(t):
    return 0.0 if abs(t) < 1e-12 else -2.0 * np.log(abs(np.sin(t / 2) / (t / 2)))


def main():
    for n in (2, 3, 4):
        B = basis_antiherm(n)
        d = len(B)
        assert all(np.allclose(b.conj().T, -b) for b in B)
        assert np.allclose([-np.trace(b @ b).real for b in B], 1.0)
        H = hess_fd(B)
        ev = np.linalg.eigvalsh(H)
        print(f"su({n}): dim={d}")
        print(f"   (1) FD Hess eigenvalues: min={ev.min():.6f} max={ev.max():.6f} spread={ev.max()-ev.min():.1e}")
        print(f"       c_N(K chart, -tr orthonormal) = {np.mean(ev):.6f}    [N/6 = {n/6:.6f}]")
        print(f"       c_N(P chart, tr T^aT^b=delta/2) = {np.mean(ev)/2:.6f}  [N/12 = {n/12:.6f}]")
        rng = np.random.default_rng(1)
        vcart = np.zeros(d, complex)
        # a unit vector along the Cartan direction diag(1,0,...,-1)/sqrt2  in K coordinates
        vcart[-1] = 1.0
        for tag, v in (("random", rng.normal(size=d) / np.linalg.norm(rng.normal(size=d))), ("Cartan", vcart)):
            v = v / np.linalg.norm(v)
            t = 1e-4
            print(f"   (2) Rayleigh({tag:6s}) = 2Phi(tv)/t^2 = {2*phi(t*v, B)/t**2:.6f}  (v.H.v = {np.real(v @ H @ v):.6f})")
        # (3) Cartan closed form in x-coordinates: Phi(i diag(x)) = sum_{k<l} g(x_k-x_l)
        x = np.array([0.21, -0.13, -0.08] + [0.0] * (n - 3))[:n]
        x = x - x.mean()
        X = 1j * np.diag(x)
        coords = np.array([-np.trace(X @ b).real for b in B])          # expand X in the K basis
        direct = sum(g(x[k] - x[l]) for k in range(n) for l in range(k + 1, n))
        print(f"   (3) Cartan closed form vs spectral: Phi(i diag x) = {phi(coords, B):.10f} vs {direct:.10f}")
        # (4) root coordinate: X = t*H_alpha with alpha(X)=t  (H = diag(1,-1,0,..)/2)
        Hc = 0.5 * np.diag([1.0, -1.0] + [0.0] * (n - 2))
        Xr = 1j * Hc
        cr = np.array([-np.trace(Xr @ b).real for b in B])
        t0 = 1e-4
        print(f"   (4) root-coordinate Phi'' = 2 Phi(t H_alpha)/t^2 = {2*phi(t0*cr, B)/t0**2:.6f}   [N/12 = {n/12:.6f}]")
        print()

    import importlib.util
    spec = importlib.util.spec_from_file_location("hc", "hc_spherical.py")
    hc = importlib.util.module_from_spec(spec); spec.loader.exec_module(hc)
    for n in (2, 3):
        B = hc.su_n_basis(n)
        d = len(B)
        norms = [round(float(-np.trace(b @ b).real), 3) for b in B]
        herm = [bool(np.allclose(b.conj().T, b)) for b in B]
        print(f"hc_spherical.su_n_basis({n}): -tr(b^2) = {norms}  hermitian? {herm}")
        X = 1e-3 * B[0]
        print(f"   ||hc.ad_matrix(1e-3 b_0)||_max = {np.max(np.abs(hc.ad_matrix(n, X))):.3e}  (should be O(1e-3))")
        # reproduce the script's own Hessian and show its eigenvalues
        Bh = hc.su_n_basis(n)

        def f(coords, n=n, Bh=Bh):
            Xm = sum(c * b for c, b in zip(coords, Bh))
            return hc.phi_flat(n, Xm)
        Hh, _ = hc.hessian(f, np.zeros(d), h=1e-3)
        evh = np.linalg.eigvalsh(Hh)
        print(f"   hc HessPhi eigenvalues = {np.round(evh, 6)}")
        print(f"   hc mean={np.mean(evh):.6f} spread={evh.max()-evh.min():.6f}")
    print()
    print("DICTIONARY: c_K = N/6 ; c_P = N/12 ; root-coordinate = N/12 ; invariant form Hess = -(1/12) B_Killing")


if __name__ == "__main__":
    main()
