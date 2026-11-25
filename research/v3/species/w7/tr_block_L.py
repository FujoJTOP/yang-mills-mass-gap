#!/usr/bin/env python3
r"""
tr_block_L.py -- TR species, test 2: does BLOCKING save the tropical squeeze?

Lead's survival hypothesis: "the two-sided log-sum-exp squeeze on rho gives a two-sided bound on
m = log(l0/l1) whose width stays non-zero as L -> oo because after blocking the number of terms is
N = e^{o(L)}".  This script tests it exactly on the P4 kernel, where everything (eigenvalues, the
max-plus power, the exact Collatz-Wielandt width) can be computed without approximation.

For the L-step (blocked) operator  T_L = K^L :
   [T]  e^{beta L rho_inf} <= rho(T_L) <= e^{beta L rho_inf} * N_eff^{(L)} ,
        N_eff^{(L)} = max_i sum_j (K^L)_ij e^{beta(u_j-u_i)} / e^{beta L rho_inf} ,
   so the PER-STEP width is  w_L = log N_eff^{(L)} / (beta L),  and the exact per-step gap is
   Delta = (1/L) log(l0(T_L)/l1(T_L)) = -log(l1(K)/l0(K)), L-INDEPENDENT.
Report w_L, Delta, w_L/Delta, and the entropy rate s_L = (1/L) log N_eff^{(L)}.

Also: the "second tropical eigenvalue" reading of the gap (rigid reading):
   delta_2 = rho_inf - max_{i not critical} W_ii ,   tpred(L) = beta L delta_2 ,
compared with the exact Delta*L.  (If delta_2 = O(h^2) the rigid reading is empty.)

Run: py -3.12 -I tr_block_L.py.  numpy only, deterministic.  Labels [T]/[V].  No QED.
"""
import numpy as np
from scipy.special import sici

EULER = 0.5772156649015328606
TWO_PI = 2.0 * np.pi


def Cin(x):
    x = np.asarray(x, float)
    _, ci = sici(x)
    return np.where(x > 0, EULER + np.log(np.where(x > 0, x, 1.0)) - ci, 0.0)


def rho1(z):
    az = np.abs(np.asarray(z, float))
    return np.where(az <= TWO_PI, np.pi * (Cin(TWO_PI) - Cin(np.minimum(az, TWO_PI))), 0.0)


def nodes(nz):
    xg, wg = np.polynomial.legendre.leggauss(nz)
    z = TWO_PI * xg
    return z, rho1(z) * (TWO_PI * wg)


def kernel(z, r1w, beta):
    s = np.sqrt(np.maximum(r1w, 0.0))
    M = (s[:, None] * np.exp(beta * np.cos(0.5 * (z[None, :] - z[:, None])))) * s[None, :]
    return 0.5 * (M + M.T)


def karp(W):
    n = W.shape[0]
    Fs = np.empty((n + 1, n))
    Fs[0] = 0.0
    for k in range(n):
        Fs[k + 1] = np.max(Fs[k][:, None] + W, axis=0)
    ks = np.arange(n)
    vals = (Fs[n][None, :] - Fs[:n, :]) / (n - ks)[:, None]
    return float(np.max(np.min(vals, axis=0)))


def tropical_vec(W, rho, iters=3000):
    u = np.zeros(W.shape[0])
    for _ in range(iters):
        un = np.max(W + u[None, :], axis=1) - rho
        un -= un.max()
        if np.max(np.abs(un - u)) < 1e-14:
            return un
        u = un
    return u


def main():
    nz, beta = 320, 2.5
    z, r1w = nodes(nz)
    K = kernel(z, r1w, beta)
    ev = np.linalg.eigvalsh(K)[::-1]
    Delta = -np.log(ev[1] / ev[0])
    W = np.log(K) / beta
    rho = karp(W)
    u = tropical_vec(W, rho)
    # critical set of the tropical problem (self-loops attaining rho)
    diag = np.diag(W)
    crit = np.abs(diag - rho) < 1e-9
    d2 = rho - np.max(diag[~crit]) if np.any(~crit) else 0.0
    print("#" * 100)
    print("# TR test 2 -- BLOCKING (T_L = K^L) and the per-step squeeze width")
    print("#" * 100)
    print(f"    nz={nz}  beta={beta}  rho_inf={rho:.6f}  exact Delta={Delta:.6f}  "
          f"l1/l0={ev[1]/ev[0]:.6f}")
    print(f"    tropical critical self-loops: {int(crit.sum())} of {nz};  "
          f"second-diagonal deficit delta_2 = {d2:.3e}  => rigid prediction beta*delta_2 = "
          f"{beta*d2:.3e}")
    print()
    print(f"    {'L':>3}{'lam0(T_L)^(1/L)':>18}{'Delta_L (per step)':>20}{'rho_inf^(L)/L':>15}"
          f"{'N_eff^(L)':>13}{'w_L':>9}{'w_L/Delta':>11}{'s_L=logN/L':>12}{'beta*L*d2':>11}")
    T = np.eye(nz)
    for L in range(1, 7):
        T = T @ K
        evT = np.linalg.eigvalsh(T)[::-1]
        D_L = -np.log(evT[1] / evT[0]) / L
        # exact Collatz-Wielandt width for the blocked operator, with the (per-step) tropical vector u
        A = (T * np.exp(beta * (u[None, :] - u[:, None]))).sum(axis=1)
        logNeff = np.log(A.max()) - beta * L * rho
        wL = logNeff / (beta * L)
        print(f"    {L:3d}{evT[0]**(1.0/L):18.6f}{D_L:20.6f}{rho:15.6f}"
              f"{np.exp(logNeff):13.4e}{wL:9.4f}{wL/D_L:11.3f}{logNeff/L:12.4f}{beta*L*d2:11.3e}")
    print()
    print("    reading: Delta_L is EXACTLY L-independent (eigenvalues are powers) and so is w_L:")
    print("    blocking multiplies numerator and denominator by L.  The per-step width is the")
    print("    entropy rate s_L/beta with s_L ~ log N_eff (POSITIVE), so N_eff^(L) = e^{L s}, not")
    print("    e^{o(L)}.  The lead's survival hypothesis fails quantitatively on this kernel.")
    print()
    print("VERDICT -> TR-TROPICAL.md.  No QED.")


if __name__ == "__main__":
    main()
