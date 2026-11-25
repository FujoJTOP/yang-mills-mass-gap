#!/usr/bin/env python3
"""audit_b1_tr_flat.py -- independent 20-line recomputation of the TR flat model.

K_ij = exp(beta cos(theta_i - theta_j)),  theta_k = 2 pi k / n  (circulant).
Tropical (max-plus) datum: W_ij = beta cos(theta_i - theta_j); every W_ii = beta is the row max
=> critical graph = the whole diagonal (n self-loops, all mean beta) => tropical gap == 0.
Exact gap: Delta = -log(lambda_1/lambda_0) from the FFT of the first row (n-independent).
Squeeze width: width = (1/beta) log sum_j exp(beta (cos th_j - 1))  ~ (log n)/beta.
"""
import numpy as np

for n, beta in ((64, 2.5), (4096, 2.5), (64, 0.5), (4096, 10.0)):
    th = 2 * np.pi * np.arange(n) / n
    row = np.exp(beta * np.cos(th))
    lam = np.sort(np.fft.fft(row).real)[::-1]
    Delta = -np.log(lam[1] / lam[0])
    width = np.log(np.exp(beta * (np.cos(th) - 1.0)).sum()) / beta
    W = beta * np.cos(th[:, None] - th[None, :])
    deficit = np.min(beta - (W - beta * np.eye(n))[~np.eye(n, dtype=bool)])   # beta - max_{i!=j} W_ij
    print(f"n={n:5d} beta={beta:5.2f}: l1/l0={lam[1]/lam[0]:.6f} Delta_exact={Delta:.6f} "
          f"Delta_trop=0 width={width:.4f} width/Delta={width/Delta:.4f} log(n)/beta={np.log(n)/beta:.4f} "
          f"deficit(beta-max_offdiag)={deficit:.2e}")

# explicit tropical check: max cycle mean = beta, attained by all n self-loops (n-fold degeneracy)
n, beta = 64, 2.5
W = beta * np.cos(2 * np.pi * (np.arange(n)[:, None] - np.arange(n)[None, :]) / n)
means = [W[i, i] for i in range(n)]
print(f"[trop] self-loop means: min={min(means)}, max={max(means)}, #ties={sum(1 for m in means if abs(m-max(means))<1e-12)}/{n}")
print(f"[trop] best 2-cycle mean  = {max((W[i,j]+W[j,i])/2 for i in range(n) for j in range(n) if i!=j):.6f}")
print(f"[trop] => tropical eigenvalue rho_inf = {max(means):.6f}, gap to any other cycle = 0")
