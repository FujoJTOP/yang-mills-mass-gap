#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
gy3_chain_m6.py -- GY3 task (b): model-level counterexample hunt.

[B] 1xL SU(2) chain (Q6-H1-STATUS.md:178): flat-chart links U_i = exp(i xi_i.sigma/2),
    Gibbs measure dmu ~ prod_i exp(beta Re tr U_i U_{i+1}^dag) dU_i (Haar link measure).
    Transfer operator T(U,V) = exp(2 beta cos th(U V^dag)) = sum_j c_j(beta) chi_j(U V^dag);
    eigenvalue on the j-isotypic: lambda_j = c_j/(2j+1)  =>  m_j = -ln(lambda_j/lambda_0)
    L-INDEPENDENT.  Carrier O = xi_3 (H1 flat-chart coordinate, per link):
    KL measure nu = sum_j w_j delta_{m_j^2}, w_j = <O P_j O> (Haar geometry, beta-indep),
    chi_L = Var(L^{-1/2} sum xi_3) = sum_j w_j h_L(m_j) (closed form).  MC cross-check.

[M6] osc(x)free (N02-FOURMODEL.md 2.3): corpus m_L = min(omega, (1/2)(2 pi/L)^2).
     Quantized completion of the free direction = free particle on ring [0,L):
     E_n = 2 pi^2 n^2/L^2 (n>=1), uniform ground state. Exact KL atom tables for carriers
     b and :b^2: (connected): <n|b|0>: |.|^2 = L^2/(4pi^2 n^2); <n|b^2|0>: |.|^2 =
     L^4(1/(4pi^2 n^2) + 1/(4pi^4 n^4)). Atoms at s_n = E_n^2 (mass^2 axis), weight
     2|<n|O|0>|^2; on the corpus's energy-axis reading m_L = E_1 = 2 pi^2/L^2.

py -3.12 -I.  Log: w7/gy3_chain_m6.log
"""
import sys, time
import numpy as np

sys.path.insert(0, 'D:/Dev/math/research/lattice-cpb')
import w4_blocked as W

out = []
def P(*a):
    s = " ".join(str(x) for x in a)
    out.append(s); print(s, flush=True)

def gauleg(n, a, b):
    x, w = np.polynomial.legendre.leggauss(n)
    return 0.5 * (b - a) * x + 0.5 * (a + b), 0.5 * (b - a) * w

NT = 600
th, wth = gauleg(NT, 0.0, np.pi)
sin2 = np.sin(th) ** 2

def chi_j(j, thv):
    return np.sin((2 * j + 1) * thv) / np.sin(thv)

def c_j(beta, j):
    """c_j = int dU e^{2 beta cos th} chi_j ; int dU = 1 => (1/2pi) int_0^pi sin^2 th (...) dth."""
    # int dU f(cos th) = (2/pi) int_0^pi sin^2 th f dth   (SU(2) class-function Haar, int dU = 1)
    return float(np.sum(wth * sin2 * np.exp(2 * beta * np.cos(th)) * chi_j(j, th)) * 2 / np.pi)

P("=" * 72)
P("[B] 1xL SU(2) chain: transfer spectrum (exact quadrature, NT=%d)" % NT)
P("    lambda_j = c_j/(2j+1);  m_j = -ln(lambda_j/lambda_0)  (L-independent)")
P("=" * 72)
P(f"  sanity beta=0: c_0 = {c_j(0.0,0.0):.12f} (want 1), c_1 = {c_j(0.0,1.0):.2e} (want 0)")
js = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0]
mchain = {}
for beta in (0.5, 1.0, 2.5, 5.0):
    l0 = c_j(beta, 0.0)
    mchain[beta] = {}
    P(f"  --- beta={beta} ---   lambda_0 = c_0 = {l0:.8f}")
    P(f"  {'j':>5} | {'c_j':>13} | {'lambda_j':>13} | {'lambda_j/lambda_0':>17} | {'m_j':>9}")
    for j in js:
        lj = c_j(beta, j) / (2 * j + 1)
        r = lj / l0
        m = -np.log(r) if r > 1e-300 else float('inf')
        mchain[beta][j] = m
        P(f"  {j:>5} | {c_j(beta,j):>13.5e} | {lj:>13.5e} | {r:>17.8e} | {m:>9.5f}")

# ---- w_j for the xi_3 carrier:  w_j = (2j+1) <O P_j O>,  P_j f = (2j+1) int chi_j(V^-1 U) f(V) dV
P("")
P("  KL layer weights w_j = <O P_j O> for O = xi_3 (Haar geometry, beta-independent):")
NC = 500
cc, wc = gauleg(NC, -1.0, 1.0)
cosT, sinT = np.cos(th), np.sin(th)
def w_j_clean(j):
    acc = 0.0
    for iU in range(NT):
        cU, sU, tU = cosT[iU], sinT[iU], th[iU]
        wU, w2U = wth[iU], sin2[iU]
        for iV in range(NT):
            cV, sV, tV = cosT[iV], sinT[iV], th[iV]
            cg = cU * cV + cc * (sU * sV)
            thg = np.arccos(np.clip(cg, -1.0, 1.0))
            int_c = float(np.sum(wc * cc * chi_j(j, thg))) / 2.0
            acc += wU * wth[iV] * w2U * sin2[iV] * tU * tV * int_c
    # dU = (2/pi) sin^2 th dth dOm/(4pi); angular: E[n3 n3' F(c)] = (1/3) E[c F(c)]
    return (2 * j + 1) * acc * (4 / np.pi ** 2) / 3.0

t0 = time.time()
wj = {}
for j in (0.5, 1.0, 1.5, 2.0, 2.5, 3.0):
    wj[j] = 4.0 * w_j_clean(j)   # corpus xi = 2 theta nhat  =>  x4
    P(f"  w_{j} = {wj[j]:.10f}   (cumulative {sum(wj.values()):.6f}, {time.time()-t0:.0f}s)")
v_haar = 4.0 * (2 / (3 * np.pi)) * float(np.sum(wth * sin2 * th ** 2))
P(f"  <xi_3^2>_Haar direct = {v_haar:.10f}   closed form (4pi^2-6)/9 = {(4*np.pi**2-6)/9:.10f}"
  f"  (x4 rescale verified)")
P(f"  Parseval: sum_j w_j (j=1/2..3) = {sum(wj.values()):.6f} vs Var = {v_haar:.6f}"
  f"  (tail j>3 carries the rest)")

def h_L(m, L):
    # ring: Var(a) = sum_{d=0}^{L-1} C(d); C(d) = w e^{-md} (+ image e^{-m(L-d)} counted by
    # the symmetric d-range; leading-layer closed form: 1 + 2(r - r^L)/(1 - r))
    if not np.isfinite(m) or m > 40:
        return 1.0
    r = np.exp(-m)
    # exact leading-layer ring: 1 + sum_{d=1}^{L-1} (e^{-md} + e^{-m(L-d)}) - (L-1) e^{-mL}
    return 1.0 + 2.0 * (r - r ** L) / (1.0 - r) - (L - 1) * r ** L

def chi_L_closed(beta, L):
    l0 = c_j(beta, 0.0)
    tot = 0.0
    for j, w in wj.items():
        lj = c_j(beta, j) / (2 * j + 1)
        m = -np.log(lj / l0)
        tot += w * h_L(m, L)
    return tot

def pred_C(beta, d):
    l0 = c_j(beta, 0.0)
    return sum(w * np.exp(-m * d) for j, w in wj.items()
               for m in [-np.log(c_j(beta, j) / (2 * j + 1) / l0)])

# ---- MC: exact chain heatbath, even/odd vectorized Gibbs
def su2_log_xi3(U):
    # corpus flat chart: U = exp(i xi.sigma/2), xi = alpha nhat, alpha = 2 theta in [0, 2 pi)
    c0 = np.clip(U[..., 0, 0].real, -1, 1)
    thv = np.arccos(c0)
    s = np.sin(thv)
    f = np.where(s < 1e-12, 1.0, thv / np.where(s < 1e-12, 1.0, s))
    return 2.0 * np.stack([U[..., 0, 1].imag * f, U[..., 0, 1].real * f, U[..., 0, 0].imag * f], axis=-1)

def run_chain(L, beta, ncfg, therm, gap, seed):
    rng = np.random.default_rng(seed)
    U = W.random_su2(rng, (L,))
    ie = np.arange(0, L, 2); io = np.arange(1, L, 2)
    def sweep(ids):
        S = W.adj(np.roll(U, -1, axis=0)[ids]) + W.adj(np.roll(U, 1, axis=0)[ids])
        A = U[ids]
        U[ids] = W.su2_heatbath(beta, A @ S, 0, 1, rng) @ A
    for _ in range(therm):
        sweep(ie); sweep(io)
    chis = np.empty(ncfg)
    Cd = np.zeros(21)
    for c in range(ncfg):
        for _ in range(gap):
            sweep(ie); sweep(io)
        xi = su2_log_xi3(U)[:, 2]
        chis[c] = (xi.sum() / np.sqrt(L)) ** 2
        for d in range(21):
            Cd[d] += float(np.mean(xi * np.roll(xi, -d)))
    return chis, Cd / ncfg

P("")
P("  MC cross-check (exact heatbath chain; chi_L = mean (L^-1/2 sum xi_3)^2):")
for beta in (0.5, 2.5):
    P(f"  --- beta={beta} ---   m(j=1/2, lightest w/ flat-chart weight) = {mchain[beta][0.5]:.5f}   "
      f"m(j=1) = {mchain[beta][1.0]:.5f}   m(j=2) = {mchain[beta][2.0]:.5f}")
    P(f"  {'L':>6} | {'chi_L MC':>12} | {'chi_L closed':>13} | {'C(1) MC':>11} | {'C(1) closed':>12} | "
      f"{'C(2) MC':>11} | {'C(2) closed':>12}")
    for L in (8, 32, 128, 512):
        t0 = time.time()
        ncfg = 12000 if L == 8 else (3000 if L <= 128 else 1500)
        chis, Cd = run_chain(L, beta, ncfg, 300, 2, seed=1234 + L + int(beta * 100))
        x = chis - chis.mean(); v = x.var(); tau = 1.0
        for lag in range(1, len(x) // 4):
            c_lag = float(np.mean(x[:-lag] * x[lag:]) / v)
            if c_lag < 0:
                break
            tau += 2 * c_lag
        se = float(np.sqrt(v * tau / len(x)))
        P(f"  {L:>6} | {chis.mean():>8.4f}±{se:.4f} | {chi_L_closed(beta, L):>13.4f} | "
          f"{Cd[1]:>11.5f} | {pred_C(beta,1):>12.5f} | {Cd[2]:>11.5f} | {pred_C(beta,2):>12.5f}"
          f"   ({time.time()-t0:.0f}s, tau={tau:.1f})")
    l0 = c_j(beta, 0.0)
    m1 = mchain[beta][1.0]
    P(f"  chi_L(inf) closed = {chi_L_closed(beta, 100000):.5f};  correlation length 1/m_1 = {1/m1:.2f}")

# ------------------------------------------------------------------ M6
P("")
P("=" * 72)
P("[M6] osc (x) free, free direction quantized = free particle on ring [0,L):")
P("     E_n = 2 pi^2 n^2 / L^2 (n >= 1); carriers b and :b^2: (connected KL)")
P("=" * 72)
P(f"{'L':>5} | {'m_L=E_1 (corpus axis)':>21} | {'m1^2 mass^2-axis = E_1^2':>25} | "
  f"{'f1^2(b)':>13} | {'f1^2(b^2)':>13} | {'chi_KL(b)':>12} | {'chi=Var(b)':>11}")
nn = np.arange(1, 200001, dtype=np.float64)
for L in (4, 8, 16, 32, 64):
    E1 = 2 * np.pi ** 2 / L ** 2
    f1_b = L ** 2 / (2 * np.pi ** 2)
    f1_b2 = L ** 4 * (1 / (2 * np.pi ** 2) + 1 / (2 * np.pi ** 4))
    chi_b = float(np.sum(L ** 4 / (4 * np.pi ** 4 * nn ** 4)))
    P(f"{L:>5} | {E1:>21.6e} | {E1**2:>25.6e} | {f1_b:>13.5e} | {f1_b2:>13.5e} | "
      f"{chi_b:>12.4e} | {L**2/12:>11.4f}")
P("")
P("  Parseval cross-checks (exact):")
for L in (4, 8):
    s1 = float(np.sum(2 * L ** 2 / (4 * np.pi ** 2 * nn ** 2)))
    s2 = float(np.sum(2 * L ** 4 * (1 / (4 * np.pi ** 2 * nn ** 2) + 1 / (4 * np.pi ** 4 * nn ** 4))))
    P(f"  L={L}: sum|<n|b|0>|^2   = {s1:.10f} vs Var(b)   = L^2/12     = {L**2/12:.10f}")
    P(f"  L={L}: sum|<n|b^2|0>|^2 = {s2:.10f} vs Var(:b^2:) = 4L^4/45    = {4*L**4/45:.10f}")
P("")
P("  Verdict M6: omega>0 does NOT exclude the model as a counterexample family -- the")
P("  free half is the 0+1d free field (acid F1 itself): m_1 (corpus energy axis) =")
P("  2 pi^2/L^2 = EXACTLY the corpus m_L (N02:73); f_1^2(b) = L^2/(2pi^2) -> inf;")
P("  b^2-carrier f_1^2 ~ L^4(1/2pi^2+1/2pi^4) -> inf. Raw M6 remains 'not statable'")
P("  (improper flat measure in b, no OS transfer matrix) -- the quantized completion")
P("  (particle on a ring) is proper and violates the theorem-form trivially.")

with open("D:/Dev/math/research/v3/w7/gy3_chain_m6.log", "w", encoding="utf-8") as fh:
    fh.write("\n".join(out) + "\n")
print("\n# log written: w7/gy3_chain_m6.log", flush=True)
