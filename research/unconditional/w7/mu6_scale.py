"""MU6-NUM: magnitude test of sigma = 2 g^2 N Lambda_MS^2.

All inputs are *quoted from verified literature values already checked in
research/v3/species/L2-GLUEBALL-LAMBDA.md* (AT2021 = arXiv:2106.00364v3):
  - r0*sqrt(sigma) = 1.160(6), r0 = 0.472(5) fm           (AT2021 text, p.18)
  - Lambda_MS/sqrt(sigma): SU(2) 0.5806(21)[570], SU(3) 0.5424(13)[185] (AT2021 Tab.15)
  - per-beta a*sqrt(sigma): AT2021 Table 13 (SU(3), beta=6/g^2), Table 12 (SU(2), beta=4/g^2)
No number below is from memory.
"""
import numpy as np

FM = 5.0677307  # GeV^-1 per fm (hbar c)

# --- 1. sqrt(sigma) and Lambda_MS -------------------------------------------
r0_sqrtsig = 1.160          # AT2021 text p.18
r0_fm = 0.472               # AT2021 text p.18
sqrtsig = r0_sqrtsig / (r0_fm * FM)      # GeV
sigma = sqrtsig**2                        # GeV^2

Lam = {"SU3": 0.5424 * sqrtsig, "SU2": 0.5806 * sqrtsig}   # GeV (Tab.15 ratios)
print(f"sqrt(sigma) = {sqrtsig:.6f} GeV = {sqrtsig*1e3:.2f} MeV   sigma = {sigma:.6f} GeV^2")
for k, v in Lam.items():
    print(f"Lambda_MS[{k}] = {v*1e3:.3f} MeV      sigma/Lambda^2 = {sigma/v**2:.4f}")

# --- 2. two-loop running of g^2, pure gauge SU(N) --------------------------
#  d g^2/d ln mu^2 = -b0 g^4 - b1 g^6 ,  b0 = 11N/(48 pi^2), b1 = 34 N^2/(768 pi^4)
def coeffs(N):
    return 11.0 * N / (48 * np.pi**2), 34.0 * N**2 / (768 * np.pi**4)

def run_two_loop(g2_0, mu0, mu, N):
    """exact 2-loop solution: t - t0 = (1/b0)[(u-u0) - (b1/b0) ln((b0 u+b1)/(b0 u0+b1))]"""
    b0, b1 = coeffs(N)
    u0, dt = 1.0 / g2_0, 2 * np.log(mu / mu0)
    u = u0 + b0 * dt                      # 1-loop start
    for _ in range(60):                   # fixed point on the log term
        u = u0 + b0 * dt + (b1 / b0) * np.log((b0 * u + b1) / (b0 * u0 + b1))
    return 1.0 / u

a_inv = 2.91   # GeV, beta = 6 (task-given anchor)
for N, key in ((3, "SU3"), (2, "SU2")):
    g2_L = run_two_loop(1.0, a_inv, Lam[key], N)
    b0, b1 = coeffs(N)
    print(f"[N={N}] b0={b0:.6f} b1={b1:.6f}  g2(a^-1=2.91GeV)=1.0  ->  "
          f"g2(Lambda_MS)={g2_L:.4f}  =>  2 g2 N = {2*g2_L*N:.3f}   "
          f"(at beta=6: 2 g2 N = {2*1.0*N:.2f})")

# --- 3. ratios -------------------------------------------------------------
print("\n--- sigma/Lambda^2  vs  2 g^2 N ---")
for N, key in ((3, "SU3"), (2, "SU2")):
    meas = sigma / Lam[key]**2
    for tag, g2 in (("g2=6/beta=1", 1.0), ("g2(Lambda_MS) 2-loop", run_two_loop(1.0, a_inv, Lam[key], N))):
        form = 2 * g2 * N
        print(f"  {key:4s} {tag:22s} formula={form:6.3f}  measured/formula = {meas/form:.4f}")

# --- 4. per-beta sigma a^2  vs  exp(-1/(2 b0 g^2)) -------------------------
tab13 = [(5.6924, .4010), (5.80, .31603), (5.8941, .2613), (5.99, .21959),
         (6.0625, .19509), (6.235, .14899), (6.338, .12948), (6.50, .10319),
         (6.60, .09024), (6.70, .07898)]                       # SU(3), beta=6/g^2
tab12 = [(2.509, .17857), (2.60, .13279), (2.65, .11342), (2.70, .09698),
         (2.75, .08375), (2.80, .07242)]                       # SU(2), beta=4/g^2

for N, tab, key in ((3, tab13, "SU3"), (2, tab12, "SU2")):
    b0, _ = coeffs(N)
    b = np.array([t[0] for t in tab]); asig = np.array([t[1] for t in tab])
    g2 = 2 * N / b                                    # beta = 2N/g^2
    sa2 = asig**2
    x = np.exp(-1.0 / (2 * b0 * g2))                  # one-loop a*Lambda exponential
    R = sa2 / x
    slope = np.polyfit(np.log(x), np.log(sa2), 1)[0]  # sigma a^2 ~ x^n
    print(f"\n--- {key}: beta in [{b.min():.4f},{b.max():.4f}], g2 = 2N/beta ---")
    print(f"  {'beta':>7} {'a sqrt(sig)':>11} {'sigma a^2':>10} {'x=e^-1/(2b0g2)':>15} {'R=sa2/x':>9}")
    for i in range(len(b)):
        print(f"  {b[i]:7.4f} {asig[i]:11.5f} {sa2[i]:10.6f} {x[i]:15.3e} {R[i]:9.2f}")
    print(f"  sigma a^2 falls x{sa2.max()/sa2.min():.2f};  x falls x{x.max()/x.min():.2f};  "
          f"R runs x{R.max()/R.min():.2f}  =>  exponent n = dln(sa2)/dln(x) = {slope:.3f}")

# --- 5. cross-check of a^-1 at beta=6 against Table 13 ---------------------
a_sig_6 = np.interp(6.0, [t[0] for t in tab13], [t[1] for t in tab13])
print(f"\na*sqrt(sigma)(beta=6, interp. Tab.13) = {a_sig_6:.5f}  =>  "
      f"a^-1 = {sqrtsig/a_sig_6:.3f} GeV   (task anchor: 2.91 GeV)")
print(f"  sigma a^2(beta=6) = {a_sig_6**2:.6f} (Tab.13)  vs  sigma/2.91^2 = {sigma/a_inv**2:.6f}")
x6 = np.exp(-1.0 / (2 * coeffs(3)[0] * 1.0))
print(f"  e^(-1/(2 b0 g2)) at beta=6 (g2=1) = {x6:.4e}")
print(f"  R(beta=6) = {a_sig_6**2/x6:.2f} (Tab.13 a)   {sigma/a_inv**2/x6:.2f} (a^-1=2.91)")
