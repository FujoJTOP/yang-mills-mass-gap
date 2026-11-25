"""
DM-FIFTH-MODE witness 1/2.

Purpose: quantify the ONE-SIDEDNESS of every finite data mode.

Object:  d_N(eps; nu) := max{ nu'([0,eps)) : nu' >= 0 on [0,S], D_N(nu') = D_N(nu) }
         = "moment/Laplace floor defect"  (blindness of the finite floor test at level eps).

Claims tested:
  (a) d_N(eps) > 0 for every finite N  -> finite data never certifies  floor >= eps   [blindness]
  (b) d_N(eps) -> 0 as N -> oo         -> the data *approach* the floor test (Hausdorff determinacy)
  (c) min over the same fiber of nu'([0,eps)) = 0  -> finite data never certifies  floor < eps either
      i.e. BOTH sides of the floor test fail at finite N  => the floor test is two-sided-undecidable.
"""
import numpy as np
from scipy.optimize import linprog

np.set_printoptions(precision=6, suppress=True)

S = 4.0
n = 401
grid = np.linspace(0.0, S, n)
ds = grid[1] - grid[0]

# ---- reference measure: YM-like KL measure -----------------------------------
# edge atom f1^2 = 0.40 at m1^2 = 0.25  +  cut on [1,4], density ~ sqrt(s-1)
m1sq, f1sq = 0.25, 0.40
w = np.zeros(n)
i0 = int(np.argmin(abs(grid - m1sq)))
w[i0] += f1sq
mask = grid > 1.0
dens = np.sqrt(np.maximum(grid[mask] - 1.0, 0.0))
dens = dens / dens.sum()
w[mask] += (1.0 - f1sq) * dens
M0 = w.sum()
print(f"reference: M0={M0:.6f}  floor={m1sq}  atom={f1sq}  ngrams={int((w>0).sum())}")

# ---- the two data maps ------------------------------------------------------
def data_moments(N):
    A = np.array([(grid / S) ** k for k in range(N + 1)])     # (N+1, n)
    return A, A @ w

def data_laplace(N, dt=0.5):
    ts = np.arange(N + 1) * dt
    A = np.exp(-np.outer(ts, grid))                          # (N+1, n)
    return A, A @ w

# ---- LP: max mass below eps, min mass below eps ------------------------------
def defect(A, b, eps):
    below = (grid < eps).astype(float)
    hi = linprog(-below, A_eq=A, b_eq=b, bounds=(0, None), method="highs")
    lo = linprog(+below, A_eq=A, b_eq=b, bounds=(0, None), method="highs")
    assert hi.status == 0 and lo.status == 0, (hi.status, lo.status)
    return -hi.fun / b[0], lo.fun / b[0]      # (max, min) of mass in [0,eps)

eps_list = [0.05, 0.125, 0.20]
print("\n=== (alpha) finite moments: mass that can hide below eps ===")
print("  N |   eps=0.05 max/min |   eps=0.125 max/min |   eps=0.20 max/min")
for N in [1, 2, 3, 4, 5, 6, 8, 10, 12, 14]:
    A, b = data_moments(N)
    row = f"{N:3d} |"
    for eps in eps_list:
        mx, mn = defect(A, b, eps)
        row += f"  {mx:.3e}/{mn:.1e} |"
    print(row)

print("\n=== (gamma) Laplace samples C(t), t in [0,T] (dt=0.5): mass below eps ===")
print("  N (T=N/2) |   eps=0.05 max/min |   eps=0.20 max/min")
for N in [1, 2, 3, 4, 5, 6, 8, 10, 12]:
    A, b = data_laplace(N)
    row = f"{N:3d} ({N*0.5:4.1f}) |"
    for eps in [0.05, 0.20]:
        mx, mn = defect(A, b, eps)
        row += f"  {mx:.3e}/{mn:.1e} |"
    print(row)

# ---- exponential rate of the defect ----------------------------------------
print("\n=== decay of the moment defect d_N(0.20) ===")
Ns, vals = [], []
for N in range(1, 15):
    A, b = data_moments(N)
    mx, _ = defect(A, b, 0.20)
    Ns.append(N); vals.append(mx)
for N, v in zip(Ns, vals):
    print(f"  N={N:3d}  d_N={v:.4e}")
ok = np.array([v > 1e-14 for v in vals])
print(f"  strictly positive for all N<=14 : {bool(ok.all())}")
if len(vals) > 4:
    r = np.polyfit(Ns[4:], np.log(np.array(vals[4:]) + 1e-300), 1)[0]
    print(f"  fitted log-slope over N=5..14  : {r:.4f}  (d_N ~ e^({r:.3f} N))")

# ---- the canonical fifth datum: indicator of G_eps --------------------------
print("\n=== 1_{G_eps} and continuous functionals along the gapless ray nu_h=(1-h)nu+h*delta_0 ===")
print("     h      floor   D8_L1(nu_h,nu)   M1        chi       1_{G_0.2}  sup  d_8(0.2)")
A8, b8 = data_moments(8)
for h in [1e-1, 1e-2, 1e-3, 1e-4, 1e-6, 1e-9]:
    wh = (1.0 - h) * w
    wh[0] += h
    floor = grid[np.argmax(wh > 0)]                    # inf supp (grid)
    bim = wh.sum()
    M1 = (grid * wh).sum() / bim
    chi = (wh[(grid > 0)] / grid[(grid > 0)]).sum() / bim
    l1 = np.abs(A8 @ wh - b8).sum() / bim
    supremum = np.abs(A8 @ wh - b8).max() / bim
    print(f"  {h:8.1e}  {floor:5.3f}   {l1:.3e}     {M1:.6f}  {chi:.6f}     {int(floor>=0.2):d}     {supremum:.2e}")
print("  -> every weak-* continuous datum moves by O(h); the floor jumps at h=0 only.")
print("  -> 1_{G_0.2} is 0 for every h>0 and 1 at h=0 : discontinuous exactly at the ray endpoint.")

# ---- two-sidedness summary --------------------------------------------------
print("\n=== SUMMARY (one-sidedness) ===")
print("  max mass below eps  > 0  for all finite N   => 'floor >= eps' unprovable from finite data")
print("  min mass below eps  = 0  for all finite N   => 'floor <  eps' unpinnable from finite data by itself")
print("  both directions fail at every finite N; only the infinite datum (or a discontinuity) decides.")
