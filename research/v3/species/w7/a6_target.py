"""a6_target.py -- correct beta targets for a startable continuum limit.

Two-loop perturbative scaling (published ASYMPTOTIC-SCALING form):
    b0 = 11N/(3 (4 pi)^2),  b1 = 34N^2/(3 (4 pi)^4)  =>  b1/(2 b0^2) = 51/121 (all N)
    ln(xi/a) = 1/(2 b0 g^2) - (51/121) ln(b0 g^2) + C
with Wilson normalization beta = 2N/g^2  ->  g^2 = 4/beta (SU(2)), 6/beta (SU(3))
    => 1/(2 b0 g^2) = 12 pi^2 beta/(11 N^2)

Anchors (all from the corpus, each traceable):
  SU(2)  beta=2.5 : xi/a = 1.0 (0.53 .. 1.30)          A4-SCALE-SEP / tm_gap CSV
  SU(3)  beta=6.0 : xi/a = 0.93(8);  beta=6.30 : 1.13(20)   lattice-cpb REPORT
  SU(3)  r0 ladder : Necco-Sommer r0/a at beta=6.0/6.2/6.4
                     -> xi/a = a^-1/m(0++), m=1.71 GeV, r0=0.5 fm
No dependencies.  Pure stdlib.
"""
import math

pi2 = math.pi ** 2
HBARC = 0.1973269804  # GeV fm


def b0(N):
    return 11.0 * N / (48.0 * pi2)


def gsq(beta, N):
    return 2.0 * N / beta


def t_of(beta, N):
    """1/(2 b0 g^2) -- the exponential argument."""
    return 12.0 * pi2 * beta / (11.0 * N * N)


def log_xi(beta, N, C):
    return t_of(beta, N) - (51.0 / 121.0) * math.log(b0(N) * gsq(beta, N)) + C


def xi(beta, N, C):
    return math.exp(log_xi(beta, N, C))


def calib(beta_a, xi_a, N):
    return math.log(xi_a) - t_of(beta_a, N) + (51.0 / 121.0) * math.log(b0(N) * gsq(beta_a, N))


def solve(beta_a, xi_a, target, N, lo=None, hi=None):
    C = calib(beta_a, xi_a, N)
    lo = beta_a if lo is None else lo
    hi = beta_a + 0.5 if hi is None else hi
    while xi(hi, N, C) < target:
        hi += 0.5
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if xi(mid, N, C) < target:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi), C


print("=" * 72)
print("[0] convention check: g^2 = 4/beta (SU(2)), 6/beta (SU(3)); 1/(2 b0 g^2)")
b2 = b0(2)
b3 = b0(3)
print(f"   SU(2) b0 = {b2:.8f}   SU(3) b0 = {b3:.8f}")
print(f"   beta=6.0 (SU3): 1/(2 b0 g^2) = {t_of(6.0,3):.6f}   [corpus: 7.18]")
print(f"   slope dln(xi/a)/dbeta one-loop SU(3) = 12 pi^2/(11 N^2) = {12*pi2/(11*9):.6f}"
      f"   [corpus 1.196]")
print(f"   slope dln(xi/a)/dbeta one-loop SU(2) = {12*pi2/(11*4):.6f}")
print(f"   two-loop prefactor exponent b1/(2b0^2) = 51/121 = {51/121:.6f}")

print("=" * 72)
print("[1] SU(2): anchor xi/a(2.5) = 1.0   (domain 0.53 .. 1.30)")
for lab, xa in (("central", 1.0), ("low", 0.53), ("high", 1.30)):
    C = calib(2.5, xa, 2)
    b3t, _ = solve(2.5, xa, 3.0, 2)
    b7t, _ = solve(2.5, xa, 7.0, 2)
    print(f"   {lab:8s} xi(2.5)={xa:5.2f} -> beta(xi=3)={b3t:.3f}  beta(xi=7)={b7t:.3f}"
          f"  g^2={4/b3t:.3f}/{4/b7t:.3f}")
print("   cross-check: measured xi/a ladder at beta=2.5 (four rows of A4-SCALE-SEP):")
print(f"      beta=2.5 xi=1.30 -> beta(xi=3)={solve(2.5,1.30,3.0,2)[0]:.3f}, "
      f"beta(xi=7)={solve(2.5,1.30,7.0,2)[0]:.3f}  (upper edge of the corpus domain)")
# deconfinement-ladder slope cross-check
bc4, bc6 = 2.2986, 2.43
slope_dec = math.log(6.0 / 4.0) / (bc6 - bc4)
print(f"   cross-check slope from deconfinement ladder (N_t=4,6): {slope_dec:.3f}/beta"
      f"  [two-loop {12*pi2/(11*4):.3f}]")
for xa in (1.0,):
    for tgt in (3.0, 7.0):
        print(f"      empirical slope route: beta(xi={tgt:g}) = "
              f"{2.5 + math.log(tgt/xa)/slope_dec:.3f}")

print("=" * 72)
print("[2] SU(3): anchor A = own lattice data  xi/a(6.0)=0.93, (6.30)=1.13")
C3 = calib(6.0, 0.93, 3)
print(f"   two-loop prediction xi(6.30) = {xi(6.30,3,C3):.3f}   [REPORT: 1.36; measured 1.13]")
for tgt in (3.0, 7.0):
    bt, _ = solve(6.0, 0.93, tgt, 3)
    print(f"   target xi/a={tgt:g}: beta={bt:.3f}  g^2={6/bt:.3f}")
print(f"   with the lower anchor 0.93-0.08: beta(xi=3)={solve(6.0,0.85,3.0,3)[0]:.3f},"
      f" beta(xi=7)={solve(6.0,0.85,7.0,3)[0]:.3f}")
print(f"   with the upper anchor 0.93+0.08: beta(xi=3)={solve(6.0,1.01,3.0,3)[0]:.3f},"
      f" beta(xi=7)={solve(6.0,1.01,7.0,3)[0]:.3f}")

print("=" * 72)
print("[3] SU(3): anchor B = r0 ladder (Necco-Sommer r0/a) + m(0++)=1.71 GeV, r0=0.5 fm")
r0_inv = HBARC / 0.5  # GeV
r0a = {6.0: 5.368, 6.2: 7.380, 6.4: 9.660}
m0 = 1.71
lad = {b: v * r0_inv / m0 for b, v in r0a.items()}
for b in sorted(lad):
    print(f"   beta={b:.1f}: r0/a={r0a[b]:6.3f}  a^-1={r0a[b]*r0_inv:5.3f} GeV"
          f"  xi/a={lad[b]:.3f}   (xi_sigma/a={(r0a[b]*r0_inv)/(0.44*math.sqrt(0.44*0+1)):.2f} n/a)")
C_B = calib(6.0, lad[6.0], 3)
print(f"   two-loop calibrated at beta=6.0 -> predicts xi(6.2)={xi(6.2,3,C_B):.3f}"
      f" (ladder {lad[6.2]:.3f}), xi(6.4)={xi(6.4,3,C_B):.3f} (ladder {lad[6.4]:.3f})")
for tgt in (3.0, 7.0):
    bt, _ = solve(6.0, lad[6.0], tgt, 3)
    print(f"   target xi/a={tgt:g}: beta={bt:.3f}  g^2={6/bt:.3f}")
# pure empirical log-linear slope through the three ladder points
import statistics as st
ys = [math.log(lad[b]) for b in (6.0, 6.2, 6.4)]
xs = [6.0, 6.2, 6.4]
mx, my = st.mean(xs), st.mean(ys)
s = sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / sum((x - mx) ** 2 for x in xs)
print(f"   empirical log-linear slope from the ladder: {s:.3f}/beta  [two-loop 1.196]")
for tgt in (3.0, 7.0):
    print(f"   empirical route: beta(xi={tgt:g}) = {6.0 + (math.log(tgt)-math.log(lad[6.0]))/s:.3f}")

print("=" * 72)
print("[4] distances (new, like-for-like): target / beta_c")
bc2, bc3 = 2.2986, 5.6924
for lab, bt, bc in (("SU(2) xi=3", solve(2.5, 1.0, 3.0, 2)[0], bc2),
                    ("SU(2) xi=7", solve(2.5, 1.0, 7.0, 2)[0], bc2),
                    ("SU(3) xi=3", solve(6.0, 0.93, 3.0, 3)[0], bc3),
                    ("SU(3) xi=7", solve(6.0, 0.93, 7.0, 3)[0], bc3)):
    print(f"   {lab}: beta={bt:.2f}  beta/beta_c={bt/bc:.2f}  "
          f"(beta+0.1 -> {bt+0.1:.2f}: ratio {(bt+0.1)/bc:.2f})")
print("   OLD target beta=2.5: SU(2) 2.5/2.2986=%.3f, SU(3) 2.5/%.4f=%.3f"
      % (2.5 / bc2, bc3, 2.5 / bc3))
print("   g^2 at old point: SU(2) %s, SU(3) %s" % (4 / 2.5, 6 / 2.5))
