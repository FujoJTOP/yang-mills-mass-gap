"""
blast_string_tension.py -- witness for BLAST-STRING-TENSION.md
Cheap, numpy only, one job.

KL model ν: atom f2 at m1^2, plus a continuum.
  chi = int s^{-1} dnu,  M = int dnu = f2 + continuum weight,  m1^2 = inf supp.
Checks:
  (1) Lemma G:  f2/m1^2 <= chi <= M/m1^2   (exact, [T]).
  (2) The atomless "bump at the floor" saturates chi -> M/m1^2 as width->0
      => the bound chi <= M/m1^2 is TIGHT at the gapless configuration.
  (3) H1 (chi <= 1/lam) is toggled by lam with IDENTICAL sigma,m1,M,chi
      => the string tension does NOT determine H1.
  (4) The derived sufficient condition:  H1  iff  lam <= 16 sigma / M  (with m1^2>=16 sigma).
"""
import numpy as np

m1sq = 1.0          # lightest mass squared = floor
M    = 1.0          # total weight

def chi_uniform(atom_f2, s_hi=6.0, n=400000):
    """atom f2 at m1sq  +  uniform continuum on [m1sq, s_hi] carrying (M-atom_f2)."""
    chi = atom_f2 / m1sq
    rest = M - atom_f2
    if rest <= 0:
        return chi
    s  = np.linspace(m1sq, s_hi, n + 1)
    sm = 0.5 * (s[:-1] + s[1:])
    w  = np.full(n, rest / n)
    chi += np.sum(w / sm)
    return chi

def chi_bump_at_floor(halfwidth, n=400000):
    """NO atom; atomless bump of half-width h at the floor, support [m1sq, m1sq+2h], total M.
       As h->0 this concentrates at the floor and saturates chi -> M/m1sq."""
    h  = halfwidth
    s  = np.linspace(m1sq, m1sq + 2*h, n + 1)
    sm = 0.5 * (s[:-1] + s[1:])
    w  = np.full(n, M / n)
    return np.sum(w / sm)

print("=" * 70)
print("(1) Lemma G:  f2/m1^2 <= chi <= M/m1^2     [T]")
print("=" * 70)
for f2 in [0.0, 0.1, 0.5, 0.9, 1.0]:
    chi = chi_uniform(f2)
    lo, hi = f2 / m1sq, M / m1sq
    ok = (lo - 1e-7 <= chi <= hi + 1e-7)
    print(f"  f2={f2:<4}: chi={chi:.5f}  lower={lo:.5f}  upper={hi:.5f}  OK={ok}")

print()
print("=" * 70)
print("(2) atomless bump at floor saturates chi -> M/m1^2  (tight at gapless)  [V-num]")
print("=" * 70)
for h in [0.5, 0.1, 0.01, 1e-3, 1e-4]:
    chi = chi_bump_at_floor(h)
    print(f"  halfwidth h={h:<7}: chi={chi:.6f}   (upper M/m1^2 = {M/m1sq:.6f})")

print()
print("=" * 70)
print("(3) H1 (chi<=1/lam) is decided by lam, NOT sigma: identical nu, two lam  [V-num]")
print("=" * 70)
chi = chi_uniform(0.5)
print(f"  fixed nu: chi={chi:.5f}  =>  H1 holds iff lam <= 1/chi = {1/chi:.5f}")
for lam in [1.0, 2.0]:
    print(f"    lam={lam}: chi<=1/lam ? {chi <= 1/lam + 1e-12}   (H1 {'holds' if chi<=1/lam+1e-12 else 'FAILS'})")
print("  => same sigma,m1,M,chi; H1 toggled by lam. sigma does not determine H1.")

print()
print("=" * 70)
print("(4) derived sufficient condition: H1 iff lam <= 16*sigma/M  (bridge m1^2>=16 sigma)  [C]")
print("=" * 70)
sigma_max = m1sq / 16.0          # bridge m1 >= 4 sqrt(sigma) saturated
print(f"  bridge m1>=4 sqrt(sigma) => sigma <= m1^2/16 = {sigma_max}")
print(f"  sufficient-for-H1 threshold: 16*sigma/M = {16*sigma_max/M:.5f}  (= m1^2/M)")
print(f"  actual H1 threshold (for this nu): 1/chi = {1/chi:.5f}")
print("  => the string-tension bound is CONSERVATIVE (16 sigma/M = 1.0 < 1/chi);")
print("     H1 follows if lam <= 16 sigma/M, a coarse-RG input sigma does not supply.")
