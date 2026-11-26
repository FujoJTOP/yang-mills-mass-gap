"""RS8 — does a THIN (1D-along-the-tube) polymer geometry let the activity be large?

Witness for research/v3/species/RS8.md.  py -3.12 -I, numpy only, deterministic, <2 s.

Four numbers only:
 A. 1D entropy: #intervals of length n through a fixed point = n   (polynomial)
    4D entropy: #connected plaquette-animals of size n through a fixed plaquette
    = 1, 20, 458, 11132, 280175  (corpus C03-BYPASS:111, exact enumeration) -> ratio ~25
 B. KP threshold for |w(gamma)| = w^{|gamma|}:  max{ w : exists a>0 with
    S(a) = sum_{gamma contains x} |w(gamma)| e^{a|gamma|} <= a }   -- in 1D and in 4D.
 C. The exact 1D radius (hard rods of length l on Z): z_c(l) = -(l-1)^{l-1}/l^l,
    versus what KP guarantees:  z_KP(l) ~ 1/(e l^2).  Ratio, and the l -> infinity limit.
 D. The YM tube as a 1D object: |P| = 4T, w_T = u^{4T}, u = u(beta=2.5) = 0.625.
"""

import math

# ---------------------------------------------------------------- A
def intervals_through_point(n):
    """Number of intervals of length n in Z containing a fixed site: exactly n."""
    return n

ANIMALS_4D = [1, 20, 458, 11132, 280175]  # closed/connected plaquette animals through p0

print("=" * 72)
print("A. ENTROPY PER POINT  (the only place the dimension enters)")
print("=" * 72)
print(f"{'n':>2} {'1D: #intervals':>16} {'4D: #animals':>14} {'ratio 4D/1D':>12}")
for n in range(1, 6):
    a1, a4 = intervals_through_point(n), ANIMALS_4D[n - 1]
    print(f"{n:>2} {a1:>16} {a4:>14} {a4 / a1:>12.2f}")
g = [ANIMALS_4D[i + 1] / ANIMALS_4D[i] for i in range(len(ANIMALS_4D) - 1)]
print("4D growth ratios a_{n+1}/a_n :", [round(x, 2) for x in g], "-> rising, >= 20")
print("1D growth ratio  (n+1)/n     :", [round((n + 1) / n, 3) for n in range(1, 5)], "-> ->1")
print()

# ---------------------------------------------------------------- B
def kp_sum_1d(w, a, nmax=4000):
    """S(a) = sum_n n * w^n * e^{a n}  = w e^a / (1 - w e^a)^2 (closed form)."""
    x = w * math.exp(a)
    if x >= 1:
        return math.inf
    return x / (1 - x) ** 2

def kp_threshold_1d(lo=1e-9, hi=1.0):
    """sup w such that exists a>0 with S(a) <= a."""
    best = 0.0
    for i in range(2001):          # scan w
        w = lo + (hi - lo) * i / 2000
        ok = False
        for j in range(1, 3001):   # optimise a
            a = 12.0 * j / 3000
            if kp_sum_1d(w, a) <= a:
                ok = True
                break
        if ok:
            best = w
        else:
            break
    return best

def kp_threshold_4d(lo=1e-12, hi=0.05, extra_ratio=25.0):
    """S(a) = sum_n a_n (w e^a)^n with a_n extrapolated past n=5 by ratio 25."""
    an = list(ANIMALS_4D)
    while len(an) < 60:
        an.append(an[-1] * extra_ratio)
    def S(w, a):
        x, s = w * math.exp(a), 0.0
        for n, c in enumerate(an, start=1):
            s += c * x ** n
            if s > 1e12:
                return math.inf
        return s
    best = 0.0
    for i in range(2001):
        w = lo + (hi - lo) * i / 2000
        if any(S(w, a) <= a for a in [0.5 * j / 200 for j in range(1, 401)]):
            best = w
        else:
            break
    return best

w1 = kp_threshold_1d()
w4 = kp_threshold_4d()
print("=" * 72)
print("B. KP THRESHOLD  (|w(gamma)| = w^{|gamma|}, per-point form)")
print("=" * 72)
print(f"1D  (polynomial entropy)      w_KP <= {w1:.6f}")
print(f"4D  (lambda ~ 25 per step)    w_KP <= {w4:.6e}")
print(f"ratio 1D/4D                   {w1 / w4:.3e}")
print("  -> the 4D threshold is smaller by the 4D entropy constant, not by a constant")
print()

# ---------------------------------------------------------------- C
def zc_rods(l):
    """Nearest singularity of the 1D hard-rod gas (rods of length l on Z):
    Xi_L = Xi_{L-1} + z Xi_{L-l};  branch point at z_c = -(l-1)^{l-1}/l^l  (l >= 2).

    l = 1 is the DEGENERATE case: rods of length 1 exclude nothing, Xi_L = (1+z)^L is a
    polynomial, the activity expansion is finite, radius = infinity.  It is the ONLY
    l with infinite radius -- i.e. 'radius = infinity' <=> bounded support of w."""
    if l == 1:
        return None                      # entire: (1+z)^L
    return -((l - 1) ** (l - 1)) / (l ** l)

def z_kp_rods(l):
    """KP for rods: rods of length l containing a fixed site <= l, weight z,
    condition l * z * e^{a l} <= a  ->  z <= a e^{-a l}/l, maximise over a: a=1/l:
    z_KP = 1/(e l^2)."""
    return 1.0 / (math.e * l * l)

print("=" * 72)
print("C. 1D EXACT RADIUS vs WHAT KP GUARANTEES  (hard rods, length l)")
print("=" * 72)
print(f"{'l':>2} {'|z_c| exact':>14} {'z_KP = 1/(e l^2)':>18} {'exact/KP':>10}")
for l in range(1, 9):
    zc, zk = zc_rods(l), z_kp_rods(l)
    if zc is None:
        print(f"{l:>2} {'inf (entire)':>14} {zk:>18.8f} {'inf':>10}   <- degenerate: w has bounded support")
    else:
        print(f"{l:>2} {abs(zc):>14.8f} {zk:>18.8f} {abs(zc) / zk:>10.2f}")
print("asymptotics (l>=2):  |z_c| ~ 1/(e l)   vs   z_KP ~ 1/(e l^2)  -> both -> 0")
print("=> the 1D radius is FINITE, shrinks with the polymer length, and is only")
print("   a factor ~l better than KP.  'l -> infinity' makes it SMALLER, not larger.")
print()

# ---------------------------------------------------------------- D
print("=" * 72)
print("D. THE YM TUBE AS A 1D OBJECT  (|P| = 4T, w_T = u^{4T})")
print("=" * 72)
beta = 2.5
u = beta / 4.0                      # SU(2) exact per-plaquette activity
print(f"beta = {beta},  u = beta/4 = {u}")
print(f"1D-along-tube test  u^4 = {u**4:.6f}  (<1  =>  geometric sum in T converges)")
print(f"rate                  4 ln(1/u) = {4 * math.log(1 / u):.6f}   (= m1 a of P3)")
print("=> along the tube the tube-series DOES converge at beta=2.5 (needs only u<1).")
print("   The failure is NOT the T-direction: it is the 4D decoration entropy,")
print(f"   whose per-point sum at beta=2.5 exceeds the KP bound by 35.5x (R4-PERCENT).")
print()
print("per-point KP sum, |X|=6 minimal closed polymer, beta=2.5 (corpus R4):")
phi = 6.809423e-02
print(f"  32 * phi = {32 * phi:.6f}  vs  best a = 1/6 = {1/6:.6f}  ->  {32 * phi / (1 / 6):.2f}x over")
print("=" * 72)
