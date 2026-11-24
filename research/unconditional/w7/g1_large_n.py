"""
g1_large_n.py -- G1 at large N: is there an IR fixed point for pure YM at N=infinity?

't Hooft-scaled 2-loop beta for SU(N) with N_f = x*N fundamentals (large-N leading):
  lambda := g^2 N,
  beta_lambda(lam,x) = -(2/3)[ (11-2x) lam^2/(16 pi^2) + (34-12x) lam^3/(16 pi^2)^2 ].
(The N factors cancel in the 't Hooft limit; derived from
 beta_g = -b0 g^3/(16pi^2) - b1 g^5/(16pi^2)^2,  b0=(11N-2Nf)/3,  b1=34N^2/3 - 12 N Nf/3 + O(1).)

Sign of beta_lambda for lam>0: controlled by (11-2x) + (34-12x) lam/(16 pi^2).
  - Pure YM  x=0: b0,b1 > 0  =>  beta_lambda < 0 for all lam>0  (no zero).
  - IR fixed point  lam* = 16 pi^2 (11-2x)/(12x-34) > 0  iff  34/12 < x < 11/2,
    i.e. the (2-loop) conformal window N_f/N in (2.833..., 5.5).
Pure YM N_f/N = 0 sits FAR BELOW the window's lower edge 34/12 ~ 2.833.

NOT a Yang-Mills statement: 2-loop perturbation theory only.
"""
import numpy as np

TWO_PI2 = 16 * np.pi**2   # 16 pi^2

def beta_lambda(lam, x):
    b0 = (11 - 2*x)      # b0*N-scaled (factor 3 and 16pi^2 in the formula)
    b1 = (34 - 12*x)
    return -(2.0/3.0) * (b0 * lam**2 / TWO_PI2 + b1 * lam**3 / TWO_PI2**2)

def zero_of_beta(x):
    """Positive fixed point of the 2-loop beta_lambda, or None."""
    # 11-2x + (34-12x) lam/(16pi^2) = 0
    denom = (34 - 12*x)
    if abs(denom) < 1e-12:
        return None
    lam = -TWO_PI2 * (11 - 2*x) / denom
    return lam if lam > 0 else None

print("="*76)
print("'t Hooft 2-loop beta_lambda; pure YM x=N_f/N=0 vs the conformal window")
print("="*76)
xs = [0.0, 2.0, 3.0, 4.0, 5.0, 6.0]
print(f"{'x=N_f/N':>8} | {'b0-sign':>7} {'b1-sign':>7} | {'fixed pt lambda*':>16} | sign(beta) for large lam")
for x in xs:
    z = zero_of_beta(x)
    b0s = '+' if (11-2*x) > 0 else ('0' if abs(11-2*x)<1e-9 else '-')
    b1s = '+' if (34-12*x) > 0 else ('0' if abs(34-12*x)<1e-9 else '-')
    zstr = f"{z:.4f}" if z is not None else "none"
    # large-lam sign (lam^3 term)
    big = np.sign(beta_lambda(1e3, x))
    print(f"{x:>8.2f} | {b0s:>7} {b1s:>7} | {zstr:>16} | {big:+.0f}")

print()
print("window (2-loop): 34/12 < x < 11/2  =  %.4f < x < %.4f" % (34/12, 11/2))
print("pure YM x=0 is %s the window (lower edge %.4f)" % ("BELOW" if 0 < 34/12 else "in", 34/12))

print()
print("="*76)
print("beta_lambda vs lam for pure YM (x=0) -- must be < 0 for ALL lam>0")
print("="*76)
lams = np.array([0.1, 1.0, 5.0, 10.0, 30.0, 100.0])
bl = beta_lambda(lams, 0.0)
for l, b in zip(lams, bl):
    print(f"  lam={l:>7.2f}  beta_lambda={b:+.6e}  {'<0 OK' if b<0 else 'NONNEG(!)'}")
print("  -> all negative:", bool(np.all(bl < 0)), " ; any zero:", zero_of_beta(0.0))

print()
print("="*76)
print("inside the window (x=4): beta crosses zero at lam* > 0 (a BZ fixed point)")
print("="*76)
z4 = zero_of_beta(4.0)
print(f"  x=4: beta_lambda(lam*)={beta_lambda(z4,4.0):+.3e} at lam*={z4:.4f}  (should be ~0)")
print(f"       beta_lambda just below: {beta_lambda(0.9*z4,4.0):+.3e}  just above: {beta_lambda(1.1*z4,4.0):+.3e}")

print()
print("VERDICT (2-loop only): pure YM (N_f/N=0) has beta_lambda<0 for all lam>0, no IR fixed point.")
print("The conformal window exists for N_f/N in (2.833, 5.5); pure YM is far below it.")
print("NOT a proof: 2-loop perturbation theory. Nonperturbative 'no fixed point at N_f=0' = large-N confinement (open).")
