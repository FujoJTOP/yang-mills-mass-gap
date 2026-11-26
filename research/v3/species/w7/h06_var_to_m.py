# H06-VARM : can  int_0^2.5 Var(S_W) dbeta = 6V u_p(2.5) = 3.04V  be turned into m ?
# ---------------------------------------------------------------------------------
# d07's claims under test:
#   (I1)  sum_{n>=1} |<n|S_W|G>|^2 = Var(S_W)                 [claimed [T]]
#   (I2)  int Var dbeta = 6V u_p(2.5) = 3.04V                 [T, FTC]
#   (I3)  Jensen: int sqrt(Var) >= sqrt(T int Var) = 2.76 sqrt(V)   <-- SIGN
# This script checks (I1) exactly in a many-body toy, (I2) numerically, (I3)'s direction,
# and converts whatever survives into a statement about Delta (= m).
import numpy as np
from scipy.integrate import quad
from scipy.special import iv
from math import comb, log

def out(s): print(s, flush=True)
B1, B2 = 0.0034, 2.5
T = B2 - B1
VO = 1.0                       # Var_O(Omega) placeholder (local O), as in d07_g.py
BUDGET = np.log(64.0/1e-6)     # 17.974  (D04 (4.1), eps=1e-6, r0=1/64)

# ---------- exact SU(2) single-plaquette moments (identical convention to d07_g.py)
def _Z(b): return quad(lambda t: np.sin(t)**2*np.exp(b*np.cos(t)), 0, np.pi, epsabs=1e-14, limit=400)[0]
def mom(b, k):
    Z = _Z(b)
    return quad(lambda t: np.cos(t)**k*np.sin(t)**2*np.exp(b*np.cos(t)), 0, np.pi, epsabs=1e-14, limit=400)[0]/Z
def u_p(b):   return mom(b, 1)                    # <cos t> = <(1/2)Re tr U_p>
def var_s(b): return mom(b, 2) - mom(b, 1)**2     # s = 1 - cos t  =>  Var(s)=Var(cos t)

out("="*84)
out("[0] SU(2) one-plaquette: quadrature vs closed form <cos> = I1'/I1 - 1/b")
b = 2.5
cf = 0.5*(iv(0,b)+iv(2,b))/iv(1,b) - 1.0/b
out(f"    quad u_p(2.5)   = {u_p(2.5):.9f}")
out(f"    closed form     = {cf:.9f}      residual {abs(u_p(2.5)-cf):.2e}")
out(f"    Var(s)(2.5)     = {var_s(2.5):.9f}   Var(s)(0)= {var_s(0.0):.9f}  (=1/4 Haar)")
Ivar1 = quad(lambda x: var_s(x), 0.0, 2.5, epsabs=1e-13, limit=400)[0]
out(f"    (I2) int_0^2.5 Var(s) dbeta = {Ivar1:.9f}   vs  u_p(2.5)-u_p(0) = {u_p(2.5)-u_p(0):.9f}")
out(f"         residual {abs(Ivar1-(u_p(2.5)-u_p(0.0))):.2e}  => FTC identity VERIFIED [T]")
out("="*84)
out("[0b] CONVENTION WARNING: 0.5072 is the INDEPENDENT(single)-plaquette value.")
out("     int_0^2.5 Var(S_W) dbeta = 6V u_p^{4D}(2.5) EXACT, but the number 3.043V comes from")
out("     u_p^{single}=0.5072.  4D MC/independent differ at O(10%) => 3.04V is a PROXY [C]/[E].")
out("="*84)

# ---------- [1] test (I1) in an exactly solvable many-body model
# N independent 2-level plaquettes s_i in {0,1}, Gibbs weight p;  S_W = sum_i s_i
def toy(N, p, beta=1.0):
    from scipy.special import gammaln
    k = np.arange(N+1, dtype=float)
    lg = gammaln(N+1)-gammaln(k+1)-gammaln(N-k+1) + k*np.log(p) + (N-k)*np.log(1-p)
    lg -= lg.max(); w = np.exp(lg); w = w/w.sum()      # exact Gibbs weights in log space
    E = k                                            # spectrum of S_W: E_k = k
    mean = (w*E).sum(); var = (w*E**2).sum() - mean**2
    E0, p0 = E[0], w[0]
    num = (w[1:]*(E[1:]**2)).sum()                   # sum_{n>=1} |<n|S_W|G>|^2  (E0=0)
    return mean, var, num, E0, p0
out("[1] TEST OF (I1)  sum_{n>=1}|<n|S_W|G>|^2 = Var(S_W)   ?   (exact, N independent 2-level)")
out(f"{'N=6V':>8}{'p':>7}{'<S_W>':>11}{'Var(S_W)':>12}{'sum_n>=1':>12}{'ratio':>9}   {'<S_W>^2-p0 E0^2':>16}")
for N, p in [(6*4**4, 0.4928), (6*6**4, 0.4928), (6*8**4, 0.4928)]:
    mean, var, num, E0, p0 = toy(N, p)
    out(f"{N:>8}{p:>7.4f}{mean:>11.4g}{var:>12.6g}{num:>12.6g}{num/var:>9.4g}   {mean**2-p0*E0**2:>16.6g}")
out("    => (I1) is FALSE unless <Omega|S_W|G> = <S_W>; gap = <S_W>^2 - p0 E0^2 = O(V^2),")
out("       Var(S_W) = O(V)  => the sum is V-times larger than Var. Exact corrected identity:")
out("       sum_{n>=1}|<n|S_W|G>|^2 = Var(S_W) + <S_W>^2 - p0 E0^2      [T]")
out("")
out("    corollary: the sharp Cauchy-Schwarz numerator is sqrt(Var+<S_W>^2-p0E0^2) ~ <S_W> = 6V(1-u_p)")
out("    i.e. O(V) -- NOT O(sqrt V).  D04's 12V is right in scaling; d07's 'sqrt V sharpening' is not.")
out("="*84)

# ---------- [2] direction of the Jensen step (I3)
xs = np.linspace(0.0, 2.5, 4001)
f = np.array([np.sqrt(var_s(x)) for x in xs])
int_sqrt = np.trapezoid(f, xs)
cs_ub = np.sqrt(2.5*Ivar1)
out("[2] DIRECTION of int sqrt(Var) vs sqrt(T int Var)  (exact single-plaquette Var(beta))")
out(f"    int_0^2.5 sqrt(Var(s)) dbeta = {int_sqrt:.9f}")
out(f"    sqrt(2.5 * int Var)          = {cs_ub:.9f}     ratio {int_sqrt/cs_ub:.6f}")
out(f"    => int sqrt(Var) <= sqrt(T int Var): the inequality runs <=, Jensen(sqrt concave) is an UPPER bound.")
out(f"       d07 (2.3) writes >= : SIGN ERROR.  d07's 2.758 sqrt(V) is a ceiling, not a floor.")
out("="*84)

# ---------- [3] what the exact int Var buys, in the two numerator variants
V = {L: L**4 for L in (4, 6, 8, 16)}
ivar4 = 6.0*u_p(B2)*np.ones(4)                       # int Var = 6V u_p(2.5)  [T] (proxy)
out("[3] INT Var CONVERTED: certificate needs  int A/Delta <= 17.974   (A = sqrt(VO*NUM))")
out("    N1 = sqrt(int Var)          (d07 numerator, Var-only)      ~ 1.744 sqrt(V)")
out("    N2 = int sqrt(Var+<S_W>^2)  (correct numerator, ~<S_W>)     ~ 2.96 V")
out(f"{'L':>4}{'V':>9}{'N1=sqrt(intVar)':>18}{'N2~<S_W>':>12}{'reqDelta_min N1':>17}{'reqDelta_min N2':>17}{'req delta N1':>14}{'req delta N2':>14}")
rows = []
for L in (4, 6, 8, 16):
    v = V[L]
    N1 = np.sqrt(3.043*v)                 # = sqrt(int Var), VO=1
    mean = 6.0*v*(1.0-u_p(B2))
    N2 = np.sqrt(v*0.25 + mean**2)        # sqrt(Var + <S_W>^2), Var=6V*Var(s)(2.5)
    d1, d2 = N1/BUDGET, N2/BUDGET                            # sup-form: Delta_min >= N/budget
    i1 = np.sqrt((1.0/B1) - 1.0/B2)                          # sqrt(int dbeta/beta^2), Delta=beta*delta
    s1, s2 = N1*i1/BUDGET, N2*i1/BUDGET                      # L^2-form: delta >= N*sqrt(int1/b^2)/budget
    rows.append((L, v, N1, N2, d1, d2, s1, s2))
    out(f"{L:>4}{v:>9}{N1:>18.4g}{N2:>12.4g}{d1:>17.4g}{d2:>17.4g}{s1:>14.4g}{s2:>14.4g}")
out("    (req delta = delta_{S_W} with Delta(beta)=beta*delta ; strong coupling gives delta=O(1) [T])")
out("    => N1 route needs delta >= %.3g L^2  ;  N2 route needs delta >= %.3g L^4   (Var_O=1)"
    % (rows[0][6]/rows[0][0]**2, rows[0][7]/rows[0][0]**4))
out("    both violate delta=O(1) already at L=4  => (b): same loop, and the TRUE numerator makes it worse")
out("="*84)

# ---------- [4] the only direction Var actually constrains Delta: from ABOVE
out("[4] CAN WE BACK-SOLVE Delta FROM int Var ?  direction test")
out("    Var(beta) >= p1*Delta^2 with p1 ~ e^{-beta*Delta}/Nstates  =>  Delta <= (1/beta)(log(Nstates)+log(Delta^2/Var))")
out("    (entropy Nstates = e^{O(V)} swamps it)   and Var-heavy tail: linear-response Var = sum p_n (E_n-<E>)^2")
out("    at beta=2.5:  Var(S_W)/(6V)^2 = %.4f  (2%% level: NOT exponentially small)" % (6*var_s(2.5)/(6.0)**2))
out("    => no exponential smallness => the low-T inversion is unavailable => Var cannot pin m from below.")
out("    The only clean inequalities that exist run  Delta <= f(Var)  (upper bounds on the gap).")
out("="*84)

# ---------- [5] four-model acid test: does int Var separate m=0 from m>0 ?
out("[5] ACID TEST: int Var per volume, with the gap")
out("    4D SU(2) beta<=2.5 : int Var / V = %.4f      gap Delta = ? (target, unknown)      -> [C]" % (6*u_p(2.5)))
out("    2D YM SU(2) exact  : int Var / V = %.4f      gap Delta = 3/(2beta)=0.6  [T] > 0" % Ivar1)
out("    free Gaussian      : int Var / V = (3/2)(1/beta1-1/beta2) = %.4g  (Delta=2pi/L -> 0, m=0)"
    % (1.5*(1.0/B1-1.0/B2)))
out("    single plaquette   : int Var / V = %.4f      no transfer => Delta undefined" % (6*Ivar1/6.0))
out("    => int Var/V is O(1) in ALL models, INCLUDING the massless free field (m=0, f1^2=0).")
out("       It does not separate m=0 from m>0 => vacuous as a carrier of m. (same mechanism as d07 (c))")
out("="*84)
out("[6] SUMMARY")
out("    (I1) FALSE off the diagonal: sum_n>=1 = Var + <S_W>^2 - p0E0^2  (extra O(V^2))")
out("    (I2) TRUE (FTC of <S_W>); number 3.043V is a single-plaquette PROXY")
out("    (I3) direction reversed: 2.758 sqrt(V) is an UPPER bound")
out("    => int Var = 3.04V appears only as sqrt(V) (or, correctly, as <S_W>=O(V)) multiplying 1/Delta")
out("    => (b): same loop, quantitatively closed; no m-lower-bound.  (c) holds as the m-blindness leg.")
