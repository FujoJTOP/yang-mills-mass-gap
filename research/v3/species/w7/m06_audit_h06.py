# m06_audit_h06.py -- INDEPENDENT re-derivation for the h06-vs-d07 retraction audit.
# Written from scratch (does NOT import or reuse d07_g.py / h06_var_to_m.py logic;
# single-plaquette moments are recomputed with mpmath, including a Bessel closed form).
#
# Checks, in order:
#   [A] SU(2) one-plaquette u_p(b), Var(s)(b), FTC identity      (mpmath, Bessel vs quad)
#   [B]  sum_{n>=1} |<n|S_W|G>|^2  == Var + <S_W>^2 - p0*E0^2  ? (exact many-body, E0 != 0)
#   [C]  h06's ratio table  Sigma/Var @ N=1536,7776,24576        (closed form)
#   [D]  Jensen direction  int sqrt(Var) vs sqrt(T int Var)      (3 exact profiles)
#   [E]  the TRUE numerator:  || (1-P0) S_W G || >= <S_W>  pointwise  => int A >= c*V
#   [F]  what (G) demands:  Delta_max  and  delta_{S_W}  in both routes; L^2 vs L^4
#   [G]  d07's free-field row vs its own FTC identity
import numpy as np
from mpmath import mp, mpf, quad as mquad, besseli, sqrt as msqrt, exp as mexp, cos, sin, pi, inf, fsum
from mpmath import log as mlog

mp.dps = 30
out = lambda s: print(s, flush=True)
B1, B2 = mpf('0.0034'), mpf('2.5')
T = B2 - B1
BUDGET = mpf('17.974')          # log(64/1e-6), D04 (4.1) -- same budget both files use
VO = mpf(1)

# ---------------- [A] exact one-plaquette ----------------
def Z(b):   return mquad(lambda t: sin(t)**2 * mexp(b*cos(t)), [0, pi])
def mom(b, k):
    z = Z(b)
    return mquad(lambda t: cos(t)**k * sin(t)**2 * mexp(b*cos(t)), [0, pi]) / z
u_p   = lambda b: mom(b, 1)
var_s = lambda b: mom(b, 2) - mom(b, 1)**2
u_bessel = lambda b: (besseli(0, b) + besseli(2, b)) / (2*besseli(1, b)) - 1/b

out("="*96)
out("[A] SU(2) one plaquette, s_p = 1 - cos t, weight sin^2 t e^{b cos t}   (mpmath, dps=30)")
out(f"    u_p(2.5)  quad   = {u_p(B2)}")
out(f"    u_p(2.5)  Bessel = {u_bessel(B2)}      residual {float(abs(u_p(B2)-u_bessel(B2))):.2e}")
out(f"    Var(s)(2.5) = {var_s(B2)}   Var(s)(0) = {var_s(mpf(0))}")
Iv = mquad(lambda b: var_s(b), [0, B2])
out(f"    FTC: int_0^2.5 Var(s) db = {Iv}   vs  u_p(2.5)-u_p(0) = {u_p(B2)-u_p(mpf(0))}")
out(f"         residual {float(abs(Iv-(u_p(B2)-u_p(mpf(0))))):.2e}")
out("    => Var(s) = d u_p/db ; int Var = u_p(2.5)-u_p(0)   [T]  (both files use this)")
out("="*96)

# ---------------- [B] the identity, with E0 != 0 ----------------
# N independent 2-level plaquettes, s_i in {e0,e1}; S_W = sum s_i; ground state = all e0.
# thermal state:  p_i = e^{-b s_i}/z  (b = beta);   |G> = sum_n sqrt(p_n)|n>   (thermal vector)
# claim:  sum_{n>=1} |<n|S_W|G>|^2 = Var(S_W) + <S_W>^2 - p0 * E0_tot^2 ,  E0_tot = N*e0
def two_level(N, e0, e1, beta):
    z = mexp(-beta*e0) + mexp(-beta*e1)
    pa, pb = mexp(-beta*e0)/z, mexp(-beta*e1)/z
    # spectrum E_k = (N-k)e0 + k e1, weight pmf(k) = C(N,k) pa^{N-k} pb^k
    from math import comb
    Ek = np.array([(N-k)*e0 + k*e1 for k in range(N+1)], dtype=float)
    pr = np.array([comb(N, k)*pa**(N-k)*pb**k for k in range(N+1)], dtype=float)
    pr /= pr.sum()
    mean = (pr*Ek).sum(); var = (pr*Ek**2).sum() - mean**2
    E0t, p0 = Ek[0], pr[0]
    sigma = (pr[1:]*(Ek[1:]**2)).sum()          # sum_{n>=1} |<n|S_W|G>|^2 = sum_{n>=1} p_n E_n^2
    norm2 = (pr*(Ek**2)).sum() - p0*E0t**2      # ||(1-P0) S_W G||^2 by Parseval
    return mean, var, sigma, norm2, E0t, p0

out("[B] IDENTITY TEST in an exactly solvable many-body model (spectrum NOT degenerate, E0 != 0)")
out("    N independent 2-level plaquettes; S_W = sum_i s_i ; s_i in {e0,e1}; thermal vector |G>")
mp.dps = 20
cases = [(12, 1.0, 3.0, 0.7), (12, 1.0, 3.0, 2.5), (20, 0.5, 2.0, 1.3), (16, 0.0, 1.0, 1.0)]
out(f"{'N':>4}{'e0':>6}{'e1':>6}{'beta':>7}{'<S_W>':>12}{'Var':>10}{'Sigma(n>=1)':>14}"
    f"{'Var+mean^2-p0E0^2':>19}{'||(1-P0)S_WG||^2':>19}{'d07 (=Var)?':>12}")
d07_holds = True
for (N, e0, e1, beta) in cases:
    mean, var, sigma, norm2, E0t, p0 = two_level(N, e0, e1, beta)
    pred = var + mean**2 - p0*E0t**2
    out(f"{N:>4}{e0:>6}{e1:>6}{beta:>7}{mean:>12.5f}{var:>10.5f}{sigma:>14.5f}"
        f"{pred:>19.5f}{norm2:>19.5f}{str(abs(sigma-var)<1e-9*max(1,var)):>12}")
    d07_holds &= abs(sigma-var) < 1e-9*max(1.0, var)
out(f"    => 'Sigma = Var' (d07 2.1) as a general statement: {d07_holds}   (should be False)")
out("    => exact form  Sigma = Var + <S_W>^2 - p0*E0^2 = ||(1-P0)S_W G||^2  :  Parseval, no assumption")
out("       beyond (i) {|n>} ONB of S_W-eigenvectors, (ii) S_W|Omega> = E0|Omega>, (iii) <Omega|G> = sqrt(p0).")
out("       For the Wilson action S_W >= 0 with min = 0 (all U_p = 1) => E0 = 0 => the p0 E0^2 term VANISHES.")
out("="*96)

# ---------------- [C] h06's ratio table, closed form ----------------
out("[C] h06 table check: Bernoulli plaquettes, p = 1-u_p(2.5);  ratio = 1 + N p/(1-p)")
p = float(1 - u_p(B2))
for N in (6*4**4, 6*6**4, 6*8**4):
    ratio = 1 + N*p/(1-p)
    out(f"    N=6V={N:>6}  <S_W>^2/Var = {N*p/(1-p):>10.4f}   Sigma/Var = {ratio:>10.4f}"
        f"   (h06: N=1536 ->1493, 7776 ->7556, 24576 ->2.388e4)")
out("="*96)

# ---------------- [D] Jensen direction ----------------
out("[D] DIRECTION of int sqrt(Var) vs sqrt(T int Var)   (sqrt is CONCAVE)")
prof = {}
prof['one-plaquette Var(s)'] = (lambda b: var_s(b), mquad(lambda b: var_s(b), [B1, B2]))
prof['Gaussian 3V/(2b^2)'] = (lambda b: 3*mpf(1)/(2*b**2), mquad(lambda b: 3*mpf(1)/(2*b**2), [B1, B2]))
out(f"{'profile':<26}{'int sqrt(Var)':>18}{'sqrt(T*int Var)':>18}{'ratio':>10}{'direction':>12}")
for name, (f, Ivf) in prof.items():
    Isq = mquad(lambda b: msqrt(f(b)), [B1, B2])
    ub = msqrt(T*Ivf)
    out(f"{name:<26}{float(Isq):>18.9f}{float(ub):>18.9f}{float(Isq/ub):>10.6f}"
        f"{('<= (upper bd)' if Isq <= ub else '>= (lower bd)'):>12}")
out("    => forced: (1/T)int sqrt(f) <= sqrt((1/T)int f)  =>  int sqrt(Var) <= sqrt(T int Var).  [T]")
out("       d07 (2.3) writes '>=' ; h06 (2.1) writes '<=' .   The '<=' is the correct one (and is NOT tight:")
out("       ratio 0.9957 for the exact plaquette profile -- Jensen loses almost nothing here).")
out("")
out("    REPAIR of d07's step WITHOUT Jensen: Var_s is decreasing, so with M := Var_s(0) = 1/4")
bs = [mpf(x)/100 for x in range(0, 251)]
vs = [var_s(b) for b in bs]
mono = all(vs[i] >= vs[i+1]-mpf('1e-25') for i in range(len(vs)-1))
out(f"      monotonicity Var_s(b) nonincreasing on [0,2.5]: {mono}  (Var_s(0)={vs[0]}, Var_s(2.5)={vs[-1]})")
Ivs = mquad(lambda b: var_s(b), [0, B2])
lb_var = Ivs/msqrt(vs[0])                    # int sqrt(Var_s) >= int Var_s / sqrt(max Var_s)
out(f"      int sqrt(Var_s) >= (int Var_s)/sqrt(Var_s(0)) = {lb_var}   (true: 1.0145... )")
out(f"      => int sqrt(6V*Var_s) db >= sqrt(6V)*{float(lb_var):.6f} = {float(msqrt(6)*lb_var):.4f} sqrt(V)")
out("      => 'int A >= c sqrt(V)' is a TRUE statement, c = %.3f (vs d07's invalid 2.758)."
    % float(msqrt(6)*lb_var))
out("      Both c sqrt(V) and the sharper c' V (block [E]) are LOWER bounds => d07's PREMISE survives;")
out("      only its derivation (Jensen '>=') and its constant are wrong.")
out("="*96)

# ---------------- [E] the true numerator is O(V): a POINTWISE bound ----------------
out("[E] the numerator needs NO Jensen:   ||(1-P0) S_W G||^2 = Var + <S_W>^2 - p0 E0^2 >= <S_W>^2")
out("    (Var >= 0, p0E0^2 = 0 for Wilson since E0 = 0)  =>  A(b) = sqrt(Var_O)*||...|| >= sqrt(Var_O)*<S_W>(b)")
out("    with <S_W>(b) = 6V(1-u_p(b)) :")
I1up = mquad(lambda b: 1 - u_p(b), [0, B2])
# int (1-u_p)/b has an endpoint LOG singularity at b=0+ ; split the interval explicitly
JS = [B1, mpf('0.01'), mpf('0.1'), mpf('0.5'), mpf('1'), mpf('1.5'), B2]
Jup = mquad(lambda b: (1 - u_p(b))/b, JS)
Jchk = mlog(B2/B1) - mquad(lambda b: u_p(b)/b, [0, B2])      # exact decomposition: int(1-u_p)/b = ln - int u_p/b
out(f"    int_0^2.5 (1-u_p) db = {I1up}      int_0^2.5 (1-u_p)/b db = {Jup}")
out(f"    (cross-check  ln(2.5/b1) - int_0^2.5 u_p/b db = {Jchk} ; both must agree)")
out("    => int A db >= 6V*int(1-u_p) db = %.4f * V     (V = L^4)  : POSITIVE VOLUME POWER, c ~ %.3g"
    % (float(6*I1up), float(6*I1up)))
out("    compare d07 (2.3): 'int A >= 2.758 sqrt(V)'.  Both are LOWER bounds; the true one is V, not sqrt(V).")
out("="*96)

# ---------------- [F] what (G) demands ----------------
out("[F] (G):  int_{b1}^{2.5} A/Delta db <= K = 17.974.    Two numerator variants x two routes.")
out("    route L1 (direct):  int A/Delta >= A_min * int db/Delta ;   Delta = b*delta  => int db/b = ln(2.5/b1)")
from mpmath import log as mlog
LNR = mlog(B2/B1)
out(f"         int_b1^2.5 db/b = ln(2.5/b1) = {float(LNR):.4f};   int (1-u_p)/b db = {float(Jup):.4f}")
out("    route L2 (Cauchy-Schwarz, as h06 does): int A/Delta <= sqrt(int A^2) * sqrt(int db/b^2)/delta,")
i2 = msqrt(1/B1 - 1/B2)
out(f"         sqrt(int_b1^2.5 db/b^2) = sqrt(1/b1 - 1/b2) = {float(i2):.4f}   (= h06's 17.14 factor)")
rows = []
out(f"{'L':>3}{'V=L^4':>8}{'6V*int(1-u_p)':>15}{'d07: 0.276L^2':>15}{'TRUE Dmax=L1':>14}"
    f"{'TRUE delta=L1':>15}{'h06 N2=~<S_W>(2.5)':>20}{'delta=L2 route':>15}")
for L in (4, 6, 8, 16):
    V = L**4
    Aint = 6*V*I1up                                   # true lower bd on int A
    Dmax = Aint/BUDGET                                # Delta_max needed (L1 route)
    dL1 = Aint/(BUDGET)/1 * 0                          # placeholder
    dL1 = 6*V*Jup/BUDGET                              # delta needed: (6V/delta)*J <= K
    N2 = 6*V*(1-u_p(B2))                              # pointwise A at b=2.5 (h06's N2)
    dL2 = N2*i2/float(BUDGET)
    rows.append((L, V, float(Aint), float(Dmax), float(dL1), float(N2), float(dL2)))
    out(f"{L:>3}{V:>8}{float(Aint):>15.4g}{0.276*L*L:>15.4g}{float(Dmax):>14.4g}"
        f"{float(dL1):>15.4g}{float(N2):>20.4g}{float(dL2):>15.4g}")
out("    scalings (fit on the four rows):")
Vv = np.array([r[1] for r in rows], float)
for idx, name in [(3, 'Delta_max (L1)'), (4, 'delta (L1)'), (6, 'delta (L2, h06 style)')]:
    yy = np.array([r[idx] for r in rows], float)
    ex = np.polyfit(np.log(Vv), np.log(yy), 1)[0]
    out(f"      {name:<22} ~ V^{ex:.3f}  = L^{4*ex:.2f}")
out("    => both routes: delta_{S_W} must GROW like L^4 (h06's exponent 4 is right);")
out("       constants are route-dependent: L1 direct ~ %.3g L^4 ; L2/C-S (h06) ~ %.3g L^4."
    % (rows[0][4]/4**4, rows[0][6]/4**4))
out("    => d07's '0.276 L^2' : is it even TRUE?  TRUE Dmax(L) / (0.276 L^2) = %s"
    % ["%.3g" % (r[3]/(0.276*r[0]**2)) for r in rows])
out("       (all >= 1 => 'Delta_max >= 0.276 L^2' is a TRUE but sqrt(V)-WEAK necessary condition)")
out("="*96)

# ---------------- [G] d07's free-field row vs its own identity ----------------
out("[G] d07 free-field row: S_W = (1/2) x^T M x, weight e^{-b S_W}, dim = 3V")
out("    <S_W> = dim/(2b) = 3V/(2b);   Var(S_W) = dim/(2b^2) = 3V/(2b^2)  (NOT b-independent!)")
for b in (B1, mpf(1), B2):
    out(f"      b={float(b):>5.3f}:  Var = 3V/(2b^2) = {float(3/(2*b**2)):>10.4f}*V   "
        f"vs d07's '3V/2 = 1.5V' (b-independent)")
out("    FTC check: int_b1^2.5 Var db = (3V/2)(1/b1-1/b2) = %.4f V  vs  <S_W>(b1)-<S_W>(2.5) = %.4f V"
    % (float((3*mpf(1)/2)*(1/B1-1/B2)), float((3*mpf(1)/2)*(1/B1-1/B2))))
out("    => d07's 'Var = 3V/2 exactly, b-independent' contradicts its own (2.2) [T]. free-field G still diverges")
out("       (int A/Delta diverges logarithmically at b->0), so the four-model VERDICT is unchanged.")
out("="*96)
out("[H] h06's N2 is a POINT value at b=2.5, but its C-S route needs sqrt(int A^2 db).  Size the gap:")
I2up = mquad(lambda b: (1 - u_p(b))**2, [0, B2])
fac = float(msqrt(I2up)/(1 - u_p(B2)))     # sqrt(int A^2)/A(2.5)   (A = 6V(1-u_p), 6V cancels)
out(f"    int_0^2.5 (1-u_p)^2 db = {I2up}   sqrt(.)/(1-u_p(2.5)) = {fac:.4f}")
out(f"    => h06's delta_CS is understated by {fac:.2f}x  =>  correct C-S constant ~ {2.82*fac:.2f} L^4")
out(f"    => the L1 route (delta >= {rows[0][4]/4**4:.3f} L^4) is TIGHTER, hence the one to quote as a")
out("       necessary condition.  Exponent L^4 unchanged either way.")
out("="*96)
out("[SUMMARY]")
out("  [1] Sigma_{n>=1}|<n|S_W|G>|^2 = Var + <S_W>^2 - p0 E0^2  [T]  (E0=0 for Wilson => Var+<S_W>^2)")
out("      d07 (2.1) 'Sigma = Var' holds iff <Omega|S_W|G> = <S_W>, i.e. ONLY on the diagonal Omega=G.")
out("  [2] int sqrt(Var) <= sqrt(T int Var)  [T]  -- d07's '>=' is a sign error (h06 right).")
out("      BUT the premise 'int A >= c sqrt(V)' is nevertheless TRUE (pointwise: A >= 6V(1-u_p)) --")
out("      so withdrawing 'proposition 3.1' and 'Delta_max >= 0.276 L^2' is an OVER-retraction:")
out("      the implication is valid, 0.276 L^2 is true (weak); the correct constant is O(L^4).")
out("  [3] true numerator O(V) => delta_{S_W} >= c L^4 (c = %.2f L1-route, %.2f L2-route)." % (rows[0][4]/4**4, rows[0][6]/4**4))
