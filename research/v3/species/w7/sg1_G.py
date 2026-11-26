# -*- coding: utf-8 -*-
"""
SG1 :  s1 := inf supp nu_w  ;  G := s1 - F  ;  four-model machine check of
      (I) G >= 0 everywhere   (II) G(M6) != 0   (III) signature   (IV) G1-covariance of F
Writes only: research/v3/species/w7/sg1_G.py | sg1_G.log     (no git commit)
py -3.12 + mpmath
"""
import mpmath as mp
mp.mp.dps = 30

L = []
def P(s=""):
    L.append(str(s)); print(s)

def hdr(t):
    P("=" * 92); P(t); P("=" * 92)

KEY_GPP, KEY_DD, KEY_F1 = "Gpp", "Delta_dyn", "f1sq"

# ---------------------------------------------------------------- 0. s1 per model
hdr("(0) s1 per model   [s1 := inf supp nu_w ; covariance s1(tau_a nu) = a + s1(nu), B05-DUAL:80]")

def E_2dym(j, g2=mp.mpf(1)):
    return g2 / 2 * j * (j + 1)

s1_2dym = E_2dym(mp.mpf(1) / 2) - E_2dym(0)
om = mp.mpf(1)
s1_m6 = om
s1_free = mp.mpf(0)

P(f"  2D YM      : E_j=(g^2/2)j(j+1), g^2=1 :  E(0)={E_2dym(0)}  E(1/2)={E_2dym(mp.mpf(1)/2)}  E(1)={E_2dym(1)}")
P(f"               s1 = 3g^2/8 = {mp.nstr(s1_2dym,12)}   (exact rational 3/8)")
P(f"  M6         : s1 = omega = {s1_m6} ;  Gpp = Delta_dyn = 2*omega = {2*om} ; f1sq = 1/(2*omega) = {mp.mpf(1)/(2*om)}")
P(f"  free m=0   : s1 = 0 (k=0 mode) ;  Gpp = 0 , Delta_dyn = 0 , f1sq = 0")
P(f"  single plaq: no transfer matrix => nu_w undefined => s1 UNDEFINED  (ATOM2:45, L2C-SD:24)")

# ---------------------------------------------------------------- 1. covariance test (i)
hdr("(1) COVARIANCE TEST (i):  F(tau_a nu) - F(nu)   vs   a      [nu = 0.5*delta_0 + 0.5*delta_1]")

nu = [(mp.mpf(0), mp.mpf('0.5')), (mp.mpf(1), mp.mpf('0.5'))]

def shift(mu, a):
    return [(e + a, w) for (e, w) in mu]

def s1_of(mu):
    return min(e for (e, w) in mu if w > 0)

def mom(mu, k):
    return mp.fsum(w * e**k for (e, w) in mu)

def mean(mu):
    return mom(mu, 1) / mom(mu, 0)

def var(mu):
    m = mean(mu); return mom(mu, 2) / mom(mu, 0) - m * m

def f1sq(mu):
    e0 = s1_of(mu); tot = mp.fsum(w for (e, w) in mu)
    return mp.fsum(w for (e, w) in mu if e == e0) / tot

cands = [("F = s1 (target)", s1_of),
         ("F = mean", mean),
         ("F = Var", var),
         ("F = f1sq", f1sq),
         ("F = m2", lambda mu: mom(mu, 2)),
         ("F = 0 (vacuous)", lambda mu: mp.mpf(0))]

P(f"  base nu : s1={mp.nstr(s1_of(nu),10)}  mean={mp.nstr(mean(nu),10)}  Var={mp.nstr(var(nu),10)}  f1sq={mp.nstr(f1sq(nu),10)}")
P("")
a_grid = [mp.mpf('0.1'), mp.mpf('0.5'), mp.mpf(1), mp.mpf(2)]
P("  candidate            " + "".join(f"{'a='+mp.nstr(a,3):>14s}" for a in a_grid) + "     verdict")
for name, F in cands:
    row, ok = [], True
    for a in a_grid:
        d = F(shift(nu, a)) - F(nu)
        row.append(mp.nstr(d, 8))
        if abs(d - a) > mp.mpf('1e-20'):
            ok = False
    P(f"  {name:20s} " + "".join(f"{v:>14s}" for v in row) + ("     (i) OK   <-- covariant" if ok else "     (i) FAIL"))
P("")
P("  => COVARIANT set (i) among all candidates = { F = s1 , F = mean }  (mean shifts by a: B05-DUAL:82 'mean is covariant').")
P("     Var / f1sq / m2 / const are G1-INVARIANT => for them G = s1 - F shifts by +a:  G(tau_a nu) = G(nu) + a")
P("     => the spec's second half 'G(tau_a nu) = G(nu)' FAILS for every one of them.  [T]")
P("     And the covariant pair splits by (ii):  mean - s1 = Var-in-money  >= 0 always (equality iff nu is one atom).")
P("     => (i) AND (ii) together are met by NO computable moment functional except F = s1 itself.")

# ---------------------------------------------------------------- 2. lower-side test (ii)
hdr("(2) LOWER-SIDE TEST (ii):  F <= s1 ?    [mean >= s1 always, equality iff nu is a single atom]")
P(f"  mean - s1 = {mp.nstr(mean(nu) - s1_of(nu),12)}  > 0   => (ii) FAILS for the only covariant moment F  (B05:141-143 NG-3)")
P(f"  Var = {mp.nstr(var(nu),12)}  vs  s1 = {mp.nstr(s1_of(nu),6)} :  Var <= s1 is model-dependent => cannot hold for all nu")

# ---------------------------------------------------------------- 3. thin atom no-go
hdr("(3) THIN-ATOM NO-GO (A3-VERDICT:61-69):  finitely-many moments cannot see s1")
s2, lam = mp.mpf(1), mp.mpf(1)
base = [(s2, mp.mpf(1))]
P(f"  base nu = delta at s2 = {s2}   (s1 = {s2},  m_k = {s2}^k)")
P(f"  nu_eps = (eps/(2*lam)) delta_eps + (1-eps/(2*lam)) delta_s2 ,  lam = {lam}")
P("  " + f"{'eps':>10s}{'s1(nu_eps)':>16s}{'w_eps':>14s}{'|dm1|':>13s}{'|dm2|':>13s}{'|dm3|':>13s}")
for eps in [mp.mpf('1e-1'), mp.mpf('1e-2'), mp.mpf('1e-3'), mp.mpf('1e-4'), mp.mpf('1e-6')]:
    w = eps / (2 * lam)
    mue = [(eps, w), (s2, 1 - w)]
    d = [abs(mom(mue, k) - mom(base, k)) for k in (1, 2, 3)]
    P("  " + f"{mp.nstr(eps,3):>10s}{mp.nstr(s1_of(mue),6):>16s}{mp.nstr(w,6):>14s}"
      + "".join(f"{mp.nstr(x,6):>13s}" for x in d))
P("")
P("  => s1 drops 1 -> eps (an O(1) change) while every moment drifts by O(eps).")
P("     F continuous in f.g. moments with F <= s1 for ALL nu  =>  F(delta_s2) <= 0 for every s2  =>  F <= 0 on the class")
P("     => ' 0 < c0 <= F <= s1 ' (spec (6)) is IMPOSSIBLE for any moment-continuous F.  [T]")

# ---- (3b) Theorem SG1-1 : (i) AND (ii) AND moment-continuity are mutually exclusive
hdr("(3b) THEOREM SG1-1:  (i) covariance  AND  (ii) F<=s1  AND  finite-moment computability = EMPTY  [T]")
P("  proof by the atom family. For any a>0 put  nu_eps^a = (eps/2lam) delta_eps + (1-eps/2lam) delta_a .")
P("  eps -> 0 : the moment vector converges to that of delta_a ,  but s1(nu_eps^a) = eps -> 0 .")
P("  (ii) for all eps  =>  F(delta_a) <= 0 .   (i) on atoms =>  F(delta_a) = F(delta_0) + a .")
P("  => F(delta_0) + a <= 0  for ALL a > 0  =>  contradiction.  qed")
P("")
P(f"  {'a':>8s}{'s1(nu_eps^a) @eps=1e-6':>26s}{'|dm1| vs delta_a':>20s}{'(ii) forces F<=':>16s}{'(i) forces F=':>16s}{'consistent?':>14s}")
for a in [mp.mpf('0.5'), mp.mpf(1), mp.mpf(2), mp.mpf(5)]:
    eps = mp.mpf('1e-6'); w = eps / (2 * lam)
    mue = [(eps, w), (a, 1 - w)]
    dm = abs(mom(mue, 1) - mom([(a, mp.mpf(1))], 1))
    P(f"  {mp.nstr(a,3):>8s}{mp.nstr(s1_of(mue),8):>26s}{mp.nstr(dm,6):>20s}{mp.nstr(eps,4):>16s}{mp.nstr(a,3):>16s}{'NO':>14s}")
P("  => no F can satisfy both sides.  The only escape is to give up finite-moment computability.")
P("     Full-moment sequence DOES determine nu (compact support => unique measure) => the unique (i)+(ii) object")
P("     reachable that way is s1 itself => G = 0.  Finite truncation is exactly what the atom family kills.")

# ---------------------------------------------------------------- 4. four-model table
hdr("(4) FOUR-MODEL TABLE:  G = s1 - F  +  checks (I) G>=0 , (II) G(M6)!=0 , (III) signature")

beta = mp.mpf('2.5')
dens = lambda r: mp.sin(r)**2 * mp.exp(-beta * (1 - mp.cos(r)))
Z = mp.quad(dens, [0, mp.pi])
Delta_dyn_plaq = mp.quad(lambda r: (beta * mp.sin(r) / r) * dens(r), [0, mp.pi]) / Z
S2 = mp.quad(lambda r: (beta * mp.cos(r)) * dens(r), [0, mp.pi]) / Z
S1 = mp.quad(lambda r: (beta * mp.sin(r)) * dens(r), [0, mp.pi]) / Z
S12 = mp.quad(lambda r: (beta * mp.sin(r))**2 * dens(r), [0, mp.pi]) / Z
Ddyn_chk = S2 - (S12 - S1**2)
c_n = lambda n: 2 * (n + 1) * mp.besseli(n + 1, beta) / beta
c0, c1 = c_n(0), c_n(1)
P(f"  single plaquette SU(2), beta={beta}  (computed here, own convention fixed explicitly):")
P(f"     conv A  E[(dS/dr)/r] over Haar*exp        = {mp.nstr(Delta_dyn_plaq,10)}    <-- does NOT match corpus (0.589671194)")
P(f"     cross-check E[b cos r] - Var[b sin r]        = {mp.nstr(Ddyn_chk,10)}")
# conv B: N02's 'conditioned at x3=0' reading: density in the 3-vector chart is J(r)*exp(-bS), J=(sin r/r)^2
dB = lambda r: (mp.sin(r)/r)**2 * mp.exp(-beta*(1-mp.cos(r))) * r      # 2D shell weight at x3=0
ZB = mp.quad(dB, [0, mp.pi])
Delta_dyn_cond = mp.quad(lambda r: (beta*mp.sin(r)/r) * dB(r), [0, mp.pi]) / ZB
P(f"     conv B  E_cond[(dS/dr)/r | x3=0] (J=(sin r/r)^2) = {mp.nstr(Delta_dyn_cond,10)}    <-- compare corpus N02: 0.589671194")
P(f"             relative gap = {mp.nstr(abs(Delta_dyn_cond-mp.mpf('0.589671194'))/mp.mpf('0.589671194'),3)}")
P(f"     c0={mp.nstr(c0,8)}  c1={mp.nstr(c1,8)}   'gap' -log(c1/c0) = {mp.nstr(-mp.log(c1/c0),8)}   (corpus ATOM2: -0.0143)")
P(f"     classical curvature d2S/dr2|_0 = beta = {beta} (Haar-normalized Gpp-like)")
P(f"     s1 = UNDEFINED (no transfer matrix)  =>  G = s1 - F is UNDEFINED on this model")

rows = [
    ("free m=0",       s1_free,  {KEY_GPP: mp.mpf(0), KEY_DD: mp.mpf(0),     KEY_F1: mp.mpf(0)}),
    ("2D YM",          s1_2dym,  {KEY_GPP: mp.mpf(0), KEY_DD: mp.mpf(0),     KEY_F1: mp.mpf(0)}),
    ("M6 (per-a)",     s1_m6,    {KEY_GPP: 2*om,      KEY_DD: 2*om,          KEY_F1: mp.mpf(1)/(2*om)}),
    ("M6 (for-all-v)", s1_m6,    {KEY_GPP: 2*om,      KEY_DD: mp.mpf(0),     KEY_F1: mp.mpf(1)/(2*om)}),
]
P("")
P(f"  {'model':16s}{'s1':>9s}{'F=Gpp':>9s}{'F=Ddyn':>9s}{'F=f1sq':>9s} |{'G=Gpp':>9s}{'G=Ddyn':>9s}{'G=f1sq':>9s} | G>=0 all ?")
for name, s1, Fv in rows:
    Gv = {k: s1 - v for k, v in Fv.items()}
    ok = all(g >= 0 for g in Gv.values())
    P(f"  {name:16s}{mp.nstr(s1,6):>9s}{mp.nstr(Fv[KEY_GPP],6):>9s}{mp.nstr(Fv[KEY_DD],6):>9s}{mp.nstr(Fv[KEY_F1],6):>9s} |"
      f"{mp.nstr(Gv[KEY_GPP],6):>9s}{mp.nstr(Gv[KEY_DD],6):>9s}{mp.nstr(Gv[KEY_F1],6):>9s} | {ok}")
P(f"  {'single plaq':16s}{'undef':>9s}{mp.nstr(beta,6):>9s}{mp.nstr(Delta_dyn_plaq,6):>9s}{'0':>9s} |{'undef':>9s}{'undef':>9s}{'undef':>9s} | undef")
P("")
P("  reference rows:")
P(f"    F=0  : G = s1 = (0, {mp.nstr(s1_2dym,6)}, {mp.nstr(s1_m6,6)}, undef)  signature (0, +, !=0, undef)")
P("           -- passes the signature TEXTUALLY, but F carries zero information => G IS the target (NS2-1: circular)")
P("    F=s1 : G = 0  = (0, 0, 0, undef)  signature (0,0,0,undef) -- passes 4/4 => falls under N09-1")
P("")
P("  => every NON-vacuous candidate gives, on M6,  G = 1 - 2 = -1 < 0   => (I) G>=0 FAILS.")
P("  => every moment-based candidate violates (i) (see (1)) and the only covariant one (mean) violates (ii) (see (2)).")

# ---------------------------------------------------------------- 5. M6 mass
hdr("(5) M6:  m6 = min(omega, (2pi/L)^2/2) -> 0   while  G(M6) = omega = 1  (candidate F=0)")
for Lv in [8, 16, 32, 64, 128, 256, 512, 1024]:
    P(f"  L={Lv:<6d} m6 = {mp.nstr(min(om, (2*mp.pi/Lv)**2/2),8)}")
P("  => m6 -> 0 while G(M6) = omega = 1:  G is not a function of m at the level of VALUES.")

# ---------------------------------------------------------------- 6. verdict
hdr("(6) VERDICT")
P("  P1 filter (G>=0 on the whole grid AND G(M6)>0): passes ONLY on the vacuous row F=0, where G = s1 = the target.")
P("  Non-vacuous rows : (I) fails on M6 (G = -1);  (i) fails for every moment F;  (ii) fails for F = mean.")
P("  (3) => 'F <= s1 with F >= c0 > 0' is IMPOSSIBLE for every moment-continuous F  => (2)+(6) of the spec cannot hold.")
P("  (b) TRIGGERED : G obtainable only from the five classes (Gpp / Delta_dyn / Var / f1sq / H1).")
P("  (c) TRIGGERED for the only surviving (vacuous) option: F = s1-type is the Krein-Birman/SSF transport itself (NS2 §3 row 3).")
P("  => (b) DEGENERATE.  No (a).  Clay 0%.")

open(r"D:\Dev\math\research\v3\species\w7\sg1_G.log", "w", encoding="utf-8").write("\n".join(L) + "\n")
print("\n[written] sg1_G.log")
