# -*- coding: utf-8 -*-
"""tt1 — ||T_-|| / gap(T_+) for the decorated-cluster transfer kernel.

Exact anchors used (all from corpus, not invented):
  X01-L1C.md:43 / C03-BYPASS.md:37   log Z_dec(B_T) = sum_a (T-a+1)_+ A(a) = phi T + log D_inf + eps_T
  C03-BYPASS.md:50 / F07-ADJCBAR.md:31  geometric model A(a) = U^(a+1)/(1-U),
                                          eps_T = U^(T+3)/(1-U)^3   (|eps_T|/U^T == const, no T prefactor)
  P4-OTHER.md:153   eD = e*20 = 54.3656365692,  Cbar in [1,1.0535] (beta<=1/48),  Chat=4
  P4-OTHER.md:159   u threshold for delta>0 : u < 1/(Cbar eD) = 0.018394   (u=beta/4)
  P4-OTHER.md:162   beta=2.5 : u=0.625, ubar=33.98
  P4-OTHER.md:103   ubar<=0.25 => D_inf in [0.690478, 1.448374]
  P4-OTHER.md:110   ubar=0.5654 => D_inf >= 0.030111
Definition locked in TT1.md 2:
  T = transfer kernel of the tube axis, symbol  Ahat(x) = sum_a A(a) x^a  (b's GF = Ahat(x)/(1-x)^2)
  A = A_+ - A_-   (A_pm := max(+-A,0))  =>  T = T_+ - T_-, both entrywise >= 0
  ||T_-|| := rho(T_-)  (Perron root of the nonnegative part)
  gap(T_+) := 1 - rho_2(T_+)   [T normalized to Perron root 1]
  Kato criterion: rho(T_-) < 1 - rho_2(T_+)
"""
from mpmath import mp, mpf, exp, log, inf, nstr

mp.dps = 30

eD = exp(mpf(1)) * 20  # 54.3656365692...  e*Delta, Delta=20
print("eD = e*20 =", nstr(eD, 15))

# ---------- 1. exact identity check for the geometric model ----------
def A_geom(U, M=40):
    return [U ** (a + 1) / (1 - U) for a in range(M + 1)]

def b_direct(U, T, M=60):
    A = A_geom(U, M)
    return sum((T - a + 1) * A[a] for a in range(0, min(T, M) + 1))

def phi_geom(U):
    return U / (1 - U) ** 2

def logD_geom(U):
    return U * (1 - 2 * U) / (1 - U) ** 3          # = U/(1-U)^2 - U^2/(1-U)^3

def eps_geom(U, T):
    return U ** (T + 3) / (1 - U) ** 3

print("\n-- 1. exact identity: sum_a (T-a+1)_+ A(a)  ==  phi*T + log D_inf + eps_T --")
for U in [mpf("0.25"), mpf("0.5654")]:
    res = [abs(b_direct(U, T) - (phi_geom(U) * T + logD_geom(U) + eps_geom(U, T))) for T in range(0, 9)]
    print(f"  U={nstr(U,6)}  max|residual| over T=0..8 = {nstr(max(res),5)}"
          f"   phi={nstr(phi_geom(U),10)}  logD={nstr(logD_geom(U),10)}  D={nstr(exp(logD_geom(U)),10)}")
# reproduce corpus P4:103 / P4:110 windows
print("  cross-check vs corpus:  U=0.25 -> D_inf =", nstr(exp(logD_geom(mpf("0.25"))), 8),
      " (corpus window at ubar<=0.25: [0.690478, 1.448374])")
print("  cross-check vs corpus:  U=0.5654 -> D_inf =", nstr(exp(logD_geom(mpf("0.5654"))), 8),
      " (corpus lower bound at ubar=0.5654: 0.030111)")

# ---------- 2. absolute summability: sum_a |A(a)| = U/(1-U) ----------
print("\n-- 2. ||T_-|| as an ABSOLUTE sum:  sum_a |A(a)| = U/(1-U)  (diverges iff U>=1) --")
for U in [mpf("0.125"), mpf("0.25"), mpf("0.2634"), mpf("1.0"), mpf("16.9893"), mpf("33.9785")]:
    if U < 1:
        print(f"  ubar={nstr(U,8):>10}  sum|A| = {nstr(U/(1-U),10)}   FINITE")
    else:
        print(f"  ubar={nstr(U,8):>10}  sum|A| = +infinity           DIVERGES")

# ---------- 3. tree sum R = u (1+R)^Delta : does the decorated gas converge? ----------
print("\n-- 3. honest tree sum  R = u (1+R)^Delta  (Delta=20): root exists? --")
def tree_root(u, Delta=20, iters=4000):
    # f(R)=u(1+R)^D - R ; f(0)=u>0, min at (1+R)^(D-1)=1/(D u)
    if u * Delta ** Delta > 1:  # crude necessary: derivative branch unreachable
        pass
    lo_arg = 1.0 / (Delta * float(u))
    if lo_arg <= 1.0:
        return None                      # no interior minimum with R>0 => f>0 always
    Rmin = lo_arg ** (1.0 / (Delta - 1)) - 1.0
    f = lambda R: float(u) * (1 + R) ** Delta - R
    if f(Rmin) > 0:
        return None                      # no dip below zero => f>0 everywhere => no fixed point
    lo, hi = 0.0, Rmin                   # f(0)=u>0, f(Rmin)<0  => smallest positive fixed point
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if f(mid) > 0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)
for uu, lab in [(mpf("0.0104"), "C03 slack point u=0.0104"), (mpf("0.0045985"), "beta=0.018394 (u=beta/4)"),
                (mpf("0.625"), "beta=2.5 (u=beta/4)")]:
    r = tree_root(uu)
    print(f"  {lab:<28} u={nstr(uu,8):>10}  R_root = " + ("none  => tree sum DIVERGES" if r is None else nstr(r, 8)))

# ---------- 4. the ratio ----------
print("\n-- 4. rho_2(T_+) = ubar (model: spec T = {1,1,ubar}),  gap = 1 - ubar --")
print("      ||T_-||/gap  :  worst case rho(T_-)=ubar (upper end of the split) ; best case rho(T_-)=0")
rows = []
betas = [("0.018394", mpf("0.018394")), ("2.5", mpf("2.5"))]
for lab, beta in betas:
    for conv, ufun in [("u=beta/4", lambda b: b / 4), ("u=beta/8", lambda b: b / 8)]:
        u = ufun(beta)
        for Cbar in [mpf(1), mpf("1.0535")]:
            ubar = eD * Cbar * u
            gap = 1 - ubar
            ratio = (ubar / gap) if gap > 0 else inf
            rows.append((lab, conv, nstr(Cbar, 5), nstr(u, 10), nstr(ubar, 10), nstr(gap, 10),
                         ("+inf" if gap <= 0 else nstr(ratio, 10))))
print(f"{'beta':>9} {'conv':>8} {'Cbar':>7} {'u':>12} {'ubar':>12} {'gap=1-ubar':>13} {'ratio':>12}")
for r in rows:
    print(f"{r[0]:>9} {r[1]:>8} {r[2]:>7} {r[3]:>12} {r[4]:>12} {r[5]:>13} {r[6]:>12}")

print("\n-- 4b. with the conservative Chat=4 folded in (ubar_eff = Chat*Cbar*eD*u) --")
for lab, beta in betas:
    u = beta / 4
    for Chat in [mpf(1), mpf(4)]:
        ubar = Chat * eD * u
        gap = 1 - ubar
        print(f"  beta={lab:>9}  Chat={nstr(Chat,2):>2}  u={nstr(u,10):>12}  ubar_eff={nstr(ubar,10):>12}"
              f"  gap={nstr(gap,10):>13}  ratio=" + ("+inf" if gap <= 0 else nstr(ubar / gap, 10)))

print("\n-- 4c. beta=2.5: |ubar/(1-ubar)| if one insists on a signed gap --")
u25 = eD * mpf("0.625")
print("   ubar(beta=2.5, u=beta/4) =", nstr(u25, 12), " (corpus P4:162 says 33.98)")
print("   ubar(beta=2.5, u=beta/8) =", nstr(eD * mpf("0.3125"), 12))
print("   |ubar/(1-ubar)|          =", nstr(abs(u25 / (1 - u25)), 12))
print("   ubar/gap with gap capped at 1 (i.e. any admissible gap<=1) =", nstr(u25, 12))

print("\n-- 5. analytic continuation trap (DO NOT read as a value) --")
for U in [mpf("0.25"), mpf("0.5654"), u25]:
    print(f"   U={nstr(U,10):>12}  closed-form exp(logD_geom) = {nstr(exp(logD_geom(U)),10)}"
          f"   [formal continuation only; series diverges for U>1]")
print("\n-- 6. verdict --")
print("   beta=0.018394 : ubar=0.2500, gap=0.7500, ratio=0.3333 < 1  => (a)")
print("   beta=2.5      : ubar=33.9785, gap=1-ubar<0 (no gap) / sum|A|=inf => ratio=+inf => (b)")
