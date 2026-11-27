"""RR9-C3 closed-form arithmetic (zero MC, deterministic).
delta_L floor slot: exponent ledger + (U3) threshold + endpoint numbers.
py -3.12 -I rr9_c3_check.py
"""
from math import log, log10, exp, sqrt

def show(tag, val):
    print(f"{tag} = {val}")

print("== [1] spectral range of truncated magnetic term, one time slice ==")
# repo convention (QQ8 (R4)): 6*L^3 plaquettes per slice; per-plaquette 1-(1/2)trU in [0,2]
for L in (4, 6, 8, 12):
    Np = 6 * L**3
    R = 2 * Np
    print(f"L={L}: Np=6L^3={Np}, R_L=2*Np={R}, ln R_L={log(R):.6f}, log10 R_L={log10(R):.6f}")

print()
print("== [2] delta_L <= R_L/D_phys ; D_phys >= 2^c(L), c(L)=E-V+1 ==")
# E in {3L^3 (spatial links only), 4L^3 (spatial+temporal)}; V=L^3
ln2 = log(2.0)
for L in (4, 6, 8, 12):
    R = 12 * L**3
    print(f"L={L}: R_L={R}")
    for name, E in (("E=3L^3", 3 * L**3), ("E=4L^3", 4 * L**3)):
        c = E - L**3 + 1
        log10d = log10(R) - c * log10(2.0)
        print(f"   {name}: c={c}, log10(delta_L) <= {log10d:.4f}, ln(delta_L) <= {log(R)-c*ln2:.4f}")

print()
print("== [3] sigma(nlev) per L^3: c(L)= (E-V) log nlev ==")
for nlev in (2, 3, 4):
    print(f"nlev={nlev}: sigma_spatial=2*ln(nlev)={2*log(nlev):.6f}, sigma_full=3*ln(nlev)={3*log(nlev):.6f}")

print()
print("== [4] what a floor costs: K_L >= R_L/delta0 ==")
for L in (4, 8, 12):
    R = 12 * L**3
    for d0 in (1.0, 1/64):
        print(f"L={L}: delta0={d0:.6g} -> K_L >= {R/d0:.6g}")

print()
print("== [5] (U3) threshold: beta_c <= 4*sqrt(3R sqrt(k)/tau), R=L^6/(delta |E01|) ==")
# best case: delta*|E01| ~ D_phys^{-1}, sqrt(k)=1, tau=1/3
tau = 1/3
for L in (4, 6, 8):
    for name, sigma in (("sigma=ln2", ln2), ("sigma=3ln2", 3*ln2)):
        # delta*|E01| >= e^{-sigma L^3} ; R <= L^6 e^{+sigma L^3}
        lnR = 6*log(L) + sigma * L**3
        lnbeta = 0.5 * (log(12.0/tau) + lnR)   # 4 sqrt(3R/tau) -> ln = ln4 + 0.5 ln(3R/tau)
        lnbeta = log(4.0 * sqrt(3.0)) - 0.5*log(tau) + 0.5*lnR
        print(f"L={L} {name}: ln R <= {lnR:.4f}, ln beta_c <= {lnbeta:.4f}, beta_c <= {exp(lnbeta):.4e}")

print()
print("== [6] commuting (factorized) endpoint: sumset of {0,2} over 6L^3 terms ==")
for L in (4, 8, 12):
    print(f"L={L}: spectrum = {{0,2,...,{12*L**3}}} -> min spacing = 2 (L-independent)")

print()
print("== [7] budget branch numbers (restated, from RR2/QQ9) ==")
print("Jacobi pair: TV(ln e)=0.13561 vs B_pair=ln(0.95769382/w0)=%.6f -> margin %.4f" % (log(0.95769382/(1/64)), log(0.95769382/(1/64))/0.13561))
print("KS 6-level:  TV=43.529345 vs B=ln(21.46752)=%.6f -> violation %.4f" % (log(21.46752), 43.529345/log(21.46752)))
print("counter-term mass per slice: pairs=10*Np=60L^3")
for L in (4, 8, 12):
    print(f"   L={L}: 60L^3 = {60*L**3}")

print()
print("== [8] exact table for the file ==")
for L in (4, 8, 12):
    R = 12 * L**3
    c2, c3 = 2*L**3, 3*L**3
    print(f"L={L}: R_L/R_1={R/12:.0f}; -log10 delta (c=2L^3) = {-(log10(R)-c2*log10(2.0)):.2f}; (c=3L^3) = {-(log10(R)-c3*log10(2.0)):.2f}")
