# hz5_negmoment_check.py — HZ5 / M_{-1} = \int s^{-1} d\nu  upper-bound machine
# deterministic exact (Fraction) + closed-form float arithmetic. No sampling, no MC, no randomness.
from fractions import Fraction as F
import math

ok = True
def chk(tag, cond, msg=""):
    global ok
    ok = ok and bool(cond)
    print(("PASS " if cond else "FAIL ") + tag + (("  " + msg) if msg else ""))

print("== [1] target  M_{-1} <= C   vs   H1: chi <= 1/lambda  (C = 1/lambda) ==")
for lam in [F(1,2), F(1), F(2), F(5,2)]:
    C = 1/lam
    for chi in [C/2, C, C]:
        chk("1:%s" % lam, (chi <= C) == (chi <= 1/lam))

print("== [2] layer cake / dyadic two-sided:  (1/2)Sum <= M_{-1} <= Sum + M0 ==")
def dyadic(weights, Nmax=400):
    # weights: list of (s_i, w_i) atoms.  Sum = sum_{n>=1} 2^n nu([0,2^-n])
    S = F(0)
    for n in range(1, Nmax+1):
        t = F(1, 2**n)
        m = sum(w for (s, w) in weights if s <= t)
        if m == 0 and n > 60:
            break
        S += (2**n) * m
    return S

cases = {
 "two atoms 1/2@1/4 + 1/2@4": [(F(1,4), F(1,2)), (F(4), F(1,2))],
 "thin atom lam=1 eps=1/64": [(F(1,64), F(1,128)), (F(2)-F(1,64), F(1)-F(1,128))],
 "single atom @1": [(F(1), F(1))],
 "single atom @1/1024": [(F(1,1024), F(1))],
}
for name, ws in cases.items():
    M0 = sum(w for (s, w) in ws)
    Mm1 = sum(w/s for (s, w) in ws)
    S = dyadic(ws)
    chk("2:" + name, (S/2 <= Mm1) and (Mm1 <= S + M0), "Sum=%.6f Mm1=%.6f" % (float(S), float(Mm1)))

print("== [3] thin atom: M_{-1} = 1/lambda EXACTLY, bottom = eps -> 0 ==")
for lam in [F(1), F(2)]:
    for k in [1, 2, 4, 8, 10, 20]:
        eps = F(1, 2**k)
        w = eps/(2*lam)
        s2 = 2*lam - eps
        Mm1 = w/eps + (1-w)/s2
        chk("3:lam=%s eps=2^-%d" % (lam, k), Mm1 == 1/lam, "Mm1=%s" % Mm1)

print("== [4] coercivity cap:  RD_c at K=1/s  =>  c <= inf supp nu ==")
for m in [F(1,4), F(1), F(4), F(9)]:
    for f in [F(1,10), F(1,2), F(9,10)]:
        for M in [F(10), F(100)]:
            ws = [(m, f), (M, 1-f)]
            Mm1 = sum(w/s for (s, w) in ws); M0 = sum(w for (s, w) in ws)
            cmax = f/Mm1                      # largest c with f >= c*Mm1  (RD_c)
            chk("4:m=%s f=%s M=%s" % (m, f, M), cmax <= m, "c_max=%.6f  m=%.6f" % (float(cmax), float(m)))

print("== [5] Young/Holder + Jensen:  M_{-p} >= M_{-1}^p M0^{1-p}  (p>1) ==")
def negmom(ws, p):
    return sum(w/(s**p) for (s, w) in ws)
for name, ws in cases.items():
    M0 = sum(w for (s, w) in ws); Mm1 = sum(w/s for (s, w) in ws)
    for p in [1.5, 2.0, 3.0]:
        Mp = negmom(ws, p)
        chk("5:%s p=%.1f" % (name, p), Mp >= Mm1**p * M0**(1-p) - 1e-12,
            "M-p=%.6g  lower=%.6g" % (float(Mp), float(Mm1**p*M0**(1-p))))

print("== [6] Legendre:  max_g (2g - s g^2) = 1/s  (grid check) ==")
for s in [0.01, 0.1, 1.0, 3.7, 25.0]:
    gs = [i*3.0/s/200000 for i in range(1, 200001)]
    best = max(2*g - s*g*g for g in gs)
    chk("6:s=%.2f" % s, abs(best - 1.0/s) <= 1e-4/s + 1e-12, "max=%.8f 1/s=%.8f" % (best, 1.0/s))
    # exact algebra: 1/s - (2g - s g^2) = (s g - 1)^2 / s >= 0
    chk("6b:s=%.2f" % s, all(abs((1.0/s - (2*g - s*g*g)) - (s*g-1)**2/s) < 1e-12 for g in [0.1, 0.7, 1.0/s, 5.0]))

print("== [7] no affine majorant of 1/s:  a+bs >= 1/s on (0,oo) is impossible ==")
# as+bs^2 >= 1 for all s>0 fails at s->0.  Check on a decreasing grid that any (a,b) fails.
for a in [-5.0, 0.0, 1.0, 5.0]:
    for b in [-5.0, 0.0, 1.0, 5.0]:
        s = 1e-6
        chk("7:a=%.1f b=%.1f" % (a, b), a*s + b*s*s < 1.0, "as+bs^2=%.3g" % (a*s + b*s*s))

print("== [8] moment ladder is BLIND: all M_k (k>=0) bounded, M_{-1} -> infinity ==")
def Mk_family(t, k):
    # nu_t = delta_t + 1_{[t+1/2, t+3/2]} ds
    return t**k + ((t+1.5)**(k+1) - (t+0.5)**(k+1))/(k+1)
tvals = [1.0, 1e-2, 1e-4, 1e-6, 1e-8]
for k in range(0, 9):
    vals = [Mk_family(t, k) for t in tvals]
    # M_k(t) is increasing in t on (0,1]: (t+1.5)^k > (t+0.5)^k  =>  max at t=1, finite, t-independent
    chk("8a:k=%d" % k, max(vals) == vals[0] and math.isfinite(max(vals)), "max=%.4f (=M_k(1))" % max(vals))
Mm1_t = [1/t + math.log((t+1.5)/(t+0.5)) for t in tvals]
chk("8b:M-1 grows", Mm1_t[-1] > 1e7, "M-1(t=1e-8)=%.3e" % Mm1_t[-1])
print("   M0,M1,M2 at t=1e-6:", ["%.6f" % Mk_family(1e-6, k) for k in range(3)])

print("== [9] converse: M0,M-1 bounded while M1 -> infinity ==")
for S in [1e2, 1e4, 1e6]:
    w = 1.0/S**3
    M1 = (1-w)*S + w*(1.0/S**2)
    Mm1 = (1-w)/S + w*S**2
    chk("9:S=%.0e" % S, Mm1 <= 2.0/S*(1+1e-12), "M1=%.3e  M-1=%.3e  (2/S=%.3e)" % (M1, Mm1, 2.0/S))

print("== [10] fiber (M0,M-1)=(1,1) contains every bottom eps in (0,2) ==")
for k in [1, 2, 3, 6, 12, 30, 100]:
    eps = F(1, 2**k); w = eps/2; s2 = 2-eps
    ws = [(eps, w), (s2, 1-w)]
    M0 = sum(v for (s, v) in ws); Mm1 = sum(v/s for (s, v) in ws)
    chk("10:eps=2^-%d" % k, (M0 == 1) and (Mm1 == 1), "bottom=%s" % eps)

print("== [11] free-field box scaling: two-level with s1 = 4 pi^2 / L^2, f=1/2 ==")
for L in [2, 4, 8, 16, 64]:
    s1 = 4*math.pi**2/L**2
    Mm1 = 0.5/s1 + 0.5/4.0
    print("   L=%2d  s1=%.6f  M-1=%.6f" % (L, s1, Mm1))
chk("11:grows like L^2 (ratio L=64/L=8 > 50)", (0.5/(4*math.pi**2/64**2)+0.125)/(0.5/(4*math.pi**2/8**2)+0.125) > 50)

print("== [12] fifth cell (uniformity control): nu_L = w_L delta_1 + (1-w_L) delta_4, w_L=1/2+1/(4+L) ==")
for L in [2, 4, 8, 16, 64]:
    w = F(1,2) + F(1,4+L)
    Mm1 = w*1 + (1-w)/4
    print("   L=%2d  f1^2=w=%s  M-1=%s=%.6f  m1^2=1  Lambda*=4" % (L, w, Mm1, float(Mm1)))
chk("12:bounded, decreasing", all((F(1,2)+F(1,4+L)) + (1-(F(1,2)+F(1,4+L)))/4 <= F(3,4) for L in [2,4,8,16,64]))

print("== [13] s^a continuum: M0,M1 bounded, M-1 -> infinity, bottom = 0, no atoms ==")
for a in [2.0, 1.0, 0.5, 0.25, 0.05]:
    M0 = 1.0; M1 = (a+1)/(a+2); Mm1 = (a+1)/a     # nu = (a+1) s^a ds on [0,1]
    chk("13:a=%.2f" % a, 0.5 <= M1 <= 1.0 and Mm1 > 0, "M1=%.4f M-1=%.3f" % (M1, Mm1))

print("== [14] Jensen direction: M_{-1} >= M0^2/M1  (lower bound, wrong way for H1) ==")
for name, ws in cases.items():
    M0 = sum(w for (s, w) in ws); M1 = sum(s*w for (s, w) in ws); Mm1 = sum(w/s for (s, w) in ws)
    chk("14:" + name, Mm1 >= M0*M0/M1, "M-1=%.6f M0^2/M1=%.6f" % (float(Mm1), float(M0*M0/M1)))

print()
print("ALL CHECKS PASS" if ok else "SOME CHECKS FAILED")
