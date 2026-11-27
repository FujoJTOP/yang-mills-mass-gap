# HZQ-EPSFLOOR 确定性算术核。无采样、无随机、无 MC、无 numpy。
# 只核本文件承重的 8 组算术；全部有理数，浮点仅用于打印。
from fractions import Fraction as F

out = []
def sec(t): out.append(f"[{t}]")

sec(1)  # w0 = r0 / (N^2/4) = 1/(4N^4)
for N in range(2, 7):
    r0 = F(1, 16 * N * N); W = F(N * N, 4); w0 = r0 / W
    out.append(f"  N={N}: r0={r0} W=N^2/4={W} w0={w0}=1/{w0.denominator} "
               f"1/(4N^4)={F(1, 4 * N ** 4)} eq={w0 == F(1, 4 * N ** 4)}")

sec(2)  # per-L reading: 2x
for N in (2, 3):
    perL = F(1, 8 * N * N); W = F(N * N, 4)
    out.append(f"  N={N}: perL={perL} /W -> {perL / W} = 1/(2N^4)={F(1, 2 * N ** 4)} "
               f"ratio={(perL / W) / F(1, 4 * N ** 4)}")

sec(3)  # margins vs w0=1/64
for x in (0.269, 0.246, 0.3427):
    out.append(f"  eps*={x}: margin vs 1/64 = {x * 64:.4f}x")

sec(4)  # Popoviciu sharpness on [0,N]
for N in (2, 3):
    mu = [(F(1, 2), F(0)), (F(1, 2), F(N))]
    m = sum(w * s for w, s in mu); v = sum(w * (s - m) ** 2 for w, s in mu)
    out.append(f"  N={N}: two-atom {{0,N}} equal -> Var={v} vs N^2/4={F(N * N, 4)} eq={v == F(N * N, 4)}")
mu = [(F(1, 3), F(0)), (F(1, 3), F(1)), (F(1, 3), F(2))]
m = sum(w * s for w, s in mu); v = sum(w * (s - m) ** 2 for w, s in mu)
out.append(f"  3-atom {{0,1,2}} equal -> Var={v} < 1 = W^2/4 (N=2) strict={v < 1}")

sec(5)  # route (a), literal readings
l1, l2 = F(1), F(9, 10)
out.append(f"  decreasing l1=1,l2=9/10: (l2-l1)/l2={ (l2 - l1) / l2 } <=0 -> vacuous")
out.append(f"  increasing l1=0,l2=Delta: (l2-l1)/l2=1 -> 'eps>=1' false (measured 0.269)")

sec(6)  # route (a), corrected form 1-l2/l1=1-e^{-m}: thin atom w=s1/(2lam)=m^2/(2lam)
for m in (F(1, 10), F(1, 2), F(9, 10)):
    for lam in (F(1), F(2)):
        w = m * m / (2 * lam); lb = m - m * m / 2   # 1-e^{-m} > m-m^2/2 > m/2 for 0<m<1
        out.append(f"  m={m} lam={lam}: w=s1/(2lam)={w} vs 1-e^-m > {lb}: violate={w < lb}")

sec(7)  # GM6-3 divided by M0 : nu = 0.1 d_1 + 0.9 Unif[4,9]
eps = F(1, 10); M0 = F(1); M1 = eps * 1 + (1 - eps) * F(13, 2); Lam = F(4); s1 = F(1)
b = (Lam * M0 - M1) / (Lam - s1)
out.append(f"  dilution eps*=1/10: M1={M1} (4M0-M1)/3={b} negative={b < 0} ; f1^2/M0={eps}")
out.append(f"  eps-coordinate: (Lam-M1/M0)/(Lam-s1)={(Lam - M1 / M0) / (Lam - s1)} same={((Lam - M1 / M0) / (Lam - s1)) == b}")

sec(8)  # beta=0 boundary: rank-one transfer => f1^2=0, M0>0
out.append("  beta=0: T=|1><1| rank 1 => no excited state => f1^2=0 ; M0=Var(O)>0 => eps*=0 < eps0")

print("\n".join(out))
print("ALL CHECKS DONE")
