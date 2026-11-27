"""
GP2-ESTIMATOR — difference estimator + resummation audit for delta_1.
Zero MC, numpy only, no npz I/O (archived numbers are quoted from
GM4-NUMERIC.md:185-188). Units: a1 = 1 throughout.

    F = 2 sinh(m2/2)/lam - 2 sinh(m1/2)
    G = 2 sinh(m1/2) - F/(lam^2-1)
    R := delta1/(D^3/24) = 24 F /(G^3 (lam^2-1)),   m_i = a_i Delta_i
"""
import numpy as np
import time

T0 = time.time()
out = []


def p(*a):
    s = " ".join(str(x) for x in a)
    out.append(s)
    print(s)


def R2(m1, m2, lam):
    w1 = 2.0 * np.sinh(m1 / 2.0)
    w2 = 2.0 * np.sinh(m2 / 2.0) / lam
    t1, t2 = 1.0, lam ** 2
    d1 = (w2 - w1) / (t2 - t1)
    Dh = w1 - d1 * t1
    return d1 / (Dh ** 3 / 24.0), d1, Dh


def R2c(m1, m2, lam):
    F = 2.0 * np.sinh(m2 / 2.0) / lam - 2.0 * np.sinh(m1 / 2.0)
    G = 2.0 * np.sinh(m1 / 2.0) - F / (lam ** 2 - 1.0)
    return 24.0 * F / (G ** 3 * (lam ** 2 - 1.0))


def m_free(m1, lam):
    a1D = 2.0 * np.sinh(m1 / 2.0)
    return 2.0 * np.arcsinh(lam * a1D / 2.0)


def RN(mh, lams):
    """exact polynomial in t=a^2 through all points (a1=1)."""
    ts = np.array([1.0] + [l ** 2 for l in lams])
    ws = np.array([2.0 * np.sinh(mh[0] / 2.0)] +
                  [2.0 * np.sinh(m / 2.0) / l for m, l in zip(mh[1:], lams)])
    c = np.linalg.solve(np.vander(ts, len(ts), increasing=True), ws)
    return c, 24.0 * c[1] / c[0] ** 3


p("=== [0] (SI-1) identity omega<->Delta is a bijection ===")
for a in (0.3, 0.03):
    for om in (0.5, 1.0):
        d = 2.0 / a * np.arcsinh(a * om / 2.0)
        p(f"    a={a:<6} omega={om:<5} roundtrip_err={2/a*np.sinh(a*d/2)-om:+.3e}")
p("    => reparametrisation adds no information ([T]); section 2.")
p("    R2 == R2_closed check:")
for (m1, m2, lam) in [(0.55, 0.16, 0.29), (1.9, 0.55, 0.2895), (1.2, 1.5, 1.25)]:
    p(f"      {m1},{m2},{lam:.4f}: {R2(m1,m2,lam)[0]-R2c(m1,m2,lam):+.2e}")

p("")
p("=== [1] ACID A: free field (delta_1=0) -> R must be 0 ===")
for (m1, lam) in [(0.55, 1 / 3.45), (1.0, 0.5), (1.9, 0.6667), (0.05, 0.25),
                  (3.0, 0.3)]:
    m2 = m_free(m1, lam)
    r, d1, Dh = R2(m1, m2, lam)
    p(f"    m1={m1:<5} lam={lam:<8.4f} m2={m2:.6f}  R={r:+.3e}  "
      f"delta1={d1:+.3e}  D={Dh:.12f}")

p("")
p("=== [2] ACID B: 2D (Delta_a == D const) -> R must be 1 ===")
for (m1, lam) in [(0.1, .6667), (0.3, .6667), (0.55, .6667), (0.55, .29),
                  (1.0, .6667), (2.0, .5)]:
    m2 = lam * m1
    r, d1, Dh = R2(m1, m2, lam)
    bias = (m1 ** 2 + m2 ** 2) / 80.0
    p(f"    m1={m1:<5} lam={lam:<7.4f} R={r:.12f}  1+(m1^2+m2^2)/80="
      f"{1+bias:.12f}  resid={r-(1+bias):+.2e}")
p("    3-point (Richardson) on 2D and on free field, same window:")
for lam2 in (0.6667, 0.5):
    mh = [0.55, 0.55 * lam2, 0.55 * lam2 ** 2]
    c, r = RN(mh, [lam2, lam2 ** 2])
    mh2 = [0.55, m_free(0.55, lam2), m_free(m_free(0.55, lam2), lam2)]
    c2, r2 = RN(mh2, [lam2, lam2 ** 2])
    p(f"      lam={lam2}: 2D R3={r:.12f} (exact 1) | free R3={r2:+.3e}")

p("")
p("=== [3] ROUND TRIP: 2-pt bias is EXACTLY dR = 24*d1/D^3 + 24*kappa*(m1^2+m2^2) ===")
D = 1.0
for (m1, lam) in [(1.9, 0.29), (1.0, 0.5), (0.55, 0.29)]:
    for (d1, kap) in [(0.0, 0.0), (1.0, 0.0), (0.0, 1e-2), (1.0, 1e-2)]:
        a1 = 2.0 * np.sinh(m1 / 2.0) / D
        a2 = lam * a1
        om1 = D + d1 * a1 ** 2 + kap * D ** 5 * a1 ** 4
        om2 = D + d1 * a2 ** 2 + kap * D ** 5 * a2 ** 4
        M1 = 2.0 * np.arcsinh(a1 * om1 / 2.0)
        M2 = 2.0 * np.arcsinh(a2 * om2 / 2.0)
        r, _, _ = R2(M1, M2, lam)
        Du = a1 * D
        ku = kap * Du ** 5
        d1u = d1 * a1 ** 3
        pred = 24 * (d1u + ku * (1 + lam ** 2)) / (Du - ku * lam ** 2) ** 3
        p(f"    m1={m1:<5} lam={lam:<6} d1={d1:<4} kappa={kap:<6} "
          f"R={r:+.9f} pred={pred:+.9f} resid={r-pred:+.1e} "
          f"LO=24k(m1^2+m2^2)={24*kap*(M1**2+M2**2):+.6f}")

p("")
p("=== [4] SYSTEMATIC COMMON MODES: E0 (anchored 1-pt) vs E1 (2-pt) ===")


def shift_eps(m1, lam, eps, d1=1 / 24.0):
    """multiplicative common mode omega -> (1+eps) omega, both points.
    Base theory is 2D-like (d1 = D^3/24, R=1) so that the shift is visible."""
    a1 = 2.0 * np.sinh(m1 / 2.0)
    a2 = lam * a1
    om1 = (1 + eps) * (D + d1 * a1 ** 2)
    om2 = (1 + eps) * (D + d1 * a2 ** 2)
    M1 = 2.0 * np.arcsinh(a1 * om1 / 2.0)
    M2 = 2.0 * np.arcsinh(a2 * om2 / 2.0)
    return R2(M1, M2, lam)[0]


for (m1, lam) in [(0.5, 0.5), (1.9, 0.6667)]:
    r0 = shift_eps(m1, lam, 0.0)
    r1 = shift_eps(m1, lam, 1e-2)
    p(f"    m1={m1:<4} eps=1%: R(0)={r0:+.9f} dR_E1={r1-r0:+.6f} "
      f"(pred -2eps*R={-2e-2*r0:+.6f}) ; dR_E0=24eps/m1^2={24e-2/m1**2:+.6f} "
      f"gain={abs((24e-2/m1**2)/(r1-r0)):.1f}x")

p("")
p("=== [5] SENSITIVITY: rho_req for 4 sigma on the 0-vs-1 discrimination ===")


def grad(m1, m2, lam, h=1e-7):
    f = lambda x, y: R2c(x, y, lam)
    return ((f(m1 + h, m2) - f(m1 - h, m2)) / (2 * h),
            (f(m1, m2 + h) - f(m1, m2 - h)) / (2 * h))


def rho_req(m1, m2, lam, target=0.25):
    g1, g2 = grad(m1, m2, lam)
    return target / np.hypot(g1 * m1, g2 * m2)


p("    m1     m2     lam      E1 rho_req   E0 rho_req   rho-gain  N-gain")
for (m1, m2) in [(1.9, 1.3), (1.9, 1.0), (1.9, 0.7), (1.9, 0.55), (1.0, 0.55),
                 (0.55, 0.29), (0.55, 0.18), (0.3, 0.1)]:
    lam = m2 / m1
    r1 = rho_req(m1, m2, lam)
    r0 = min(m1, m2) ** 2 / 96.0
    p(f"    {m1:<6} {m2:<6} {lam:<8.4f} {r1:<12.3e} {r0:<12.3e} "
      f"{r1/r0:<9.2f} {(r1/r0)**2:<7.1f}")

p("")
p("=== [6] TWO-POINT LEGITIMACY: need |kappa| <= 1/(192(m1^2+m2^2)) ===")
p("    kappa = delta2/D^5 (dimensionless); power counting kappa ~ mhat^-5")
for (m1, m2) in [(1.9, 1.3), (1.9, 1.0), (1.9, 0.55), (1.0, 0.55), (0.55, 0.18)]:
    bd = 1.0 / (192.0 * (m1 ** 2 + m2 ** 2))
    pc = min(m1, m2) ** -5
    p(f"    m1={m1:<5} m2={m2:<5} bound={bd:.3e}  pc~{pc:.3e}  "
      f"OVER by {pc/bd:8.1f}x")
p("    2D exception kappa=1/1920: needs m1^2+m2^2 <= 10")
for (m1, m2) in [(1.9, 0.55), (1.0, 0.55), (2.5, 1.5), (3.0, 2.0)]:
    ok = (1 / 1920.0) * (m1 ** 2 + m2 ** 2) <= 1.0 / 192.0
    p(f"      m1={m1:<5} m2={m2:<5} 2D |kappa|(m1^2+m2^2)="
      f"{(m1**2+m2**2)/1920:.4f} <= 1/192={1/192:.4f}? {'PASS' if ok else 'FAIL'}")

p("")
p("=== [7] PADE [2/2] in t: Taylor coefficient is exact by construction ===")


def pade22(ts, ws):
    A = np.zeros((5, 5))
    b = np.zeros(5)
    for i, (t, w) in enumerate(zip(ts, ws)):
        A[i] = [1.0, t, t ** 2, -w * t, -w * t ** 2]
        b[i] = w
    x = np.linalg.solve(A, b)
    p0, p1, p2, q1, q2 = x
    a1 = p1 - p0 * q1
    poles = np.roots([q2, q1, 1.0]) if abs(q2) > 1e-14 else np.array([])
    return a1, poles, np.linalg.cond(A)


d1true = 0.10
p("    model omega = D + a1 t + c t^2/(1+t/tc)  (cutoff pole at tc)")
for tc_rel in (100.0, 10.0, 3.0, 1.0):
    ts = np.linspace(0.05, 1.0, 5)
    ws = 1.0 + d1true * ts + 0.5 * ts ** 2 / (1 + ts / tc_rel)
    a1, poles, cond = pade22(ts, ws)
    rng = np.random.default_rng(0)
    acc = [pade22(ts, ws + 1e-3 * np.std(ws) * rng.standard_normal(5))[0]
           for _ in range(400)]
    amp = np.std(acc) / 1e-3
    pt = np.abs(poles.real).min() if len(poles) else np.inf
    p(f"    tc/tmax={tc_rel:<6} a1={a1:.9f} (true {d1true}) pole|t|={pt:8.3f} "
      f"cond={cond:.2e} noise-amp={amp:.2e} in-window={'YES' if pt<=1 else 'no'}")
p("    => a1 is exact at every tc (0 dof). The pole is an OUTPUT, not a")
p("       control. Noise amplification ~1e2 is the whole cost.")

p("")
p("=== [8] ARCHIVED GW7 TRIPLE (GM4-NUMERIC.md:185-188) ===")
trip = [(2.50, 8, 12, 3.576e-3, 1.9247, 0.1522),
        (2.59, 10, 15, 2.849e-3, 2.2537, 0.4092),
        (2.66, 12, 18, 2.386e-3, 2.4078, 0.3817)]
p("    pair lam_vol R          sigma_R   sigma from 0   sigma from 1")
for i in range(3):
    for j in range(i + 1, 3):
        lam = trip[i][1] / trip[j][1]
        m1, s1 = trip[i][4], trip[i][5]
        m2, s2 = trip[j][4], trip[j][5]
        r = R2c(m1, m2, lam)
        g1, g2 = grad(m1, m2, lam)
        sr = np.hypot(g1 * s1, g2 * s2)
        p(f"    {i+1}-{j+1} {lam:.5f} {r:+.6f}  {sr:.6f}   {abs(r)/sr:6.2f}  "
          f"{abs(r-1)/sr:6.2f}")
p("    cutoff null R_null(m0) at lam=8/12:")
for m0 in (1.6, 1.9, 2.195, 2.4, 2.6):
    p(f"      m0={m0:<6} R_null={R2c(m0, m0, 8 / 12):+.6f}")
p("    gate 1 (ratio): m2/m1 must equal lam on a constant-physics line")
for i in range(3):
    for j in range(i + 1, 3):
        p(f"      {i+1}-{j+1}: m2/m1={trip[j][4]/trip[i][4]:.6f} "
          f"lam={trip[i][1]/trip[j][1]:.6f} ratio="
          f"{trip[j][4]/trip[i][4]/(trip[i][1]/trip[j][1]):.4f}")
p("    gate 2 (log slope): p = d ln mhat/d ln a, a=1/L")
for i in range(3):
    for j in range(i + 1, 3):
        m1, s1 = trip[i][4], trip[i][5]
        m2, s2 = trip[j][4], trip[j][5]
        p_ = np.log(m2 / m1) / np.log(trip[i][1] / trip[j][1])
        sp = np.hypot(s1 / m1, s2 / m2) / abs(np.log(trip[i][1] / trip[j][1]))
        p(f"      {i+1}-{j+1}: p={p_:+.4f} +- {sp:.4f}  "
          f"({abs(p_)/sp:.2f}s from 0, {abs(p_-1)/sp:.2f}s from 1)")
p("    what 1800 s on the L_s=12 point buys (N=65000, rho 18%->7.7%):")
for s3 in (0.18, 0.18 * np.sqrt(12000 / 65000)):
    m1, s1 = trip[0][4], trip[0][5]
    m3 = trip[2][4]
    p_ = np.log(m3 / m1) / np.log(8 / 12)
    sp = np.hypot(s1 / m1, s3 / m3) / abs(np.log(8 / 12))
    p(f"      sigma3/m3={s3:.4f}: p={p_:+.4f} +- {sp:.4f}  "
      f"({abs(p_-1)/sp:.2f}s from 1=mass)")

p("")
p("=== [9] COST MODEL (throughput anchored L_s=8 <-> mhat=1.9) ===")
C8, RHO_REF = 4.8e-3, 0.079


p("    m1     m2     lam     rho_req    expensive point             TOTAL")
for (m1, m2) in [(1.9, 1.3), (1.9, 1.0), (1.9, 0.7), (1.9, 0.55), (1.9, 0.4)]:
    lam = m2 / m1
    rq = rho_req(m1, m2, lam)
    tot = 0.0
    det = []
    for m in (m1, m2):
        c = C8 * (1.9 / m) ** 4
        N = 12000 * (RHO_REF / rq) ** 2
        tot += N * c
        det.append(f"m={m:.2f}: N={N:.2e} {c:.2e}s/cfg {N*c:.1e}s")
    p(f"    {m1:<6} {m2:<6} {lam:<7.4f} {rq:<10.2e} " + " | ".join(det) +
      f"  TOT={tot:.1e}s={tot/3600:.1f}h")
p("    30 min = 1800 s. At L_s=8 that is 3.75e5 cfg -> rho=1.41% (one point).")

p("")
p(f"=== WALL CLOCK: {time.time()-T0:.2f} s ===")
with open(__file__.replace(".py", ".log"), "w", encoding="utf-8") as f:
    f.write("\n".join(out) + "\n")
