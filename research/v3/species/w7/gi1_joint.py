"""GI1-JOINT witness — is there a CROSS-SECTOR sum rule that puts charged and neutral
states in one budget?  Minimal, deterministic, numpy only, single process, < 1 s.

Sections
  S1  sector purity of the KL measure of a local operator; cross terms vanish identically
  S2  "neutral = total - charged" fixes the neutral WEIGHT, never the EDGE
  S3  joint variance is exactly additive (Cov = 0); a joint upper bound needs a charged LOWER bound
  S4  one linear datum (the sum rule) cannot pin the edge  (discrete + the gapless anomaly witness)
  S5  2D: nu_{trF2} is the empty measure while the neutral gap is exact and positive
  S6  additive budget cannot see a min over sectors
"""
import numpy as np

np.set_printoptions(precision=12, suppress=False)
out = []


def say(s=""):
    out.append(s)
    print(s)


# ---------------------------------------------------------------- S1
say("== S1  sector purity / no cross terms ==============================")
# 4-dim toy: |0> vacuum (q=0), |1> neutral excited (q=0), |2>,|3> charged pair (q=+1,-1)
E = np.array([0.0, 1.7, 2.3, 2.3])
q = np.array([0, 0, 1, -1])
a, b = 1.0, 0.5
O0 = np.zeros((4, 4)); O0[0, 1] = O0[1, 0] = a          # N-ality 0
Oq = np.zeros((4, 4)); Oq[0, 2] = Oq[2, 0] = b          # N-ality 1 (intertwiner H_0 -> H_1)

def kl_weights(O):
    """|<n|O|0>|^2 for n = 1,2,3 (vacuum n=0 subtracted); returns s = E_n^2."""
    v = O[:, 0]
    return E[1:] ** 2, np.abs(v[1:]) ** 2

s0, w0 = kl_weights(O0)
sq, wq = kl_weights(Oq)
st, wt = kl_weights(O0 + Oq)
say(f"  nu_O0 support (s=m^2, weight): " + str([(float(x), float(y)) for x, y in zip(s0, w0)]))
say(f"  nu_Oq support (s=m^2, weight): " + str([(float(x), float(y)) for x, y in zip(sq, wq)]))
say(f"  nu_(O0+Oq)                   : " + str([(float(x), float(y)) for x, y in zip(st, wt)]))
say(f"  additivity  max|w_tot - (w0+wq)| = {np.max(np.abs(wt - (w0 + wq))):.3e}   [exact 0]")
# cross terms <n|O0|0><0|Oq|n> vanish term by term because the supports are disjoint in n
cross = np.abs(O0[:, 0] * Oq[:, 0])
say(f"  max_n |<n|O0|0>*<n|Oq|0>| = {np.max(cross):.3e}   (charged states get weight 0 from O0)")
say(f"  <0|O0|2> = {O0[2,0]:.1f} , <0|Oq|1> = {Oq[1,0]:.1f}   (grading: mixed matrix elements = 0)")
say(f"  full-H spectrum = {E[1:]}  includes the charged 2.3; nu_O0 sees only {s0[0]:.2f}")
say()

# ---------------------------------------------------------------- S2
say("== S2  total - charged fixes WEIGHT, not EDGE =====================")
for mn in (1.0, 0.3):
    En = np.array([0.0, mn, 2.0, 2.0])
    s0, w0 = En[1:] ** 2, np.array([1.0, 0.0, 0.0])
    sq, wq = En[1:] ** 2, np.array([0.0, 0.25, 0.0])
    tot, ch, ne = w0.sum() + wq.sum(), wq.sum(), w0.sum()
    edge_n = np.sqrt(s0[w0 > 0][0])
    say(f"  m_n={mn:.2f}: total weight={tot:.2f}  charged={ch:.2f}  neutral={tot-ch:.2f}"
        f"  inf supp nu_0={edge_n:.2f}")
say("  -> identical (total, charged, neutral) budget; different neutral edge.")
say()

# ---------------------------------------------------------------- S3
say("== S3  joint variance is exactly additive (Cov = 0) ===============")
p = np.array([0.4, 0.2, 0.2, 0.2])                     # any probability; p2=p3 => <xiq>=0 (grading)
xi0 = np.array([0.0, 1.0, 0.0, 0.0])                   # supported in the neutral sector
xiq = np.array([0.0, 0.0, 1.0, -1.0])                  # supported in the charged sectors
m0 = p @ xi0; mq = p @ xiq
cov = p @ (xi0 * xiq) - m0 * mq
v0 = p @ (xi0 - m0) ** 2; vq = p @ (xiq - mq) ** 2; vj = p @ ((xi0 + xiq) - (m0 + mq)) ** 2
say(f"  Cov(xi0,xiq) = {cov:.3e}   [exact 0: disjoint supports]")
say(f"  Var(xi0)={v0:.6f}  Var(xiq)={vq:.6f}  Var(xi0+xiq)={vj:.6f}"
    f"  diff={abs(vj-(v0+vq)):.3e}")
say("  joint upper bound V on Var(xi0+xiq) alone gives Var(xi0) <= V (the marginal bound, no gain).")
for v_lo in (0.0, 0.4):
    V = 1.0
    say(f"  if Var(xiq) >= {v_lo:.1f}:  Var(xi0) <= V - Var(xiq) <= {V - v_lo:.1f}"
        + ("   (needs Var(xiq)=0: no info)" if v_lo == 0 else "   (needs a CHARGED lower bound)"))
say()

# ---------------------------------------------------------------- S4
say("== S4  one linear datum cannot pin the edge =======================")
Lam = 1.0
# (a) discrete: M_{-2} = 1 for both; edge in s is 1 vs eps=0.01 (mass 1 vs 0.1)
eps, w = 0.01, 1e-3
S2 = (1 - w) / (Lam - w / eps)            # keep M_{-2} = w/eps + (1-w)/S2 = 1
M2_B = w / eps + (1 - w) / S2
say(f"  nu_A = delta(s=1): M_-2 = 1.0, edge (in s) = 1.0, mass edge = 1.0")
say(f"  nu_B = {w} delta(s={eps}) + {1-w:.4f} delta(s={S2:.6f}): "
    f"M_-2 = {M2_B:.12f}, edge (in s) = {eps}, mass edge = {eps**0.5:.2f}")
say(f"  difference in M_-2 = {abs(M2_B-1.0):.3e}   -> same sum-rule datum, different edge")
# (b) continuous anomaly witness: dnu = s^2 e^{-s/Lam} ds
S = 200.0
num = np.trapezoid(np.exp(-np.linspace(0, S, 200001) / Lam), np.linspace(0, S, 200001))
say(f"  gapless witness dnu = s^2 exp(-s/Lam) ds:  M_-2 = int_0^inf exp(-s/Lam) ds = {Lam}")
say(f"    truncated integral (S=200Lam) = {num:.10f}   (converges to Lam),  inf supp = 0, atom = 0")
say("  -> a POSITIVE total budget coexists with gap = 0.")
say()

# ---------------------------------------------------------------- S5
say("== S5  2D: nu_{trF2} empty, neutral gap exact and positive =======")
C2 = {0.0: 0.0, 0.5: 0.75, 1.0: 2.0}          # SU(2): C2(R) = R(R+1), R = 0, 1/2, 1
weights = [C2[R] * (1.0 if R == 0 else 0.0) for R in (0.0, 0.5, 1.0)]
say(f"  eigenvalues of trF^2 on characters R=0,1/2,1: {[C2[R] for R in (0.0,0.5,1.0)]}")
say(f"  <0|trF^2|R> = C2(R) delta_{{0R}}  ->  nu weights {weights}  (empty measure)")
say(f"  int rho w = 0 for every w;  neutral gap = (2N^2/(N^2-1)) * m_ch = 2.6667 m_ch > 0 (exact)")
say()

# ---------------------------------------------------------------- S6
say("== S6  additive budget cannot see a min over sectors ==============")
for name, e0, e1 in (("A", 1.0, 0.5), ("B", 0.4, 1.7)):
    tot = 1.0 + 1.0
    ch = 1.0
    say(f"  family {name}: sector-0 edge={e0}, sector-1 edge={e1}; total weight={tot}, "
        f"charged weight={ch}, neutral weight={tot-ch}; neutral edge={e0}, full gap={min(e0,e1)}")
say("  -> same additive data (total 2, charged 1), different neutral edge; the target is a MIN.")
say()
say("ALL CHECKS DONE")
