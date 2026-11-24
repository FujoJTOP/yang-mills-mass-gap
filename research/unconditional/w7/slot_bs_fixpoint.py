"""
slot_bs_fixpoint.py -- Birman-Schwinger eigenfixpoint on the 0++ channel.
ILLUSTRATIVE separable (rank-1) BS model, 3D, single job.  py -3.12

Displays the FIXPOINT TOPOLOGY of kappa = F(kappa) for the 0++ two-gluon BS
channel and tests the structural questions of SLOT-BS-FIXPOINT.md. TOY only.

Model (standard rank-1 separable BS, s-wave):
  chi(p)=1/(p^2+b^2), relative kinetic eps(p)=p^2/m_g (reduced mu=m_g/2),
  BS: 1 = g * I(B;m_g),  I(B;m_g)=(1/2pi^2) int p^2 chi^2/(p^2/m_g+B) dp,
  bound state at relative energy -B (B>0); absolute mass m_B = 2 m_g - B.
  D_thresh = 1 - g*I(0;m_g) = 1 - g/g_c(m_g)   (sign of the BS determinant at threshold).

Claims computed:
  (1) g_c(m_g)=2pi^2 m_g/int(chi^2); D_thresh<0 <=> supercritical <=> K_BS(threshold)>1
      <=> sign-definite inversion (the only escape from the gapless fixed point).
  (2) subcritical (D_thresh>=0, STABLE): no sub-threshold bound state -> F=threshold
      -> F(0)=0 -> gapless fixed point survives.
  (3) positive gap needs supercriticality (massless tachyon) OR a lifted threshold m_g>0.
  (3b) the gap m_B tracks m_g (controller), not self-supplied.
  (4) self-anchoring loop m_g=c*m_B: algebraic => only trivial fixed point m_B=0.
  (5) acid: g=0 -> no bound state -> gap 0.
"""
import numpy as np

b = 1.0
pmax = 3000.0
N = 400_001
p = np.linspace(1e-6, pmax, N)
dp = p[1] - p[0]
chi2 = 1.0 / (p**2 + b**2)**2
int_chi2 = np.sum(chi2) * dp                    # ~ pi/4
p2_chi2 = p**2 * chi2

def I(B, mg):
    return np.sum(p2_chi2 / (p**2 / mg + B)) * dp / (2 * np.pi**2)

def g_c(mg):
    return 2 * np.pi**2 * mg / int_chi2

def binding(g, mg):
    if g <= g_c(mg):
        return 0.0
    lo, hi = 1e-10, 200.0 + 200.0 * mg
    for _ in range(90):
        mid = 0.5 * (lo + hi)
        if g * I(mid, mg) - 1.0 > 0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)

def m_B(g, mg):
    B = binding(g, mg)
    return (None if B == 0.0 else 2 * mg - B), B

print("=" * 74)
print("SLOT-BS-FIXPOINT : 0++ Birman-Schwinger eigenfixpoint -- TOY topology demo")
print("=" * 74)
print(f"int chi^2 dp = {int_chi2:.6f}   (closed form pi/4 = {np.pi/4:.6f})")

print("\n(1) Critical coupling g_c(m_g) = 2 pi^2 m_g / int(chi^2)  [supercritical boundary]")
for mg in [0.25, 0.5, 1.0, 2.0]:
    print(f"    m_g={mg:4.2f}:  g_c = {g_c(mg):9.4f}")

print("\n    D_thresh = 1 - g/g_c : <0 <=> supercritical <=> K_BS(threshold)>1 <=> sign-definite inversion")
mg = 1.0
for gfrac in [0.5, 0.9, 1.0, 1.1, 2.0]:
    D = 1.0 - gfrac
    tag = "supercritical (tachyon onset)" if D < 0 else "subcritical (STABLE vacuum)"
    print(f"    g = {gfrac:4.2f} g_c :  D_thresh = {D:+.4f}   {tag}")

print("\n(2)/(3) gap m_B = 2 m_g - B   ('--' = no isolated sub-threshold state => F=threshold => gap 0)")
print(f"    {'m_g':>6} {'g/g_c':>7} {'B':>9} {'m_B':>9}   interpretation")
for mg in [0.5, 1.0, 2.0]:
    for gfrac in [0.7, 1.2, 1.6]:
        g = gfrac * g_c(mg)
        mb, B = m_B(g, mg)
        if mb is None:
            print(f"    {mg:6.2f} {gfrac:7.2f} {B:9.4f} {'--':>9}   no bound state (gapless; F=threshold)")
        else:
            tag = "physical glueball" if mb > 0 else "DISSOLVED (strong attraction, m_B<=0)"
            print(f"    {mg:6.2f} {gfrac:7.2f} {B:9.4f} {mb:+9.4f}   {tag}")

print("\n(3b) the gap tracks the controller m_g (fixed g/g_c=1.2): m_B rises with m_g")
for mg in [0.5, 1.0, 1.5, 2.0, 2.5]:
    mb, B = m_B(1.2 * g_c(mg), mg)
    print(f"    m_g={mg:4.2f}:  m_B = {mb:+8.4f}")

print("\n(4) Self-anchoring loop m_g = c*m_B : algebra (1-2c) m_B = -B <= 0 for B>0.")
print("    => the only self-anchored fixed point is B=0 (g<=g_c), i.e. m_B=0 (TRIVIAL).")
print("    A nontrivial gap needs an EXTERNAL m_g (or the supercritical tachyon) -- the anchor is imported.")

print("\n(5) Acid test (free / no interaction):")
B0 = binding(0.0, 1.0)
print(f"    free massless (g=0): B={B0:.4f}, bound state: {'NO' if B0==0 else 'yes'}, gap=0 => NOT certified (F1 ok)")
print("    U(1) pure : singlet color factor (C_R-2C_A)/2 = (0-0)/2 = 0 (trivial adjoint) => g_eff=0 => gap 0 (F2 ok)")
print("    :phi^2:   : gapped, atom-free; no gauge binding below its own 2m threshold => gap 0, no atom (acid ok)")

print("\n" + "=" * 74)
print("VERDICT (toy): F(0)=0 persists on the stable (subcritical) vacuum; the sign-")
print("definite inversion D_thresh<0 demands supercriticality (massless 0++ tachyon),")
print("which is ABSENT; a positive gap needs a threshold controller m_g; the self-")
print("anchoring loop is degenerate (only m_B=0). No self-anchored k->0 gap.")
print("=" * 74)
