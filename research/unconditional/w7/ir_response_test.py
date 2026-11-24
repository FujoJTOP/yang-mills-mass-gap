#!/usr/bin/env python3.12
"""
ir_response_test.py  --  decisive, honest test of the candidate "6th species"
    N_IR^(n) := IR-non-analytic part of the connected n-th-order causal nonlinear response.

Claimed dichotomy property:  gapped => analytic at omega=0 => N=0 ;  gapless => non-analytic => N!=0.
We test:
  (1) LINEAR response:  analytic at 0  <=>  gap  (Kallen-Lehmann).  => reformulation (circular). [T]
  (2) NONLINEAR (n=3):  free field is GAUSSIAN => all connected n>=3 responses == 0
                        => N^(3)(free massless) = 0 = the "gapped" value  => MISLABELS a gapless
                        field as gapped  => FAILS the free-field acid test.  [T + V-here]
  (3) ANOMALY: convergent-moment (sum-rule) functionals are satisfied by BOTH gapped and a
               gapless witness dnu = s^2 e^{-s/Lambda} ds  => blind / non-discriminating. [V-here]

Single job, modest memory. No external data.  numpy only.
"""
import numpy as np

np.random.seed(12345)
OUT = {}
def hr(t): print("\n" + "=" * 78 + "\n" + t + "\n" + "=" * 78)

# ---------------------------------------------------------------------------
hr("PART 1 [T] -- Wick theorem: a GAUSSIAN (free) field has ALL connected n>=3 responses = 0")
print("""
For a free scalar phi (any mass m), the measure is Gaussian. The connected generating
functional W = log Z is quadratic. Hence every connected correlator of order n>=3 vanishes
IDENTICALLY (not approximately):
    chi^(n)_connected  ==  0     for all n >= 3,   for ALL m (gapped OR gapless).
Therefore   N_IR^(n>=3)  ==  0   EXACTLY, for the free MASSIVE and the free MASSLESS scalar.
=> a GAPLESS free field gives N^(3) = 0 = the GAPPED value.  The detector is invalid. [T]

Numerical confirmation that the cumulant ESTIMATOR reproduces this, and that the SAME gapless
field shows its gaplessness only in the LINEAR (2-pt) channel:""")

# --- free Gaussian field on a 3D periodic lattice, sampled by FFT diagonalisation --------
def sample_free_field_3d(L, m, Nconf, nbatch=200):
    """Return arrays phiconf[Nconf, L, L, L] sampled from the free Gaussian measure
       cov G(k) = 1/(m^2 + sum_i 2(1-cos k_i)).  Real field via Hermitian-symmetric FFT."""
    k = 2 * np.pi * np.fft.fftfreq(L)                 # FFT order (matches fftn)
    KX, KY, KZ = np.meshgrid(k, k, k, indexing='ij')
    lam = m * m + 2 * (1 - np.cos(KX)) + 2 * (1 - np.cos(KY)) + 2 * (1 - np.cos(KZ))
    lam = np.maximum(lam, 1e-12)
    sqrtG = 1.0 / np.sqrt(lam)
    out = np.empty((Nconf, L, L, L))
    filled = 0
    while filled < Nconf:
        nb = min(nbatch, Nconf - filled)
        w = np.random.standard_normal((nb, L, L, L))    # white noise (real, unit var)
        wf = np.fft.fftn(w, axes=(1, 2, 3))
        phif = wf * sqrtG[None, :, :, :]                # colour the noise (diagonal in k)
        phi = np.fft.ifftn(phif, axes=(1, 2, 3)).real   # Hermitian-symmetric -> real field
        out[filled:filled + nb] = phi
        filled += nb
    return out

print(f"{'m':>6} {'xi=1/m':>8} {'<phi^2> (2pt,: grows => gapless)':>34} "
      f"{'connected 3pt: max|z|':>22} {'connected 4pt cumul (z)':>24}")
res_p2 = []
for m in [1.0, 0.5, 0.2, 0.05]:
    L = 8
    Nc = 20000
    phi = sample_free_field_3d(L, m, Nc)
    # 2-point (variance) -- grows as m -> 0  =>  longer correlation length  =>  gapless regime
    var = float(np.mean(phi * phi))
    # connected 3-point at a few DISTINCT triples -> should be exactly 0 (Gaussian)
    zs3 = []
    for (a, b, c) in [((0,0,0),(1,0,0),(0,1,0)), ((0,0,0),(2,0,0),(0,0,1)),
                      ((0,0,0),(0,0,0),(1,1,1)), ((3,1,0),(0,2,0),(1,0,2))]:
        prod = phi[:, a[0], a[1], a[2]] * phi[:, b[0], b[1], b[2]] * phi[:, c[0], c[1], c[2]]
        mu, se = float(np.mean(prod)), float(np.std(prod) / np.sqrt(Nc))
        # subtract the (vanishing) lower connected pieces; <phi>=0 exactly for the free measure
        zs3.append(abs(mu) / max(se, 1e-30))
    # connected 4-point cumulant at a 4-tuple:  <x1x2x3x4> - sum_{pairings} <xixj><xkxl>
    def G(i, j):
        return float(np.mean(phi[:, i[0], i[1], i[2]] * phi[:, j[0], j[1], j[2]]))
    p = [(0,0,0),(1,0,0),(0,1,0),(1,1,0)]
    k4 = (G(p[0],p[1])*G(p[2],p[3]) + G(p[0],p[2])*G(p[1],p[3]) + G(p[0],p[3])*G(p[1],p[2]))
    prod4 = phi[:,0,0,0]*phi[:,1,0,0]*phi[:,0,1,0]*phi[:,1,1,0]
    mu4, se4 = float(np.mean(prod4)), float(np.std(prod4)/np.sqrt(Nc))
    z4 = abs(mu4 - k4) / max(se4, 1e-30)
    print(f"{m:6.3f} {1.0/m:8.1f} {var:34.5f} {max(zs3):22.2f} {z4:24.2f}")
    res_p2.append((m, var, max(zs3), z4))
OUT["part2_free_gaussian"] = res_p2
print("""
Reading: <phi^2> grows without bound as m->0 (correlation length xi=1/m -> infinity: the field
is GAPLESS in the IR), yet every connected n>=3 cumulant stays at |z| <~ few (consistent with 0).
=> the SAME gapless field that the LINEAR 2-pt sees as gapless gives N^(3)=0. [T confirmed V-here]
   The gaplessness is visible ONLY in the linear (KL) channel -- which is the gap reformulated.""")

# ---------------------------------------------------------------------------
hr("PART 2 [C/V-here] -- frequency-space 3rd-order response: gapped analytic, gapless singular, FREE = 0")
print("""
Leading (tree) connected cubic response for a phi^3 coupling g, internal propagators D(w)=1/(w^2-m^2+i0),
with w1+w2+w3=0.   chi^(3)(w) = g * D(w1) D(w2) D(w3).   We set (w1,w2,w3)=(w,w,-2w):
""")
print(f"{'m (gap)':>8} {'g':>7} {'|chi^(3)(w)| @w=1e-3':>22} {'|chi^(3)| @w=1e-2':>20} "
      f"{'|chi^3|@1e-3 / @1e-2':>20}   verdict")
for (m, g) in [(1.0, 1.0), (1.0, 0.1), (1.0, 0.0), (0.0, 1.0), (0.0, 0.1), (0.0, 0.0)]:
    def D(w, m):
        return 1.0 / (w*w - m*m + 1e-12j)
    def chi(w):
        return g * D(w, m) * D(w, m) * D(-2*w, m)
    a = abs(chi(1e-3)); b = abs(chi(1e-2))
    ratio = a / max(b, 1e-300)
    if g == 0.0:
        verd = "FREE: identically 0 (Gaussian) -- for BOTH m=1 and m=0"
    elif m > 0:
        verd = "gapped: smooth, Taylor-analytic near w=0 (radius ~ m)"
    else:
        verd = "gapless: diverges as w->0 (non-analytic IR)"
    print(f"{m:8.2f} {g:7.2f} {a:22.4e} {b:20.4e} {ratio:20.3e}   {verd}")
# quantitative: IR exponent of |chi^(3)| for the gapless interacting case
ws = np.logspace(-4, -1, 40)
amp = np.array([abs(1.0 * D(w,0.0)*D(w,0.0)*D(-2*w,0.0)) for w in ws])
slope = np.polyfit(np.log(ws), np.log(amp), 1)[0]
print(f"\n  gapless (m=0,g=1): d log|chi^(3)| / d log w = {slope:.3f}  (expected 6: three 1/w^2 factors)")
OUT["part3_tree_chi3"] = {"gapless_IR_slope": float(slope)}
print("""
Decisive: as the coupling g -> 0 the response vanishes IDENTICALLY (free), while the field stays
gapless (m=0).  So N^(3) is proportional to g and -> 0 in the free gapless limit:
    N^(3)(free massless) = 0 = N^(3)(any gapped theory).
The nonlinear functional tracks INTERACTION STRENGTH (non-Gaussianity), NOT the gap.  [C]""")

# ---------------------------------------------------------------------------
hr("PART 3 [T/V-here] -- LINEAR response (KL): 'analytic at 0 <=> gap' via the even-moment criterion")
print("""
Kallen-Lehmann:  G_R(w)=int_0^inf rho(s)/(w^2 - s^2 + i0) ds.  Expand near w=0:
    G_R(w) = - sum_{n>=0} w^{2n}  M_n ,      M_n := int_0^inf rho(s) s^{-2n-2} ds .
Analytic at w=0  <=>  the Taylor radius is >0  <=>  ALL moments M_n finite.
  * gapped (rho = delta(s-m))      : M_n = m^{-2n-2}   finite for all n  => analytic.  [T]
  * gapless witness rho = s^2 e^{-s/Lambda} ds : M_0 = Lambda finite, but M_1 = int e^{-s/Lambda}s^{-2}ds = +inf
                                     => Taylor series does NOT exist => non-analytic at w=0.  [T]
We confirm M_1 divergence numerically with an IR cutoff eps:""")
print(f"{'measure':>34} {'M_0':>14} {'M_1(eps=1e-2)':>16} {'M_1(eps=1e-4)':>16} {'M_1(eps=1e-6)':>16}")
Lam = 1.0; mgap = 0.7
M0_g = 1.0/mgap**2; M1_g = 1.0/mgap**4      # gapped delta(s-mgap)
print(f"{'gapped delta(s-0.7): exact':>34} {M0_g:14.5f} {M1_g:16.5f} {M1_g:16.5f} {M1_g:16.5f}")
row = [0.0, 0.0, 0.0]
M0_w = Lam
for idx, eps in enumerate([1e-2, 1e-4, 1e-6]):
    s = np.linspace(eps, 60*Lam, 2_000_000)
    rho = s*s*np.exp(-s/Lam)                            # = dnu/ds for the gapless witness
    row[idx] = float(np.trapezoid(rho * s**(-4), s))    # M_1 = int rho s^{-4} = int e^{-s/Lam} s^{-2} ds
print(f"{'gapless witness s^2 e^{-s/Lam}':>34} {M0_w:14.5f} {row[0]:16.4f} {row[1]:16.4f} {row[2]:16.4f}")
OUT["part4_kl_moments"] = {"gapped_M0": M0_g, "gapless_M0": M0_w,
                           "gapless_M1_eps": [1e-2, 1e-4, 1e-6], "gapless_M1_vals": row}
print(f"""
Reading: M_0 is FINITE and positive for BOTH measures (this is the anomaly-type sum-rule moment --
see Part 4: it is BLIND to the gap).  M_1 is finite for gapped, but DIVERGES ~ 1/eps for the
gapless witness (values above grow as eps shrinks).  The gap is encoded EXACTLY in the divergence
of the higher KL moments = the IR non-analyticity of the LINEAR response.  [T]
=> "analytic at w=0 <=> gap" is CORRECT -- but it is the KL spectral representation restated.
   To DECIDE it you must know whether all M_n are finite, i.e. the spectrum at 0, i.e. the gap.
   The linear functional is a REFORMULATION of the gap (circular as a tool).""")

# ---------------------------------------------------------------------------
hr("PART 4 [V-here] -- the anomaly/low-energy sum rule is BLIND to the gap (non-discriminating)")
print(f"""
The trace-anomaly / dilaton sum rule pins the CONVERGENT moment M_0 = int rho(s)/s^2 ds.
For the SAME two measures:
   gapped delta(s-0.7):     M_0 = 1/m^2 = {M0_g:.4f}
   gapless witness:         M_0 = Lambda = {M0_w:.4f}
BOTH finite, BOTH positive, same functional form -- the anomaly moment cannot tell them apart.
This is exactly the witness dnu = s^2 e^{{-s/Lambda}} ds: positive, gapless (support down to s=0),
yet saturates the convergent sum rule.  So a functional that is DEFINED/pinned by the anomaly
moment inherits this blindness => NON-DISCRIMINATING (it returns the same value gapped vs gapless).
The gap lives in the NON-convergent higher moments (Part 3), i.e. precisely in the IR
non-analyticity -- which the anomaly does NOT fix.  Hence:
   * N defined BY the anomaly moment  => blind / non-discriminating        [horn (b)]
   * N defined AS the IR non-analytic part => discriminating in principle, but computing it
       = knowing the higher moments are finite = knowing the gap = circular  [horn (a)]""")

# ---------------------------------------------------------------------------
hr("PART 5 [V-here] -- nonzero connected 3-pt is about INTERACTION+SKEW, not about being gapless")
print("""
To nail the orthogonality: a GAPPED but NON-GAUSSIAN theory generically has chi^(3)_connected != 0
(equal-time), while a GAPLESS but GAUSSIAN theory has it == 0.  So the equal-time cubic correlator
does NOT track the gap either.  1D lattice phi^4 (gapped, m>0) Metropolis, pointwise 3rd cumulant:
""")
def metropolis_phi4_1d(L, m, lam, h, delta, nsweep, nequil, nthin=2):
    phi = np.zeros(L)
    V = lambda x: 0.5*m*m*x*x + lam*x**4 + h*x
    def dS(i, nv):
        ip = (i+1) % L; im = (i-1) % L
        # kinetic 1/2 (dphi)^2 term (a=1)
        old_k = 0.5*((phi[i]-phi[im])**2 + (phi[ip]-phi[i])**2)
        o = phi[i]; phi[i] = nv
        new_k = 0.5*((phi[i]-phi[im])**2 + (phi[ip]-phi[i])**2)
        dv = V(nv) - V(o)
        phi[i] = o
        return (new_k - old_k) + dv
    samples = []
    for sw in range(nsweep):
        for i in range(L):
            nv = phi[i] + delta*np.random.standard_normal()
            d = dS(i, nv)
            if d <= 0 or np.random.random() < np.exp(-d):
                phi[i] = nv
        if sw > nequil and (sw % nthin == 0):
            samples.append(phi.copy())
    return np.array(samples)

def third_cumulant(x):
    # kappa_3 = E[(x-mu)^3]  (central third moment)
    mu = x.mean(); d = x - mu
    return float(np.mean(d**3)), float(np.std(d**3)/np.sqrt(x.size))

L = 96
for (m, lam, h, tag) in [(1.0, 1.0, 0.0, "Z2-symmetric (h=0): symmetry => 3pt=0"),
                         (1.0, 1.0, 0.4, "explicit breaking (h=0.4): non-Gaussian skew => 3pt!=0"),
                         (1.0, 0.0, 0.0, "free Gaussian (lam=0,h=0): 3pt=0")]:
    smp = metropolis_phi4_1d(L, m, lam, h, delta=1.2, nsweep=6000, nequil=1500, nthin=3)
    x = smp.reshape(-1)
    k3, se = third_cumulant(x)
    print(f"  m={m}, lam={lam}, h={h}:  kappa_3 = {k3:+.4f} +- {se:.4f}   [{tag}]")
OUT["part5_phi4"] = {"note": "nonzero connected 3pt appears for gapped non-Gaussian skewed; zero for free and for symmetric"}
print("""
=> a GAPPED interacting (broken-symmetry) theory has chi^(3)_connected != 0.  The cubic correlator,
   like the frequency-space IR function, is sensitive to INTERACTION + ORDER PARAMETER + SYMMETRY,
   not to the existence of a gap.  [V-here]""")

# ---------------------------------------------------------------------------
hr("SUMMARY TABLE -- the dichotomy, model by model")
print("""
model                      | gap? | Gaussian? | linear N (analytic@0) | nonlinear N^(3) IR | acid?
---------------------------+------+-----------+------------------------+--------------------+---------------------------
 free massive scalar       | yes  | yes       | analytic  (N=0)        | 0  (chi^3==0)      | ok
 free MASSLESS scalar      | NO   | yes       | NON-analytic (N!=0)    | 0  (chi^3==0)      | *** FAILS ***
 interacting GAPPED        | yes  | no        | analytic  (N=0)        | 0 (support>gap)    | ok
 interacting GAPLESS/crit. | NO   | no        | NON-analytic           | !=0 (propto g)     | ok, but propto g
---------------------------+------+-----------+------------------------+--------------------+---------------------------
 KEY ROW = free massless: GAPLESS yet N^(3)=0 = the gapped value  => false negative (Gaussianity).
 The linear column is correct but is the KL spectrum (reformulation); the nonlinear column is
 nonzero only when gapless AND interacting AND the operator couples with non-analytic IR structure.
""")
OUT["verdict"] = {
    "free_gapless_N3": 0.0,
    "nonlinear_parity": "N^(3) propto interaction strength g; -> 0 as g->0 while spectrum stays gapless",
    "linear": "faithful (analytic@0 <=> gap) but == Kallen-Lehmann spectrum => reformulation/circular",
    "anomaly_moment": "blind (M_0 finite & positive for gapped AND gapless witness) => non-discriminating if anomaly-pinned",
    "overall": "(a) REFORMULATION, with the nonlinear horn additionally broken by Gaussianity (free-vacuous, in the mislabeling direction); partially (b) if anomaly-derived. NOT (c).",
}

import json, os
p = os.path.join(os.path.dirname(os.path.abspath(__file__)), "ir_response_test_results.json")
with open(p, "w") as f:
    json.dump(OUT, f, indent=2)
print(f"[wrote {p}]")
print("\nDONE. No QED. The free-field Gaussianity breaks the nonlinear discriminator (N^(3)(free massless)=0);")
print("the linear one is the KL spectrum (reformulation); the anomaly moment is blind. Verdict: (a)/(b), not (c).")
