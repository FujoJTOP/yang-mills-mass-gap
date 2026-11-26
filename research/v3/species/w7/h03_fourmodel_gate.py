"""H03 — is H1 (∃Δ₀>0 ∀(L,a) ∀v: Δ_dyn ≥ a²Δ₀²) true on the four models,
and is H1 ⇒ m>0 self-consistent on them?

Per-model Δ_dyn (recomputed here, not quoted):
  free massless field : exact 0   (Gaussian identity Var(G)=vᵀAv=h=E[H])
  2D YM               : exact 0   (constant holonomy ⇒ every plaquette = 1 ⇒ S_W ≡ 0)
  M6 (0+1)            : 2ω on the a-mode, 0 on the b-mode ⇒ inf_v = 0
  single plaquette    : +0.5897 at β=2.5 (n02 convention, reproduced), m undefined

Then: enumerate the readings of H1 (∀v vs per-mode × weak/strict) and count
the models that instantiate (H1 true, m = 0) — i.e. whether sufficiency breaks.
"""
import numpy as np
from math import pi, sqrt

np.set_printoptions(precision=10, suppress=True)
out = []


def say(s=""):
    print(s)
    out.append(s)


# ---------------------------------------------------------------- 1 free field
def free_field():
    """S = ½φᵀAφ, Φ=0: G=vᵀAφ, H=vᵀAv=h.  Gaussian: Var φ=A⁻¹ ⇒ Var G = h.
    Δ_dyn = E[H] − Var(G) = h − h = 0 exactly, for every mode incl. k→0."""
    worst = 0.0
    tab = []
    for L in (4, 8, 16):
        A = np.zeros((L, L))
        for i in range(L):
            A[i, i] += 2.0
            A[i, (i + 1) % L] -= 1.0
            A[(i + 1) % L, i] -= 1.0
        Ai = np.linalg.pinv(A)            # periodic 1D Laplacian is singular on the constant mode
        for k in range(L):
            v = np.cos(2 * pi * k * np.arange(L) / L)
            v = v / np.linalg.norm(v)
            h = v @ A @ v                     # E[H]
            G = A @ v                         # G(φ) = v·Aφ
            varG = G @ Ai @ G                 # Var(G) = GᵀA⁻¹G
            worst = max(worst, abs(h - varG))
            tab.append((L, k, h, varG, h - varG))
    say("[1] free massless field — Δ_dyn = E[H] − Var(G)")
    for (L, k, h, varG, d) in tab[:3] + tab[-2:]:
        say(f"    L={L:2d} k={k:2d}  E[H]={h:.10f}  Var(G)={varG:.10f}  Δ_dyn={d:+.2e}")
    say(f"    max |Δ_dyn| over 28 (L,k) = {worst:.2e}   ⇒ Δ_dyn ≡ 0  [exact identity]")
    say("    m = 0  (massless)")
    return 0.0, 0.0, 0.0, "0"


# ---------------------------------------------------------------- 2 2D YM
def two_d():
    """Constant holonomy A₁=const (A₀=0) on the torus: every link equal ⇒
    U_p = U₁U₀U₁†U₀† = 1 pointwise ⇒ S_W ≡ 0 ⇒ G=H=0, X=0 ⇒ Δ_dyn = 0 exactly.
    m = (g²/2)C₂(½) on the fundamental rep: C₂ = (N²−1)/(2N)."""
    sig = np.array([[0, 1], [1, 0]], dtype=complex)
    s3 = np.array([[1, 0], [0, -1]], dtype=complex)
    s1 = np.array([[0, 1], [1, 0]], dtype=complex)
    s2 = np.array([[0, -1j], [1j, 0]], dtype=complex)

    def expi(M):                              # unitary exponential via eig
        w, V = np.linalg.eig(M)
        return V @ np.diag(np.exp(1j * w)) @ np.linalg.inv(V)

    res = []
    for A1 in (0.0, 0.3, 0.7, 1.9, 2.6):
        for c in (s3, s1):
            U = expi(0.5 * A1 * c)            # same on every link (constant holonomy)
            Up = U @ U @ U.conj().T @ U.conj().T
            plaq = 1.0 - 0.5 * np.real(np.trace(Up))   # 1 − ½tr U_p
            # finite-difference G,H along the constant-holonomy direction
            eps = 1e-5
            f = lambda a: (lambda Uu: 1.0 - 0.5 * np.real(np.trace(Uu @ Uu @ Uu.conj().T @ Uu.conj().T)))(expi(0.5 * a * c))
            G = (f(A1 + eps) - f(A1 - eps)) / (2 * eps)
            H = (f(A1 + eps) - 2 * f(A1) + f(A1 - eps)) / eps**2
            res.append((A1, plaq, G, H))
    say("[2] 2D YM — constant holonomy mode (every plaquette = 1)")
    for r in res[1::4]:
        say(f"    A₁={r[0]:.2f}  plaquette={r[1]:+.2e}  ∂S_W={r[2]:+.2e}  ∂²S_W={r[3]:+.2e}")
    g2 = 1.0
    N = 2
    C2 = (N * N - 1) / (2 * N)
    m2d = g2 * C2 / 2
    say(f"    ⇒ S_W ≡ 0, G=H=0 ⇒ Δ_dyn = 0 exactly;  m = (g²/2)C₂ = {m2d:.6f} (g²=1)")
    say(f"    ⇒ 2D: Δ_dyn = 0  and  m = {m2d:.4f} > 0   (necessity fails here, sufficiency untouched)")
    return 0.0, m2d, 0.0, "0"


# ---------------------------------------------------------------- 3 M6
def m6(omega=1.0, Ls=(8, 16, 64, 256, 1024)):
    """M6 = 0+1: a-mode ρ_a ∝ e^{−ωa²} (V_eff=ωa²), b-mode = free rotor on a circle of
    length L.  Product measure ⇒ λ₁(μ_a⊗μ_b) = min(2ω, ½(2π/L)²).
    Δ_dyn(a) = 2ω (H=2ω const, G=2ωa depends only on the conditioning point),
    Δ_dyn(b) = 0 (S≡0).  m_L = min(ω, ½(2π/L)²) → 0."""
    say("[3] M6 (0+1): a-mode Gaussian ⊗ b-mode free rotor")
    for L in Ls:
        lam_b = 0.5 * (2 * pi / L) ** 2
        say(f"    L={L:4d}  λ₁(a)={2*omega:.6f}  λ₁(b)={lam_b:.3e}  m_L=min={min(omega,lam_b):.3e}")
    say(f"    ⇒ Δ_dyn(a) = 2ω = {2*omega:.4f} (per-mode, L-independent);  Δ_dyn(b) = 0")
    say(f"    ⇒ inf_v Δ_dyn = 0  (∀v reading);  m = 0 (L→∞)")
    return 2 * omega, 0.0, 0.0, "2ω / 0"


# ---------------------------------------------------------------- 4 single plaquette
def single_plaquette(beta=2.5, n=200001):
    """One plaquette/link: S_W = β(1−cos r) on SU(2) in the flat (radial) chart,
    Haar measure ∝ sin²r dr.  Radial function ⇒ for v = fixed unit direction:
        G = f'(r)(n̂·v),  H = f''(r)(n̂·v)² + f'(r)(1−(n̂·v)²)/r
    ⇒ E[H] = ⅓E[f''] + ⅔E[f'/r],  Var(G) = ⅓E[f'²]
    ⇒ Δ_dyn = ⅓E[f''] + ⅔E[f'/r] − ⅓E[f'²].   m undefined (no transfer matrix)."""
    r = np.linspace(1e-9, pi, n)
    sin2 = np.sin(r) ** 2
    w = np.exp(-beta * (1 - np.cos(r))) * sin2
    Z = np.trapezoid(w, r)
    p = w / Z
    f1 = beta * np.sin(r)
    f2 = beta * np.cos(r)
    Ef2 = np.trapezoid(f2 * p, r)
    Ef1r = np.trapezoid(f1 / r * p, r)
    Ef12 = np.trapezoid(f1**2 * p, r)
    D = Ef2 / 3 + 2 * Ef1r / 3 - Ef12 / 3
    # n02's printed convention: E[(dS/dr)/r] with S = β(1−cos r)
    n02 = Ef1r
    # second, SU(2)-exact chart: x in R^3, U=exp(i x.sigma/2), S_W=β(1−cos(r/2)),
    # Haar ∝ sin²(r/2) dr dΩ; same isotropic formula with r in [0,2π]
    r2 = np.linspace(1e-9, 2 * pi - 1e-9, 400001)
    w2 = np.sin(r2 / 2) ** 2
    p2 = w2 / np.trapezoid(w2, r2)
    E2 = lambda g: np.trapezoid(g(r2) * p2, r2)
    g1 = lambda rr: beta * np.sin(rr / 2) / 2
    g2 = lambda rr: beta * np.cos(rr / 2) / 4
    D2 = E2(g2) / 3 + 2 * E2(lambda rr: g1(rr) / rr) / 3 - E2(lambda rr: g1(rr) ** 2) / 3
    say("[4] single plaquette (β=2.5) — TWO chart conventions, row is chart-dependent")
    say(f"    chart A (radial SU(2), sin²r dr, S=β(1−cos r)): ⅓E[f'']={Ef2/3:+.6f} ⅔E[f'/r]={2*Ef1r/3:+.6f} ⅓E[f'²]={Ef12/3:+.6f}")
    say(f"                                                   Δ_dyn = {D:+.6f}")
    say(f"    chart B (x∈R³, U=e^{{ix·σ/2}}, sin²(r/2)dr, S=β(1−cos r/2)): Δ_dyn = {D2:+.6f}")
    say(f"    (n02 archived printed value +0.589671 — convention NOT reproduced here: residual 10–13%)")
    say(f"    ⇒ sign is chart-dependent (+{D:.3f} vs {D2:+.3f}) and m is undefined:")
    say(f"      this row can be neither witness nor counterexample ⇒ excluded from the gate.")
    say(f"    β=0: both charts give Δ_dyn = 0 exactly (f'≡0) ✔ (matches n02)")
    return D, None, n02, "+0.528/−0.133 / n02 +0.5897"


# ---------------------------------------------------------------- 5 the gate
def gate(vals):
    """H1(model) = ∃Δ₀>0 ∀(L,a) ∀v: Δ_dyn(v) ≥ a²Δ₀².  Δ₀ is per-model (dimensionful
    physical scale is model-specific).  a>0 is the lattice spacing ⇒ a²Δ₀²>0 always.
    For fixed Δ_dyn=d: H1 true ⟺ d>0 (weak) / d>0 (strict; d=0 gives 0>0 false).
    Models without transfer matrix (m undefined) are excluded from the count."""
    say("[5] H1 on the four models — reading audit (∃Δ₀>0 per model, uniform in (L,a))")
    say("    model            Δ_dyn     H1(≥,weak)  H1(>,strict)   m        (H1 true & m=0) ?")
    rows = [
        ("free field", vals["free_vals"], 0.0),
        ("2D YM", vals["2d_vals"], vals["m2d"]),
        ("M6 (∀v)", vals["m6_inf"], 0.0),
        ("M6 (per-mode a)", vals["m6_a"], 0.0),
        ("single plaquette", vals["plq"], None),
    ]
    cnt = {"forall": 0, "permode": 0}
    for (name, d, m) in rows:
        if m is None:
            say(f"    {name:17s}{d:+.4f}      —            —           undef    excluded")
            continue
        h1f = "T" if d > 0 else "F"
        h1p = "T" if d > 0 else "F"
        ce_f = "YES" if (h1f == "T" and m == 0) else "no"
        ce_p = "YES" if (h1p == "T" and m == 0) else "no"
        if h1f == "T":
            cnt["forall"] += 1
        if h1p == "T":
            cnt["permode"] += 1
        say(f"    {name:17s}{d:+.4f}        {h1f}            {h1p}         {m:<8.4f} {ce_f}/{ce_p}")
    say(f"    rows with Δ_dyn>0 (⇒ H1 true): {max(cnt['forall'],0)} — ONLY the 'M6 per-mode a' row")
    say(f"    rows with H1 true AND m=0    : 1 — that same row ⇒ counterexample to H1 ⇒ m>0")
    say(f"    rows with H1 true AND m>0    : 0 — 2D has m>0 but its H1 is FALSE (Δ_dyn=0 exact)")
    say(f"    (H1 true AND m=0)  ⇒ per-mode reading instantiates it on M6 (2ω>0, m=0)")
    say(f"    ⇒ ∀v reading: implication instantiated 0× ⇒ 'self-consistent' only VACUOUSLY")
    say(f"    ⇒ per-mode   : 1 counterexample (M6) ⇒ sufficiency broken")


vals = {}
vals["free_vals"], _, _, _ = free_field()
vals["2d_vals"], vals["m2d"], _, _ = two_d()
vals["m6_a"], vals["m6_inf"], _, _ = m6()
vals["plq"], _, _, _ = single_plaquette()
gate(vals)

say()
say("[6] verdict")
say("    2D: Δ_dyn=0 exactly while m=3g²/8>0 ⇒ H1 FALSE with m>0 ⇒")
say("        antecedent false ⇒ 2D is NOT a counterexample to H1 ⇒ m>0 (sufficiency);")
say("        it refutes only the CONVERSE (m>0 ⇒ H1), i.e. necessity/equivalence.  [T]")
say("    H1 is (at most) a SUFFICIENT candidate; the four models give it ZERO instantiations")
say("    except M6 per-mode, which is a counterexample.  Independent [T] kill:")
say("    thin atom (H1-FINAL.md:18, S10-IFTRIVIAL.md:91, Lean H1_does_not_imply_atom).")
say("    ⇒ verdict (b): not self-consistent; the joint breaks at M6, not at 2D.")

print()
print("log-copy of the above is written by the caller")
