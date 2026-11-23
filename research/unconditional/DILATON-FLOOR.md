# DILATON-FLOOR — can the dilaton decay constant `f_G` (the non-anomaly residue) be bounded below?

> 2026-11-23 · Target: a rigorous lower bound `f_G ≥ c > 0` — equivalently the trace decay constant `⟨0|T^μ_μ|G⟩ = f_G m_G²` (INEQ-R7's `F_T`), the non-anomaly residue the anomaly cannot supply (G2B) — from **(a)** reflection positivity + the exact `H`, **(b)** 't Hooft anomaly matching, **(c)** K-complexity. **判词 (verdict): NO-GO, 2/5.** None supplies a lower bound; the one identity that exists says `f_G` *is* the residue rescaled by `m_G⁴/(β/2g)²`. So `f_G ≥ c` ⟺ `f₁² ≥ r₀` ⟺ H*/I2, the isolation/upper-gap rung. **The dilaton decay constant is the unprotected order parameter in dilaton clothing.** No `5/5 ACHIEVED`.

## 1 · The exact dictionary
By the trace identity `θ^μ_μ ≡ (β/2g)O_S`, `O_S=¼trF²` (INEQ-T4), the dilaton is normalized `⟨0|θ^μ_μ|G⟩=f_G m_G²`, and the same identity gives `f₁≡⟨0|O_S|G⟩ = F_T/(β/2g)`. Hence, identically,
  **`f₁² = f_G² m_G⁴ / (β/2g)²`**,
so a floor `f_G ≥ c` is *literally* the residue floor `f₁² ≥ r₀` with `r₀ = c²m_G⁴/(β/2g)²`. The soft-scaling `f₁²∝m⁴` is exactly the structural blindness G2B found. The question is whether *any* protected structure fixes the *magnitude* `f_G`; the algebra says it is as unprotected as `f₁`.

## 2 · Reflection positivity + exact H — upper only
`OS`/KL positivity gives `ν_S ≥ 0`; `f_G²` is a modulus-square. Every positivity functional (Cauchy–Schwarz, Hadamard/determinant, OS inequalities) is Schur-concave and bounds a matrix element from **above** (GAP-G2-edge; EDGE-f1sq §sign; GOLDS-POS-MASS). Concretely `ν_S ≥ f₁²δ_{m₁²}` with any fixed total caps the atom's share, `f₁² ≤ M₀` — an **upper** bound; and for dim-4 `O_S` even `M₀` diverges, leaving only the anomaly log coefficient, still an aggregate cap (GAP-G2-edge). As `m_G→0` the cap `f_G² ≤ (total)/m_G⁴` *weakens*, so it cannot forbid `f_G→0`.

The one positivity→atom bridge is Goldstone, needing a conserved charge pinning the atom at `m=0`, weight `‖Q|0⟩‖²` (GOLDS-POS-MASS). YM has **no** conserved dilatation charge: the anomaly `θ^μ_μ=(β/2g)O_S≠0` explicitly breaks it (`β≠0`), and the Lorentz-scalar condensate `⟨O_S⟩` does not break Poincaré (no degenerate vacua), so Goldstone's lemma does not apply. There is no `[H,Q_D]=0` anchor for a *positive-mass* `f_G`. **RP + exact H give `f_G ≤ C`, never `f_G ≥ c`.**

## 3 · 't Hooft anomaly matching — aggregate / line-sector
Both relevant anomalies fix aggregates. **(i) Trace/scale matching** fixes the UV coefficient (central charge / the `β²` log of `⟨T^μ_μT^ν_ν⟩`), a *moment* of `ν_S` — an upper bound by the moment lemma; and at the pole the Ward `f₁²⟨G|θ|G⟩=(β/2g)f₁²g_GGG` makes `f₁²` **cancel** (identity, not inequality), pinning the *forward/diagonal* `g_GGG`, never the off-diagonal `f_G` (INEQ-T4). **(ii) 1-form `Z_N^{(1)}` matching** lives in the line/defect Hilbert space; a local `trF²` cannot carry a 1-form charge, so it is provably blind to the `0⁺⁺` residue (PROOF-1FORM-ANOMALY). An anomaly forces a state to *exist*, not to be *heavy* or *coupled* (GAP-IN11 §4). Matching buys existence + contact-term normalization — the diagonal/aggregate side; `f_G` is off-diagonal/edge.

## 4 · K-complexity — a moment functional, not a 4D theorem
Krylov (K-)complexity is built from Lanczos coefficients `b_n` of the Liouvillian; the `b_n` are the continued-fraction data of the moments `M_k=∫s^k dν_S` (Parker et al. 2019). It is therefore a **moment functional**, and the moment lemma is decisive: countably many positive moments do not determine `inf supp ν` or its atom (`δ₀+δ₂` and `2δ₁` share `M₀,M₁`; GAP-G2-edge). K-complexity is moment-blind to the edge — the same wrong-sign class as every positivity functional. Independently, in 4D continuum it needs a UV cutoff (Lanczos on an unbounded Liouvillian) and is cutoff-dependent, so it is not even a theorem here. **No.**

## 5 · Honest answer — it is the unprotected order parameter
`⟨0|θ|G⟩=f_G m_G²` is a *definition*; `f_G` is the scale-symmetry-breaking order parameter, the analogue of `f_π` in chiral-scale EFT (Crewther–Tunstall). `f_π` is fixed by the conserved axial current; scale symmetry has no such current because the anomaly is precisely what makes the dilaton massive *and* leaves `f_G` unprotected — and the anomaly is an RG-*covariant aggregate* ("constancy, not edge"; CONTINUUM-VARIANCE-DEFECT). The INEQ-T4 flip condition is verbatim this target: "a theorem forcing the trace decay constant `⟨0|T^μ_μ|G⟩` to saturate the protected forward charge — a one-state-dominance / spectral-saturation statement. None known; tantamount to the IR solution." So `f_G ≥ c` = `f₁² ≥ r₀` = H*/I2 = the gap one rung up. **Boundable neither below by these tools nor independently of the residue — it is the residue.**

## 5-point verdict (target `f_G ≥ c`)
1. **Nonperturbative — PASS.** The candidates (anomaly, soft theorem, Krylov) are model-independent.
2. **Controls SIGN / lower-bounds the atom — FAIL.** Positivity Schur-concave; anomaly gives `f_G² ≤ (aggregate)/m_G⁴` (upper, weakening as `m_G→0`); K-complexity moment-blind.
3. **Non-circular — FAIL.** `f_G ≥ c` ⟺ `f₁² ≥ r₀` ⟺ H*/I2 ⟺ the gap.
4. **Non-SUSY + fails free field honestly — PASS (vacuous).** Free massless and even *gapped* free `:φ²:` have `f_G=0`; tools correctly silent, no spurious atom.
5. **Uniform `λ` over `[1/48, 2.5]`, `a→0` — FAIL.** No uniform lower bound; `f_G` can fall with `m_G`.

**2/5 — NOT 5/5.**

**Sharpest line.** A floor on the dilaton decay constant and a floor on the `trF²` residue are the *same inequality* under the trace identity. The anomaly fixes only the aggregate `f_G²m_G⁴≤∫ρ_θ` (upper, *weakening* as `m_G→0`); 't Hooft matching fixes diagonal/aggregate or line-sector data; K-complexity is a spectral-moment functional — all blind to the edge atom. With no conserved dilatation charge, the one positivity→atom bridge is unavailable by *symmetry*, not computation. Proving `f_G ≥ c` is proving the isolation H*, the gap shifted one rung. **The dilaton floor is not a new door; it is the same unprotected order parameter.**

*Refs.* INEQ-T4 (pole Ward homogeneous in `f₁²`; forward≠off-diagonal; flip condition); INEQ-R7 (`f₁=F_T/(β/2g)`, per-state collapse); GAP-G2-edge (moment lemma); GOLDS-POS-MASS (only bridge pins `m=0`; no `Q_D`); PROOF-1FORM-ANOMALY (local `trF²` blind to `Z_N^{(1)}`); GAP-IN11 §4; G2B (soft-scaling `f₁²∝m⁴`); EDGE-f1sq-attack (Schur-convex wall); PROOF-F1SQ (gapped free `:φ²:` ⟹ `f₁²=0`); CONTINUUM-VARIANCE-DEFECT (anomaly = RG-covariant aggregate); G2-ASSEMBLE-100 (`H*=I2=H1`). External: 't Hooft 1979; Gaiotto–Kapustin–Seiberg–Willett arXiv:1412.5148; Jaffe, Phys. Rev. 158 (1967) 1454; Osterwalder–Schrader CMP 31 (1973) 83; Parker–Cao–Avdoshkin–Scaffidi–Altman, Phys. Rev. X 9 (2019) 041017, arXiv:1812.08657; Crewther–Tunstall, Phys. Rev. D 91 (2015) 034035, arXiv:1411.3972.
