# G1-UNIFORM-ZA — marginal gauge-covariant BPHZ / uniform Z_a: existence reduces to one analytic estimate, the dim-4 class closes only order-by-order, the sign is BPHZ-invisible and is (★)=H1

> **判词 (verdict).** Three-way split; only (A) closes. **(A)** The *combinatorial/perturbative* core is a theorem: BRST (ghost-#0, dim 4) leaves `{Tr F², Tr F F̃}` (θ a total derivative, drops), and asymptotic freedom absorbs the marginal log into `g_j²`, so `Σ_j g_j⁴=ζ(4)/(16β₀⁴)<∞` — the "`Σ_j g_j²~log log(1/a)` wall" *is* the running coupling, not a wall. **(B)** The **nonperturbative uniform Z_a** (`L3 = sup_a Z_a⁻¹‖H_a‖_D<∞` on a common core = `UVS(L)`) does **not** close: it equals **(★)∧(★dec)**, a one-sided *upper* bound on force variance, and every proven tool gives the opposite (lower) edge. The dim-4 class is closed order-by-order, **not** nonperturbatively. **(C)** The **sign** is **H1** and is *BPHZ-invisible*: even a complete uniform Z_a yields fixed-volume existence, never `G″(0)>0`. G1 thus advances from "no known idea" to "fixed-volume existence conditional on one UV analytic estimate, ⊥ gap"; it does not close, and the residue is hard analysis (B), not a no-go, plus the gap itself (C).

## 1. What "uniform Z_a" is
Common core `D` = smooth smeared gauge-invariant field-strength polynomials, graph norm `‖ψ‖_D=‖ψ‖+‖Fψ‖`; `Z_a` = UV-block fluctuation variance of `δ_a tr F²` (A2 §L3). `L3`: `sup_a‖Ĥ_a|_D‖<∞`, `Ĥ_a=Z_a⁻¹ι_aH_aι_a⁻¹`. By Trotter–Kato this is (i) `Z_a` matches scheme C [L2, proven]; (ii) a block Poincaré `Λ_UV>0` [PROVE-H2, proven but `e^{−2β_jN_p}`]; (iii) `Λ_UV` uniform over the `J~log(1/a)` tower [open] — the whole content.

## 2. What closes (and why no no-go blocks it)
L1 (scheme C) is a finite reparametrisation — definitional. L2 is a genuine theorem: Reisz lattice BPHZ + Symanzik gives order-by-order lattice↔RT matching with the *same* counterterm structure; the dim-4 class is finite and θ drops. The scale-sum closes: `g_j²=1/(2β₀j)+O(j⁻²)`, `Σ g_j⁴` converges, `δg_J=O(J^{−3/2}log J)→0` (BREAK-IN-1). `UVS(L)` then yields compactness, a diagonal subsequence, Minlos, and OS0/1/3 (GL3/GL4; RP weakly closed, free). Lemma B gives `W≥δ_W`; the free box gap `c_L=4π²/L²` is `a`-uniform. No no-go blocks this: the corpus no-gos (translation-rigidity, index-rigidity, nested cone) all target the *sign/zero mode*, not the UV bound. Templates — BBS d=4 marginal, GN₂, KT — make (B) plausible, though it is hard analysis at exactly d=4, where RS's degree condition and Balaban's uniformity stop.

## 3. What does not close — and why the class is not "closed nonperturbatively"
The missing estimate is one inequality per shell, uniform in `j,U′,a`:
`Var_{μ(·|U′)}(⟨Pᵀ∇S_j,v⟩) ≤ E_{μ(·|U′)}⟨PᵀH_jP v,v⟩ + W_j(v) − δ`  **(★)**,
plus band decay `|Cov(s_B,s_B′)|≤Ce^{−κ|B−B′|/a}` (★dec, sketched only). The proven direction gives `H⪰0` on the typical region (cone-empty) — a *lower* edge — but `G″=E[H]−Var+W` and `Var≥0` **subtracts**, so `H⪰0` says nothing about `G″`. One-loop exactness gives `E[H]−Var=0`; the only positive term is `W≈1/Var(a)`. Hence `G″≥δ ⟺ Var(a)≤1/δ`: an *upper* variance bound no positivity supplies. The `a⁻⁴` coefficient `c(g)=Σc_ℓg^{2ℓ}` is an asymptotic IR-renormalon series whose Borel ambiguity is the gluon condensate `⟨Tr F²⟩`, Stokes constant undetermined on ℝ⁴ (CONTINUUM-ASYMPTOTIC-BEYOND). An exact nonperturbative `Z_O` exists as a running-coupling ratio, but multiplies an unknown base and cannot manufacture a uniform bound. "Class closed nonperturbatively" is precisely L3 — open.

## 4. Existence vs sign — the exact split
**Existence (may close).** Fixed `L`, `a→0`. Residue = L3 alone, a UV one-derivative statement at weak coupling `g_j→0`, ⊥ the gap and *strictly milder* than the two-sided Hessian sign (c). If (★)∧(★dec) held, G1-EXISTENCE becomes unconditional: subsequential `μ_L` with OS0/1/3. It does **not** presuppose the gap. **Sign (H1).** At physical coupling `g≈1.5` (β=2.5), typical-field/high-momentum cell, `k→0`: `G″(0)≥λ ⟺ Var_μ(a_{k→0})≤1/λ ⟺ m_eff²(0)≥λ` = H1 = the gap (W6-BE: 21–42% negative directions, `−1.12β`, 90–97% in the UV block; `PᵀHP` min≈−0.3). Ward `G″·Var(a)=1` makes sign and variance bound one number. BPHZ constrains counterterm structure, not Hessian sign. Even a complete L3 gives no clustering, uniqueness, gap, or `L→∞` (IN-6/G2).

## 5. Ledger
| item | status |
|---|---|
| L1 scheme C | closed |
| L2 all-orders match; dim-4 `{Tr F², Tr F F̃}` | closed (order-by-order) |
| scale-sum `Σg⁴` | closed given L2 |
| `UVS(L)⇒μ_L`, OS0/1/3 | closed implication; RP free |
| nonperturbative `Z_a` = L3 = (★)∧(★dec) | **OPEN** |
| class "closed nonperturbatively" | **OPEN** = L3 |
| sign / H1 / gap | not addressed, BPHZ-invisible |
| `L→∞`, uniqueness, clustering | not addressed (G2/IN-6) |

**Sharpest line.** Marginal BPHZ is the right *existence* language and the wrong *sign* language: it closes the operator class and the log, then stops where the estimate becomes one-sided — one upper bound on variance. That bound is hard analysis for the UV tower (existence, no no-go) and H1 verbatim at the physical zero mode (the gap). It never presupposes the gap; it also never yields it.

*Refs.* PROOF-G1-BPHZ §1–5; W3-marginal-BPHZ (0)–(4); W2-regularity-RG (1)–(4); CONTINUUM-FREDHOLM §3–5; CONTINUUM-TIGHTNESS [1]–[3]; CONTINUUM-ASYMPTOTIC-BEYOND §1–4; A2-scheme-match §L1–L3; L3-scale-uniform; PROOF-L3 §1–5; BREAK-IN-1; BREAK-IN-7(b); G1-EXISTENCE-THEOREM; G1-fixed-volume; W4-polchinski-target §1–4; W4-blocked-positivity; ABC-CONSTRUCTIVE §3; ABC2-SYNTHESIS §3–4; PROOF-IR-H1; k0-TRANSLATION-RIGIDITY; P1-L0; PROVE-H2; W6-BE-at-2.5; Bauerschmidt–Brydges–Slade LNM 2242; FMRS CMP 103 (1986) 67; Fröhlich–Spencer CMP 81 (1981) 527; Balaban CMP 109 (1987) 249, 116 (1988) 1, 122 (1989) 175/355; Magnen–Rivasseau–Sénéor CMP 155 (1993) 59; Aizenman–Duminil-Copin arXiv:1912.07973; Hairer Invent. Math. 198 (2014) 269.
