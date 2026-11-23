# INEQ-T5 — sum-rule / positivity literature hunt: does ANY proven inequality LOWER-bound the atom `f₁²`?

> 2026-11-23 · **判词：对 `f₁²=ν({m₁²})` 没有真正的例外。** 所有"由正性/解析性导出"的下界，要么下界**质量**（Weingarten/Nussinov：`m_Γ≥m_π`），要么下界**矩**（FESR Cauchy–Schwarz：`F₋₂≥F₀³/F₁²`），要么（唯一的例外）把原子钉在**零质量**（Goldstone）。对 `m₁²>0` 的单个原子，正性只给平凡的 `f₁²≥0`；peer-reviewed 文献明文：正的下界 "requires extra input"。无捏造。

**Target.** `ν≥0` (OS/reflection positivity), `supp ν⊂[m₁²,∞)`, `f₁²=ν({m₁²})`. Question: any proven `f₁²≥c>0` from positivity/analyticity alone, without assuming the answer?

**Classification — what each actually proves, the object it bounds, and the direction on the atom:**

| Tool (real refs) | Actually proves | Object | On `f₁²` |
|---|---|---|---|
| Weingarten PRL **51**,1830 ('83); Witten PRL **51**,2351 ('83); Vafa–Witten PRL **53**,535 ('84); review Nussinov–Lampert Phys.Rept. **362**,193 ('02) | `m_Γ≥m_π` (flavor-nonsinglet), baryon`≥((n−2)/(n−3))m_π`; engine = RP + γ₅-Hermiticity + Cauchy–Schwarz | **MASS**, anchored to `m_π≠0` | none; glueball sector still C–S ⇒ **upper** on coupling |
| Källén Helv.Phys.Acta **25** ('52); Lehmann Nuovo Cim. **11** ('54); spectral positivity | `ρ(s)≥0`, `ρ=Σf_i²δ+cont` | density | trivial `f₁²≥0` only; `f₁²≤M_k` is **upper** |
| FESR positive-weight (Hölder/Cauchy–Schwarz; cf. SVZ Nucl.Phys. **B147** '79) | `F₋₂≥F₀³/F₁²≥F₀²/F₂≥…`, `F_k=∫sᵏρ` | **MOMENTS** | none on atom — aggregate-to-aggregate |
| Once-subtracted dispersion, proven `Im Π≥0` sign | `Π(Q²)=∫ρ/(s+Q²)`, `Q²>0` | spacelike correlator | atom`≤Π` **upper**; a lower bound impossible (continuum adds `≥0`) |
| Oehme–Zimmermann Commun.Math.Phys. **97** ('85) superconvergence | `∫ρ σ²dσ²=λ/λ*` (gauge / indefinite metric) | moment equality | constrains moments; no atom sign |
| Bourrely–Caprini–Lellouch PRD **79**,013008; Caprini PoS CD12,006 | `ΣB_mn a_m a_n≤1`, `\|F\|≤B·exp` | form-factor coeffs | **UPPER** |
| SMEFT "beyond positivity" (arXiv:1903.08664) | lightest residue vs. cutoff under UV assumptions | EFT coupling | not YM positivity alone |

**Decisive confirmation — in Yang–Mills itself.** Caron-Huot–Pokraka–Zahraee, *Two-point sum-rules in 3D Yang–Mills*, JHEP **01** (2024) 195 (arXiv:2309.04472): with `ρ(z)≥0` alone the lightest glueball residue obeys only the *trivial* `g₁²≥0` — "positivity alone forces the minimal residue to zero" (the continuum absorbs everything: the free-field edge) — and a genuine lower bound on the lightest coupling "require[s] extra input" (there: non-rigorous Borel/resonance saturation). 3D, but the KL-positivity conclusion is dimension-independent. This *is* the campaign's sign wall, stated in the peer-reviewed sum-rule literature.

**The two genuine "lower-bound-from-positivity" families, both off-target.**
1. *RP ⇒ mass orderings.* Weingarten's Cauchy–Schwarz is the only rigorous positive lower bound in lattice SU(N) — but it bounds a **mass**, needs a nonzero anchor (`m_π`; pure glue has none), and has the form `(cross)²≤(same)(same)` ⇒ upper on couplings.
2. *Positivity ⇒ moment inequalities.* FESR Hölder genuinely lower-bounds a negative-weight moment by positive ones — but it is **saturated by a single delta** (`M_{-2}M_2=M_0²` iff ν is an atom), so it is *easiest* to satisfy for concentrated measures and says nothing forcing an atom; and isolating `f₁²` needs the continuum bounded **above** (resonance dominance/duality) — an assumption, not a theorem.

**The one genuine exception in method — and why it misses.** The Goldstone theorem (Goldstone–Salam–Weinberg '62; rigorous in Wightman/Strocchi, *Symmetry Breaking*) is the *sole* positivity+analyticity consequence forcing a specific atom to nonzero weight: a broken continuous symmetry ⇒ `ρ⊃Zδ(m²)`, `Z>0` from Hilbert-space positivity. But it pins the atom at **zero** mass, needs a spontaneously broken *continuous global* symmetry, and carries no scale — pure SU(N) YM has none (center/1-form is discrete, adjoint-blind). It cannot produce `ν({m₁²})`, `m₁²>0`.

**Sharpest line.** Positivity is convex: every OS/Cauchy–Schwarz/Hölder/unitarity inequality has the sign `(atom)≤(aggregate)` or `(cross)²≤(same)(same)`. The reverse — `total≤C·atom` (SCHUR-CONVEX coercion) — is absent from the entire literature above; forcing it needs the continuum bounded above, which is the edge theorem. Free field (`g→0`) admits `f₁²=0`, so any true bound must inject nonperturbative IR dynamics = the gap renamed.

*Refs.* EDGE-f1sq-attack; OUT-transplant; B1-nonperturbative-residue; SCHUR-CONVEX; GAP-IN11; Weingarten PRL 51,1830; Witten PRL 51,2351; Vafa–Witten PRL 53,535; Nussinov–Lampert Phys.Rept.362,193; Källén 1952; Lehmann 1954; SVZ Nucl.Phys.B147,1979; Oehme–Zimmermann CMP97,1985; Bourrely–Caprini–Lellouch PRD79,013008; Caprini PoS CD12,006; Caron-Huot–Pokraka–Zahraee JHEP01(2024)195 [arXiv:2309.04472]; Goldstone–Salam–Weinberg 1962; Strocchi, *Symmetry Breaking*.
