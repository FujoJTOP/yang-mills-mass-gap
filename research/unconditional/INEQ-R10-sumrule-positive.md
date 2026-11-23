# INEQ-R10 — a sum rule with a PROVEN-POSITIVE integrand that `f₁²` DOMINATES: does one exist?

> 2026-11-23 · 判词：**不存在，且是定理级的不存在。** 任何"被证为正的被积函数"的求和规则，对 `f₁²` 只能是**上界**；若强收敛（superconvergence）且被积函数定号，规则**自毁**（被积函数≡0）。四族真实求和规则逐一归类，**每一族都落入"齐次"（homogeneous，`f₁²` 抵消）或"聚合"（aggregate，矩/总量）**；第三类——点隔离 coercivity（`total ≤ C·atom`）——在文献中从未出现。无捏造。

**Target.** `ν≥0` (OS/RP), `supp⊂[m₁²,∞)`, `f₁²=ν({m₁²})`. Want an exact relation with a sign-definite weight `w≥0` in which `f₁²` **dominates** a proven-positive integral ⇒ `f₁² ≥ …`.

**The sign lemma (this is the whole answer).** Every sum rule is `∫₀^∞ ρ(s)w(s)ds = C`, with `ρ≥0` (unitarity / Källén–Lehmann, 1952/54).
- *Positive weight `w≥0`.* Then `f₁² w(m₁²) ≤ ∫ρw = C`. With `C` a known positive number (perturbative/anomaly coefficient) this reads `f₁² ≤ C/w(m₁²)` — an **UPPER** bound. A **lower** bound needs `C` split as `f₁² + (rest)` and `rest` bounded **above**; that `rest≤` *is* resonance dominance, not a theorem. Same for once-subtracted dispersion `Π(Q²)=∫ρ/(s+Q²) ≥ f₁²/(m₁²+Q²)` (upper) and the Weingarten "QCD inequality" `(cross)²≤(same)(same)` (upper on coupling; INEQ-T5).
- *Sign-indefinite weight* (Weinberg difference `∫(ρ_V−ρ_A)w=0`). Integrand not proven-positive ⇒ no inequality; the relation is **homogeneous/aggregate** (a moment).
- *Sign-definite superconvergence* `∫ρw=0`, `w≥0` ⇒ `ρw≡0` — **vacuous**. So "proven-positive integrand" ∧ "superconvergence" is impossible unless `C>0`, which returns to the upper-bound case.

This is exactly the peer-reviewed YM statement: Caron-Huot–Pokraka–Zahraee, JHEP **01** (2024) 195 [arXiv:2309.04472]: with `ρ≥0` alone the lightest glueball residue obeys only `g₁²≥0`; "positivity alone forces the minimal residue to zero," a genuine lower bound "requires extra input."

**The four real families — each wears one of the two masks.**
1. *Positive-weight FESR / SVZ / KL moments* (SVZ Nucl.Phys.**B147**,385 ('79); Logunov–Soloviev–Tavkhelidze '67). `∫sᵏρ=M_k` (k even) ⇒ upper via moment lemma; Hölder `F₋₂≥F₀³/F₁²` is aggregate↔aggregate and saturated by a single delta (INEQ-T5). **Aggregate.**
2. *Weinberg sum rules* (Weinberg PRL **18**,507 ('67)). `∫(ρ_V−ρ_A)=0`, `∫s(ρ_V−ρ_A)=0`: sign-indefinite; lightest V may dominate but the rest has no definite sign; isolating `f₁²` needs the A-continuum above ⇒ duality *assumption*. **Homogeneous/aggregate.**
3. *Trace-anomaly / Chanowitz–Ellis scalar-glueball rule* (Chanowitz–Ellis PRD **7**,2490 ('73); `1/m_{f₀}²·|⟨f₀|(β/4α_s)G²|0⟩|² ≃ −4⟨(β/4α_s)G²⟩`). The atom multiplies `1/m²` against the **gluon condensate** — an aggregate, renormalon-ambiguous, convention-sign RHS that is *not* a proven-positive 4D number; and the rule already *assumes* single-resonance saturation. The atom does **not** dominate: the contact/anomaly term and continuum do (INEQ-T4). **Aggregate (and assumed).**
4. *Oehme–Zimmermann superconvergence* (CMP **97**,569 ('85)): gauge / indefinite-metric moment equality; no atom sign.

**What WOULD work — and it is the gap.** A rule `T ≤ C·f₁²`, `T=∫φ dν`, `φ>0`, `C` known: the Schur-convex **coercivity** (SCHUR-CONVEX; EDGE-f₁²). It is the third, absent class (*point-isolating*). It never appears because `f₁²=ν({m₁²})` is a measure of a **singleton** = upper-semicontinuous (INEQ-T2), while every exact sum rule pairs `ν` with a (limit of a) continuous weight and returns an **infimum** of such pairings ⇒ upper bounds. Constructing it is edge dominance = the gap.

**Only proven lower bound, and it is not a sum rule.** residue–level-gap `f₁² ≥ (ΛM₀−M₁)/(Λ−m₁²)`, `Λ=E₂` (GAP-G2-edge; INEQ-T1/T2) — consumes the second level, the one non-moment input.

*Refs.* INEQ-T1/T2/T4/T5;EDGE-f1sq;SCHUR-CONVEX;GAP-IN11;GAP-G2-edge; Weinberg PRL 18,507;SVZ Nucl.Phys.B147,385;Chanowitz–Ellis PRD 7,2490;Oehme–Zimmermann CMP 97,569;Källén 1952;Lehmann 1954;Caron-Huot–Pokraka–Zahraee JHEP 01(2024)195 [arXiv:2309.04472].
