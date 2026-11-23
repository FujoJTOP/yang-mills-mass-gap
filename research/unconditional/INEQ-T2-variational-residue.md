# INEQ-T2-variational-residue — a dual/min–max characterization of the edge atom `f₁²`: lower bound, or collapse?

> 2026-11-23 · The residue-specific dual, constructed head-on. **Verdict: the exact dual exists and is genuinely different in form from the gap's Courant–Fischer min–max — but it is an INFIMUM (it yields upper bounds); its conjugate sup degenerates to 0. No residue-dual is a lower bound from OS/protected data. The only proven residue lower bound consumes `Λ=E₂` and `m₁²`.** Everything marked (theorem) is operator/measure theory, not YM-specific. No fabricated proof.

**Setup.** `ν≥0` on `[m₁²,∞)`, `f₁²=ν({m₁²})` the bottom atom, Stieltjes `F(z)=∫dν/(t−z)`, correlator `C(t)=∫e^{−ts}dν`.

**The dual pair (constructed).** Primal = positive measures; dual = bounded continuous tests `φ`, pairing `⟨φ,ν⟩=∫φdν`. Bump functions decrease to `1_{m₁²}` (a singleton is closed), so
`f₁² = inf{ ⟨φ,ν⟩ : φ ≥ 1_{m₁²} }`. **(theorem)** Each admissible test obeys `⟨φ,ν⟩ ≥ f₁²`: an *upper* bound. The conjugate lower dual `sup{⟨ψ,ν⟩ : ψ ≤ 1_{m₁²}}` **equals `0`** — continuity forces any such `ψ≡0`, because `1_{m₁²}` is the indicator of a *point*, hence usc, and its largest continuous minorant is `0`. So the naive lower dual certifies nothing below `f₁²>0`.

**Semicontinuity obstruction (the residue-level sign wall).** For `F` closed, portmanteau gives `μ(F) ≥ limsup μₙ(F)`; with `F={m₁²}`, the atom `ν({m₁²})` is **upper semicontinuous** — an infimum of continuous functionals, never a supremum. A lower-bound dual needs lsc, which the atom is not (smeared `δₙ→δ` yet atoms `0→1`). This is the measure-level twin of EDGE-f₁²'s sign obstruction: *every* natural dual is an infimum.

**Correlator dual — the physically meaningful one.** By the spectral theorem
`e^{t m₁²}C(t) = f₁² + ∫_{rest} e^{−t(s−m₁²)}dν`, which **decreases to `f₁²` from above**; hence `f₁² = inf_t e^{t m₁²}C(t)` — again an infimum of finite-resolution (smeared) correlators, each `≥ f₁²` → upper bounds. Two exact collapses:
- **location-collapse.** The formula needs the exponent `m₁²`: with `μ>m₁²` the limit is `0`, with `μ<m₁²` it diverges. The dual is informative *only at the true gap* — `m₁²` is baked in.
- **residual-collapse.** To convert the infimum into a *lower* bound one must control `∫_{rest}e^{−t(s−m₁²)}dν ≤ ε` uniformly in `t`, i.e. bound spectral weight above `m₁²` relative to the atom — precisely the `E₂` / two-scale (second-threshold) input.

**Conic / feasibility dual (Markov–Krein).** Separating hyperplanes certify support-infeasibility via `p≥0` on `supp`, `∫p dν<0`. But `f₁²≥r₀` is not a linear/moment constraint — the atom is usc, invisible to finite linear functionals (linear-data lemma) — so no separating hyperplane exists; and the free field is feasible (`ν=δ₀`), so the dual cone is **empty**. **(theorem)**

**Consistency with the known lower bound.** The *only* proven residue lower bound, the residue–level-gap lemma `f₁² ≥ (ΛM₀−M₁)/(Λ−m₁²)`, itself consumes `Λ=E₂` and `m₁²`. The residue-dual cannot manufacture this without the second level.

**Sharpest line.** The residue is one coefficient, and its exact dual *is* a distinct object from the gap's Courant–Fischer min–max — but it is an **infimum** (upper bounds), its conjugate sup vanishes, and portmanteau makes that unavoidable; even its exponent presupposes `m₁²`. A gap-free residue lower dual would need a Schur-convex coercivity against a usc functional — the same reversed sign. The residue-dual is not an escape: it is the `E₂`+gap obstruction, confirmed at the coefficient level.

*Refs.* EDGE-f1sq-attack (sign obstruction); LEVER-nonconvex-edge §1 (min–max direction), §0B (linear-data lemma); RESIDUE-CONDITIONAL (residue–level-gap lemma, H*); Portmanteau theorem (atom is usc on closed sets).
