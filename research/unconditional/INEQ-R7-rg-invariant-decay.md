# INEQ-R7 — Does an RG-invariant decay constant exist that controls f₁²?

**Verdict.** Yes, an RG-invariant decay constant exists — and it is `f₁` in disguise. The one exactly
protected combination is the trace anomaly `T^μ_μ=(β/2g)O_S`, `O_S≡¼trF²` (INEQ-T4). Define the glueball
decay constant `F_T≡⟨0|T^μ_μ|G⟩` (dim 3). Since `T^μ_μ` is the conserved dilatation-current divergence,
`F_T` is **RG-invariant**, and the operator identity gives the *exact* relation

  `f₁(μ) = F_T / (β/2g)(μ)`  ⟹  `f₁² = F_T² / (β/2g)²`.

The anomalous dimension of `f₁` is **not unknown**: it is precisely the running of `β/2g`. So the
"scheme-independent norm" *equals* `f₁²` up to a known factor; nonzero-`F_T` and nonzero-`f₁` are the
*same* statement. It cannot **independently** lower-bound `f₁²`. The only quantity provably nonzero
without dynamics — the UV anomaly/`⟨T^μ_μT^ν_ν⟩` coefficient (free central charge, exactly computable) — is a
*short-distance/aggregate* datum; by the moment lemma it **upper**-bounds the edge (and for dim-4 `O_S` the
total moment diverges, only the log coefficient is scheme-independent, GAP-G2-edge). Wrong sign.

**Construct where control holds.** At `β<1/48` (SZZ; λ>96) the convergent cluster expansion gives a uniform
gap; `F_T` is a convergent series with leading term `O(1/a³)≠0`. A *residue* bound `f₁²≥r₀` is a stronger
extraction than the cited decay-rate theorem — plausible inside the convergence radius, **not** a verbatim
corollary (flagged, no overclaim). In the AF UV the anomaly coefficient is the clean, rigorously-nonzero
RG-invariant.

**Transport — and the collapse.** RG-invariance ⇒ `F_T` is *constant along the line of constant physics*,
**per state**. Single confining phase (crossover, no transition) ⇒ continuity across `β<1/48→2.5`. Then the
exact relation *determines* `f₁²(β)=F_T²/(β/2g)²>0` — **for the state that carries `F_T`.** Here it breaks:
`F_T` is a **per-state** matrix element; `f₁²=ν({m₁²})` is an **edge** label (lightest state of the `O_S`
measure). RG-invariance tracks *its state*, not the edge. The 0⁺⁺ spectrum is unprotected; an excitation may
dive below the protected state across the crossover (no no-level-crossing theorem). `F_T` then stays nonzero
yet is **no longer the edge residue** — the lightest state's `⟨0|O_S|lightest⟩`, hence `f₁²`, can shrink to
zero. This is verbatim INEQ-T4's one-state-dominance flip.

**Bottom line.** RG-invariance buys *constancy*, not *nonvanishing-at-the-edge*. The decay constant of the
*lightest* glueball presupposes its spectral ordering — the IR gap. Transporting a nonzero value to `β≈2.5`
**is** a no-level-crossing / spectral-saturation theorem, equivalent to `H*(r₀)` (RESIDUE-CONDITIONAL).
**Exists: yes. Provably nonzero at the crossover: no — collapse to the gap.**

**Flip condition.** Prove the protected trace-charged state never leaves the 0⁺⁺ edge (a spectral-ordering
identity), or replace `f₁²` of the lightest by a channel-blind protected coercion — neither exists.
