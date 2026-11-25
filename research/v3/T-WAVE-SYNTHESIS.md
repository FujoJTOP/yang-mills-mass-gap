# T-WAVE-SYNTHESIS — no standard expansion reaches the binding criterion at β=2.5; the corrections: the carrier exists, binding is a mechanism not an axis, and isolation needs a bounded ratio

> **Correction (U wave, 2026-11-24).** Two formulations below are superseded. (i) The uniform residue
> floor `W: w₁≥w₀` is a **weight** datum, not isolation: `W⇒A` only, and `W` is incomparable with
> `I` (U6). The phrase ``w₁≥w₀ = isolation = H*'' is a weight-vs-support conflation. (ii) The gap has
> two independent sufficient routes, `W∧H1` (squeeze) and `I∧(ρ²≤C²)` (window), not one. The
> corrected picture is in `U-WAVE-SYNTHESIS.md` §2.

> **VERDICT (no QED, Clay 0%).** Six routes attacked the binding/isolation criterion left by the R
> wave. None proves it. The wave's value is a **correction package** plus a **sharp no-go**:
>
> 1. **The carrier exists** (T4). The total-normalized zero mode of a gauge-invariant density at a
>    fixed physical scale has `Var(q)=χ=O(1)`, not `V^{-1}`. The corpus's `Var∝1/V` was the
>    mean-normalized mode — a bookkeeping artifact, not a no-go. What remains is UV renormalization
>    and the ceiling direction, not carrier existence.
> 2. **Binding is not an independent axis** (T6). On a Källén–Lehmann measure
>    `I ⇒ B_spec ⇒ A` (both strict), so spectral binding is implied by isolation; dynamical binding
>    is the atom's mechanism. The independent axes are `A` (atom visibility) and `I` (isolation);
>    the gap `G` is orthogonal and needs H1/ceiling plus an absolute scale.
> 3. **No standard expansion reaches β=2.5** (T1, T2, T5). The polymer activity is `1.014`
>    (`u⁴>1`), the Ursell ratio `4.6×10³`, the block-spin help `1.03×`; the tube overlap has only
>    upper bounds; the center-vortex route is neutral-blind. The missing factor is `40–55×`,
>    equivalent to the gap.
> 4. **Isolation needs a bounded ratio** (T3). The window theorem: with `Δ=(ρ²−1)μ₁`, isolation plus
>    `1+c≤ρ²≤C²` gives `μ₁≥δ/(C²−1)`; the double well has `ρ²→∞` and no gap. The upper bound is
>    load-bearing.
>
> The wave does not reduce Clay. It fixes the **statements** of the four predicates (`A,W,I,G`) and
> shows that the missing objects are a genuinely nonperturbative weight floor and an independent
> isolation bound.

---

## §0 · The six verdicts

| Route | Verdict | Result |
|---|---|---|
| **T1** polymer RG at β=2.5 | `[T]`/`[NG]`/`[O]` | Exact SU(2) activity `u=c₁/c₀=0.0104` at the P3 wall but **1.014** at β=2.5, so `u⁴=1.06>1` and the four-wall tube does not decay (`m₁a=−0.057<0`). Trivial-character fraction falls `0.98→0.17` (dense gas); Ursell ratio `3.1→4.6×10³`. A one-loop block flow gives `u_eff=0.988` over `ξ/a≈1.3`, only `1.03×` help, because there is no scale separation. Missing: `u_eff≤1/λ_a`, equivalently `|w(P)|≤e^{−κ|P|}` with `κ≥ln λ_a`; missed `40–55×`. |
| **T2** flux-tube overlap | `[T]`/`[NG]`/`[O]` | Minimal polymer: straight four-wall tube, weight `a₀u^{4T}`, `a₀>0`. Scheme-free overlap is the atom fraction `w₁=f₁²/M₀`; RP gives only log-convexity (holds at `w₁=0`), PF is vacuum-anchored, the transfer matrix gives only a decreasing upper bound. The object does not collapse with tube length (`w₁` is `T`-independent; roughening is a mass); the estimator does. Missing: uniform `w₁≥w₀>0`, equivalently leakage `Σ_{n≥2}w_n≤1−w₀`. This is the **weight** datum `W`, not isolation; `W` and `I` are incomparable (U6). |
| **T3** order-2 TP₂ | `[T]`/`[NG]`/`[O]` | Isolation `ρ²=μ₂/μ₁` is the second eigenvalue of the 0⁺⁺ transfer generator; exact for heat kernels: 2D YM `ρ²=2N²/(N²−1)∈(2,8/3]` (SU(2)=8/3), rotor 4, Mehler 2, box 4. The order-2 certificate is strict total positivity (TP₂ = Sturm–Liouville heat kernel), forcing `μ₁<μ₂`; PF only gives `μ₁>0`. TP₂ cannot be stated on the orderless quotient `𝒜/𝒢` `[NG]`. **Window theorem:** `Δ=(ρ²−1)μ₁`; isolation plus `1+c≤ρ²≤C²` gives `μ₁≥δ/(C²−1)`. Double well: `ρ²→∞`, `Δ=O(1)`, gap→0. The upper half is tractable; the lower half at β=2.5 needs a binding-induced TP₂ order (unbuilt `O(u)` wall). |
| **T4** gauge-invariant carrier | `[T]`/`[V]` | **A non-vacuous carrier exists.** The total-normalized zero mode `q=(1/√V)Σ(O−⟨O⟩)` of a gauge-invariant density at a fixed physical scale has `Var(q)=Σ_z C(z)=χ=O(1)`, independent of `V`; not `V^{-1}`. Witness: free massive scalar `Var=1/m²` for every `L` in `d=1,2,4` (slope 0); `:φ²:` `χ=4.109/0.689` `L`-independent; one-plaquette `Var=1/4`. The corpus's `Var∝1/V` is the **mean-normalized** mode `(1/V)ΣO` — a normalization artifact. Residual: bare `tr F²` has `χ~a^{-4}`, so renormalize at a fixed scale; `χ` carries only the ceiling `m₁²≤M₀/χ`. |
| **T5** center-vortex binding | `[T]`/`[NG]`/`[O]` | Tomboulis–Yaffe gives `σ>0` in the charged (N-ality-1) channel; `tr F²` is neutral, every vortex center phase it sees is `1` (witness: neutral phase 1.000000, fundamental ω). The vortex condensate is provably neutral-blind. `σ>0` does not imply an atom: 2D pure YM confines (`σ_fund=g²C₂/2>0`) yet `f₁²(tr F²)=0`. The `m_{0⁺⁺}∝√σ` bridge is measured, not proved. Positive dual-surface measures exist for Z(N)/U(1) (FK/FKG) but are sign-killed for SU(N) in 4D (Weingarten). Missing: a sign-cancelling positive regrouping in the N-ality-zero sector — the non-abelian random-surface sign problem. |
| **T6** adversarial separation | `[T]`/`[V]` | `supp ν_ac⊆supp(ν−f₁²δ)`, so `Λ★≤s_th`, hence `I⇒B_spec⇒A` (both strict). The ``isolation without binding'' cell is empty; the (binding, isolation) square is a triangle. Dynamical binding is the atom's mechanism, not an axis. `B∧I` still do not give the gap: the double well keeps `A` and large `λ₂−λ₁` but the transfer gap Kramers-collapses (`λ₁~e^{−S}`) and `ρ→∞`. True axes: `A` and `I`; `G` is orthogonal. Correction: 2D `ρ²=8/3` is the holonomy's isolation, not the `tr F²` atom (which is 0). |

---

## §1 · The four predicates

The R wave's ``atom + isolation'' pair is now replaced by four predicates:

```
   A  (atom visibility)   f₁² > 0
   W  (weight floor)      w₁ = f₁²/M₀ ≥ w₀ > 0
   I  (isolation)         Λ★ = inf supp(ν − f₁²δ_{m₁²}) > m₁²
   G  (gap)               m₁² = inf supp ν > 0
```

with the strict implications

```
   I  ⟹  B_spec  ⟹  A          (B_spec: the atom lies strictly below its continuum)
   W  ⟹  A
   W  and  I  are INCOMPARABLE  (threshold atom: W∧¬I; thin isolated atom: I∧¬W)
   G  ⊥  (A, W, I)  without a side condition
```

The gap has **two independent sufficient routes**:

```
   W ∧ H1        ⟹  G      (squeeze: m₁² ≥ f₁²/χ ≥ w₀λ)
   I ∧ (ρ²≤C²)   ⟹  G      (window theorem: Δ=(ρ²−1)μ₁, μ₁ ≥ δ/(C²−1))
```

The double well is the witness that neither route is automatic: it keeps `W` and `I` (uniform
`w₁≈0.97`, large `Δ`) while `λ₁→0`, `χ→∞`, `ρ²→∞`, so `G` collapses. So both ``binding alone'' and
``isolation alone'' fail; the gap needs `W∧H1` or `I∧(ρ²≤C²)`.

---

## §2 · Why the standard binding routes fail

The failures are not quantitative; they are structural.

1. **No scale separation at β=2.5.** The exact activity is `1.014`, the trivial-character fraction
   `0.17`, the Ursell ratio `4.6×10³`; the block-spin flow helps by `1.03×`. The missing factor is
   `40–55×`, which is the gap. A bare or single-layer polymer expansion cannot cross.
2. **The overlap has only upper bounds.** Reflection positivity gives log-convexity, which is
   satisfied at `w₁=0`; Perron–Frobenius is vacuum-anchored; the transfer matrix gives a decreasing
   upper bound. No lower bound survives.
3. **The TP₂ order is unbuilt.** The heat-kernel models give exact ratios, but the lower half at
   β=2.5 requires a binding-induced TP₂ order; RP supplies only the `e^{−tL}` form, the polymer
   expansion only the constant.
4. **The center-vortex route is neutral-blind.** `σ>0` lives in the charged channel; `tr F²` is
   N-ality zero and sees phase `1`. The 2D counterexample (`σ>0`, `f₁²=0`) shows `σ` is not enough.
   The missing object is a sign-cancelling positive regrouping — the non-abelian random-surface sign
   problem.

---

## §3 · The carrier correction

Q6's carrier hazard was half genuine and half artifact. The artifact: `Var∝1/V` is the
**mean-normalized** mode `(1/V)ΣO`, not the **total-normalized** mode `(1/√V)ΣO` that H1 uses. For the
latter,

```
   Var(q) = (1/V) Σ_{x,y} C(x−y) = Σ_z C(z) = χ_O = O(1) ,
```

independent of `V`. The carrier exists; the witness checks a free massive scalar in `d=1,2,4`
(`Var=1/m²`, slope `0`), the composite `:φ²:` (`χ=4.109/0.689`), and the one-plaquette model
(`Var=1/4`). The genuine residual is UV: the bare `tr F²` has `χ~a^{-4}` and must be renormalized at
a fixed physical scale. The ceiling direction remains: `χ` is a residue-weighted harmonic mean, so
`m₁²≤M₀/χ`, never the edge.

---

## §4 · What survives from the R wave

- The **exact anchor** P3 at β<1/48 and the flux-tube mechanism.
- The **isolation parent** statement `I⇒A`, now sharpened to the three-predicate picture.
- The **resolution bound** `sup|E−g|≥f₁²/2`; every atom-reader is circular.
- The **adversarial classification**: existence, residue floor, and isolation have different
  statuses; `H1` plus a fixed gap does not force an atom.
- The **2D correction**: `ρ²=8/3` is the holonomy's isolation, not the `tr F²` atom (which is 0).

---

## §5 · What the wave did not win, and the next object

- No proof of `A`, `W`, `I`, or `G` at β=2.5.
- No standard expansion reaches the binding criterion; the missing factor is `40–55×`.
- No reduction of Clay: the gap needs `W∧H1` or `I∧(ρ²≤C²)`, and both inputs are open.

**The next objects are a pair.** (i) A uniform lower bound `W: w₁≥w₀>0` on the atom fraction of the
minimal four-wall tube at β=2.5 (leakage `Σ_{n≥2}w_n≤1−w₀`). (ii) The isolation bound
`I` with the window `ρ²≤C²`. They are **independent** (U6); neither forces the other, and neither
alone forces the gap. No variance, positivity, spectral, Perron, RG, RMT, duality, center-vortex,
Wilson-loop, entanglement, or Harnack tool currently supplies either.

**One line.** *The carrier exists; binding is the atom's mechanism, not an axis; the predicates are
A, W, I, G, with `I⇒A`, `W⇒A`, `W⊥I`, and the gap via `W∧H1` or `I∧(ρ²≤C²)`; no standard expansion
reaches β=2.5, missing 40–55×. No proof, no QED, Clay 0%.*

*Refs:* `T1-POLYMER-RG.md`, `T2-TUBE-OVERLAP.md`, `T3-TP2-ISOLATION.md`, `T4-CARRIER.md`,
`T5-CENTER-VORTEX.md`, `T6-BINDING-ADVERSARIAL.md`; witnesses `w7/t1_*.py` … `w7/t6_*.py`;
`R-WAVE-SYNTHESIS.md`, `Q-WAVE-SYNTHESIS.md`.
