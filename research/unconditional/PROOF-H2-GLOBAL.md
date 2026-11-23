# PROOF-H2-GLOBAL — block-local ⇒ global by an inter-block contraction at β≈2.5: attempt and exact failure

**Verdict.** The global LSI at β≈2.5 does **not** close by an inter-block contraction, and no better norm, small-volume/boundary trick, or two-scale Stein decomposition closes it either: all three routes collapse onto one inequality — the **block-Dobrushin (Zegarlinski) contraction** on the *future-boundary mixture* of block conditionals, whose only rigorous bound is `C(β)≈15≫1` at β=2.5. It is **not** Savvidy (intra-block, absorbed by `C_HS=2βN_p(B)` in the proven block gap); it is the inter-block coupling. Two facts must stay separate: (i) even granting `C(β)<1`, one proves only the **global fiber LSI** (H2-global, a *UV* statement at fixed `U′`), and (ii) the **system** gap additionally needs the coarse **IR** Poincaré, which fails at the k→0 mode (`G″(0)≈+0.043−0.208=−0.165<0`). The contraction is thus both unprovable at β=2.5 *and* insufficient for the Clay object. No bound is fabricated.

## The exact missing inequality
Order the fiber into blocks `B_1,…,B_n`, `G_k=σ(ξ_{B_1..k})`. The proven block gap is LSI for the **direct** conditional `μ(ξ_{B_k}|ξ_{∂B_k},U′)`, constant `1/Λ_UV`, uniform in `ξ_∂`. But the Lu–Yau/Nachtergaele variance peeling uses the **past-only** conditional
`μ(·|G_{k-1})=∫ μ(·|ξ_{∂B_k})\,dP(ξ_{∂B_k\past}\,|\,G_{k-1})`,
a **mixture** over the integrated-out future boundary. LSI is *not* closed under mixtures (two `N(±a,1)` mixed have LSI `≈a²`, not `1`), so `1/Λ_UV` does **not** transfer. What is missing is mixture-stability of the block conditional:
```
Var_{μ(·|G_{k-1})}(f) ≤ (1/Λ_UV)/(1−C(β)) · E|∇_{B_k} f|²,  uniformly in k,U′,V,
  C(β) := sup_k Σ_{j>k, j∼k} δ_{kj},   δ_{kj} := sup_{ξ_past} osc_{ξ_{B_j}} log μ(ξ_{B_k}|ξ_{∂B_k}).
```
Zegarlinski's sufficient condition is `C(β)<1`. **This is the exact missing contraction inequality**; it is true as a reduction, and its hypothesis is what fails.

## Why it fails, and each route
**Magnitude.** `δ_{kj}≤2β·N_straddle(B_k,B_j)`: only straddling plaquettes move, each by `≤2β=5`. A shared face has `N_straddle≥1` and a block has up to 4 future neighbours, so `C(β)=sup_kΣδ_{kj}≳15≫1`. A **single** straddling plaquette already gives `δ≥2β=5>1` — the one-coupling term alone violates `C<1`.

- **(i) Better norm / renormalized lattice.** This `C` *is* the block-Dobrushin coefficient. χ²/entropy or Wasserstein norms change the per-plaquette influence to `O(β²)` or `(1/Λ_UV)β`, still `≥1` per shared plaquette at β=2.5. Coarsening lowers `C∝β_eff·N_straddle`, but one factor-2 step leaves `β_eff≈β` (g²≈1.6, slow running); `β_eff<1` is the flow to strong coupling = area law = **the gap** (IN-9). Larger blocks worsen it (`N_straddle∼b²`). Circular.
- **(ii) Small-volume + boundary.** The free-box `a`-uniform gap gives the single-block conditional uniform in *every* boundary — that *is* `Λ_UV`, already in hand, and it does **not** compose: the torus is a nested *mixture* of boundary-pinned blocks, whose constant multiplies by `(1−C)^{-1}` per scale — finite only if `C<1`, divergent with cycle length if `C>1`. Boundary-uniformity is necessary, not sufficient.
- **(iii) Two-scale Stein/FP.** REPAIR-2SCALE absorbs the score into the *global* fiber generator `L_⊥^{-1}` using `Λ_UV`, reducing the theorem to `H1` (coarse SLC) `+` `H2-global` (global fiber LSI). "Inter-block fibre coupling = UV high band" is correct physics, but proving its global LSI is exactly `H2-global=C(β)<1` — the same mixture wall; and even granting it, `m=min(λ/2,…)≤0` since `H1` fails at the soft mode. Stein relocates, but supplies neither piece.

**Sharpest line.** The block gap is UV; the system gap is UV⊗IR, and the block→global contraction is the wrong door to either — unprovable at β=2.5 (`C≈15`), and even if granted it yields only the high-band fiber LSI, leaving the k→0 coarse mass, which *is* the mass gap and is disfavoured (`G″(0)≈−0.165`), untouched. The last step is not in the coupling; it is the IR soft mode in disguise.

*Refs.* PROVE-H2-uvgap (`Λ_UV`); H2-globalization; H2-b-uniform; H2-TIEBREAKER; REPAIR-2SCALE; SCORE-HYPOTHESIS; W6-BE-at-2.5; W6-induced-mass; CONDITIONAL-LATTICE-GAP; VAR-FINAL; P4 §4–5 (L\*); Zegarlinski, CMP 1989; Guionnet–Zegarlinski; Lu–Yau, 1993; Nachtergaele, CMP **175** (1996); Chafaï (LSI–mixture); Shen–Zhu–Zhu, arXiv:2204.12737 (β<1/48); IN-9 (`SPEC-RES`).
