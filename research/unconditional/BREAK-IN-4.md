# BREAK-IN-4 — can Proposition B be a theorem in 4D?

**Verdict: no.** Prop B (`(H1)(H2)(H3) ⇒ z(q)≤Ce^{−cq²}`) stays a *conditional* lemma, and both gaps close as **negatives**: (i) `H1⟹H2` is **false** — an explicit counterexample, not a missing estimate; (ii) even granting H2, **BL-5** makes the `q²` gain identically zero at the crossover. Prop B is therefore vacuous exactly where IN-4 needs it; neither a proof of `H1⟹H2` nor a bound on `c₀=T₂/m²` is reachable.

**1. `H1⟹H2` is refuted, not open.**
- **H1 is not even a theorem.** BL-1: `∏_{p∈∂c}U_p=1` is the *flux* cocycle `dF=0`; in non-abelian YM it constrains the flux sum, not the *support* (a large isolated plaquette is an open-support defect). H1 is an **assumption**, a theorem only abelian-quantized (BL-1′).
- **Granting H1, H2 still fails.** BL-3: the slice `x₄=0` (all spatial plaquettes) is a closed 2-chain (`∂S=0`) supported on a **codim-1 slab**, `|S|=3L³`, with per-plaquette cost `β(1−Re tr U_p/N)=O(1/g²)` **independent of q**. This legitimate gauge configuration auto-satisfies `dF=0`; yet "closed 2-chain" ⇏ "codim-2 membrane" ⇏ "tension `σ(q)=c₀q²`."
- **No salvageable fragment.** No gauge-invariant functional measures support codimension; any bridge ("surface defects energy-dominate," "energy selects d=2") is either equivalent to H2 (circular) or needs structure pure YM lacks. BL-4: `dF=0` fixes entropy (`γ_all→γ_cyc`, q-independent) and leaves `σ(q)` untouched.

**2. `c₀=T₂/m²` has no rigorous bound.** `T₂` (membrane tension) and `m` are both nonperturbative dynamics; `m>0` **is** the Clay gap. The only rigorous window (strong coupling, `β<β*≪2.5`) has `m~e^{−c/β}`, disjoint from the crossover. `c₀` is a *ratio*, not universal. Numerically `c₀` is a Gaussian-tail artifact (blocked field Gaussian; `c₀=2.4–10` across thresholds — `lattice-cpb` §9.4). The only provable inequality, `c₀>0 ⟸ (T₂>0 ∧ m<∞)`, is vacuous.

**3. BL-5 kills the mechanism at `q≈1`.** Theorem: `q^α|_{q=1}=1 ∀α`; at `q_c≈0.835`, `q²−q≈−0.14`, negligible against `κq=0.42` and `γ_ent≳1`. The operative criterion is the **absolute** `σ(1)=c₀`: single block `c₀≥0.5+ln2≈1.19`; series `c₀≥0.5+γ_ent+ln3≈2.6`. The exponent is irrelevant at the point of use.

**4. What closes.** BL-1, BL-3, BL-4, BL-5 are rigorous, and conditional Prop B holds. IN-4 stays → G2. The 4D port is *additionally* blocked by NP-7 (the center vortex carries charged-sector charge; the gap is neutral).

**EXACT obstruction.** `H1⟹H2` fails at BL-3: a closed 2-chain of codim-1 support with q-independent cost is in the ensemble, and nothing gauge-invariant selects codim-2. `c₀` is a ratio whose denominator `m` is the target, controlled only in a strong-coupling window disjoint from the crossover. Even if both fell, BL-5 makes the `q²` improvement exactly zero at `q=1` (and ≈0.14 at `q_c≈0.84`). Prop B cannot become a 4D theorem; the true requirement is the dynamical absolute `σ(1)≳1.2–2.6`, not the exponent.
