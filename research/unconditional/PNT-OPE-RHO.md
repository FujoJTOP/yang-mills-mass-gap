# PNT-OPE-RHO — bootstrap on `ρ(s)`: **no edge-sensitive constraint beyond finite moments reads the KL atom; moment-blindness is intrinsic (`[T]`)**

> **VERDICT (no QED).** Bounding the `0⁺⁺` atom directly from the spectral function `ρ(s)` — using the
> exact normalizations plus **any** of the candidate "edge-sensitive" data (UV/AF tail, threshold matching,
> renormalon-free UV condition, a second-operator inequality) — **cannot force an atom**. The reason is
> structural and proof-level: the support floor is a **topological (weakly-discontinuous)** property of the
> measure, while every candidate datum (moments, tail, threshold, operator bounds) is a **weakly-continuous
> integral** functional; the two live in different topologies, so no integral data pins the floor. The only
> functional that *reads* the atom is a **concentration** functional (non-measurable in QFT) or the **pole
> residue of the full Stieltjes transform `Π`** (= A3, read-not-force). **Moment-blindness is intrinsic;
> the exact missing step is the pole residue = `H1`.** `Clay = 0%`.

---

## §0 · Labels

**[T]** 真定理 · **[C]** 条件定理 · **[O]** 精确缺失步 · **[NG]** no-go for the named class · **[V-num]** numeric witness (schematic). Acid: `≤0`/`=0` on free massless **and discriminatingly** on `:φ²:` (`f₁²=0`).

---

## §1 · Setup and the four edge-sensitive candidates

`ρ(s) ≥ 0` on `s ≥ 0` (`s = m²`); KL atom = `f₁² δ(s − m₁²)`; gap ⇔ `f₁²>0` with `m₁² = inf supp ρ > 0`.
The renormalon-free constraints on `ρ` are:
- **moments** `M_j = ∫ s^j ρ ds`: `M₋₁ = χ = G″(0)⁻¹` (Ward), `M₀ = W` (norm), `M₋₂ = A` (trace-anomaly scale), and higher `M_j` (OPE coefficients);
- **UV/AF tail** `ρ(s) ~ s (ln s)^{−b}` as `s→∞` (`[T]` leading behavior; renormalon-free);
- **threshold** `ρ(s)=0` for `s<Δ_thr` (the bootstrap gap assumption);
- **second-operator inequality** `ρ(s) ≤ B(s)` (a cross-correlator bound).

`SLOT-BOOTSTRAP` already proved the **moments alone** are floor-blind (`Δ_min(K)=0 ∀K`). The question here:
does **any** of the other three read the atom?

---

## §2 · Theorem: intrinsic moment-blindness (why no integral datum can pin the floor)

**Statement `[T]`.** Let `𝒟` be any collection of constraints each of the form `∫ ρ φ_i = c_i` (moments,
tail integrals, any fixed test-functionals) together with pointwise `0 ≤ ρ ≤ B` (an upper spectral bound),
finitely many in number. Then for every `ε>0` there is `ρ_ε ≥ 0` with `inf supp ρ_ε < ε` satisfying all of
`𝒟` whenever the true `ρ` satisfies `𝒟` (i.e. the feasible set is floor-unbounded below).

**Proof idea `[T]`.** `supp ρ` is **not continuous in the weak topology**: `ρ ↦ inf supp ρ` is
**upper-semicontinuous**, not continuous, and is not a functional of the weak limit. Every constraint in
`𝒟` is `∫ ρ φ_i` — a **weakly-continuous** linear functional of `ρ`. A weakly-continuous functional cannot
detect a change of **support** that preserves all integrals (two measures with equal integrals can have
arbitrarily different support floors — a smoothing/shifting of mass within `[ε, ·]`). Concretely: to match
`M₋₂ = ∫s⁻²ρ` with floor `Δ→0`, concentrate a **continuous** bump near `Δ` where `s⁻² ~ Δ⁻²` is huge; the
required bump weight `w(Δ) = C·Δ² → 0` (this is `SLOT-BOOTSTRAP`'s `w₀*(Δ) ≈ CΔ²`), so it never conflicts
with the other moments or the tail. ⇒ the floor is a free direction orthogonal to `𝒟`.

**Consequence.** No number of moments, no UV tail, no threshold, no operator **upper** bound can certify
`Δ_min>0`. Forcing an atom requires a **lower** bound on a spectral functional or a **concentration**
constraint (a bound on `∫ φ(log ρ)`, an entropy/point-mass functional) — none is measurable in QFT.

---

## §3 · The four candidates, individually refuted

- **UV/AF tail `[NG]`.** The tail `ρ ~ s(ln s)^{−b}` is a statement at `s→∞`; the atom sits at `s→0`
  (the lightest glueball). The two **decouple**: `M₋₂` is dominated by the lowest `s` and receives
  negligible weight from the tail, while the tail is fixed by high `s`. Witness (Case 3): fixing the
  **weight and first moment above `s=20`**, a floor-`Δ`, atom-free measure matches all four moments **and**
  the tail to `≤2.4×10⁻¹³` for `Δ` down to `3×10⁻³`. The renormalon-free UV condition does **not** help.
- **Threshold `Δ_thr` `[NG]` (asymmetry).** Imposing `ρ=0` below `Δ_thr` and matching moments gives the
  **largest feasible** `Δ_thr = Δ_max` — an **upper** bound on the floor (witness: `Δ_max = 1.006` vs true
  floor `1.000`). It can **never** certify a **lower** bound `Δ_min>0`: feasibility at `Δ_thr` says only
  "the floor could be ≤ `Δ_thr`." This is the positivity/lsc asymmetry (`INF supp` is usc). The bootstrap
  threshold direction is **structurally** an upper-bound machine.
- **Renormalon-free UV condition `[NG]`.** Same as the tail: a UV datum. It constrains `s→∞`; it is silent
  on the small-`s` support floor. (It removes the `u=2` renormalon from the *data*, which is valuable — but
  it does not add IR edge-sensitivity.)
- **Second-operator inequality `ρ ≤ B` `[NG]`.** An **upper** bound on `ρ` can only **exclude** atoms
  (cap them), never **force** one. Forcing needs `ρ` to be **unable** to spread out — a lower bound or a
  total-mass/concentration condition — which is exactly `f₁²>0` = the gap.

---

## §4 · What WOULD read the atom (and why it is `H1`)

The atom weight is the **residue of the pole of the Stieltjes transform** `Π(Q²)=∫ρ(s)/(s+Q²)ds` at
`s = −m₁²`. The atom is visible **iff** one knows `Π` as an analytic function and locates its pole — i.e.
`f₁² = −Res_{s=−m₁²} Π`. That is **A3** (the Fredholm determinant / pole read), explicitly **read-not-force**:
it requires the **full** nonperturbative `Π`, which is the theory itself. The two non-moment functionals
that parametrize the atom —
- the **pole residue** of `Π` (`[T]`: only analytic continuation reads it; = `H1`),
- a **concentration functional** `∫φ(log ρ)` (`[T]`: non-measurable in QFT; no operator computes it),
— are the only candidates, and both are blocked. **There is no third, edge-sensitive, renormalon-free,
finite, measurable functional.**

---

## §5 · Witness

`w7/pnt_ope_rho_witness.py` (numpy+scipy, `py -3.12 -I`, `<1 s`). True measure: atom `f₁²=0.5` at `m₁²=1`
+ continuum `ρ~s` on `[4,300]`. Moments `[M₋₂,M₋₁,M₀,M₁]=[0.500,0.513,2.5,400.63]`.
- **Case 1 (moments only):** atom-free, nonneg competitor matches to `≤1.9×10⁻¹³` for floor
  `Δ = 1.0, 0.3, 0.1, 0.03, 0.01, 0.003`. Floor is free.
- **Case 3 (moments + fixed UV tail above s=20):** same, `≤2.4×10⁻¹³` for `Δ` down to `3×10⁻³`. Tail
  does not bite.
- **Case 2 (threshold):** `Δ_max = 1.006` (true floor `1.000`) — only an upper bound; the search for
  `Δ_min>0` is structurally empty.
- **Concentration:** the gapless match at floor `10⁻³` has residual `2×10⁻¹³`; the atom is gone.

(Case 2's `Δ_max` is *slightly* above the true floor purely from the finite continuum grid; it is an upper
bound as required, and tightens to `1.000` with a finer grid.)

---

## §6 · Honest bottom line

- **`[T]` Intrinsic moment-blindness.** No integral datum (moments, UV tail, threshold matching,
  renormalon-free UV condition, operator upper bounds) can lower-bound the KL support floor; `Δ_min = 0`
  in every such class. Forcing an atom needs a **concentration** functional or the **pole residue of `Π`**.
- **`[NG]`** the finite/integral-data class (extends `SLOT-BOOTSTRAP` to the UV tail and thresholds).
- **`[O]` exact missing step:** the pole residue `f₁² = −Res Π` — the full nonperturbative `Π` — which is
  `A3 = H1` verbatim. Running the bootstrap *directly on `ρ`* changes nothing: the obstruction is not the
  choice of object (`Π` vs `ρ` vs moments) but the **finite-data topologies**.
- **The one quantitatively new statement:** the UV/AF tail is a **large-`s`** datum and is provably
  **decoupled** from the small-`s` support floor — so "use the renormalon-free UV behavior" cannot build the
  strut (the renormalon-free data is IR-blind, exactly as the OPE data is UV-ambiguous).
- **Clay = 0%; no QED.**

*Refs:* `ROOF-OPE.md`, `SLOT-BOOTSTRAP.md`, `INVENT-SYNTHESIS.md`; A3 (`PROOF-F1SQ.md`), `PROFILE.md`.
Witness: `w7/pnt_ope_rho_witness.py`.
