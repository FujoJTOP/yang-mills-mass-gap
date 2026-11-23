# STATUS — handoff-proof (q₀ ≥ q₁)

> 2026-11-22 · 4D pure Yang–Mills · full report: `REPORT.md`
> **Step 0 confirmed by audit** (`research/handoff-audit/16-input-audit.md`): the 16 inputs do NOT contain the handoff. It is a **new input/estimate**, not a hidden one. The interface keeps `cert : ℝ → ℝ` (`gluing` §2.7 is now corrected to this shape); the constant-`gapBar` caveat is an interface rule, not an open issue.

## 0. Verified Lean artifacts (P1′)

**One-command check:** `bash lean/check-handoff.sh` — compiles all three against the pinned Mathlib cache; final line `ALL THREE COMPILE`. The coordinator ran it: exit 0. The script now also enforces the axiom audit: `HandoffMargin.lean` carries `#guard_msgs` blocks freezing `[propext, Classical.choice, Quot.sound]` for both `handoff_iff_margin` and `required_lambda`, so a leaking `sorry` or extra axiom fails the same command.

The three files it checks:

- `lean/YMStatement/Handoff.lean` — `handoff_gap` (with `cert : ℝ → ℝ`) and `mass_gap_of_certificate`.
- `lean/YMStatement/RGHandoff.lean` — one-loop flow, scheme covariance, `q` monotonicity, `gap_lower_bound`.
- `research/handoff-proof/HandoffMargin.lean` (**P1′**) — `handoff_iff_margin` and the proved inverse `required_lambda`. Axiom audit: only `propext`, `Classical.choice`, `Quot.sound`; no `sorryAx`.

(`research/lean-infra/YangMillsInfra/RG/Handoff.lean` is an earlier source copy, not part of `check-handoff.sh`.)

**P1′ only formalizes the pure-number equivalence `q₀≥q₁ ⇔ M≥0` and the inverse form. It does not prove `H3`; `CPB-W` remains the single missing estimate.**

## 0.5 Primary deliverable: RB2 (2-loop running bound)

**Goal:** `|β(λ) + Aλ² + Bλ³| ≤ C₃λ⁴` on `0<λ≤7.3` (or the one-sided sign `β ≤ −Aλ²−Bλ³`), with `A=1/C_RG`, `B=51/(121C_RG²)`, `C_RG=24π²/11`. Then `t_up ≤ t_up^{2L}` and `q₀ ≥ (m/Λ)exp(−t_up^{2L})`.

**Proved (algebra):** `t_up^{2L}(λ) = C_RG/λ − (51/121) ln(1 + 121 C_RG/(51λ))`. At `λ=6`: `2.639`, `q₀=0.571` (`m/Λ=8`); `q₀=1` at `λ≈7.28` (`g²≈2.43`).

**Exact missing nonperturbative input (`NPB`):** in a fixed scheme (block-spin effective coupling), a bound `|β(λ) − β_{≤2}(λ)| ≤ C₃λ⁴` (or the sign), with explicit `C₃`, uniform in UV cutoff and volume, for `0<λ≤7.3`. Not available: Balaban-type RG stops at `λ~1`; the 3-loop coefficient is perturbative and scheme-dependent; there is no rigorous nonperturbative beta function at these scales.

**Two falsification tests:**
1. A 3-loop/lattice `t_up` at `g²≈2` above `2.639` beyond the allowed `Δ`.
2. A zero or inflection of the exact `β` below `λ=7.3`.

`RB2` is the quantitative core of `CPB-W`; `CPB-W` remains the single missing estimate. No other routes.

## 0.6 One-step class preservation (CPB-W first block)

Delivered `CPBW-one-step.md` (U(1), fixed mass, fixed volume, first block; **C² norm on the ball**) and `CPBW-skeleton.md` (9 lemmas with (i) statement, (ii) source, (iii) whether an unsolved item is CPB-W/NPB itself).

- **Norm:** `C²` on the ball (`sup + Lip + Hessian`), not Lip alone. The conditional expectation preserves sup (no contraction) and contracts Lip (upper bound); the **Hessian lower bound** is needed for non-degeneracy. With an upper bound alone the Gaussian part can be eaten and `C_W` silently collapses to the empty set.
- **U(1) scope (must read):** U(1) is free, `F=dA` makes `dF=0` automatic; there is **no Bianchi cancellation**. The U(1) lemma tests only the Gaussian/blocking step (covariance splitting, fluctuation log, C² contraction/non-degeneracy). It does **not** test the signed-gas defect mechanism or the `e^{−cq²}` Bianchi scaling; it must not be read as “CPB-S verified”.
- **U(1), conditional success:** `‖V'−V∘A_∥‖_{C²} ≤ C₁ε²L(q)+C₂ε³L(q)²`, `L(q)=sup_xΣ_y|G_⊥(x,y)|² ≤ C₃ln(1/q)+C₄`; under `εL^{1/2}≤1/2` and the non-degeneracy `ε<c₀m²/C`, also `S_{k+1}''≥(c₀m²−Cε)Id>0`, `m_{k+1}=m(1+O(ε))`. Gives `δ≤C(g²/16π²)ln(1/q)+O(λ³)`. U(1)'s one-loop coefficient vanishes, so the actual `δ` is smaller.
- **Nonabelian, failure (named):** skeleton lemma 6, the **large-field C² bound** `‖V‖_{C²,R}≤ε` with `ε<c₀m²/C`, `R~g/q`. The commutator gives `‖V''‖≲R³/g²+R²/g`; at `q~1, g²~2, R~1.4` this is `≈2.8 > c₀`. Balaban's large-field bounds (CMP 122 (1989) 175–354; UV stability CMP 102 (1985) 255–275) assume `g` small (UV) and do not cover `q~1`. **Nonabelian one-step class preservation is unproven.** Lemma 6 is CPB-W's class-preservation core, not NPB.
- **Important distinction:** failure of `εL^{1/2}≤1/2` only says **this estimate (C²-on-ball) cannot prove class preservation**; it does **not** say class preservation is false. The large-field region may still be preserved; the bound is too crude. CPB-W is not refuted.
- **Gaussian-weighted class (alternative, `CPBW-one-step.md` §5b):** use `μ_G`-weighted local L², `C_W^G(ε)={S=S^G+V: sup_x∫|V_x|²dμ_G≤ε²}`. Then `εL^{1/2}` drops from `~2` to `O(0.2–1)`: central `0.24` (passes), conservative `0.96` (marginal), depending on `c_F=⟨F²⟩/g²` at `q~1`.
- **(a)/(b) judgment:** the needed interpolation is a **local Bernstein/cumulant bound** `|log E[e^{-V}]+E[V]−½Var(V)|≤C₃E|V−EV|³+C₄E|V−EV|⁴`. It uses **finite moments, not `‖V‖_∞`**, so it is (b), not (a).
- **Correction on `t_*`:** `M(t)=E[e^{t|V|}]` converges only for **`t<1`**, not `t~1/(c₄σ⁴)`. `V=S−S_G`'s large-field leading term is the negative **quadratic** `−S_G~−F²/(2g²)` (the quartic is only the first small-field correction), so `|V|~S_G` and `e^{(t−1)S_G−tS}` converges iff `t<1`. Use `t_*=1/2` and a **finite-order Taylor expansion with Lagrange remainder** (third cumulant finite on the compact group), not an infinite cumulant series. The radius 1 is independent of `c_F`; `c_F` enters only the constants (`c_F^{k/2}`), an O(1) factor at `c_F=1,2`. **New number:** the remainder raises `εL^{1/2}` from `0.24` to about `0.4–0.7` — **marginal**, no longer a comfortable 2× margin.
- **Verdict after correction:** (b) is **not closed**; `t_*<1` imposes no constraint on `c_F` (any `c_F` allows `t_*=1/2`). The wall remains a **local cumulant bound + computing `c_F`**, now with a finite-order remainder.
- **Local cumulant bound (one-page: `CPBW-cumulant.md`):** target is the third **connected** cumulant `|κ₃^ξ(V_x,V_y,V_z)|≤C_cum e^{−q(|x−y|+|y−z|+|z−x|)/2}`, with the local density `c₃=sup_xΣ_{y,z}|κ₃|` uniform in volume and UV cutoff. Proof route: Mayer/cluster expansion; the input is a **local activity bound** (compact group, `S∈[0,2/g²]`), **not the global `sup|V|`**, so logically (b) survives. But at `g²~2` the standard Kotecký–Preiss bound typically gives `C_cum≥2`.
- **Pass/fail (`εL^{1/2}=0.24+C_cum R₃`, `R₃=0.16–0.46`):** `C_cum=1` → `0.40–0.70` (marginal; only the low end passes); `C_cum=2` → `0.56–1.16` (fails); `C_cum=3` → `0.72–1.62` (fails). The pass condition is `C_cum≲0.6–1.6`. **If the cluster expansion only gives `C_cum≥2`, the Gaussian-weighted route is numerically dead.**
- **Bianchi/sign-cancellation lemma (one-page: `CPBW-bianchi-lemma.md`):** the `q²` form `|w(C)|≤e^{−c₀q|C|−c₁q²|C|}` is **neutral at `q≈1`** (`q²≈q`), so it cannot move `C_cum` from 2 to 1.6 or `z(1)` to 1/2. The only hope is the **signed cancellation factor** `κ_c<1` (fixed-boundary sum of signed closed surfaces; `κ_c` is distinct from the KP constant `κ=1/2`): `|Σ_{∂C=∂C₀}w(C)|≤e^{−c₀q|C₀|}κ_c^{|C₀|}`. CPB-W needs `κ_c≤0.8`; CPB-S needs stronger; the same `κ_c` serves both. **Judgment: `q²` alone is insufficient; the signed cancellation is untested. If `κ_c≈1`, both CPB-W and CPB-S close.** Falsification: measure `z_signed/z_abs` per block at `q≈1` (not just `S_def/|C|`).
- **The wall moves from “C²-on-ball” to “local cumulant bound + computing `c_F`”; it does not disappear.**
- **RB2/NPB untouched.** CPB-W is not claimed.

## 1. Routes tried and ruled out

| Route | Tried? | Why ruled out (as a standalone proof of `H3`) | Retained as |
|---|---|---|---|
| **(a) Balaban multiscale RG** | yes (framework) | Expansion parameter is the coupling; by construction it stops at `λ~1`; it controls UV but does **not** generate the mass or track `m/Λ`. At `λ~10` the large-field regions are not suppressed. | The class-preservation language for **`CPB-W`** (weak reach). |
| **(b) 2026 sigma-model template** | yes (method) | The all-coupling signed gas is a theorem only in 2D O(N); the 4D YM port (gauge constraint, Bianchi instead of an O(N) hard constraint) is unproven. | The signed-gas / effective-density method for **`CPB-S`** (strong reach). |
| **(c) CNS strong-coupling expansion** | evaluated | It can raise the strong-side reach `q₁`, but does not touch the weak side; alone it cannot establish the overlap at the physical boundary `g²≈2`. | A concrete strong-side increment; target `β≈2.6–3.3` for `q₁≈1`. |
| **(d) Stochastic quantization / Bakry–Émery** | evaluated | At weak coupling the Wilson Hessian has negative directions from the commutator term (the same anti-screening that makes the theory asymptotically free); a positive curvature at the crossover is essentially equivalent to the mass gap. | A **diagnostic** (`κ_eff`, lowest Hessian eigenvalue) for `CPB-S`. |
| **(e) Small-`h^∨` priority** | evaluated | A reduction/constant table, not a proof. | Group checks: SU(2), SO(3), SU(3), G₂. |

**Invented object:** the **crossover activity** `z(q) := sup_x Σ_{C∋x}|w_q(C)| e^{κq|C|}` (signed-gas polymer activity with the mass gap as the small parameter) and the **two-sided reach estimates** `CPB-W` / `CPB-S`.

## 2. Current best statement of `q₀ ≥ q₁`

With `B := 11h^∨/(24π²)`, `1/g²(μ) = B ln(μ/Λ)`, `λ* := h^∨g*²`, and `m` the **physical** gap (`q(μ)=m/μ` exactly):

```
q₀ = (m/Λ) exp(−24π²/(11λ*))
H3: q₀ ≥ q₁  ⟺  M := ln(m/Λ) − 24π²/(11λ*) − ln q₁ ≥ 0
    ⟺  t* := 24π²/(11λ*) ≤ t_S := ln(m/Λ) + ln(1/q₁)
```

**Explicit numbers (one-loop, `m/Λ=8`):**

| `λ*` | `q₀` | verdict for `q₁=1` |
|---:|---:|---|
| 3 (`g²=1`, SU(3)) | 6.1e−3 | **fails by ~2 orders** |
| 6 (`g²=2`) | 0.22 | fails; needs `q₁≤0.22` |
| 9 (`g²=3`) | 0.73 | fails; needs `q₁≤0.73` |
| 10.4 (`g²=3.45`, `α_s≈0.275`) | 1.00 | passes |

**Required reach (one-loop):** `q₁=1 → λ*≈10.4`; `q₁=0.5 → 7.77`; `q₁=0.2 → 5.84`. **Group:** SU(2) needs `g²` 1.5× SU(3) at the same `λ*`.

**Physical calibration (`numerics/REPORT.md`; NOT a proof):** `g²=1` is not the handoff — it sits at `μ≈158Λ`. The physical handoff is `g²≈1.8–2.3` (`μ=Δ=(4–8)Λ`), where the 3-loop/running-mass model gives `t_up≈2.0` and `q₀≈0.4–1` for `m/Λ=4–8` — the same order as `q₁~1`, **not** a 20× gap. Caveats: `m/Λ=(4–8)` is a literature input with circularity; the lattice `β≈6` overlap is O(1) evidence (30–50% scaling violation); the numerical report does not prove the overlap.

**One-loop vs 3-loop matters at `g²≈2`:** one-loop `t_up=3.59` gives `q₀=0.22`; 3-loop/numerical `t_up≈2.0` gives `q₀≈1.1`. So `CPB-W` must also give a **≥2-loop bound on `t_up`**.

**2-loop bound (RB2):** see §0.5 for the statement and the named gap `NPB`. Values: `t_up^{2L}(6)=2.639`, `q₀=0.571` (`m/Λ=8`); `q₀=1` at `λ≈7.28` (`g²≈2.43`); `t_up^{2L}(9)=1.592`, `q₀=1.628`. At `λ=6`, `q₁=1` needs `m/Λ≥14.0`. This is a genuine 2-loop value, not an inversion of a 3-loop number.

**Normalization (not an error):** this report uses `T(f)=1/2`: `b₀=11h^∨/(48π²)`, `B=2b₀=11h^∨/(24π²)`, `t_up=24π²/(11λ*)`; SU(3) `g²=1` gives `7.178`, matching numerics. `group-uniformity`'s `12π²` is `T(adj)=2h^∨` with `g²` half; after conversion the formulas agree. **Mixing** them flips `q₀` by `e^{3.59}≈36` at `g²=1`, so always state the normalization.

**Conditional reduction (Theorem 2):** if `CPB-W` gives weak reach `λ_c` (and `t_up(λ_c)`) and `CPB-S` gives strong reach `q_c`, then
```
H3  ⟺  m/Λ ≥ q_c · exp(t_up(λ_c)).
```
One-loop `t_up=24π²/(11λ_c)`: ≈6.0 for `q_c=1, λ_c=12`. Numerical `t_up≈2.0` at `g²≈1.8`: ≈7.4 for `q_c=1`, ≈3.7 for `q_c=0.5`. Neither estimate alone gives `H3`; they leave an explicit mass bound.

## 3. Single missing estimate

**`CPB-W` (bottleneck):** the weak-side (massive background-field) class is preserved up to
```
λ_c ≈ 7.3   (2-loop; SU(3) g²≈2.43, α_s≈0.19;  3-loop/numerical ≈6;  one-loop-conservative ≈12),
relative error ≤ 1/2,  PLUS the 2-loop bound RB2 on t_up.
```
Its small parameter is `(g²/(16π²))·ln(1/q) ≈ 0.01–0.03`, **not** `g²≈2`; the mass gap supplies the IR cutoff. This is the one estimate that decides the required `m/Λ`.

**Conversion warning:** do **not** invert the one-loop formula with a 3-loop `t_up`. One-loop `t_up=24π²/(11λ)` at `t=2.0` gives `λ≈10.76`; a result of 21 is the `48π²` slip. The 3-loop `t(λ)` is a **different function** (≈2.0 at `λ=6`, smaller at `λ=12`). The requirement is a **2-loop bound at `λ≈6`**, not a larger `λ_c`; keep `λ_c≈12` as a one-loop statement only.

Secondary (new object): **`CPB-S`** — the signed-gas activity satisfies `z(q_c)≤1/2` at `q_c≈1` (equivalently, effective Hessian positive on the gauge-invariant subspace).

**`CPB-S` convention (corrected):** `z(q)=sup_x Σ_{C∋x}|w_q(C)|e^{κq|C|}`, `κ=1/2`; assume `|w_q(C)|≤e^{−c₀q|C|}` with `c₀>κ`. The combined per-block weight is `e^{−(c₀−κ)q|C|}` (0.61 at `q=1`, 0.37 at `q=2` for `c₀=1`). Then `z` is non-increasing in `q`, worst at `q_c`; check once. Earlier the report wrote `|w(C)|` without the `q`-dependence, which made `z` look increasing; the `+` is the standard KP factor, the error was omitting the `q`-dependence and `c₀>κ`.

**Bianchi scaling (Proposition B, conditional):** if the defects are closed 2-surfaces (`dF=0`) with area-extensive action and tension `σ(q)~q²`, then `z(q)≤Ce^{−cq²}`; if only a domain wall (tension `~q`), `z~1/q` and the `q²` claim fails. The missing input is the coefficient `c₀` in `σ≥c₀q²`. Consequences (2-loop, `m/Λ=8`): `c₀=1, 0.5, 0.25` → `q_c=0.83, 1.18, 1.67` → `λ*=6.80, 7.78, 9.10`; the `q`-linear branch (`q_c≥2`) gives `λ*=10.02`.

## 4. What I need next

1. **A computation (highest value):** SU(3) lattice at **`β≈2.6–3.3` (`g²≈1.8–2.3`, physical handoff)**; conservative cross-check at `β≈2.2` (`g²≈2.7`) and `β≈1.7` (`g²≈3.5`). Measure `ξ/a`, connected-4-point ratio `R₄` of block variables, defect density `ρ`, and the lowest Hessian eigenvalue `κ_eff`. Pass: `z`-proxy small / `κ_eff>0`; fail: `R₄≫1` or `κ_eff<0`. **`ξ/a` is numerical evidence, not a certificate**: the certificate is the rigorous inequality `q₁^σ≤0.143`; evidence needs a 2σ lower bound and `L/ξ≥4` (L=12 saturates at `ξ/a≈3–4`).
2. **A number from the literature:** the exact `β` threshold in Cao–Nissim–Sheffield (does it reach `β≈2.6–3.3`, i.e. `q₁≈1`?), and an independent lattice `m/Λ` (glueball mass in units of `Λ`, not the circular `Δ/Λ∈[4,8]` input).
3. **A definition:** the signed-gas norm `z(q)` for the **equivariant** block spin (face-average with parallel transport, `gauge-invariance` session) and the exact weak class `C_W(ε₀)` at `λ_c`.
4. **An audit:** confirm the `24π²` convention with the main session; confirm `gluing` §2.7 is changed to `cert : ℝ → ℝ`.

## 5. Blockers / confidence

- **Blockers:** `CPB-W` and `CPB-S` are conjectures; no 4D constructive YM exists at `α_s~0.3`; `m/Λ` has no independent rigorous lower bound.
- **High confidence:** the equivalence `H3 ⇔ M≥0` (Lean-verified, P1′: compiles exit 0, axiom-clean); the normalization conventions; the number table; `cert : ℝ → ℝ` necessity.
- **Medium:** Step 0 (audit-level, not Lean-rebuilt).
- **Low:** `CPB-W`/`CPB-S` truth; `H3` under standard constants.
- **What would change the verdict:** `m/Λ≥20` (then `q₁=1` needs `λ*≈7.4`); strong coupling to `q₁≤0.2` (needs `λ*≈5.8`); a proof of `CPB-W`; the Bianchi `e^{−cq²}` bound lowering `q_c` to ~1.2 (c₀≈0.5); a stronger `RG1–RG5` found in a Lean-level recheck.
