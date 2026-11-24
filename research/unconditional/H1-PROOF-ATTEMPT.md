# H1-PROOF-ATTEMPT — the strongest rigorous chain to H1 we can actually prove now, the single open cell, and the exact missing lemma

> 2026-11-23 · Direct proof attempt at H1 (not a new species). H1: `Var_μ(a_{k→0}) ≤ 1/λ`, `(L,a)`-uniform; equivalently `G″(k→0)=E[H]−Var(∇S|0)+W_geo ≥ λ` ⟺ `m_eff²(k→0)≥λ`.
> **判词 (verdict, up front).** The (field-strength × momentum) plane is closed by banked **真定理** everywhere **except one cell**, and that cell is **H1 verbatim** with zero slack. Restated as a region table: **6 of 9 cells are closed unconditionally at the physical coupling** (low-field whole row, high-field whole row), a 7th is closed only in the `g→0` limit (AF-UV tower, **inactive at β=2.5, g≈1.5**), and the **typical-field column — its high-momentum fiber negatives and its `k→0` target mode are the *same number* by Ward/Schur — is the one open object = H1.** Hence **region-wise ≈ 90% of the plane is rigorously signed, but 0% of H1's actual assertion is closed**, because the target coarse mode (`k→0` at typical field) sits *inside* the open cell and none of the banked closures contain it. The exact single missing lemma is the conditioned-force-variance bound `(◆)` below; it is H1 renamed, not a weaker input. **Scope caveat (REDUCTION-AUDIT): H1 is a susceptibility bound, necessary but not sufficient for the Clay gap — the gap additionally needs the residue/ordering datum `H*` (`I1∧I2`). So "single missing lemma for H1" ≠ "single missing lemma for Clay."** Clay stays **0%**.

---

## §1 · The rigorous chain we CAN prove now (banked theorems only)

Fix compact simple `SU(N)`, `N≥2`; Wilson action on `T⁴_L` with spacing `a`; `μ = Z⁻¹ e^{−S_W}dU`, Haar per link; flat chart `X=log U ∈ su(N)`; `V := S_W + Φ`, `Φ := −log J` (Haar Jacobian in exponential coords); mode `a = ⟨v,X⟩`, `f = ∂_n V = ⟨v,∇V⟩`, `h = V_vv`, conditioned on `a=0`.

### §1.0 The exact identity layer (真定理, all β, nonperturbative)

Differentiating the marginal `p(a)=e^{−V_eff(a)}` twice gives, exactly,
```
   G″(0)[v] = E_μ[V_vv | a=0] − Var_μ(⟨v,∇V⟩ | a=0) = V_eff″(0)[v] = m_eff²(k→0)[v].      (G2.0)
```
With `Δ_dyn = E[H_vv]−Var(G_v)`, `W_geo = E[Φ″]−Var(Φ′)`, `X = −2Cov(G_v,Φ′)` (all conditioned), this splits as `G″ = Δ_dyn + W_geo + X` (PROOF-JACOBIAN-SAVES §1). Three further exact statements, all proved for *every* `β`, no convexity, no positivity:

- **Ward (ABC-NEW-INEQ §1):** `E_μ[f·a]=1`; if `p(a)` is Gaussian then `G″·Var_μ(a)=1`, and conversely. Off-Gaussian the identity is `G″·Var(a)=1` up to the non-Gaussianity slack.
- **Schwinger–Dyson / Stein (ABC-NEW-INEQ, FORGE5-CONE-VARIANCE, SUB1):** unconditioned, `E_μ[H_n]=E_μ[f²]`, hence `Var_μ(f) ≤ E_μ[H_n]` — but this is the **unconditioned** variance (equipartition at `E f=0`); it is *gap-blind* (free massless saturates it with `Var(a)=∞`).
- **SCBL identity (INVENT-VARIANCE-INEQ §1):** `G″(0)=S(v)−D(v)` with `S := E_μ[h_n−⟨K⁻¹c,c⟩|0]`, `D := Var_μ(f|0)−E_μ[⟨K⁻¹c,c⟩|0]`, `K=Π₀ Hess V Π₀`, `c=Π₀ Hess V v`. Then `D≤0 ⇔ K⪰0`, `D=0 ⇔` Gaussian, `D>0 ⇔` fiber Savvidy. `G″>0 ⟺ D(v)<S(v)`.

These are the exact dictionaries. None of them supplies a sign by itself.

### §1.1 The region table (field-strength `φ := 2g|F|` × momentum `k`)

Two thresholds: the **shallow-band top** `c_Φ=N/12` (SU2 `0.167`, SU3 `0.25`; Lemma B) and the **Savvidy cone gate** `μ²` (where `−D²+2ad F_B` develops negatives). "CLOSED" = a banked theorem yields the sign with the stated quantifier.

| Region | Status | Theorem / source | Constant · quantifier | Verify / fine print |
|---|---|---|---|---|
| **Low field** `φ ≤ N/12`, **any `k` incl `k→0`** | ✅ **CLOSED / 真定理** | `Φ″⪰N/12` ⟹ `K^V = K^S+Φ″ ⪰ −2g\|F\|+N/12 ≥ 0` ⟹ fibre Brascamp–Lieb fires ⟹ `D≤0` ⟹ `G″≥S` | all `L`, all `k`; explicit `G″≥S≥N/12` at `β=0` (`χ≤12/N`) | ABC-CONSTRUCTIVE §2; ABC2-SYNTHESIS §2; D3 §2. *Fine print:* the corpus banks `G″≥S>0`; the explicit uniform `λ` on the whole band for `β>0` is inherited from the `β=0` floor `N/12`, not independently re-optimised. |
| **High-field cone** `φ > μ²`, **any `k` incl `k→0`** | ✅ **CLOSED / 真定理** | WSTC: the cone's variance *mass* is exponentially starved | `Var_{μ₀}(f·1_{C_κ}) ≤ (cβκ_max)² e^{−cβ(κ²−κ_typ²)N_□} → 0`, **`L`-uniform** | ABC-NEW-INEQ §2. *Verify:* right sign; the large field `f~β‖F‖` is multiplied by its own `e^{−β‖F‖²V}` weight, so the variance mass is elsewhere. |
| **AF-UV tower** `g_j→0` | ✅ **CLOSED / 真定理** | AF value-tail; `Σ_j g_j⁴ = ζ(4)/(16β₀⁴)<∞` | destabilising `≤ 2g_j\|F\| e^{−cμ_j⁴/g_j⁴} ≤ μ_j² e^{−cμ_j⁴/g_j⁴}→0` | ABC-CONSTRUCTIVE §2. *Verify:* genuine for `g_j→0`, **inactive at β=2.5** (`g≈1.5`); closes exactly the shells the continuum removes. |
| **Typical field** `N/12 < φ ≲ μ²`, **`k→0` target** | ❌ **OPEN = H1** | — | `PᵀHP` min `≈ −0.3` on the block-uniform mode | ABC2-SYNTHESIS §3; MARGINAL; OTHER-H1-POINCARE |
| **Typical field**, **high-momentum fibre negatives** | ❌ **OPEN = H1 (same object)** | — | SU2 `21.7%` / SU3 `42%` directions negative, `−1.12β` (`−g\|B\|`), 90–97% in the fluctuation block | W6-BE-at-2.5; PROOF-SAVVIDY-EXACT-SIGN |
| **Typical field**, **AF-UV limit strip** | ⚠ **条件 (limit only)** | AF at `g_j→0` | not active at `g≈1.5` | ABC2-SYNTHESIS §2 |

Two whole rows (low field, high field) are closed **including the `k→0` mode**, by named theorems, uniformly in `L`. The third (typical-field) column is **not** closed at any momentum; its two faces are tied by Schur+Ward into a single number:

> `PᵀHP` min `≈−0.3` at `k→0` (the target mode) and the `21.7–42%` high-momentum negatives at `k~k_max` are **one object**: the high-momentum fibre negatives integrate (via the Schur complement `K⁻¹`) into `Var(G|0)`, which by the exact identity `G″(0)=E[H|0]−Var(G|0)+W` fixes `G″(0)`. **One cell, not two.** (ABC2-SYNTHESIS §3)

### §1.2 The global theorems (not cells — they hold everywhere)

- **`W_geo ≥ 0`** (Prékopa, PROOF-JACOBIAN-SAVES §2): `Φ` convex (`g(d)=−2log(sin(d/2)/(d/2))`, `g″≥1/6` on `|d|<2π`) ⟹ `e^{−Φ}dX` log-concave ⟹ every linear marginal `a` log-concave ⟹ `W_geo=E[Φ″]−Var(Φ′)≥0` for **every mode incl `k→0`, every `L`.**
- **`E[Φ″] ≥ N/12`** (Lemma B, PROOF-JACOBIAN-SAVES §2): pointwise, under *any* measure, any `β`, any `k`.
- **β=0 / strong-coupling anchors:** at `β=0`, `V_eff″ = W_geo ≥ N/12` pointwise ⟹ `χ=Var(a) ≤ 12/N`, all `L`, all modes (**真定理**, D3 §2; measured SU2 `0.248` vs `E[Φ″]−Var(Φ′)=0.246`). At `β<1/48`, **Lemma SC** pushes SZZ's LSI through the `O(1)`-Lipschitz block map ⟹ coarse Poincaré with `λ=C_blk⁻²λ_μ` (**真定理**, PROOF-COARSE-POINCARE §(c); SZZ, CMP **400** (2023) 805, arXiv:2204.12737). P3 gives the residue `f₁²≥r₀>0` unconditionally at `β<1/48` (**真定理**, PROOF-P3-UNCONDITIONAL).

### §1.3 What the chain composes to

`(G2.0) + (§1.2 global) + (low-field row) + (high-field row) + (AF-UV tower)` close **every cell of the plane except the typical-field column**. On that column none of the four handles is nonzero: `2g|F| > N/12` kills the Brascamp–Lieb/shallow argument (§1.1 row 1); `g≈1.5` kills the AF argument (§1.1 row 3); `κ≈κ_typ` sends the WSTC exponent `e^{−cβ(κ²−κ_typ²)N_□}→1` to `O(1)` (§1.1 row 2); and the `N/12` floor is `~10×` short of `2g|F|` (SUB1). The unconditioned SD identity is exactly *equipartition* there (`E[H]−Var≈0`), so it carries no sign; and the universal Cauchy–Schwarz bound on the cross term, `|X| ≤ 2√(Var G·Var Φ′) ≈ 0.40`, **exceeds** `W_geo ≈ 0.32`, giving `G″ ≳ −0.08 < 0` — short by `O(λ)`, *one variance bound short* (ABC2-CONSTRUCTIVE §4). The conditioning defect `G″(0)|_{open} = [E(V″|0)−E V″] − [Var(V′|0)−Var V′]` **is** H1.

---

## §2 · The single open cell, stated precisely

> **OPEN CELL (the one missing region).** *Typical field × the fibre negative modes that feed the `k→0` coarse mode.* Concretely: `μ`-typical configurations with `N/12 < 2g|F| ≲ μ²` (outside the shallow band, at/below the Savvidy cone gate), carrying the high-momentum Savvidy negative spectrum (`λ_min=−g|B|`, `90–97%` of negatives in the fluctuation block, `21–42%` of directions), whose Schur integration into the fibre variance `Var(G|0)` determines the `k→0` target curvature `G″(0)` via `G″(0)=E[H|0]−Var(G|0)+W`.

Why it is **not empty** (ACid-verified, no vacuity): W6-BE directly measures negative measure there at `β=2.5`; the cell carries the **`PᵀHP` min `≈−0.3`** of the block-uniform target mode itself. Why it is **not larger**: every other cell is closed by §1.1's theorems. Why it is **exactly H1, zero slack**: the three exact dictionaries all collapse onto it —

```
   G″(0)|_{open} ≥ λ
     ⟺ [E(V″|0) − Var(V′|0)] ≥ λ                                   (exact identity G2.0)
     ⟺ D(v) < S(v) − λ                                             (SCBL, INVENT §1)
     ⟺ Var_μ(a_{k→0}) ≤ 1/λ     (Ward; iff marginal Gaussian)       (= H1)
     ⟺ Var_ν(φ) ≤ (1/λ)∫|∇φ|²dν  (coarse Poincaré, OTHER-H1 §1)
```

No step in this chain loses a factor: high-momentum `Var(G|0)` = lower bound `G″(0)` = upper bound `Var(a)` = coarse Poincaré = `m_eff²(0)≥λ` (ABC2-SYNTHESIS §4). The cell **is** the gap in measure language; it is not a residue of an unfinished decomposition that a fifth tool could still shrink.

---

## §3 · The sharpest finish-lemma, with exact hypotheses

Two logically equivalent statements; the first is the exact target, the second is the sharpest *localised sufficient form* (obstruction concentrated in one measurable scalar).

### §3.1 The finish lemma `(◆)` — exact form

> **FINISH LEMMA `(◆)`.** There exists `λ = λ(N,β_g)>0`, **independent of `(L,a,U′)`**, such that for `μ`-a.e. fibre configuration `U′`, and for every coarse transverse mode `v ∈ ran P` with `k` below the lowest non-zero lattice momentum (`k→0`) and `‖v‖=1`:
> ```
>    Var_μ(⟨v,∇V⟩ | a=0)  ≤  E_μ[V_vv | a=0] − λ.                       (◆)
> ```
> By `(G2.0)` this is **identically** `G″(0)[v] ≥ λ`; by Ward it is **identically** `Var_μ(a_{k→0}) ≤ 1/λ` (with equality of the two forms iff the coarse marginal is Gaussian). **This is H1.**

Note `V = S_W + Φ`, so `⟨v,∇V⟩ = G_v + Φ′` and `V_vv = H_vv + Φ″`: `(◆)` correctly carries the *full* potential, so no sign-uncontrolled cross term `X` is left over. (The frequently quoted `Var_μ(⟨v,∇S⟩|0) ≤ W_geo` is *not* a rigorous standalone sufficient condition, because it drops the sign-unknown `X=−2Cov(G,Φ′)`; measured `X≈+0.22` helps, but its sign is not proved.)

### §3.2 The sharpest localised form — one scalar competition

With `K=Π₀ Hess V Π₀` (fibre Hessian), `c=Π₀ Hess V v`, `S(v)=E_μ[V_vv−⟨K⁻¹c,c⟩|0]`, `D(v)=Var_μ(∂_nV|0)−E_μ[⟨K⁻¹c,c⟩|0]`, and `G″=S−D` exactly, the finish lemma is the **single-scalar competition**
> **FINISH LEMMA `(◆′)`.** There is `ε>0` (and a curvature floor `σ>0`), uniform on the open cell, such that
> ```
>    κ(v) := D(v)/S(v)  ≤  1 − ε,      S(v) ≥ σ > 0.
> ```
> Then `G″(0) = S(1−κ) ≥ εS ≥ εσ > 0`. Measured `κ = D/S ≈ 2%` at `β=2.5` (`G″≈1/Var(a)` to `2–5%`: PROVE-Gpp `0.274` vs `0.268`; `0.485` vs `0.45`; INVENT §5 `D/S≈2%`); at `β=0`, `S=S^Φ≥N/12`, `D≤0`, so `κ≤0` and `(◆′)` is a theorem with `σ=N/12` (D3).

`(◆′)` isolates the obstruction: `D>0` **is** the indefinite fibre (`K` Savvidy-negative, `90–97%` of negatives in the fibre), and `D≤(1−ε)S` is the exact statement that the fibre Brascamp–Lieb defect stays a bounded fraction of the Gaussian–Schur curvature. It does **not** evade the wall: `D` is the conditioned force variance, so bounding it **is** H1.

### §3.3 Exact hypotheses (for either form)

1. `SU(N)` compact simple, `N≥2`; Wilson action `S_W=βΣ_□(1−N⁻¹Re tr U_□)` on `T⁴_L`, spacing `a`, volume `V=L⁴a⁴`; product Haar `dU`.
2. **Coupling at the physical crossover**, `g=O(1)` (SU2 `β≈2.5`, `g≈1.3–1.5`; SU3 `β≈5.7`). **Crucially not** `g→0` (AF-UV inactive) and **not** `β<1/48` (strong-coupling inactive).
3. **Typical field:** `N/12 < 2g|F| ≲ μ²` — outside the shallow band, at/below the Savvidy cone gate (the measure-dominant cell).
4. **Mode:** coarse block-averaged `v∈ran P`, transverse to the gauge orbit, `‖v‖=1`, `k→0`; conditioning on the symmetric slice `a=⟨v,X⟩=0`.
5. **Uniformity:** `λ` independent of `L` (volume), `a` (spacing / continuum direction), and `U′` (fibre conditioning).
6. `K=Π₀ Hess V Π₀` invertible `μ`-a.e. (fails only on the measure-zero `det K=0` set — a transverse massless mode, the fibre edge).
7. **NOT assumed:** log-concavity of `S_W`; any FKG/positive-association inequality; small `β`; the gap itself; any residue/atom (`f₁²`) input.

The hypothesis is **exactly** H1 — not a strictly weaker sufficient input. Every banked route (SCBL, WSTC, AF, Prékopa, coarse-Poincaré, induced-mass, stochastic-localization) lands on `(◆)` with zero slack (ABC2-SYNTHESIS §4; FORGE5-SL-IR-SINGLE).

---

## §4 · Honest grade

### §4.1 What is rigorously closed (真定理), and how much

| Account | Value | Basis |
|---|---|---|
| `(field × momentum)` **cells closed unconditionally at the physical coupling** | **6 / 9 ≈ 67%** | low-field row (3/3) + high-field row (3/3), §1.1 |
| Cells closed counting the `g→0` limit strip as limit-closed | **8 / 9 ≈ 89%** — the corpus's "~90% region coverage" | + typical×AF-UV (limit only) |
| Open region | **one object** (typical-field column; its two faces collapse by Schur+Ward) | §2 |
| **Fraction of H1's assertion `Var(a_{k→0})≤1/λ` closed** | **0%** | the `k→0` target mode lies *at typical field*, inside the open cell; the closed `k→0` cells are low-field (narrow) and high-field (measure-starved), and **do not contain the physical coarse mode** |

The distinction is the whole honest point: **region-coverage ≠ probability mass ≠ H1.** The closed regions are either the narrow shallow band (`φ≤N/12`) or the measure-starved high-field cone; the measure-dominant typical field (where `21–42%` of directions are negative and the target mode lives) is exactly the open object. "~90%" counts *plane position*, not where the measure is.

### §4.2 The single missing step (精确缺失的一步)

> **Missing lemma:** the conditioned-force-variance bound `(◆)` (equivalently `(◆′)`) on the typical-field × high-momentum cell, uniformly in `(L,a,U′)`. In every vocabulary it is H1: `Var_μ(⟨v,∇V⟩|0) ≤ E_μ[V_vv|0] − λ` ⟺ `G″(0)≥λ` ⟺ `Var(a_{k→0})≤1/λ` ⟺ `D(v) ≤ (1−ε)S(v)`.

Why it is genuinely *the* single missing lemma (not several): G2-ASSEMBLE-100 classifies **26 inputs, every one either a proven theorem or the same species, zero third category**; the residue/isolation arm (`I2`), the FRG face (`E1`), the coupling/Doeblin face (`D2`), the information/Fisher face (`E2`), and the negative-association face (`G2.6`) are all the *same* one-sided variance upper bound. There is no distinct second open item. **But** — see §4.3 — it is the single missing lemma *for H1*, and H1 is not Clay.

### §4.3 Scope caveat: H1 ≠ Clay (REDUCTION-AUDIT)

H1 (`Var(a_{k→0})≤1/λ`) is a **susceptibility/screening** bound (`a` is the zero-momentum mode, `Var(a)=χ=∫s⁻¹dν`), a residue-weighted *harmonic* mean of masses. It is **necessary** for the gap (gap ⇒ `χ≤M₀/λ`) but **not sufficient**: a gapped generalized free field and a massive `:φ²:` both have `f₁²=0`, so `gap ⇏ atom`; conversely an isolated atom at `inf supp ν` ⇒ gap, but only if that channel carries the lightest state. The sufficient condition is the strictly stronger **`H* = I1∧I2`** (residue `f₁²≥r₀` **and** isolation). The corpus *identifies* `H*=H1` structurally (via the residue–level-gap lemma `I1⟺I2` and Ward `G″·Var(a)=1`), but that identification is a **conjecture, not a theorem** (G2-ASSEMBLE-100 §7). So: the single missing lemma for **H1** is `(◆)`; the single missing lemma for **Clay** is `(◆)` **plus** the `H*=H1` identification (or `H*` constructed independently). Proving `(◆)` would close H1 and (on the corpus's conjecture) transport to the gap — it would **not** by itself be a proof of Clay.

---

## §5 · Scorecard (honest five-point, for the assembled H1 chain)

1. **Nonperturbative — PASS.** Exact identities `(G2.0)`/Ward/SD/SCBL, Prékopa, Brascamp–Lieb, WSTC large-deviation mass, exact Savvidy spectrum `−D²+2ad F`; no loop expansion, no `1/N`.
2. **Sign-controlling lower bound `G″(0)≥λ` — PARTIAL.** All four closed regions give the sign (genuine lower bounds); the typical-field cell does **not** — `(◆)` open.
3. **Non-circular — PASS on every closed region; FAIL at the open cell (by construction).** The closed-region arguments assume no gap; the open cell *is* H1, zero slack.
4. **Non-SUSY + fails free field (acid test) — PASS.** Free massless: `F≡0 ⇒ Φ″=0, W_geo=0, G″=0, Var(a)=∞, f₁²=0`; the closed-region inequalities return `0≥0`, no spurious `λ`. `β=0` SU(2) gives `0.248` — the discriminator.
5. **Uniform `m_eff²(k→0)≥λ` — FAIL at the open cell.** Exactly H1.

---

## §6 · Sharpest line

The strongest rigorous chain closes **two full rows of the (field × momentum) plane** — the shallow band by `Φ″⪰N/12`+Brascamp–Lieb (incl. `k→0`), the high-field cone by WSTC starvation, the AF-UV tower by `g_j→0` — and the remaining **typical-field column is one cell, not many**, because the high-momentum Savvidy negatives and the `k→0` target mode are literally the same number through the Schur complement and Ward. That cell is not a gap in the argument: it is **H1 written in closed form**, `Var_μ(⟨v,∇V⟩|0) ≤ E_μ[V_vv|0] − λ`, and its sharpest face is the single measured scalar `D/S ≈ 2% < 1` — a fibre Brascamp–Lieb defect that only overflows because the Savvidy fibre is indefinite (`90–97%` of negatives there). So the honest verdict is: **the region closes; the sign at the one physical cell does not, and that cell is the gap.** And, because H1 is a susceptibility bound, even closing it is a necessary, not a sufficient, step toward Clay.

---

*Refs (internal, each load-bearing item cited where used).* PROOF-SKELETON §0–1; REDUCTION-AUDIT (H1≠Clay, susceptibility vs residue); PROOF-JACOBIAN-SAVES (identity, `W_geo≥0`, `Φ″⩾N/12`, `(◆)`, β=0 `0.248/0.246`); PROOF-COARSE-POINCARE (coarse Poincaré exact, `PᵀHP` min −0.3, Lemma SC, `H_ff` indefinite); INVENT-VARIANCE-INEQ (SCBL `G″=S−D`, `D/S≈2%`, `S^Φ≥N/12`, acid test); ABC-NEW-INEQ (Ward `G″Var(a)=1`, WSTC `≤(cβκ_max)²e^{−cβ(κ²−κ_typ²)N_□}`, negatives high-momentum); ABC-CONSTRUCTIVE (shallow band, AF-UV, `(★)=H1`, ~10× shortfall); ABC2-SYNTHESIS (region plane, one open cell, zero slack); ABC2-NEW-INEQ (STC shells, `g_J=1.27`, `δ_J=2.5–3.8`, transport zero-sum); ABC2-CONSTRUCTIVE (open-cell identity, `|X|≤0.40 > W_geo≈0.32`); OTHER-H1-POINCARE (coarse Dirichlet exact, two-route exhaustion); PROOF-IR-H1 (orthogonal split, fast≠slow, rename); SUB1-COV-BOUND (`D≤cW_geo` false ~10×, conditioning is the gap); FORGE5-CONE-VARIANCE (SD `E[H]=E[f²]` unconditioned, gap-blind); FORGE5-SL-IR-SINGLE (`k→0` residual = H1, martingale `A_0` is the answer); D3-SUSCEPTIBILITY-VARIATIONAL (β=0 `χ≤12/N`, three exact faces, acid clean); PROVE-Gpp-positive (`G″≈1/Var` 2%, Lemmas A/B, `(★)`); k0-CONSTRAINT-EFFPOT (`V_eff″=G″`, β=0 theorem); W6-BE-at-2.5 (`−1.12β`, 21.7%/42% negative, 90–97% fibre); PROOF-SAVVIDY-EXACT-SIGN (`−D²+2ad F`, `λ_min=−g|B|`, `N_neg∝(gB)²`); PROOF-P3-UNCONDITIONAL (`f₁²≥r₀`, `β<1/48`, no `I2`); G2-ASSEMBLE-100 (§3 26-row ledger, §7 `H*=H1` conjecture, no third open item).
*Refs (external, as checked in the sources).* Brascamp–Lieb, Adv. Math. **22** (1976) 366; Prékopa, Acta Sci. Math. **34** (1973) 335; Bakry–Émery 1985; Holley–Stroock, Ann. Probab. **15** (1987) 1541; Shen–Zhu–Zhu, Comm. Math. Phys. **400** (2023) 805, arXiv:2204.12737; Savvidy, Phys. Lett. B **71** (1977) 133; Nielsen–Olesen, Nucl. Phys. B **144** (1978) 376; Leutwyler, Nucl. Phys. B **179** (1981) 129; Osterwalder–Seiler, Ann. Phys. **110** (1978) 440; Kotecký–Preiss, Comm. Math. Phys. **103** (1986) 19; Efron–Stein, Ann. Statist. **9** (1981) 586; Ledoux (LSI pushforward); von Renesse–Sturm, CPAM **58** (2005) 923; Cornwall, Phys. Rev. D **26** (1982) 1453.
