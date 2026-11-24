# SM5-TCI-HEAD — the tachyonic-compensation inequality attacked directly on the Savvidy soft mode at β=2.5: the rank-one/measure structure, the three attempted proofs and where each dies, and the exact logical strength (the strong form **is** the gap and is **false**; the true forms give only a curvature floor)

> 2026-11-24 · Attack surface **SM5**, the single most load-bearing open inequality of the UV face:
> the tachyonic-compensation inequality (TCI), on the **soft (k→0) Savvidy mode** at the head β=2.5.
> Reads: `docs/paper-v2/ch05-theorems.tex` (§Feshbach–Bogoliubov, Lemma NS, `eq:split`,
> `eq:vardecomp`, `eq:tci`); `research/v3/S3-UV-FESHBACH.md`; `H1-RG-DESCENT.md`; `H2-DEFECT-BOUND.md`;
> `H6-HEAD-SPECTRAL.md`; `H4-HEAD-COUNTEREXAMPLE.md`; `STAR2-CLASSIFICATION.md`;
> `unconditional/FORGE-BLIND.md` (§2–4); `unconditional/D-BBS-INTEGRATE-NEG.md`;
> `unconditional/SAVVIDY-DTERM.md`; `unconditional/SAVVIDY-UNIFORM-SOFTENING.md`;
> `unconditional/SAVVIDY-MULTICELL-RG.md`; `unconditional/SAVVIDY-NEGMODE-EXACT.md`.
> Witness: `research/v3/w7/sm5_tci_soft.py` (`py -3.12 -I`, numpy only, deterministic; output quoted §5).
>
> **Labels (charter §2).** **[T]** true theorem/identity · **[C]** conditional/model · **[NG]** no-go for a
> named class · **[O]** exact missing step / open · **[V]** numerical witness here · **[std]** standard
> textbook tool. No QED. Clay 0%.

---

## 0 · Verdict (no QED)

**TCI has, on the soft mode, three *inequivalent* readings — and the answer to "can it be proved?" flips
between them. That trichotomy *is* the result.**

1. **TCI(strong) — essential-infimum-uniform spectral form** (the form the two-block lemma needs to yield a
   *spectral* gap): `ess inf_{a∈supp ρ_v, ‖φ‖=1, φ∈ℋ₋} Σ(φ; a) ≥ μ > 0`. **This form is FALSE `[NG]`** at the
   Savvidy barrier. The soft mode is an exact double well, `V_eff(y)=½m²y²+¼κy⁴`, `m²=−c<0`, `c=g|B|`: the
   well curvature is `−2m²=2c>0` (where TCI holds with margin `2c`), but the **barrier** curvature is
   `m²=−c<0`, and the barrier is in the support of the marginal. Thus `ess inf G″ = −c < 0` — TCI(strong)
   fails, exactly as the STAR2 lattice counterexample `G″=A−Bβ²<0` (`β>4.307`) fails it. **`[NG]`.**

2. **TCI(typical) — pointwise at the coarse minimiser** `a₀*=argmax ρ_v`: `G″(a₀*)=−ρ_v''/ρ_v ≥ 0`. **This
   form is TRUE `[T]`** (the coarse minimiser is interior because the Haar marginal `ρ_v→0` at the chart
   edge), with the soft-mode margin `= −2m² = 2c = 2g|B|` minus the interaction Schur loss. But it is a
   **pointwise curvature floor**, NOT uniform in configuration, and **NOT the mass gap** (curvature ≠ atom:
   at β=0 `G″≥N/12` yet `f₁²=0`). **`[T]`**, gives a curvature floor only.

3. **TCI(avg) — measure-weighted (Fisher)** `E_{a~ρ_v}[G″]=I(ρ_v)=E[(∂_a log ρ_v)²] ≥ 0`: **TRUE `[T]`**
   (an integration-by-parts identity), margin exactly the coarse Fisher information. Its **uniform
   strictness `I(ρ_v) ≥ δ>0` is the coarse spectral gap = H1 `[O]`**, and the measured margin *shrinks with
   volume* (`+0.28` at L=4 → `+0.049` at L=8), the kinematic-artifact signature. **`[T]` identity,
   `[O]` margin.**

**What TCI would give if proved (item 3 of the task).** In its **strong** form TCI ⇒ (two-block lemma) a
**uniform curvature floor** `λ_min(𝔸) ≥ min{μ/(1+2‖B‖²/γ²), γ/2}>0`, and — because a *uniform* floor on
`G″` is exactly what H6-B (`λ₁(ν) ≥ ess inf G″`) turns into a Poincaré constant — the **mass gap** for the
soft mode (Ward: `G″·Var_ν(a)=1`). So **TCI(strong) is strictly stronger than the gap** (`ess-inf-floor ⇒
gap`; the converse fails — the double well has `λ₁ ~ e^{−ΔV}>0` while `ess inf G″=−c<0`). Precisely
*because* it is that strong, **it is false at the head.** The forms that are actually **true** (typical,
averaged) are **strictly weaker than the gap** and deliver only a **curvature floor / sign** — not an atom.
The object that is both true and gap-bearing is neither: it is the *tail-controlled* floor
`ess inf_{a∈bulk} G″(a) ≥ μ` with the Savvidy tail (marginal measure `≤ e^{−c/g²}`) excised = `H6 §8.1` = **H1 `[O]`**.

**The single failing term (item 2).** It is the **rank-one Savvidy softening** `D_Sav = E_μ[|⟨φ₋|c⟩|²]/(g|B|)`
at the head. Its three "faces" are the same object: the **O(1) overlap** `|⟨φ₋|c⟩|²` because the soft mode
sits *inside* the Landau band (`k~√(gB)`, no scale separation); the **IR-marginal** log-log accumulation
`Σ = C_RG Var(φ²) ln(ln(k_UV/Λ)/ln(k_IR/Λ))` across cells; and the **Landau pole** `k_IR→Λ` (`g²→∞`) where
the accumulation is decided and one-loop tools vanish. At β=2.5 this term is `O(S)` (`D/S≈0.97` measured),
so TCI becomes an `O(1)` competition the frozen-Hessian bound cannot resolve.

**Acid test (item 4).** Free massless field: `F=0 ⇒ cF=0`, no negative block, no quartic `Q=0`, no
coupling `K=0` ⇒ TCI reads `0 ≥ μ>0`, which **FAILS** — and the certificate correctly returns `0`. The
averaged form gives flat `ρ` ⇒ `I=0`. Both acid tests pass discriminatively (no false positivity). **TCI
must fail on the free field and does.**

---

## 1 · TCI on the soft mode: exact statement, rank-one structure, measure weight, term-by-term signs

### 1.1 The negative block and the rank-one projector

The soft mode is the lowest-Landau-level (LLL) rung of the quadratic Savvidy operator on transverse 1-forms,
`Q_Sav = −D² + 2 ad F_B`, in the spin–anti-aligned sector. Its eigenvalue is the tachyonic one,
```
λ₋ = − g|B| =: −cF ,     cF := g|B| =: c ,          (1.1)
```
with Landau degeneracy `N_neg = c_N (g|B|)²V/8π²` per cell (`c_{SU2}=1, c_{SU3}=3`), so the block is
**extensive** — dim `ℋ₋ ∝ V F` — even though the **negative eigenvalue is a single rank-one object per
mode.** Choose one soft wavepacket `φ₋` (`‖φ₋‖=1`), e.g. the `k→0` block-uniform combination; the
rank-one projection is
```
P₋ = φ₋ φ₋ᵀ ,     Q_Sav = K₊ − |λ₋| P₋ ,     K₊ ⪰ 0,  P₋ φ₋ = φ₋ ,  K₊ φ₋ = 0 .      (1.2)
```
(K₊ = positive part; φ₋ lies in its kernel.) This is the exact rank-one structure of `SAVVIDY-DTERM` §1:
`K⁻¹ = K₊⁺ − P₋/|λ₋|`, verified to machine zero in the witness (§5·B).

### 1.2 The three TCI terms and their signs

TCI is the demand that the **Schur complement on the negative block** be positive (`FORGE-BLIND` §2–3,
`ch05 eq:tci`):
```
TCI(soft):   ⟨φ₋, Q φ₋⟩ + ⟨φ₋, K φ₋⟩ − cF  ≥  μ > 0 ,      cF = g|B| .        (1.3)
```
The three pieces, with their signs and where each comes from:

| term | symbol | exact content on the soft mode | **sign** |
|---|---|---|---|
| **tachyonic** | `−cF` | bare quadratic Landau energy `λ₋=−g|B|` of the LLL; the instability | **negative** `[T]` |
| **quartic 2nd variation** | `⟨φ₋,Qφ₋⟩` | 2nd variation of the quartic (interaction) part evaluated **at the stabilized saddle** `φ₀²=c/λ₄`: `Q=3λ₄φ₀²=3c=3g|B|` | **positive** `[T]` |
| **Feshbach / Schur lift** | `⟨φ₋,Kφ₋⟩` | the Schur correction `−⟨φ₋,BC⁻¹Bᵀφ₋⟩` from eliminating the **gapped positive complement** `C⪰γ~g|B|`; convention **(A)** (S3 §2.4) ⇒ **hurts** the small block | **≤ 0** `[T]` |

So `Q − cF = 3c − c = 2c = 2g|B| > 0` **at the well**, and the full TCI reads
```
TCI(soft) = 2c + ⟨φ₋,Kφ₋⟩  ≥  μ ,     ⟨φ₋,Kφ₋⟩ ≤ 0 .                       (1.4)
```
Two structural facts pin the signs (both in the corpus, both `[T]`):

- **`K = 0` at strict one loop.** For the pure quadratic Savvidy operator, the LLL is an *exact* eigenspace
  of `Q_Sav` — the inter-Landau-level off-diagonal is `B=0` — so `K=−BC⁻¹Bᵀ=0` one-loop. The Schur lift is
  **interaction-generated** (`O(g²)`, from the quartic's coupling of LLL to higher rungs), and in convention
  (A) it is a **loss**. (The dual "eliminate the *negative* block" Lemma NS, `D-BBS-INTEGRATE-NEG` (b), is
  convention **(B)** where `Σ=A+B|C|⁻¹Bᵀ⪰A` *helps*; TCI is **not** that tool — conflating (A),(B) is the
  named category error of S3 §2.4.)
- **The quartic rescue is nonperturbative, not in the frozen Hessian.** `Q=3c` exists only at the *true
  saddle* `±φ₀`, `|φ₀|²=c/λ₄>0`. The **frozen one-loop Hessian at the perturbative vacuum** has `Q=0`,
  `K=0`, so `TCI = −cF = −c < 0`: *the frozen Hessian fails TCI outright.* TCI can hold only as a statement
  about the **nonperturbative effective potential at its minimum** (`FORGE-BLIND` §3: "convex at its minimum,
  pointwise and nonperturbatively"). This is exactly why the three readings differ (§0).

### 1.3 The measure weight

TCI (1.3) is a statement about the **frozen Hessian / effective potential**, and carries **no measure**. The
object the mass gap actually needs is the **measure-carried** marginal curvature (`H2`, `H4`, `eq:split`),
```
G″_soft(a₀) = E_{μ(·|U′)}[ ∂_v²V | a₀ ] − Var_{μ(·|U′)}( ∂_vV | a₀ ) = S − D ,                 (1.5)
```
`V = S_W + Φ`, `Φ=−log J` the Haar potential, `a=⟨v,X⟩` conditioned at `a₀`. On the soft mode the **measure
weight enters through the rank-one D-term** (`SAVVIDY-DTERM` §1, `H2` §1): with `K⁻¹=K₊⁺−P₋/|λ₋|`,
```
D = D₊ + D_Sav ,   D₊ := Var(f|0)−E⟨K₊⁺c,c⟩ ≤ 0  (BL, log-concave K₊⪰0) ,
                   D_Sav := E_μ[ |⟨φ₋|c⟩|² ] / (g|B|)  ≥ 0 .                                 (1.6)
```
The **measure weight** is the marginal `ρ_v(a)` of the soft amplitude `a`, equivalently the law of
`|⟨φ₋|c⟩|²` averaged over `μ(·|0)` — and `D_Sav` is the radial `φ₋`-component of `Var(∂_vV)` pushed through
the rank-one projector. The measure-weighted TCI is therefore
```
TCI(measure):   G″_soft = S − D  ≥  μ    ⟺   D_Sav < S − δ + |D₊|    (exact, since D₊≤0) .       (1.7)
```
`D_Sav` is the whole game: it is the measure-weighted softening that the quartic overcompensation `2c` must
beat. At the head `D_Sav/D ≈ 1` (`D/S≈0.83–1.10`, mean `0.97`, `H2 §1`): the rank-one term **saturates** `S`.

---

## 2 · The three attempted proofs, and the exact failing term

### 2.1 Route 1 — one-octave support (`H2 §2.1`): fails on the O(1) overlap

`φ₋` is LLL-localized in momentum at `|p|≲√(gB)`. By Cauchy–Schwarz with `‖φ₋‖=1`,
```
D_Sav ≤ E_μ[ ∫_{|p|≲√(gB)} |c(p)|² dp ] / (g|B|) .       [T]                                 (2.1)
```
The **soft mode sits inside the band** (`k ~ √(gB)`), and `c=Π₀ Hess V v` has scale `~gB`, so the right
side is `O((gB)²/(gB)) = O(gB) = O(S)` — **no small parameter.** The overlap itself is `O(1)`: the witness
(§5·C, reproducing `H2`) gives `overlap² = x²/(1+x²)`, `x=k/√(gB)`, i.e. `0.0099` at `x=0.1` (deep IR) but
`0.500` at `x=1` (**the head**). So the one-octave bound is `(L,a)`**-uniform** but **not `<1`** — it does
not close. **Exact failing term: the `O(1)` overlap `|⟨φ₋|c⟩|²` at `k~√(gB)`.**

### 2.2 Route 2 — two-loop sign protection (`H2 §2.3`, `SAVVIDY-NEGMODE-EXACT`): real but too lossy at the head

The exact (nonperturbative, true-saddle `±v`) susceptibility is bounded and `c`-independent,
`Var(φ²)→2/λ₄` as `c=g|B|→∞`, so
```
D_Sav ≤ 2 g² / λ₄ .        [T]  (λ₄-set; 2-loop sign-protected: soft ≥ 0)                    (2.2)
```
This is a genuine uniform upper bound — but at β=2.5 (`g²=4/β=1.6`, `S≈1.4`) and **natural `λ₄~g²`**:
`2g²/λ₄ ≈ 2 > S ≈ 1.4` (§5·D). The bound **overshoots the very thing it must beat.** `.5`, `2g²/λ₄=4`;
`.1`, `=2`; only at `λ₄/g²≳1.5` does it dip below `S`. The `bound − measured ≈ 0.6` is the part of `D` the
rank-one ceiling cannot see. **Exact failing term: the marginal `O(1)` inequality `2g²/λ₄ < S_bare` = the
softened sub-criticality.**

### 2.3 Route 3 — Feshbach decomposition (`S3`, `D-BBS`): exact, but the sign of the Schur term is the whole story

Eliminate the gapped positive complement `C⪰γ~g|B|`:
```
Σ = A − B C⁻¹ Bᵀ ,   A = A_S + quartic ,   convention (A):  Σ ⪯ A  (coupling HURTS) .       (2.3)
```
- **One loop (frozen Hessian):** `B=0` (LLL exact eigenspace) ⇒ `Σ = A = −cF = −c < 0`. **TCI fails.**
- **True saddle (add the nonperturbative quartic `Q=3c`):** `Σ = 2c + K`, `K ≤ 0`, `|K| = O(g²)` from the
  interaction-induced LLL↔higher-rung coupling divided by `γ~g|B|` (the exact power is unproven — the
  graph-angle series `Σ_j g_j⁴` of `S3 §2.2`/`H1 §2` is the relevant counting). So `TCI = 2c − |K|`, an **`O(1)`
  competition** between the quartic overcompensation and the Feshbach/measure loss.
- **The barrier (the essential infimum):** on the same double well, *at the barrier* `y=0`, the Feshbach
  decomposition gives `Σ = m² = −c < 0`. The two-block lemma, applied to the Hessian over the **support** of
  the marginal, must include the barrier ⇒ `ess inf Σ = −c < 0` ⇒ **TCI(strong) fails.**

So the Feshbach route is **exact and sign-transparent**, and it proves that `TCI` at the well is `2c − |K|`
and at the barrier is `−c`. It closes nothing uniformly: the **barrier is in the support** and its curvature
is negative. **Exact failing term: the barrier / essential infimum — equivalently the interaction Schur loss
`|K|` growing into the `2c` margin.**

### 2.4 Route 4 — the measure-weighted form (`H2` Theorem M, `H3`): true as an identity, uniform margin = the gap

Averaging the marginal curvature over its own marginal kills the barrier problem **in expectation**:
```
E_{a~ρ_v}[G″(a)] = −∫ ρ_v (log ρ_v)'' da = I(ρ_v) := E_{ρ_v}[(∂_a log ρ_v)²] ≥ 0 .     [T]       (2.4)
```
This is exact (integration by parts; `∫ρ_v''=0` for a normalisable density), **automatically positive**, and
**`=0` on the free massless field** (flat `ρ_v`), so the acid test passes discriminatively (§5·E). But the
margin is `I(ρ_v) = 1/Var_ρ(a)` (Gaussian) — **the coarse mass-squared itself** — and the measured `I`
**shrinks with volume** (`+0.28` at L=4 → `+0.049` at L=8). So:
```
TCI(measure-weight):   E_{ρ_v}[G″] = I(ρ_v) ≥ 0  [T]  ;   uniform I(ρ_v) ≥ δ>0  = H1  [O] .   (2.5)
```

### 2.5 The exact failing term, stated once

All four routes die on the **same object**, in four costumes. The **exact failing term** is the rank-one
Savvidy softening `D_Sav = E_μ[|⟨φ₋|c⟩|²]/(g|B|)`:
- **O(1) overlap** — the soft mode is inside the Landau band (`k~√(gB)`), so `|⟨φ₋|c⟩|²=O(1)`, no scale
  separation (Route 1; `H2 §2.1`);
- **IR-marginal softening** — accumulated across cells as `Σ = C_RG Var(φ²) ln(ln(k_UV/Λ)/ln(k_IR/Λ))`,
  dominated by the last few cells near `Λ` (`SAVVIDY-MULTICELL-RG`, `H2 §2.4`);
- **Landau pole** — `k_IR→Λ⁺`, `g²(k)→∞`, where the accumulation is *decided* and perturbation theory dies.

These are not three failures but one: `D_Sav` at `k~√(gB)`, `O(S)`, whose accumulation is IR-marginal and
decided at the Landau pole. **The single inequality left is `(★)`-softened: `D_Sav < S_bare − δ` at the
crossover, an `O(1)` comparison = H1.**

### 2.6 The β-dependent version (what a β-dependent TCI can and cannot do)

TCI is **provable at both ends of the β axis and open only at the head**:
- **Strong coupling β<1/48:** Shen–Zhu–Zhu give a volume-uniform LSI `K_S=N(½−24|β|)>0` (`[T]`, `C1`), so
  the marginal is log-concave and TCI holds with the full strong-coupling constant (§5·F). `[T]`.
- **Deep UV `g_j→0`:** the negative set has marginal measure `≤ e^{−c/g_j²}→0` (Lemma B), the Savvidy block
  shrinks (`N_neg ∝ (g_jB)²V→0`), and TCI is certifiable with `δ_W ≳ c_L=4π²/L²` (`[C]`, cone-emptiness).
- **Head β≈2.5:** `g=O(1)`, `N_neg` large, `D_Sav=O(S)`, the `O(1)` competition. `[O]`.

So a **β-dependent TCI `μ(β)`** is a genuine theorem at `β<1/48` and UV-certifiable at `g_j→0`; the head
`β≈2.5` is exactly the disjoint intermediate regime where no method banks a constant (`STAR2 §5`,
`S6 §4`). **This is the tower's head, and it is H1.** No β-dependent TCI bridges it.

---

## 3 · What TCI would give if proved — precise logical strength

The task asks: *curvature floor or mass gap?* The answer is **both, in different readings — and the strong
reading that gives the gap is precisely the one that is false.**

**(a) Strong form `TCI(ess-inf, uniform)` ⇒ mass gap, and is strictly stronger than the gap.**
If `∃μ>0: ess inf_{a∈supp ρ_v,‖φ‖=1,φ∈ℋ₋} Σ(φ;a) ≥ μ` uniformly in volume/cutoff, then
1. two-block lemma (`ch05 eq:twoblock-bound`): `𝔸 ≥ min{μ/(1+2‖B‖²/γ²), γ/2}·I > 0` — a **uniform
   curvature floor** of the full Hessian;
2. H6-B (`λ₁(ν) ≥ ess inf_a G″(a)`, `[T]`, sharp in 1-D) turns the uniform floor into the **coarse Poincaré
   constant** `λ₁(ν) ≥ μ`; Ward (`G″·Var_ν(a)=1`) makes this the **coarse mass-squared**;
3. hence the **mass gap** for the soft mode, uniformly. So `TCI(strong) ⇒ gap`.

The implication is **strict**: the Savvidy double well has `λ₁(ν) ~ e^{−ΔV} > 0` while `ess inf G″ = m²<0`
(H6 §3, `[T]`+`[V]`). So `TCI(strong)` is **sufficient but not necessary** for the gap — it is a *curvature
floor that, by being uniform, is logically at least as strong as the gap.*

**(b) But TCI(strong) is FALSE `[NG]`.** The exact double well (§1.2, §5·A) and the STAR2 lattice
counterexample (`G″=A−Bβ²<0`, `A=0.20511530`, `B=0.01105929`, `β*=4.30661`) both violate `ess inf G″≥μ`:
the barrier is a genuine negative-curvature point in the support. **TCI in its gap-bearing form cannot be
proved at the head — it is refuted.** This is not a tool gap; it is a correct no-go for the *uniform
essential-infimum* statement.

**(c) The true weakenings deliver only a curvature floor.**
- `TCI(typical)` at the coarse minimiser: `G″(a₀*)≥0` `[T]` (H4, interior minimiser). This is a **pointwise
  curvature floor**, not uniform; and curvature ≠ atom: at β=0 `G″≥N/12` yet `f₁²=0`. **Strictly weaker
  than the gap; does not give it.**
- `TCI(avg)` Fisher: `E[G″]=I(ρ_v)≥0` `[T]`, uniform margin `= gap` `[O]`. A **sign**, whose uniform strength
  is the open gap.

**(d) The object that is both true and gap-bearing (the actual target).** Neither `ess inf` over the full
support (false) nor a single typical point (insufficient), but the **tail-controlled** floor
```
(◆)   ess inf_{a ∈ bulk of ρ_v} G″(a) ≥ μ > 0 ,   with the Savvidy tail (measure ≤ e^{−c/g²}) excised  =  H6 §8.1  =  H1  [O] .
```
The Savvidy negativity sits on a set of marginal measure `≤ e^{−c/g²}` (Lemma B); if the *bulk* essential
infimum is positive and the tail is measure-smaller than the gap scale, the transfer-matrix bound `λ₁ ≥ μ`
survives. Proving `(◆)` is the campaign's H1, and is strictly **between** `TCI(strong)` (false) and
`TCI(typical)` (true, insufficient).

**One-line logical strength.** *TCI(strong) is a uniform curvature floor that is exactly the mass gap (via
H6-B) and is therefore false at the Savvidy barrier; the true forms (typical-point, Fisher-averaged) are
curvature-floor/sign statements that are strictly weaker than the gap; the gap needs the tail-controlled
floor (◆), which no TCI-local argument supplies.*

---

## 4 · Acid test, ledger, and reading

### 4.1 Acid test — free field must fail TCI (and does)

| object | free massless field (`F=0`, `Φ≡0`) | passes? |
|---|---|---|
| TCI `Q+K−cF ≥ μ` | `cF=0`, no LLL (`N_neg=0`), `Q=0`, `K=0` ⇒ `0 ≥ μ>0` **false** | **YES** — returns `0`, gapless `[T]` |
| one-octave `D_Sav` | `B=0` ⇒ no negative mode ⇒ `D_Sav=0` | YES (vacuous, not false-positive) `[T]` |
| two-loop `2g²/λ₄` | `g²→0`; the true object is `D_Sav=0` | YES (as a *bound*, vacuous off-cell) `[T]` |
| Fisher `E[G″]=I(ρ)` | flat `ρ` ⇒ `I=0` ⇒ `E[G″]=0` | **YES, discriminatively** `[T]`+`[V]` |
| `β=0` Haar contrast | `G″≥N/12>0` but `f₁²=0` — curvature, not atom | YES (shows floor≠atom) `[T]` |

The free field is reached by `cF→0`; every version of TCI correctly **fails** (returns `0`). No tool
manufactures false positivity. The charter acid test is passed.

### 4.2 Ledger

| item | statement | label |
|---|---|---|
| rank-one structure `Q_Sav=K₊−|λ₋|P₋`, `K⁻¹=K₊⁺−P₋/|λ₋|` | exact (1.2); witness split error `0` | **[T]** |
| TCI terms/signs: `−cF<0`, `Q=3c>0` (at saddle), `K=−BC⁻¹Bᵀ≤0` (conv. A) | §1.2; `Q−cF=2c` verified | **[T]** |
| double well: `G″(well)=2c>0`, `G″(barrier)=m²=−c<0` | §5·A; normal form | **[T]**+**[V]** |
| `TCI(strong)` = ess-inf-uniform ⇒ two-block floor ⇒ (H6-B) mass gap | §3(a) | **[T]** (implication) |
| `TCI(strong)` itself | **FALSE** at the Savvidy barrier | **[NG]** |
| `TCI(typical)` at coarse minimiser `G″(a₀*)≥0` | interior minimiser; curvature floor only | **[T]** (H4) |
| `TCI(avg)` = Fisher identity `E[G″]=I(ρ_v)≥0` | IBP; free-field `=0` | **[T]** |
| uniform `I(ρ_v)≥δ>0` (= gap) | measured shrinking with L | **[O]** |
| one-octave `D_Sav≤E∫_oct|c|²/(gB)=O(S)` | `(L,a)`-uniform, **not `<1`** at head | **[T]** (lossy) |
| two-loop `D_Sav≤2g²/λ₄`, natural `λ₄~g²` ⇒ `≈2>S≈1.4` | exact, **too lossy** | **[T]** (lossy) |
| Feshbach: one-loop `Σ=−c<0`; true-saddle `Σ=2c−|K|`; barrier `Σ=−c` | §2.3 | **[T]** |
| exact failing term = `D_Sav` (O(1) overlap ⊗ IR-marginal ⊗ Landau pole) | §2.5 | **[T]** structure / **[O]** endpoint |
| `β<1/48`: SZZ LSI ⇒ TCI holds strongly | `K_S=N(½−24|β|)>0` | **[T]** |
| UV `g_j→0`: bad-set measure `≤e^{−c/g²}`, TCI certifiable `δ_W≳c_L` | cone-emptiness | **[C]** |
| **head β≈2.5**: TCI margin | `O(1)` competition, no uniform constant | **[O]** |
| TCI in *any* form bridges the head / gives the gap | **false**: strong form `[NG]`; true forms curvature-only | **[NG]** (of that claim) |
| **(◆) tail-controlled bulk log-concavity** | the true gap-target = H1 | **[O]** |
| **Clay** | unbuilt | **0%** |

### 4.3 Reading

TCI on the soft mode is **not one inequality but a trichotomy**, and the trichotomy is the result. The
*strong* form (essential-infimum-uniform) is exactly a **curvature floor**, and a uniform curvature floor is
**logically the mass gap** (H6-B + Ward) — which is why it is **false** at the Savvidy barrier (exact double
well; STAR2 lattice counterexample). The *true* forms (typical-point at the coarse minimiser, Fisher-averaged)
are **strictly weaker** and give only a **curvature floor / sign**, never an atom. The single **failing term**
is the rank-one Savvidy softening `D_Sav`, whose three faces (O(1) overlap at `k~√(gB)`, IR-marginal
accumulation, Landau pole) are the same object; at β=2.5 it is `O(S)` and TCI degenerates to an `O(1)`
competition. **What TCI would give if proved in the strong form *is* the mass gap — and that is precisely
why it cannot be proved: the strong form is refuted, and the provable weakenings stop at a curvature floor.**
The gap needs the tail-controlled floor `(◆)` (bulk log-concavity, Savvidy tail `e^{−c/g²}` excised) — H1,
still open. **No QED. Clay 0%.**

---

## 5 · Witness — `research/v3/w7/sm5_tci_soft.py` (`py -3.12 -I`, deterministic)

Quoted output (abridged):

```
(A) DOUBLE-WELL: TCI at the well vs at the barrier  (c = g|B|)
  c=gB=1.00 lam4=1.00 : y0^2= 1.000  Q=3lam4 y0^2= 3.000  -cF=-1.00  TCI(well)=Q-cF= 2.000 (=2c)  G''(barrier)=-1.000
  c=gB=1.60 lam4=1.60 : y0^2= 1.000  Q=3lam4 y0^2= 4.800  -cF=-1.60  TCI(well)=Q-cF= 3.200 (=2c)  G''(barrier)=-1.600
  c=gB=0.50 lam4=2.00 : y0^2= 0.250  Q=3lam4 y0^2= 1.500  -cF=-0.50  TCI(well)=Q-cF= 1.000 (=2c)  G''(barrier)=-0.500
  => TCI holds at the WELL with margin 2c>0 ; the ESSENTIAL INFIMUM over the support is m^2=-c<0. Uniform/ess-inf TCI FAILS.

(B) EXACT RANK-ONE SPLIT  K = K+ - |lam_-| P_- ,  D = D+ + D_Sav
  <K^-1 c,c> = 3.9509005729 ;  <K+^+ c,c> - <phi,c>^2/|lam_-| = 3.9509005729 ;  split error = 0.00e+00
  D+ = Var(f|0) - E<K+^+ c,c>  <= 0 (BL) ;  D_Sav = E<phi,c>^2/|lam_-|  >= 0

(C) ONE-OCTAVE OVERLAP  overlap^2 = x^2/(1+x^2),  x = k/sqrt(gB)
     x=0.05 -> 0.0025 (deep IR) ; x=0.10 -> 0.0099 ; x=0.30 -> 0.0826 ; x=0.50 -> 0.2000
     x=1.00 -> 0.5000 (HEAD, O(1)) ; x=2.00 -> 0.8000
  => suppression only deep IR ; at the head k~sqrt(gB) overlap is O(1).

(D) TWO-LOOP SIGN-PROTECTED BOUND  D_Sav <= 2 g^2/lam4   vs   S
   beta=2.5 -> g^2=1.600 ; S_bare=1.40
     lam4/g^2=0.5 -> 4.000 NO ; =1.0 -> 2.000 NO ; =1.5 -> 1.333 YES ; =2.0 -> 1.000 YES
  => for NATURAL lam4~g^2 the exact bound ~2 > S~1.4 : too lossy at the head.

(E) FISHER IDENTITY  E_rho[G''] = I(rho) = 1/sigma^2 ; free -> 0
     Gaussian var 4.00 -> 0.2500 ; 1.00 -> 1.0000 ; 0.25 -> 4.0000 ; free (flat) -> 0.0000
     volume-shrinking margin (H2 measured): L=4 -> +0.28 ; L=8 -> +0.049
  => averaged TCI is TRUE [T]; its uniform margin is the coarse gap [O].

(F) BETA-DEPENDENT MARGIN mu(beta)
   beta=0.00000 K_S(N=2)=1.0000 LSI>0 [T] ; beta=0.02083 (1/48) K_S=0.0000 ; beta=0.5 K_S=-23 ; beta=2.5 K_S=-119
  => TCI provable (strong) for beta<1/48 ; UV certifiable (measure of bad set <= e^{-c/g^2}) ; HEAD beta~2.5 is the O(1) competition.
```

Readings: **(A)** the double well is the whole obstruction — TCI holds at the well (`2c`), fails at the
barrier (`−c`), so the essential infimum is negative; **(B)** the rank-one split `K⁻¹=K₊⁺−P₋/|λ₋|` is exact
to machine zero and `D₊≤0, D_Sav≥0`; **(C)** the overlap is `O(1)` at the head, suppressed only deep IR; **(D)**
the two-loop bound overshoots `S` for natural `λ₄~g²`; **(E)** the Fisher identity is exact and vanishes on
the free field, with a volume-shrinking margin; **(F)** TCI is proved at strong coupling, open at the head.
**Not a proof — a witness for the trichotomy, the rank-one algebra, and the precise location of the failure.**

---

## 6 · References (internal, this session)

- `docs/paper-v2/ch05-theorems.tex` — `eq:tci` (TCI), `eq:twoblock`/`eq:twoblock-bound` (two-block lemma),
  Lemma NS (convention B), `eq:split`/`eq:vardecomp` (variance decomposition), `sec:feshbach`.
- `research/v3/S3-UV-FESHBACH.md` — Lemma S (`G″=E[Σ]−D`), §2.4 conventions (A)/(B), §3 TCI vs `(★)`.
- `research/v3/H1-RG-DESCENT.md` — head soft mode `γ_j→0`, `η_j≥1`, `G″≤0` at β≈2.4–2.5; curvature≠atom.
- `research/v3/H2-DEFECT-BOUND.md` — rank-one split, three routes, Theorem M (Fisher), `D/S≈0.97`.
- `research/v3/H6-HEAD-SPECTRAL.md` — `λ₁(ν)≥ess inf G″` (H6-B), typical-point vs essential-infimum,
  Savvidy double well, `[NG]` spectral class, `(◆)` the missing object.
- `research/v3/H4-HEAD-COUNTEREXAMPLE.md` — `G″(a₀*)≥0` at the interior coarse minimiser `[T]`.
- `research/v3/STAR2-CLASSIFICATION.md` — lattice barrier counterexample, head = H1.
- `research/unconditional/FORGE-BLIND.md` §2–4 — TCI statement, two-block lemma, free-field acid test.
- `research/unconditional/D-BBS-INTEGRATE-NEG.md` — Lemma NS (convention B), `Σ⪰A` for `C≺0`.
- `research/unconditional/SAVVIDY-DTERM.md` — `D=D₊+D_Sav`, `D_Sav=E|⟨φ₋|c⟩|²/(g|B|)`, `D<S⟺G″>0`.
- `research/unconditional/SAVVIDY-UNIFORM-SOFTENING.md` (one-octave `(L,a)`-uniform, `O(1)`),
  `SAVVIDY-MULTICELL-RG.md` (IR-marginal log-log, Landau pole), `SAVVIDY-NEGMODE-EXACT.md` (`Var(φ²)≤2/λ₄`),
  `SAVVIDY-2LOOP.md` (sign protection).
- Standard (used as textbook tools, not re-fetched): Schur complement / Haynsworth inertia; Prékopa (1971);
  Brascamp–Lieb (1976); Bakry–Émery (1985); Savvidy (1977), Nielsen–Olesen (1978); Kramers metastability.
- Witness `research/v3/w7/sm5_tci_soft.py` (deterministic, `py -3.12 -I`, numpy 2.4.6; exit 0).

*No third-party sources were downloaded this session; all external results are standard and cited from the
corpus files, not re-fetched.*
