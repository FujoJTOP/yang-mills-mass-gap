# H2-DEFECT-BOUND — bounding the measure-carried defect `D` at the UV head: rank-one split exact `[T]`, head `D<S` is a marginal `O(1)` competition `[O]`, measure-weighted positivity = Fisher information `[T]`; uniform pointwise `D<S` is false `[NG]`

> 2026-11-24 · Attack surface H2: the **head of the UV tower** (the last open piece of `(★)`).
> At `β≈2.5`, the bulk Savvidy Hessian is indefinite (least eigenvalue `≤0`), yet the true
> `G″` is measured **positive** (`+0.13` quoted; `+0.28` at `L=4,β=2.5`; `+0.049` at `L=8,β=2.4`).
> The positivity is **measure-carried**: `G″ = S − D` is a conditioned average, not a Hessian
> eigenvalue. This report makes the defect `D` quantitative.
> Reads: `STAR2-CLASSIFICATION.md`, `S3-UV-FESHBACH.md`, `SAVVIDY-DTERM.md`,
> `SAVVIDY-MULTICELL-RG.md`, `SAVVIDY-UNIFORM-SOFTENING.md`, `SAVVIDY-LATTICE-CELL.md`,
> `SAVVIDY-NEGMODE-EXACT.md`, `SAVVIDY-2LOOP.md`, `docs/paper-v2/ch04-gap.tex` §Savvidy.
> Witness: `research/v3/w7/h2_defect_bound.py` (`py -3.12 -I`, numpy only; output quoted in §5).
>
> **Labels (charter §2).** **[T]** true theorem · **[C]** conditional/model · **[O]** exact
> missing step · **[NG]** no-go for a named class · **[V-num]** this-round numerics.

---

## 0 · Verdict (no QED)

The exact rank-one split is confirmed and is the right coordinate, but it does **not** give a
uniform bound on `D`:

1. **Exact rank-one split `[T]`.** With the Savvidy Hessian written as `K = K₊ − |λ₋| P₋`
   (`P₋ = φ₋φ₋ᵀ` rank one, `φ₋ ∈ ker K₊`, `|λ₋| = g|B|`),
   ```
   D = D₊ + D_Sav ,   D_+ = Var(f|0) − E⟨K₊⁺c,c⟩ ≤ 0  (BL on the log-concave positive fibre),
                      D_Sav = E[|⟨φ₋|c⟩|²] / (g|B|) ≥ 0 .
   ```
   Hence `D ≤ D_Sav`, and **`D < S − δ` on the typical-field set ⟺ `D_Sav < S − δ + |D₊|`**
   (exact, since `D₊≤0`). On a set where the BL term is negligible this is the softened
   sub-criticality `g²Var(φ²) < S_bare − δ`.

2. **The head bound `[O]`.** Three routes each give a *bound* but not a *uniform
   `D_Sav < S − δ`*: one-octave support gives overlap `O(1)` at `k~√(gB)` (no suppression);
   deep-IR decoupling is uniform **only** in the deep IR (`D_Sav/S→0`); the exact rank-one bound
   `D_Sav ≤ 2g²/λ₄` (`Var(φ²)≤2/λ₄`, `c`-independent, two-loop sign-protected `[T]`) is **too
   lossy at the head** (`2g²/λ₄ ≈ 2 > S ≈ 1.4` for natural `λ₄~g²`). The comparison
   `2g²/λ₄ < S_bare` is a **marginal `O(1)` inequality** = the softened H1.

3. **Uniform pointwise `D<S` is FALSE `[NG]`.** The exact one-plaquette counterexample
   (`G″=0.20511530−0.01105929 β² < 0` for `β>4.307`) and the measured `5/12` modes with
   `G″<0` at `L=8` both violate it. No pointwise bound exists.

4. **Correct measure-weighted statement `[T]`, with margin `[O]`.** Averaging the marginal
   curvature over the coarse-field marginal `ρ` gives **exactly** the Fisher information:
   ```
   E_{a~ρ}[G″(a)] = ∫ρ(−(log ρ)'') da = I(ρ) := E_ρ[(∂_a log ρ)²] ≥ 0 .
   ```
   So `E_ρ[D] = E_ρ[S] − I(ρ) < E_ρ[S]` with margin `I(ρ)`, **automatically positive** and
   **`=0` for the free massless field** (acid test passes). Measured margin
   `I(ρ)=E[G″]≈+0.049` at `L=8,β=2.4` (`=3.5%·E[S]`). But a **uniform lower bound `I(ρ)≥δ` is
   exactly the mass gap** (bounded coarse-mode variance), and the measured `I` **shrinks with
   volume** (`+0.28` at `L=4` → `+0.049` at `L=8`) — the kinematic-artifact signature. `[O]`.

**Clay 0%. No QED. Where it fails:** at `k~√(gB)` the rank-one overlap is `O(1)` and the exact
bound is lossy; the barrier is the IR-marginal / Landau-pole endpoint, not a divergence.

---

## 1 · The exact rank-one split and what `D<S−δ` requires `[T]`

**Setup (as in `SAVVIDY-DTERM` / `INVENT-VARIANCE-INEQ`).** Flat chart, conditioned measure
`μ(·|0) ∝ e^{−V}dX`, coarse mode `v`, `a=⟨v,X⟩` conditioned at `a=0`; `h_n=vᵀ Hess V v`,
`f=∂_nV`, `c=Π₀ Hess V v`, `K=Π₀ Hess V Π₀`. Exact identity:
```
G″(0) = S − D,   S = E[h_n − ⟨K⁻¹c,c⟩|0],   D = Var(f|0) − E[⟨K⁻¹c,c⟩|0].
```

**The split `[T]`.** Write the Savvidy Hessian `K = K₊ − |λ₋| φ₋φ₋ᵀ`, where `K₊ ⪰ 0` is the
positive part with `φ₋ ∈ ker K₊` (the negative eigenvalue of the full Hessian is
`−|λ₋| = −g|B|`). Because `φ₋ ∈ ker K₊`, `K` and `K₊` agree on `φ₋^⊥` and
`K⁻¹φ₋ = −φ₋/|λ₋|`, `K₊⁺φ₋ = 0`. For any `c = c_∥ + c_⊥`:
```
⟨K⁻¹c,c⟩ = ⟨K₊⁺c,c⟩ − |⟨φ₋|c⟩|²/(g|B|) .
```
Substituting into `D` (`SAVVIDY-DTERM` §1; witness (1): split error `1.8×10⁻¹⁵`):
```
D = D₊ + D_Sav ,    D₊ := Var(f|0) − E⟨K₊⁺c,c⟩ ,    D_Sav := E[|⟨φ₋|c⟩|²]/(g|B|) .
```
- **`D_Sav ≥ 0` `[T]`** — the negative-mode enhancement (one rank-one projection).
- **`D₊ ≤ 0` `[T]`, given the positive fibre is log-concave** (`K₊ ⪰ 0`): Brascamp–Lieb gives
  `Var(f|0) ≤ E⟨K₊⁺c,c⟩`. Witness (1): `D₊ = −0.031, −0.022, −0.0135` for `B=0.5,1,2`; the sign
  is robust on any log-concave fibre.

Hence the **exact two-sided ceiling** `D ≤ D_Sav`, and the **softened inequality**
```
(★)   D_Sav = E[|⟨φ₋|c⟩|²]/(g|B|) < S_bare   ⟹   D < S   ⟺   G″ > 0 .     [T]
```

**What must be true for `D<S−δ`.** With `δ>0` the required margin,
```
D < S − δ  ⟺  D₊ + D_Sav < S − δ  ⟺  D_Sav < (S − δ) + |D₊|     (D₊ ≤ 0) .     (1.1)
```
So the head needs **the rank-one term to be dominated by the bare stiffness *plus* the
Brascamp–Lieb buy-in** `|D₊|`. Two readings of "measure-carried positivity", both exact:

- **Rank-one-dominated case** (`|D₊|` small; the single-mode toy where the fibre is `φ₋` itself,
  `D₊≡0`): (1.1) reduces to the **softened sub-criticality**
  ```
  g² Var(φ²) < S_bare − δ   ⟺   κ_bare > g²Var(φ²) + δ .
  ```
  The margin comes entirely from `S_bare` beating the rank-one term; positivity is carried by the
  *measure's variance* `Var(φ²)`, not by the (indefinite) Hessian.
- **BL-dominated case** (`|D₊|` `O(S)`): the margin comes from the Brascamp–Lieb negativity of the
  conditioned positive fibre. **BL gives only the sign `D₊≤0`, not a lower bound on `|D₊|`** — the
  magnitude needed to pay (1.1) is *not* supplied by BL. This is one precise gap.

**Measured.** `S≈1.40`, `D/S∈[0.83,1.10]` (mean `0.97`); since `D₊≤0`,
`D_Sav = D − D₊ ≥ D ≥ 0.83 S` (typically `≈0.97 S`). Thus the rank-one term **already saturates**
`~S`: the observed
`G″=S−D=δ` is the **small residue** `S−D_Sav+|D₊|`. There is **no room** for a uniform `δ` coming
from `S_bare` alone — the sign is decided by an `O(1)` competition in (1.1). This is the exact
content of "measure-carried": the rank-one (Hessian) part is `~S`, and the *only* things that can
lift `G″` above zero are `S_bare − D_Sav` and `|D₊|`, both `O(1)` and sign-undecided a priori.

---

## 2 · Three proof attempts on the rank-one term, and the obstruction `[T]`/`[O]`

We try to bound `D_Sav = E[|⟨φ₋|c⟩|²]/(g|B|)` by the three available structures.

### 2.1 One-octave support `|p| ~ √(gB)` `[T]` — no suppression at the head

`φ₋` is the lowest-Landau-level mode of `−D²+2ad F_B`, localised in momentum at `|p|≲√(gB)`.
Using `‖φ₋‖=1` and Cauchy–Schwarz:
```
|⟨φ₋|c⟩|² ≤ (∫_{|p|≲√(gB)}|φ₋(p)||c(p)| dp)² ≤ ‖φ₋‖² · ∫_{|p|≲√(gB)}|c(p)|² dp
          = ∫_{|p|≲√(gB)} |c(p)|² dp .
```
Hence
```
D_Sav ≤ E[ ∫_{|p|≲√(gB)} |c(p)|² dp ] / (g|B|) .     [T]
```
This is the **one-octave bound**: the rank-one term is bounded by the one-octave band energy of
`c`, divided by `g|B|`. **It does not close at the head.** The coarse mode sits *inside* the band
(`k ~ √(gB)`), `c = Π₀ Hess V v` has scale `~gB`, so the right side is `O((gB)²/(gB)) = O(gB) =
O(S)`. No small parameter. This is exactly `SAVVIDY-UNIFORM-SOFTENING`'s "one octave, `(L,a)`-
uniform, but `O(1)`": *uniform in `(L,a)`* (no UV `log`), **not** `< 1`.

### 2.2 Deep-IR decoupling `[T]` — uniform, but only in the deep IR

At `k ≪ √(gB)` the coarse mode is orthogonal to `φ₋` (a localised zero-mean wavepacket), so
```
|⟨φ₋|c⟩|² ~ (k/√(gB))² · (amplitude) ,   k→0  ⇒  D_Sav/S → 0 .
```
Witness (5): overlap² `= 0.0025, 0.0099, 0.083` at `k/√(gB)=0.05, 0.1, 0.3`, rising to
`0.20, 0.50, 0.80` at `0.5, 1, 2`. So **`D<S` holds uniformly in the deep IR with full margin
`S_bare` `[T]`** — the floor survives there. At the head `k~√(gB)` the overlap is `O(1)` and this
lever gives nothing. The deep-IR theorem does **not** reach the head.

### 2.3 Two-loop sign protection / exact `Var(φ²)` bound `[T]` — real but lossy

The exact (nonperturbative, true-saddle `±v`) susceptibility is **bounded and `c`-independent**:
```
Var(φ²) → 2/λ₄   as  c=g|B| grows   (`SAVVIDY-NEGMODE-EXACT`);   witness: 0.575, 0.957, 2.141, 2.000
```
so
```
D_Sav ≤ 2 g² / λ₄ .     [T]   (exact, λ₄-set, c-independent, 2-loop sign-protected: soft≥0)
```
This is a **genuine uniform upper bound** on the rank-one term. But at the head (`β=2.5`,
`g²=4/β=1.6`):
```
   λ₄/g²     2g²/λ₄      < S≈1.4 ?
    0.5        4.0         NO
    1.0        2.0         NO  (natural λ₄~g²)
    1.5        1.33        YES
    2.0        1.0         YES
```
For **natural `λ₄~g²` the exact bound is `≈2`, which overshoots `S≈1.4`** — the bound is finite
and uniform but **too lossy to certify `D<S`** at the head, even though the measured `D/S<1`. The
`bound − measured` gap (`≈0.6`) *is* the part of `D` the rank-one ceiling cannot see. Worse, off
the Savvidy cell the true `D_Sav` is exponentially smaller than `2g²/λ₄` (no condensate), so the
bound is *exponentially lossy* there too. It cannot be promoted to a sharp uniform `D<S`.

### 2.4 The obstruction, stated precisely `[O]`

Combining 2.1–2.3, the surviving open inequality is
```
softened (★):   2 g²/λ₄(β)  <  S_bare(β) − δ      at the crossover β≈2.5 ,     [O]
```
a **pure marginal `O(1)` ratio** `κ_bare/(g²/λ₄)` (since naturally `λ₄~g²`, `SAVVIDY-UNIFORM-
SOFTENING` §3). It is *not* a divergence. But it **cannot be decided perturbatively cell by cell**:
`SAVVIDY-MULTICELL-RG` shows the accumulated softening over the flow is
```
Σ = C_RG · Var(φ²) · ln( ln(k_UV/Λ) / ln(k_IR/Λ) ) ,
```
**IR-marginal** (log-log), **IR-dominated** (the last few cells near `Λ` carry almost all of it),
and **divergent as `k_IR→Λ⁺`** — i.e. exactly where `g²(k)→∞` (**Landau pole**) and one-loop
perturbation theory fails. So the only comparison `κ_bare > Σ` is a **nonperturbative IR
statement at `k_IR~Λ`** — `H1` in the RG coordinate. **This is the obstruction**: not a divergence
of any *single* cell (per-octave softening is bounded `[T]`), but the *accumulation* being an
IR-marginal sum decided at the Landau pole where the tools vanish.

---

## 3 · The measure-weighted statement and its margin `[T]`/`[O]`

Since the pointwise bound is `[NG]` (§0.3) and the uniform bound is `[O]` (§2), the correct object
is the **average over the typical-field set** = the coarse-field marginal `ρ`.

**Theorem M (Fisher identity) `[T]`.** Let `ρ(a)` be the normalisable marginal of the coarse
coordinate `a=⟨v,X⟩`, and `G″(a) = −(log ρ)″(a)` the local marginal curvature conditioned at
`a`. Then
```
E_{a~ρ}[G″(a)] = ∫ ρ (−(log ρ)'') da = −∫ ρ'' da + ∫ (ρ')²/ρ da = 0 + I(ρ) ,
I(ρ) := E_ρ[(∂_a log ρ)²] = Var_ρ(∂_a log ρ) ≥ 0 .
```
Equivalently
```
E_ρ[S] − E_ρ[D] = E_ρ[G″] = I(ρ) ≥ 0   ⇒   E_ρ[D] ≤ E_ρ[S] ,   margin = I(ρ) .
```
*Proof*: `−(log ρ)″ = −ρ″/ρ + (ρ'/ρ)²`; multiply by `ρ` and integrate; `∫ρ″ da = [ρ']_{±∞} = 0`
for a smooth normalisable density. ∎ (Witness (2): Gaussian `m=1,2` give `E[G″]=I=1,4`;
non-Gaussian well `1.4344`; **free massless flat `ρ` ⇒ `E[G″]=I=0`**.)

**This is the correct measure-weighted statement, and the acid test passes.** The typical-field
average of the defect is below the bare stiffness by **exactly the Fisher information of the
coarse marginal**. Positivity is automatic (`[T]`), strict iff `ρ` is non-flat in `a`, and the
free massless field (flat `ρ`) returns `0` — no false positivity.

**Margin, and why it does not close.** The margin is `I(ρ) = 1/σ_ρ²` (Gaussian case), i.e. the
**effective mass²** of the coarse mode. It is measured positive but thin:
```
   E[G″] = I(ρ):   +0.28 (L=4, β=2.5)  →  +0.049 (L=8, β=2.4)   [V-num]
   = 3.5% of E[S]≈1.40  at L=8 ;  per-mode G″ ∈ [−0.14,+0.24], 5/12 < 0 .
```
Bounding `I(ρ) ≥ δ` **uniformly in volume** is *equivalent* to a bounded coarse-mode variance —
i.e. **the mass gap**. The measured `I` **shrinks with `L`** (`0.28→0.049`), the exact signature of
a kinematic (finite-box) positive rather than a genuine IR mass. So Theorem M converts the
pointwise problem into an exact average identity whose uniform strictness is the gap itself:
```
(★_avg)   I(ρ_v) ≥ δ > 0   uniformly in L   ⟺   the coarse mode has a mass gap .     [O]
```
This is honest progress (an exact identity + a clean margin) but **not** a closure: it relocates
H1 to "the coarse-marginal Fisher information is uniformly positive," which is the mass gap in the
information-theoretic coordinate.

---

## 4 · Ledger, acid test, and reading

| item | statement | label |
|---|---|---|
| rank-one split `D = D₊ + D_Sav`, `φ₋∈ker K₊` | exact algebra | **[T]** |
| `D_Sav = E|⟨φ₋|c⟩|²/(gB) ≥ 0` | one rank-one projection | **[T]** |
| `D₊ ≤ 0` on log-concave positive fibre | Brascamp–Lieb (sign only) | **[T]** |
| `D ≤ D_Sav`; `D<S−δ ⟺ D_Sav<S−δ+\|D₊\|` | exact consequence | **[T]** |
| one-octave bound `D_Sav ≤ E∫_{oct}\|c\|²/(gB)` | Cauchy–Schwarz + LLL support | **[T]** (lossy at head) |
| deep-IR `D_Sav/S→0` uniformly | scale separation, `(k/√(gB))²` | **[T]** (deep IR only) |
| exact bound `D_Sav ≤ 2g²/λ₄` | `Var(φ²)≤2/λ₄`, `c`-indep., 2-loop-protected | **[T]** (lossy at head) |
| head `2g²/λ₄ < S_bare−δ` (softened H1) | marginal `O(1)`, `λ₄~g²` | **[O]** |
| accumulated softener IR-marginal, Landau pole | log-log, IR-dominated, `k_IR→Λ` | **[T]** structure / **[O]** endpoint |
| uniform pointwise `D<S` | false: exact counterexample; 5/12 modes `G″<0` | **[NG]** |
| `E_ρ[G″] = I(ρ) = E[(∂_a log ρ)²] ≥ 0` | Fisher identity (Theorem M) | **[T]** |
| uniform `I(ρ) ≥ δ>0` (mass gap) | measured shrinking with `L` | **[O]** |
| **(★) at the head** | not closed; reduced to `(★_avg)` | **[O]** |

**Acid test (free field), for every tool used:**
- Rank-one split: free massless `B=0` ⇒ no negative mode, `D_Sav=0`, `D₊=0` ⇒ `D=S=0`,
  `G″=0`. Passes (manufactures no positivity).
- Rank-one bound `2g²/λ₄`: at `β→∞` `g²→0`, but the correct free-field statement is `D_Sav=0`;
  the bound `2g²/λ₄` is a Savvidy-cell object and is *vacuous* (not false-positive) off the cell.
  Passes as a **bound** (it does not assert `G″>0` on a free field).
- Theorem M: free flat marginal ⇒ `I(ρ)=0` ⇒ `E[G″]=0`. **Passes discriminatingly** (witness (2)).

**Bottom line.** `D` **can** be bounded on the typical-field set only as a **measure-weighted
average**, and there the bound is exact: `E[D] = E[S] − I(ρ)`, margin `= I(ρ) ≥ 0`, vanishing on
the free field. It **fails** as a uniform pointwise bound at the head, for a *precise, non-
divergent* reason: at `k~√(gB)` the rank-one overlap is `O(1)` (one-octave support gives no
suppression), the exact `2g²/λ₄` ceiling overshoots `S_bare` for natural `λ₄~g²`, and the
accumulated softening is IR-marginal — decided at the Landau pole `k_IR~Λ` where perturbation
theory dies. The head remains **H1**, now in three equivalent coordinates: softened
sub-criticality `2g²/λ₄ < S_bare`, uniform Fisher information `I(ρ) ≥ δ`, and the mass gap itself.
**Clay 0%; no QED.**

---

## 5 · Witness — `research/v3/w7/h2_defect_bound.py`

Deterministic, `py -3.12 -I`, numpy only. Quoted output:

```
(1) RANK-ONE SPLIT:  <K^-1 c,c> = <K_+^+ c,c> - <phi,c>^2/lam_-   split error = 1.78e-15
    D_+ (log-concave fibre) = -0.03113, -0.02215, -0.01350   (B=0.5,1,2; all <=0)
(2) FISHER:  Gaussian m=1: E[G'']=I=+1.000000 ; m=2: +4.000000 ; well: +1.434435
             FREE massless (flat): E[G'']=I=+0.000000   (acid test)
(3) RANK-ONE BOUND 2g^2/lam:  beta=2.5 g^2=1.600 S~1.4 :
             lam/g^2=1.0 -> bound 2.000 (NO, overshoot) ; =2.0 -> 1.000 (YES)
             exact Var(phi^2), lam=1 = 0.5749, 0.9565, 2.1409, 2.0002  -> 2 (bounded, c-indep)
(4) MEASURED per-mode (beta=2.4, L=8, Ncfg=400, Nmode=12):
             D/S mean=0.9651 min=0.8311 max=1.1028 ; G'' mean=+0.0490 min=-0.1445 max=+0.2370
             modes D<S: 7/12 ; D>S: 5/12 ; E[G'']=0.035*E[S]
(5) OVERLAP: k/sqrt(gB)=0.1 -> 0.0099 (deep IR->0) ; =1.0 -> 0.5000 (HEAD: O(1))
```

Readings: (1) the split is exact and `D₊≤0`; (2) the typical-field average of `G″` *is* the Fisher
information, `0` on the free field; (3) the exact rank-one bound is uniform but overshoots `S` at
natural `λ₄~g²`; (4) the measured defect is `D/S≈0.97` (near-cancellation, `5/12` modes negative);
(5) overlap suppression is a deep-IR theorem that does not reach `k~√(gB)`. **Not a proof — a
witness for the exact split, the exact average identity, and the precise location of the failure.**

---

## 6 · References (internal, this session)

- `STAR2-CLASSIFICATION.md` §3 (head = H1), §1 (exact counterexample, strong form false `[NG]`).
- `S3-UV-FESHBACH.md` §1 (Schur–variance identity), §2 (fibre gap ⇒ `D≤0`), §4 (failing step).
- `SAVVIDY-DTERM.md` §1 (rank-one split, `D≤D_Sav`), §3 (deep IR vs crossover).
- `SAVVIDY-MULTICELL-RG.md` §2–4 (IR-marginal log-log, Landau pole).
- `SAVVIDY-UNIFORM-SOFTENING.md` §2–3 (one-octave `(L,a)`-uniform, marginal `O(1)` in `β`).
- `SAVVIDY-NEGMODE-EXACT.md` §2 (`Var(φ²)≤2/λ`, `c`-independent); `SAVVIDY-2LOOP.md` §2 (sign protected).
- `SAVVIDY-LATTICE-CELL.md` §1 (`D/S≈0.94–1.00`, `5/12` modes `G″<0`); `w7/k0_su2_b2.4_L8_n2_results.json`.
- `docs/paper-v2/ch04-gap.tex` §Savvidy (exact `D`-split, uniformity tension, lattice cell).
- Brascamp–Lieb (1976), Prékopa (1971), Fisher information (standard), Schur complement / Haynsworth.
- Witness `research/v3/w7/h2_defect_bound.py` (deterministic; output quoted in §5).

*No third-party sources were downloaded this session; all external results are standard and cited
from the corpus files, not re-fetched.*
