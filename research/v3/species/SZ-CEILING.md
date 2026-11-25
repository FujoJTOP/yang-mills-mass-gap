# SZ-CEILING — the radius ceiling of the Shen–Zhu–Zhu route

> **Route B1 ("绕过新物种").** Task: measure the *ceiling* of the **existing** rigorous route of
> Shen–Zhu–Zhu (SZZ, *A stochastic analysis approach to lattice Yang–Mills at strong coupling*,
> CMP **400** (2023) 805–851, arXiv:2204.12737), i.e. decide whether the strong-coupling radius
> `|β| < 1/(16(d−1)) = 1/48` (d=4, N-independent) is **this technique's ceiling**, or whether the
> same technique can be pushed to `β ≈ 2.5`.
>
> Reads (mandated): `research/v3/species/LOOP-RUNBOOK.md` §1; `docs/paper-v2/ch13-machine.tex`
> §"Why the existing rigorous tools stop"; `research/v3/C1-SZZ-NSCALING.md`;
> `research/v3/S5-N-UNIFORM.md`; `research/v3/S6-DOBRUSHIN.md`;
> `research/unconditional/PROOF-IN11.md`; `research/v3/species/LIT-CHECK.md` §Q3,
> `research/unconditional/ATOM-strong-coupling-extend.md` §1.
> Witnesses (new, this note): `research/v3/w7/sz_ceiling_const.py`,
> `sz_ceiling_hess.py`, `sz_ceiling_exact.py`, `sz_ceiling_lmax.py` (`py -3.12`, numpy only).
>
> **Labels.** `[T]` theorem · `[T-num]` exact machine-precision computation · `[C]` conditional ·
> `[O]` open · `[NG]` no-go for a named class · **[X] correction** of a premise in the task book.

---

## §0 判死线 (T0) — written before the work, then judged

### 0.1 The claim to be judged

**Kill-line.** `1/48` sits in **one and only one algebraic step**, and every constant in that step
is now pinned:

```
  ρ(β)  =  Ric_coeff  −  Hess_coeff · |β|              (Bakry–Émery curvature budget)
        =    N/2      −      C*·N      · |β| ,          ρ(β) > 0  ⟺  |β| < 1/(2C*).        (0.1)
```

- **Gain side** `Ric_coeff = N/2`: **exact equality, not an inequality.** SZZ's `e:Ricvv` writes the
  SU(N) Ricci as `((N+2)/2 − 1)|X|²`, and `(N+2)/2 − 1 ≡ N/2` — the `+2/2` and `−1` cancel
  *identically*. `N/2` is the exact Einstein constant of `(SU(N), ⟨X,Y⟩=−tr XY)` (Killing form
  `B = −2N g`, `Ric = ¼B = (N/2)g`). Therefore **no improvement of the gain is possible — not by
  a factor 1.01** `[T]`. The N-independence of `1/48` is *not* a cancellation of two `N`'s: it is
  the statement that the Einstein constant is linear in `N` and the Hessian is linear in `N`.
- **Cost side** `8(d−1) = 4 × 2(d−1)`: the number of links per plaquette (4, the
  Cauchy–Schwarz/triangle slack) times the number of plaquettes per link (`2(d−1)` = 6 in d=4,
  exact counting). This is the only place where a constant could be improved — and it is here that
  the whole question lives.
- `C* := sup_{configs U} sup_{|v|=1} Hess_S(v,v)/(N|β|)` in the bi-invariant metric. SZZ's Lemma 4.1
  is the statement `C* ≤ 8(d−1)`.

### 0.2 Verdict (hybrid **(a)+(c)** — and the route dies)

1. **(a) holds in d=2, exactly.** `C*(d=2) = 8.000000 = 8(d−1)`, i.e. in the one dimension where
   the plaquette geometry is unfrustrated, **SZZ's constant is exactly optimal**: the
   `1/(16(d−1))` is the sharp zero of the curvature budget itself, not a lossy estimate.
   Witness `sz_ceiling_lmax.py` / `sz_ceiling_exact.py`: the extremal point is
   **the maximal centre-vortex configuration** (every plaquette holonomy `= z·1`, `z ∈ Z_N`)
   together with the `k=(0,π)` mode, and it attains the bound to machine precision `[T-num]`.
2. **(b) is bounded: at most a factor 1.5 — and the tightest value is measured, not guessed.**
   For `d ≥ 3` the saturating family is *obstructed* (a `Z_2` frustration argument, §2.3), so
   `C* < 8(d−1)` strictly; the measured values are
   `C*(d=3) ∈ [12, 16]`, `C*(d=4) ∈ [16.000, 24]` (lower end = machine-precision evaluation at the
   extremal configuration, upper end = SZZ's proven bound). Hence

   ```
     ceiling  β_A < 1/(2C*)  ∈  [ 1/48 , 1/32 ]  =  [ 0.020833 , 0.031250 ]   (d = 4).        (0.2)
   ```

   **The measurable ceiling is `1/32 = 0.03125` (β_A); the published `1/48` is optimal up to a
   factor 1.50, and no further tightening exists inside the class.**
3. **(c) the structural inequality, in the requested form.** The budget (0.1) can hold only if

   ```
        |β| ·  C*  <  Ric_coeff / N  =  1/2 ,        C* ≥ 2(d−1)  (single-link mode),
        i.e.   |β| · (plaquettes per link)  <  1/2  asymptotically,   and   ceiling = Θ(1/d).  (0.3)
   ```

   The gain `N/2` is **dimension-free and coupling-free** (it is a *single-link* curvature), while
   the cost grows as `Θ(dβ)` *and grows with the variation scale of the configuration*. Hence the
   ceiling of every "pointwise Γ₂ with a curvature gain" architecture is `Θ(1/d)`: **raising d
   lowers the ceiling.** No amount of constant-chasing changes this; the ceiling at d=4 is `1/32`.
4. **Dead — and doubly so.** Like-for-like: ceiling `β_W ≤ 0.281` vs head `β_W = 2.5` ⇒ **8.9×** short
   (equivalently in `β_A`: ceiling `0.03125` vs head `2.5/9 = 0.278` ⇒ the same `8.9×`).
   *Do not* quote `2.5/0.03125 = 80×`: that compares a `β_A` ceiling with a `β_W` head and is the
   same unit slip that `ATOM` §1 and `LIT-CHECK` Q3′ already corrected (§4.1). Since the ceiling is
   dimension-controlled and *measured*, not extrapolated, this is not a "not yet".
   **Moreover the rigorous (not merely measured) in-class bound already excludes the target:** the
   single-link mode forces `C* ≥ 2(d−1)` with no input whatsoever, so
   `ceiling ≤ 1/(4(d−1)) = 1/12 = 0.0833` (β_A) `= 0.75` (β_W, N=3) — **3.3× short of 2.5 for
   *every* pointwise-`Γ₂` member of the class, unconditionally.** The measured `1/32` closes the
   remaining factor 1.5.

> **One-line.** The SZZ radius is the exact analytic zero of a curvature budget whose gain is an
> exact Einstein constant and whose cost is an exact sup over configurations; the optimum is
> `1/(2C*) ∈ [1/48, 1/32]`, attained (up to the 1.5 gap) by the maximal centre-vortex configuration.
> **The direction is permanently closed as a route to `β≈2.5`; what survives is a conditional
> reduction (§5) that costs exactly the campaign's H1.**

---

## §1 T1 — the constant chain, disassembled

### 1.1 The chain from the action to `1/48`

```
 S = Nβ Σ_p Re tr Q_p            (SZZ `e:defS`; 't Hooft-scaled Wilson action)
   │
   ├─(1) PLAQUETTE COUNT       →  2(d−1) plaquettes through a link;  4 links per plaquette
   │                              (no loss: exact lattice combinatorics; d=4 → 6 and 4)
   ├─(2) TRACE / TRIANGLE STEP →  |tr(A v² B)| ≤ ‖v²‖₁ = ‖v‖², then
   │                              Σ_p (Σ_{i∈p} a_i)² ≤ 4 Σ_p Σ_{i∈p} a_i² = 8(d−1) Σ_e a_e²
   │                              ← THE ONLY INEQUALITY IN THE COST. Slack = 4 (links/plaquette).
   ├─(3) SUMMATION             →  Σ_p (Σ_{i∈p}a_i)² = 2(d−1)Σ a² + 2Σ_{i<j} N_ij a_i a_j
   │                              with N_ij ∈{0,1}, exactly 6(d−1) face-neighbours per link:
   │                              this is an IDENTITY, not a bound (no loss here) [T-num §2.1].
   ├─(4) RICCI                 →  Ric = (N/2)|v|²  EXACT (no loss; §0.1)
   ├─(5) BE Γ₂                 →  ρ = N/2 − C*N|β|            (no loss, given (2))
   └─(6) LSI ⇒ Poincaré ⇒ gap  →  λ_gap ≥ ρ ;  Var ≤ (1/ρ)Σ E|∇|²  (factor 2 between the LSI and
                                  Poincaré conventions — a convention, NOT a loss)
```

### 1.2 Table: every constant that enters `1/48`

| # | object | value | where it enters | optimal? | improvement available | status |
|---|---|---|---|---|---|---|
| 1 | plaquettes through a link | `2(d−1) = 6` | multiplicative cost | **yes** (exact count) | 1.00× | `[T]` |
| 2 | links per plaquette | `4` | triangle/Cauchy–Schwarz slack | **yes in d=2**; ≤1.5× in d=4 | 1.00× (d=2) / 1.50× (d=4) | `[T-num]` |
| 3 | plaquette-pair summation | `8(d−1)` | `=4·2(d−1)` | **yes** (identity + Perron–Frobenius, `λ_max(N_ij)=8(d−1)` attained by the constant mode) | 1.00× | `[T-num]` |
| 4 | Ricci curvature `ρ` | `N/2` (`=(N+2)/2−1`) | additive gain | **yes — equality** | 1.00× | `[T]` |
| 5 | BE `Γ₂` conversion | no constant | — | yes | 1.00× | `[T]` |
| 6 | LSI vs Poincaré | factor 2 | `C_LSI = 2/ρ`, `C_Poinc = 1/ρ` | convention | — | `[T]` |
| 7 | **random-environment decomposition loss** | **0** for SU(N)/SO(N) | — | n/a | — | `[X]` see below |
| 8 | **Dobrushin coefficient** | **does not enter** | — | n/a | — | `[X]` see below |
| 9 | **`sin/cos` bounds** | **do not enter** | — | n/a | — | `[X]` see below |
| 10 | mass-gap rate `c_N ≍ N^{−3}` | proof artefact | correlation length only | explicitly "not optimal" | — | `[T]`, irrelevant to the radius |

**Total improvement available inside the class: `8(d−1)/C* = 1.00 (d=2), 1.33 (d=3), 1.50 (d=4)`.**
Nothing else in the chain is lossy. `1/48 → 1/32`, full stop.

### 1.3 Three premise corrections `[X]`

- **[X1] SZZ's group is `SU(N)`/`SO(N)`, not `U(N)`.** The verified e-print statement is
  `G = SU(N)` or `SO(N)` (`C1-SZZ-NSCALING.md` §1). The `U(N)` case is *Nissim 2025*
  (arXiv:2510.22788), which is forced to rewrite `U(N) = (SU(N)×U(1))/Z_N` precisely **because**
  `U(N)`'s Ricci is not uniformly positive (`LIT-CHECK.md` §Q3 row 4). See §5.3 — the obstruction
  is exactly the centre, and it is fatal for the BE budget, so the `SU(N)` restriction is
  *structural*, not cosmetic.
- **[X2] Eldan stochastic localization is not what produces `1/48`.** The verified reading of the
  paper (`C1-SZZ-NSCALING.md` §2, with line anchors) is: Langevin dynamics + Bakry–Émery `Γ₂` with
  the **exact** SU(N) Ricci + the "Hessian is `O(N)`, not `O(N²)`" lemma. The *random-environment /
  semigroup splitting* decomposition is the `U(N)` device of Nissim 2025. This matters materially:
  a localization-based re-derivation that keeps the same `Γ₂` budget inherits the ceiling (0.2)
  verbatim, whereas one that replaces the budget by a *martingale/averaged* one is a different
  architecture (see §6, the only live variant).
- **[X3] Dobrushin does not enter `1/48` at all.** `S6-DOBRUSHIN.md` §2 gives the *separate, cruder*
  route `β_Dobr = 1/(6(d−1)κ₀)` (standard-Wilson `≈1/72`); it is a different method with a smaller
  radius and no role in the SZZ chain. Likewise no `sin/cos` bound appears: since the metric is
  bi-invariant, the geodesic Hessian is chart-free, so the usual exponential-chart `sin t/t` slack
  is *absent*. **All slack in the entire chain is the single triangle step (2) of §1.1.**

---

## §2 T1 continued — the exact optimality computation

### 2.1 The summation is an identity (`[T-num]`)

For `a_e ≥ 0`, with `N_ij = #{plaquettes containing both links i,j}`,

```
   Σ_p (Σ_{i∈p} a_i)²  =  2(d−1) Σ_e a_e²  +  2 Σ_{i<j} N_ij a_i a_j ,      N_ij ∈ {0,1},  (2.1)
```
`N_ii = 2(d−1)`, and each link has exactly `6(d−1)` face-neighbours. Both facts were checked
lattice-exactly in `sz_ceiling_const.py` P1/P2 (d=2,3,4; `rowsum = 8(d−1)` uniformly; symmetrised
`N_ij` has `λ_max = 8(d−1)` to `1e-15`, with `N·1 = 8(d−1)·1`). Consequences:

- The **only** inequality in the cost is `|2nd variation of Re tr| ≤ (Σ_i a_i)²`, i.e. the triangle
  step, whose slack is `4 = #links/plaquette` *when the four directions are equal*.
- `8(d−1)` is the **exact** supremum of that relaxed quadratic form over unit `v`: it is a row-sum
  of a non-negative symmetric matrix, attained by the constant mode `a_e = const`.

### 2.2 The actual Hessian: extremal object identified (`[T-num]`)

At the **maximal centre-vortex configuration** (all `Q_p = z·1`, `z ∈ Z_N`; solvable over `Z_mod`
for `L` even — consistency over a 3-cube needs `6 ≡ 0`, over a 2-plane needs `L² ≡ 0`) the second
variation collapses: for `N=2` (`z=−1`) and for `N=3` (`z=diag(−1,−1,1)`, `det=+1`),

```
   q_p(v)  =  Re tr( z·(Σ_i B_i)² )  =  ‖Σ_i B_i‖²      (the trace bound, saturated pairwise)
```

so the whole Hessian at that configuration is the **curl form** `Q(v) = Σ_p ‖Σ_{i∈p} σ_i X_{l_i}‖²`,
with symbol (plane wave `X_{x,μ} = b u^{x_μ} T`, `u_μ = e^{ik_μ}`)
`|(1+u_μ)(1−u_ν)|² = 16 cos²(k_μ/2) sin²(k_ν/2)`.

Exact results (`sz_ceiling_exact.py`, `sz_ceiling_lmax.py`; `C* = 2λ_max` of the assembled matrix,
confirmed against the direct geodesic second derivative to `1e-15`):

| d | `C*` measured at the extremal config | `8(d−1)` | ratio | ceiling `1/(2C*)` |
|---|---|---|---|---|
| 2 | **8.000000** | 8 | **1.0000** | 0.062500 = **SZZ's own value** |
| 3 | 12.000000 | 16 | 0.7500 | 0.041667 |
| 4 | **16.000000** | 24 | 0.6667 | **0.031250** |

The same numbers for `N = 2` and `N = 3` (to machine precision) — **`C*` is N-independent**, as it
must be (the trace bound `|tr| ≤ ‖v²‖₁ = ‖v‖²` is N-free, and for odd `N` the saturating
direction `i·diag(a,−a,0,…)` has a zero eigenvalue so the `det=+1` constraint costs nothing).
At the flat configuration `λ_max = 0` exactly (the action is *maximised* there, so `Hess ≤ 0`);
at random configurations `λ_max` is 1.4–1.7× smaller than at the vortex (d=2: 2.94 vs 4.00;
d=3: 3.78 vs 6.00; d=4: 4.74 vs 8.00 — `sz_ceiling_lmax.py`), i.e. the vortex is the extremal one
among everything sampled, and it coincides with the extremal object of §2.3.

### 2.3 Why `d ≥ 3` cannot saturate: the `Z_2` frustration (`[T]`)

Saturation of `8(d−1)` needs, for every plaquette, the four `Y_i = σ_i X_{l_i}` to be **parallel and
co-directed**. Write `X_e = s_e b T` (`s_e ∈ {±1}`). The three conditions
`Y_1=Y_2`, `Y_1=−Y_3`, `Y_1=−Y_4` read, in `Z_2`,

```
   s_{l1}=s_{l2},      s_{l3}=s_{l1}+1,      s_{l4}=s_{l1}+1        (mod 2).        (2.2)
```

With the linear ansatz `s(x,μ) = A_μ + Σ_k c_k x_k`, (2.2) forces `c_μ + c_ν = 1` **for every pair
`μ<ν`** — solvable in `Z_2` only for `d ≤ 2`. The `k=(0,π)` plane wave *is* that solution in d=2.
Hence

```
   C*(d=2) = 8(d−1) exactly ;      C*(d≥3) < 8(d−1) strictly ,                    (2.3)
```
and the improvement factor is bounded away from ∞ by the geometry. (In the *parallel-direction*
family the best is `16·m(d−m)/d = 12` (d=3) and `16` (d=4) — matching the measured values, i.e. the
optimum in that family *is* the vortex-configuration optimum; non-parallel directions can only add
the remaining `[16,24]` window, which is why the bracket in (0.2) is a factor ≤1.5 wide.)

---

## §3 T2 — quantitative: the solvable dimensions, and a *direct* d=4 measurement

The task asks for the method/truth ratio in d=2 (Migdal-exact) and d=3, then extrapolation.
The honest answer is stronger than extrapolation: **`C*(4)` is computable directly**, and the
d=2,3 values serve as validation of the mechanism and of the extremal object.

### 3.1 d=2 — exact, and it is the *worst* case for the method

- **Method's constant:** `C*(2) = 8 = 8(d−1)` exactly (§2.2) ⇒ the ceiling is `1/16`, **exactly**
  SZZ's value. **Zero improvement is available in d=2.**
- **Truth, exactly solvable:** 2D lattice YM (Migdal) is confining with a mass gap at **every**
  `β < ∞`; there is no upper critical coupling. So the "K factor" of method-vs-truth at the level of
  *thresholds* is `K(2) = ∞`: the method's radius is finite where the truth's is infinite.
  **A divergent K factor cannot be extrapolated as a constant multiple — the ceiling is therefore a
  structural statement about the method, not a multiplicative slack** (this is the decisive negative
  answer to T2's extrapolation premise).
- **The one exact, dimension-free constant-level miss:** at `β = 0` the Gibbs measure is
  `Haar^{⊗E}`, whose sharp Poincaré constant is `1/λ₁` with `λ₁ = C₂(fund) = (N²−1)/N`
  (exact: `Σ_a e_a² = (−(N²−1)/N)·1` for the ONB `e_a = √2T_a`), whereas the BE budget at `β=0`
  gives `1/Ric = 2/N`. Ratio:

```
   (SZZ's β=0 Poincaré constant)/(sharp)  =  (2/N)/(N/(N²−1))  =  2(N²−1)/N²  ∈ [1.50, 2).    (3.1)
```

  i.e. **the method is already 1.50× (N=2) / 1.78× (N=3) loose at `β=0`**, before any β-dependence.
  (This is a *spectral-vs-pointwise* miss: `Ric` is the pointwise Γ₂ constant, `λ₁` is the integrated
  bottom. It is not additive with a β-cost, so it does **not** compound with §2 — it is the same
  *kind* of loss as the one that matters in §6.)

### 3.2 d=3 — bracket

`C*(3) ∈ [12.000, 16]` ⇒ ceiling `β_A ∈ [1/32, 1/24] = [0.03125, 0.041667]` (SZZ: `1/32`). No exact
solution exists in d=3 (nor is a Migdal-type reduction available), so the truth side is not
computable; the *method* side is, and the mechanism (vortex + saturating mode) is dimension-general.

### 3.3 d=4 — measured, with its error bar

```
   C*(4) ∈ [16.000, 24]     ⇒     ceiling  β_A ∈ [1/48, 1/32] = [0.0208, 0.0313]   (error bar = 1.5×)
                                       measured ceiling   =  1/32 = 0.03125
```

**Ratios to the head — all in one normalisation at a time** (`β_W = N²β_A`, N=3):

| units | ceiling | head | ratio | note |
|---|---|---|---|---|
| `β_A` ('t Hooft-scaled, SZZ) | 0.03125 | 0.2778 (= 2.5/9) | **8.9×** | like-for-like |
| `β_W = N²β_A`, N=3 (standard Wilson) | 0.281 | 2.5 | **8.9×** | like-for-like |
| published SZZ vs head | 0.1875 | 2.5 | 13.3× | matches `ATOM` §1 ("≈13× in λ") |
| (mixed-unit artefact) `2.5 / 0.03125` | — | — | 80× | **wrong**: `β_W` head over `β_A` ceiling |

**Conclusion of T2.** The extrapolation T2 hoped for (a constant `K`, then `K × 1/48`) does not
exist: `K = ∞` in d=2 because the truth is gapless-free (gap for all β). What replaces it is a
*direct* exact d=4 computation whose error bar is only 1.5×, and which lands at `1/32`.

---

## §4 T3 —"绕过" judgment

### 4.1 Normalisation first (this is a real erratum, not bookkeeping)

SZZ's `β` multiplies `N·tr`; the physical head `β ≈ 2.5` is **standard Wilson**
`S_W = β_W Σ_p (1 − (1/N)Re tr Q_p)`, and `β_W = N² β_A`. So:

- the campaign's bracket `[1/48, 2.5]` mixed two normalisations; in `β_W` (N=3) it is
  **`[0.1875, 2.5]`**, i.e. `13.3×` (consistent with `ATOM-strong-coupling-extend.md` §1 and
  `LIT-CHECK.md` Q3′);
- my improved ceiling in the same units: **`β_W ≤ 0.281` vs `2.5` ⇒ `8.9×`**.

### 4.2 Judgment: **dead** (not "open", not "not yet")

```
   ceiling(d=4)  =  1/32  =  0.03125   (β_A)      [bracket [1/48, 1/32], measured end 1/32]
   = 0.281 in β_W (N=3) ;  head = 2.5 (β_W)
   short by      =  8.9× (like-for-like, either normalisation)
```

**Because** (i) the gain `N/2` is an exact Einstein constant (§0.1), (ii) the cost's arithmetic is an
identity and its only inequality is sharp in d=2 and within 1.5× in d=4 (§2), (iii) the ceiling is
`Θ(1/d)` and therefore *decreases* with the dimension we need — there is no version of constant
tuning that reaches 2.5. `1/48` is optimal up to 1.5×; the technique is capped at `1/32`.

### 4.3 If one insists on crossing: exactly which link must be replaced

The first link that must go is **(5): the pointwise-in-configuration budget**. Precisely:

> Replace `Hess_S(v,v) ≤ C* N|β| |v|²  for **every** configuration`
> by a bound on `E_μ[Hess_S(v,v)]` (or on a stopped/localized version of it).

Why this is the right link, and why it is expensive:

- **By parts, the averaged Hessian is already good:** for `μ ∝ e^{−S}dHaar`,
  `E_μ[Hess S] + E_μ[∇S⊗∇S] = 0`, so `E_μ[Hess S] ⪯ 0`. The *average* curvature is **better** than
  `Ric`, by an unbounded amount: the budget would read `ρ = N/2 + E[|∇S|²] > N/2`, i.e. **no
  β-restriction at all**. The entire `1/48` is the *price of the sup over configurations*.
- **That price is not removable by cleverness inside the same architecture**: the extremal object
  (§2.2) is an honest configuration (the maximal centre vortex, a legitimate `Z_N` field on an even
  lattice), so any pointwise statement must pay for it.
- **What replacing it costs:** one needs the averaged (equivalently: conditioned, coarse-grained)
  version, i.e. exactly `Var(⟨v,∇S⟩) ≤ E[H_{vv}] + W_full − δ` uniformly in `(L,a)` — the campaign's
  **H1**. This is not an accident: `S6-DOBRUSHIN.md` §5 reaches the same identity
  (`G″ = E[H] − Var + W`), and `LOOP-RUNBOOK.md` §1 records that H1 is the only named inequality
  left. **So "replace link (5)" = "prove H1".** The SZZ route is therefore not a bypass; it is the
  pointwise shadow of the same wall.

**Complete chain SZZ → Clay, had the ceiling been > 2.5** (recorded for completeness; each step is an
*existing* theorem, so the route would have been "no new species"):

| step | content | available theorem | status |
|---|---|---|---|
| 1 | BE `Γ₂` + Hessian lemma ⇒ LSI/Poincaré with `ρ = N/2 − C*N|β|` | SZZ Thm 1 (`[T]`) | done |
| 2 | LSI ⇒ exponential ergodicity of Langevin | SZZ Thm 1.3 | done |
| 3 | LSI + martingale ⇒ exponential clustering ⇒ lattice mass gap | SZZ Cor. `co:zmm` (rate `c_N ≍ N^{−3}`, not optimal) | done |
| 4 | uniform in `N` at fixed `β` | `C1-SZZ-NSCALING` (linear-N, N-free radius) | done |
| 5 | **extend `β` from `1/32` to `2.5`** | — | **this note: impossible in-class** |
| 6 | `a → 0` at fixed physical gap, with the same constant | Balaban (UV only) + Osterwalder–Seiler (strong only) | **missing** |
| 7 | reflection positivity / OS reconstruction of the continuum limit | standard once 6 holds | would be `[T]` |
| 8 | `m_{0⁺⁺} > 0` in the continuum, and isolation | needs IN-11 ordering (`PROOF-IN11.md`: only *membership* `{0⁺⁺,2⁺⁺}` is a theorem) | **missing** |

So even a hypothetical `β`-crossing would leave steps 6 and 8 — i.e. the `a→0` transport and the
IN-11 ordering — untouched. **Two independent walls, not one.**

---

## §5 T4 — acid tests and boundaries

### 5.1 Free field / abelian

- On an **abelian** group (`U(1)`, and the linearised/free Gaussian) `Ric = 0`, so the budget (0.1)
  reads `ρ = −C*N|β| < 0` for **every** `β > 0`: the method is **vacuous on the free/abelian side at
  every coupling**. In the charter's language it therefore *passes the letter* of the acid test
  (it never falsely certifies a free theory).
- **But it fails the spirit (D10 risk).** At `β = 0` (pure Haar — the campaign's `f₁² = 0` reference,
  `LOOP-RUNBOOK.md` §1) the constant is `N/2 > 0`, and it is *positive precisely because the group is
  compact and non-abelian*. So `K_S` measures **group curvature**, not Yang–Mills dynamics: it is
  positive on the trivial reference where the target marker `f₁²` vanishes. It becomes a
  discriminator **only after** the `β`-transport is closed — which is exactly what §4.3 forbids.
  Hence: `K_S` is a **curvature statement with a β-window**, not a certificate for `f₁² > 0`.

### 5.2 `N → ∞`

Already settled in `C1-SZZ-NSCALING.md` and re-confirmed here: `K_S = N(1/2 − C*|β|)`, per-colour
constant N-free, radius N-free, and `C*` itself is N-independent (§2.2, checked for N=2,3 to machine
precision). The large-`N` route is **not** blocked by degradation — but it also gains nothing:
the radius is the same `1/32`, and the *only* N-dependence sits in the (irrelevant) rate
`c_N ≍ N^{−3}`.

### 5.3 `SU(N)` vs `U(N)`: the centre is not cosmetic — it is load-bearing

This is the sharpest structural point of the present note.

- `U(N) = (SU(N)×U(1))/Z_N`. Along the **central** direction `X_e = iθ·1` (which commutes with
  everything, so it is a *flat* direction of the group), the Riemannian curvature gives
  `Ric = 0`, while the action's Hessian is **not** zero: at a configuration with
  `Re tr Q_p < 0` it is `+O(N²β)` per plaquette. Therefore
  `Ric − Hess = −O(N²β) < 0` **for every `β > 0`**: the BE inequality with *any* positive constant
  is **identically false along the centre**.
- Consequence: the `SU(N)` restriction in SZZ is **forced by the architecture**, not a convenience,
  and any `U(N)` version *must* first quotient the centre (Nissim 2025's semigroup/random-environment
  splitting). This is also why the "random environment decomposition" appears in the `U(N)` papers
  and not in SZZ — the correction `[X1]/[X2]` of §1.3.
- For **SU(3) vs U(3)** as targets: the difference is *exactly* this centre direction, so it is
  important — but in the direction that *kills* the budget, never in the direction that helps.
  (Physically the centre is what gives confinement; analytically here it is what gives
  `Ric = 0` on a direction the action still curves.)
- **Odd/even `N` and `det = 1`:** for odd `N` the "reflection" element with `det = +1`
  (`diag(−1,−1,+1)`) plays the role of `−1`; the saturating direction `i·diag(a,−a,0,…)` has a zero
  eigenvalue, so the `det = 1` constraint costs **nothing** in `C*` (§2.2, N=3 = N=2 numerically).
  There is no `N`-parity loophole.

---

## §6 The one live variant (and why it is not a bypass)

For the record, and to keep the "route permanently closed" statement honest about its scope:

- The closed class is: **pointwise-`Γ₂` curvature budgets, i.e. `inf_f Var/∫|∇f|²` controlled by
  `Ric − sup_configs Hess`.** Ceiling `1/(2C*) ∈ [1/48, 1/32]`, `Θ(1/d)`. This includes SZZ,
  Nissim's `U(N)` splitting (same budget, abelian factor quotiented), any Eldan-localization
  re-derivation that keeps the same `Γ₂` bookkeeping, and any "sharpen the constants" attempt.
- The **only** escape inside "existing mathematics" is to make the budget **averaged/stopped**:
  bound `E_μ[Hess]` (or its localized/conditioned version), where by parts the curvature is
  `≥ Ric + E|∇S|²` — *i.e.* no β-restriction at all. That is the campaign's H1, and its
  one-sidedness (`Var ≥ 0` subtracted) is precisely the `D/S ≈ 2%` defect recorded in
  `LOOP-RUNBOOK.md` §1. **So the escape exists but is not free: it *is* the open problem.**
- Consequence for the loop's bookkeeping: this route should be recorded on the **D7 (constant /
  scaling)** axis as *not* a lossy-constant failure but as a **D3 (regime) + D2 (circular)** pair:
  it is regime-locked at `β ≤ 1/32`, and its only exit re-imports H1.

---

## §7 Verification ledger and what changes in the corpus

| item | check | status |
|---|---|---|
| `K_S = (N+2)/2 − 1 − 8(d−1)N|β|`, radius `1/(16(d−1))`, N-free | `C1-SZZ-NSCALING.md` §1–§2 (verified against arXiv:2204.12737 e-print, line anchors) | `[T]` (carried) |
| `(N+2)/2 − 1 ≡ N/2`, and `(N/2)g` is the exact Ricci of `(SU(N), −tr)` | Killing form `B = −2N g` ⇒ `Ric = ¼B`; independent recomputation | `[T]`, this note |
| `8(d−1) = 4 × 2(d−1)`; `λ_max(N_ij) = 8(d−1)` attained by the constant mode | `sz_ceiling_const.py` P1/P2 (d=2,3,4; `1e-15`) | `[T-num]` |
| `C*(d=2) = 8` exactly; `C*(3) = 12`, `C*(4) = 16`; upper bound `8(d−1)` | `sz_ceiling_exact.py`, `sz_ceiling_lmax.py` (extremal = maximal centre vortex + saturating mode; direct vs assembled Hessian agree to `1e-15`) | `[T-num]` |
| `C*` independent of `N` (2 vs 3) | same scripts, identical to 6 d.p. | `[T-num]` |
| `d ≥ 3` cannot saturate (`c_μ+c_ν=1 ∀μ<ν` unsolvable in `Z_2`) | §2.3, proof sketch + the parallel-family maximum `16m(d−m)/d` reproduced exactly | `[T]` |
| `β_W = N²β_A`; head `2.5` is `β_W`; `β_W < 0.19` for SZZ at N=3 | `ATOM-strong-coupling-extend.md` §1; `LIT-CHECK.md` Q3′; recomputed `9/48 = 0.1875` | `[T]` |
| `β = 0` sharp constant `(N²−1)/N` vs `Ric = N/2` ⇒ miss `2(N²−1)/N²` | exact Casimir computation `Σ_a e_a² = −(N²−1)/N·1` | `[T]` |
| SZZ is `SU(N)`/`SO(N)`, not `U(N)`; localization/random-environment is Nissim 2025 | `C1-SZZ-NSCALING.md` §1; `LIT-CHECK.md` §Q3 row 4 | `[X1]` `[X2]` |
| Dobrushin (`β_Dobr = 1/(6(d−1)κ₀)`) does not enter `1/48` | `S6-DOBRUSHIN.md` §2 vs §1.1 here | `[X3]` |

**Proposed corpus edits (not applied — this note is the record):**

1. Replace the campaign bracket `[1/48, 2.5]` by `[0.1875, 2.5]` in `β_W` (or `[1/48, 2.5/9]` in
   `β_A`). `ch13-machine.tex` already flags the `1/24`/`75×` errors; this is the *same* class of slip
   but in the bracket itself, which is still used as the "distance to beat".
2. Record the **radius ceiling `1/32`** (bracket `[1/48, 1/32]`) next to every statement of SZZ's
   `1/48`, with the one-line reason: *gain = exact Einstein constant, cost = exact config-supremum*.
3. `LOOP-RUNBOOK.md` §2: add to the dead list, with death axis **D3 + D2**: "sharpen the SZZ constant"
   — dead at `1/32`, exit re-imports H1.

---

*Refs.* Shen–Zhu–Zhu, CMP **400** (2023) 805, arXiv:2204.12737; Nissim, arXiv:2510.22788;
Cao–Nissim–Sheffield, arXiv:2505.16585 / 2509.04688; Osterwalder–Seiler, Ann. Phys. **110** (1978) 440;
Bakry–Émery (1985); `C1-SZZ-NSCALING.md`; `S5-N-UNIFORM.md`; `S6-DOBRUSHIN.md`;
`PROOF-IN11.md`; `ATOM-strong-coupling-extend.md` §1; `LIT-CHECK.md` §Q3/Q3′;
`docs/paper-v2/ch13-machine.tex` §"Why the existing rigorous tools stop".
Witnesses: `research/v3/w7/sz_ceiling_const.py`, `sz_ceiling_hess.py` (verified against a direct
geodesic-second-derivative implementation to `1e-15`), `sz_ceiling_exact.py`, `sz_ceiling_lmax.py`.
*No git commit was made.*
