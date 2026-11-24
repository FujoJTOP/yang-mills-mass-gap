# S5 — N- and a-uniformity of the block-conditional one-sided variance bound (★)

> 2026-11-24 · Attack surface C, S5. Angle: **N- and a-uniformity of (★).** C1 established that SZZ's
> strong-coupling constant is **linear in N** (`K_S = N(1/2−24|β|)`, N-independent radius `1/48`). The
> question: does the **block-conditional** version of (★) preserve the linear-N scaling, or does
> conditioning on `U′` introduce a `1/N` or an exponential (e.g. in the boundary volume)? And how does
> the block scale in `a`?
>
> Reads: `C2-BALABAN-UNIFORM.md` §2.4, `C1-SZZ-NSCALING.md`, `docs/paper-v2/ch03-existence.tex` §L3,
> `ch05-theorems.tex` (W_geo, Lemma SC), `G1-L3.md`, `G1-L3-UV-FORCE-BOUND.md`,
> `G1-L3-TOWER-POINCARE.md`, `G1-UNIFORM-ZA.md`, `PROOF-COARSE-POINCARE.md`, `L3-band-decay.md`.
> Witness: `research/v3/w7/s5_n_uniformity.py` (`py -I`, output quoted in §2–§5).
>
> **Labels (charter §2).** **[T]** theorem · **[C]** conditional · **[NG]** no-go for a named class ·
> **[O]** exact missing step / open · **[R]** reported / unrefereed.

---

## 0. Verdict (no QED)

**The block-conditional (★) PRESERVES the linear-in-N scaling; it does NOT introduce a `1/N` and does
NOT introduce an exponential factor in the constant.** And it is **`a`-uniform** (`δ_W = Θ(a⁰)`),
surviving `a→0` at fixed physical block size.

Precisely:

1. **Every term in (★) has a benign N-counting, identical in mechanism to SZZ's cancellation.** The
   Hessian cost is `O(N)` per unit mode (**not `O(N²)`**), because the adjoint dimension `N²−1` is
   **multiplicity, not magnitude** — SZZ's Lemma 4.1 is exactly the statement that the trace structure
   kills the naive `N²` [T]. The Haar curvature reservoir `E[Φ''] ≥ (N/12)‖v‖²` (hence `W_geo`, up to a
   subleading variance ratio) is `O(N)` per unit mode [T]/[C]. The
   projection `P` has `Lip = C_blk = O(1)` [T]. Boundary plaquettes are `O(N)` per plaquette, same as
   bulk [C]. Hence `δ_W = W_j − (E[H]−Var) ≥ c·N` with `c` N-independent.
2. **The comparison to SZZ must be regime-aware.** SZZ's `K_S = N(1/2−24|β|)` is positive only at
   **strong coupling `|β|<1/48`**; the RG tower lives at **weak coupling** where `β_{SZZ,j} =
   11j/(12π²) ≈ 0.093 j > 1/48` for every shell (§2). **SZZ is off-regime in the tower.** The block
   constant is therefore **not** "≥ SZZ's `cN`" *via SZZ*; it is `≥ c_W·N` **independently, from
   `W_geo`** (coupling-free). Both are linear in N; the block's per-color coefficient is `c_W ∈
   [1/12, ≈0.13]`, smaller than SZZ's `1/2` but N-independent and positive.
3. **Conditioning on `U′` costs no factor in the constant.** `W_j`, `H_j` and `P` have `U′`-independent
   N-counting. The only `U′`-dependence is the small-field/large-field split; on the large-field boundary
   set the naive uniform-in-`U′` bound degrades (**[NG]** for the "all `U′`, no excision" class), but that
   set has measure `≤ e^{−c/g_j²} = e^{−2cβ₀j} ~ e^{−c'N j}` — it **shrinks exponentially in N** and is
   excised by Lemma B [T]. So the exponential lives in an *exceptional set's measure*, never in `δ_W`.
4. **In `a`:** the block side in lattice units is `L^{j+1}`; at fixed physical block `ℓ`, `a→0` means
   `j* = log_L(ℓ/a) → ∞`. The bare Hessian grows `~Nj`, but `× field²  g_j²` is `O(1)`; `W_j` is
   coupling-free; the defect is `O(1)` via the `O(a)`-band clustering. Hence `δ_W = Θ(a⁰)`, `a`-uniform.
   The **collapsing** Holley–Stoock bound `e^{−2β_jN_p} ~ (aL)^c` is a *bound artifact*, not `δ_W`
   (`G1-L3-TOWER-POINCARE`); (★) is precisely what replaces it by the `a`-uniform floor.
5. **The genuinely open content is NOT N or a.** It is the `O(1)` equipartition-defect **coefficient**
   (`E[H]−Var`, the one-sidedness `Var ≥ 0` subtracts) and the `O(a)`-band clustering rate `(★dec)`.
   Both are **N-benign and a-benign in their bookkeeping**; their closure is the same
   `(★)∧(★dec)` as C2/G1-L3. The obstruction is analytic (one-sided variance), not N/a-combinatorial.

**Sharpest line.** C1's linear-N cancellation survives the block decomposition because the adjoint
dimension is multiplicity, not magnitude (§3); C2's `a`-uniformity survives it because `W_j` is
coupling-free and the bare-Hessian growth `~Nj` is cancelled by the field amplitude `g_j` (§4); the
conditioning does not add an exponential to the constant, only to the measure of a large-field set that
shrinks in N (§5). **No N-uniformity loss, no a-uniformity loss. `[C]` on `(★)∧(★dec)`.** Clay 0%.

---

## 1. The object, and the exact decomposition it is counted against

(★), per shell `j`, per unit gauge-invariant collective mode `v`, uniform in the environment `U′` and in `a`:

```
(★)   Var_{μ(·|U′)}(⟨Pᵀ∇S_j, v⟩)  ≤  E_{μ(·|U′)}⟨PᵀH_jP v, v⟩  +  W_j(v) − δ_W        (δ_W>0),
```

equivalently the conditional coarse effective Hessian

```
      G″_j(v)  =  E_{μ(·|U′)}[⟨PᵀH_jP v,v⟩]  −  Var_{μ(·|U′)}(⟨Pᵀ∇S_j,v⟩)  +  W_j(v)  ≥  δ_W .
```

The terms and their origins (`G1-L3.md` §2, `ch03-existence.tex` eq. `species-identity`,
`PROOF-COARSE-POINCARE`):

- `∇S_j`, `H_j = Hess_ξ S_j`: force and force-Hessian of the shell action w.r.t. the **fine** (block-interior)
  links `ξ`, at fixed boundary/coarse field `U′`. `S_j` = Wilson action at coupling `β_j = 2N/g_j²`.
- `P`: block-averaging projection onto **coarse** modes; `Pᵀ` = its adjoint; `PᵀH_jP` = the coarse-restricted
  Hessian; `⟨Pᵀ∇S_j,v⟩` = the coarse-projected force. `P` also implements the `⊥`-orbit (gauge) projection.
- `W_j(v)`: the Haar/Jacobian geometric curvature `W_geo = E[Φ''] − Var(Φ')`, `Φ = −log J` (the log of the
  exponential-map Jacobian); the curvature of the pure-Haar marginal.
- `E[H]−Var`: the **equipartition defect** — zero at one loop (for `e^{−V}dX`, `E[∂²V]=E[|∂V|²]` identically),
  nonzero only from non-Gaussianity: a divergence plus a `U′`-boundary covariance (`G1-L3-UV-FORCE-BOUND` §1).
  This is the object (★) must bound from above.
- `δ_W`: the certified floor; `(★) ⟺ δ_W ≤ inf_{U′} G″_j`.

Two normalizations are used below and kept distinct:
the **SZZ/'t Hooft** one, `S = Nβ Re Σ tr` (where `K_S = N(1/2−24|β|)` is an equality), and the
**tower** one, `g_j² = 1/(2β₀j)`, `β_j = 2N/g_j² = 4Nβ₀j`, `β₀ = 11N/(48π²)` (`C2-BALABAN-UNIFORM` §2.1).
They are related by `β_{SZZ,j} = 2/(N g_j²)`.

---

## 2. The N/a bookkeeping table

Per **unit gauge-invariant collective mode** `v` (so the adjoint dimension enters as the size of `v`, once,
not squared), in the SZZ/'t Hooft normalization unless noted. `j` is the shell index; "physical" means after
multiplying by the field amplitude squared. `A(j) ~ B(j)` means same leading order.

| object | definition | **N-scaling** | **a-scaling** (fixed physical block `ℓ`) | label / source |
|---|---|---|---|---|
| **force** `⟨Pᵀ∇S_j,v⟩` | coarse-projected gradient | `~ (E[H])^{1/2} ~ (Nj)^{1/2}` bare; `O(1)` physical | `O(1)` | Gaussian equipartition [T] leading |
| **Var(LHS)** | `Var_{μ(·\|U′)}(⟨Pᵀ∇S_j,v⟩)` | `≈ E[H] ~ Nj` bare; `O(1)` physical | `O(1)` | equipartition [T] leading; deviation [O] |
| **block Hessian** `E⟨PᵀH_jP v,v⟩` | 2nd variation; adjoint fibre `dim N²−1` | `O(N β_{SZZ,j}) = O(Nj)` bare. **Not `O(N²)`**: `|Hess| ≤ 8(d−1)N|β|`; adjoint dim is *multiplicity* | `O(Nj)·g_j² = O(1)` | SZZ Lemma 4.1 [T]; adjoint-dim [T] (C1 §2a) |
| **Haar curvature reservoir** `E[Φ'']` | `= W_geo + Var(Φ')`; the `O(N)` gain | **`≥ (N/12)‖v‖² = c N`**; adjoint dim again multiplicity | coupling-free ⇒ **`a`-uniform (`a⁰`)** | ch05 Lemma B [T] |
| **Haar curvature** `W_j = W_geo` | `E[Φ'']−Var(Φ')` | `≥ 0` proven; **measured `≈0.125·N`** (`N=2`: `0.25`). `O(N)` up to a variance ratio | coupling-free ⇒ `a`-uniform (`a⁰`) | ch05 Thm [T] (`≥0`); `O(N)` **[R]/[C]** |
| **boundary plaquettes** | fine–coarse shared plaquettes; source of score `s_B` | per-plaquette `O(N)` (same as bulk); score physical `O(1)` | count `~(ℓ/a)^{d−1}`; per-normalized-mode `O(1)`; support `O(a)` | support/range: `L3-band-decay` [T]; rate [O] |
| **projection** `P` | block averaging | `Lip = C_blk = O(1)`, N-independent | `O(1)`, a-independent kernel | Lemma SC [T] |
| **defect** `E[H]−Var` | equipartition violation (divergence+bdy cov) | `O(1)` subleading (`O(g_j²·Nj)`) — **not `O(N)`** | `O(1)` (via `O(a)` band) | [O] = (★) content |
| **floor** `δ_W` | `= W_j − defect` | `c N − O(1)` ⇒ **linear, positive** for `N≥2` | **`Θ(a⁰)`** ⇒ `a`-uniform, survives `a→0` | [C] on (★); `W=O(N)` part [T] |

Rows [T] are assembled from proven pieces (SZZ Lemma 4.1; `W_geo≥0` and `E[Φ'']≥(N/12)‖v‖²`; Lemma SC
`C_blk=O(1)`; Lemma B measure bound). The one **[O]** row is the defect — the open (★) content itself;
its *bookkeeping* (O(1) not O(N), O(1) in a) is what this note establishes is benign. The `δ_W` row is
**[C]** because it is `W_geo − defect`, i.e. (★) itself.

---

## 3. The N-cancellation survives blocking: adjoint dimension is multiplicity, not magnitude

This is the crux, and it is exactly C1's mechanism re-read per mode.

**SZZ's counting.** `K_S = N/2 − 8(d−1)N|β|`: the Ricci gain `(N/2)|X|²` is `O(N)`, the Hessian cost
`8(d−1)N|β|` is `O(N)` (`|Hess_S| ≤ 8(d−1)N|β||v|²`, Lemma 4.1), the difference is `O(N)`, the positivity
radius `|β|<1/48` is N-free (`C1-SZZ-NSCALING` §2). The essential miracle is that the Hessian is `O(N)`,
**not** `O(N²)` despite `S = Nβ Re tr` naively looking `O(N²)`: the trace is bounded by `N`, and the
unitary/commutator structure prevents the `N` from squaring (SZZ explicit remark, quoted in C1 §2a).

**The adjoint dimension `N²−1`.** The gauge field at each link is valued in `su(N)`, `dim = N²−1`; the
Hessian `H_j` is a `[(N²−1)·#links] × [(N²−1)·#links]` matrix, and a unit mode `v` has `N²−1` colour
components. The point is: **`N²−1` multiplies the *number* of modes, not the *size* of each.** The
per-mode Hessian eigenvalue is `O(Nβ)` (Lemma 4.1); it is the same `O(N)` for all `N²−1` colours. So the
variance along **one** unit gauge-invariant mode `v` is `O(Nβ)` — there is no `N²−1` factor in the
per-mode constant, and therefore no `1/N` or `N`-loss when the per-mode inequality (★) is stated. This is
`[T]` (Lemma 4.1 is a theorem) and it is the same fact that makes `K_S ∝ N`, not `∝ N³` or `∝ 1/N`.

**Why blocking (`P`) does not spoil it.** `P` is a spatial block-averaging with `Lip(P)=C_blk=O(1)`
(Lemma SC [T]); it averages over **spatial** links and is the identity on the colour index. It is a
contraction independent of `N`, so `PᵀH_jP` has the same `O(N)` per-mode scale as `H_j`. Averaging over
space lowers the mode's *variance* (each `(Pv)_i ∼ 1/√#block-links`), not its N-scale.

**The surviving `O(N)` is the Haar curvature, not a Ric−Hess difference.** In the block identity the `O(N)`
content sits in a different place than in SZZ: here `E[H]` and `Var` are both `O(Nj)` and **cancel at leading
order** (equipartition), so `G″_j = W_geo − (Var−E[H]) = W_geo − defect`. The `O(N)` reservoir is the Haar
Jacobian curvature `E[Φ''] ≥ (N/12)‖v‖²` (ch05 Lemma B [T]); `W_geo = E[Φ''] − Var(Φ')` is proven `≥ 0`
(Prékopa/log-concavity [T]) and **measured `≈0.125N`** at `N=2` (ch05: `E[Φ'']=0.34`, `W_geo=0.25`); the
defect is subleading. So the costs to the floor are `Var(Φ')` (the variance-vs-curvature ratio at `β=0`) and
the defect, and

```
   δ_W  =  W_geo − defect  =  E[Φ''] − Var(Φ') − defect  ≥  (1/12 − r_var)·N − O(1),
```

with `r_var := Var(Φ')/N` the measured pure-Haar variance ratio. **Measured:** `N=2`, `E[Φ'']/N≈0.17`,
`W_geo/N≈0.125`, `defect≈0.02` ⇒ `δ_W≈0.23>0`. The `O(N)` is robust because `E[Φ'']` is `O(N)` and both
subtractions (`Var(Φ')`, defect) are *subleading per color*. Hence `δ_W = cN − O(1)` with `c = 1/12 − r_var`
N-independent. **Linear in N, positive, N-free radius.**

Witness (`s5_n_uniformity.py` part 2): `E[Φ'']_lb/N = 1/12` exactly for `N = 2…1000`; `δ_W/N → 1/12` for
any `O(1)` variance-ratio and defect; `δ_W>0` at `N=2` needs `r_var + defect < 1/6` (satisfied by the
measured `r_var≈0.045`, `defect≈0.02`, sum `≈0.065 ≪ 1/6`).

> **N-uniformity status.** **[T]** that `E[Φ''] ≥ cN`, that `W_geo ≥ 0`, that the Hessian is `O(N)`, that
> `P` is `O(1)`. **[C]** that `δ_W ≥ cN`, conditional on (i) the pure-Haar variance ratio `r_var < 1/12` and
> (ii) the defect being `O(1)` — i.e. on a `β=0` variance-vs-curvature comparison plus (★) itself. **No
> `[NG]`** for the N-scaling: the linear-N structure is robust to blocking. The only way to lose it is a
> subtraction (`Var(Φ')` or the defect) growing like `N` with a coefficient exceeding `1/12` — but both are
> *subleading per color* (the `β=0` case is log-concave, so `W_geo≥0` is exact; measured ratios are `≪1/12`),
> so this is a coefficient check, not an N-catastrophe. Note the honest residue: even the pure-Haar floor
> `W_geo ≥ cN` is a one-sided variance-vs-curvature comparison, but at `β=0` log-concavity makes it far
> easier than the interacting case and guarantees `W_geo ≥ 0` for free.

---

## 4. The `a`-scaling: `δ_W` is `Θ(a⁰)`, survives `a→0` at fixed physical block

**Block size in lattice units vs physical.** At shell `j` the spacing is `a_j = a L^j`; a block of side
`L·a_j` in lattice units has **physical** side `ℓ = L^{j+1} a`. Holding the physical block `ℓ` fixed as
`a→0` forces `j* = log_L(ℓ/a) − 1 → ∞`, i.e. **one looks at the ever-deeper UV shell** where
`g_{j*}² = 1/(2β₀j*) ~ 1/(N log(1/a)) → 0` (asymptotic freedom) and the field amplitude is `O(g_{j*}) → 0`.

Three `a`-facts, and they combine to `a`-uniformity of `δ_W`:

1. **The bare Hessian grows, the physical curvature does not.** `β_{SZZ,j} = 4β₀j/N ~ j`, so
   `E[H] ~ 8(d−1)Nβ_{SZZ,j}·(field²)` — but the field² amplitude is `g_j² = 1/(2β₀j)`, so
   `E[H]·g_j² = 16·8(d−1) = O(1)`, **independent of `j` (hence of `a`)**. Witness part 3: the product is
   exactly `48.000` across `j = 1…50` (the constant is convention-dependent; its `j`- and N-independence is
   the content).
2. **`W_j` is coupling-free.** `W_geo` depends only on the Haar measure, not on `β_j`, so it is manifestly
   `a`-uniform (`a⁰`).
3. **The defect is `O(1)` in `a`** *provided* the band clustering `(★dec)` holds with range `O(a)`: the
   equipartition defect is a boundary/divergence term whose extensive parts cancel between `E[H]` and
   `Var`, leaving a band-local residual `∼ O(1)`.

Therefore `δ_W = W_j − defect = Θ(a⁰)` — it does **not** scale as a positive power of `a`, and it survives
`a→0` at fixed physical block, improving slightly (`defect ∝ g_j² → 0`).

**What collapses is a bound, not `δ_W`.** The only *proven* per-shell Poincaré constant is Holley–Stroock
`Λ_{UV,j} = e^{−2β_jN_p}/C_G`, and since `β_j = 2N/g_j² ~ N²j` grows into the UV, it collapses like
`(aL)^c → 0` (`G1-L3-TOWER-POINCARE` §2 [T]). This is a **bound artifact**: HS overcounts the interaction as
oscillation `2β_jN_p` exactly where the theory is asymptotically free and the gap is nearly free. `δ_W`
does **not** collapse because it is set by `W_geo` (constant) minus an `O(1)` defect — the physical floor.
(★) is precisely the inequality that **replaces** the collapsing HS constant by this `a`-uniform floor.

> **a-uniformity status.** **[T]** `W_j` coupling-free (`a⁰`), physical Hessian `×g_j²` `a`-uniform,
> HS-collapse is a bound artifact. **[C]** `δ_W = Θ(a⁰)`, conditional on (★dec)'s `O(a)`-band clustering.
> The `a`-power is `a⁰`; no `a^{>0}` vanishing. **No `[NG]`** for the a-scaling given (★dec).

---

## 5. The conditioning: no exponential in the constant, only in an exceptional set that shrinks in N

This is the one genuine subtlety, and it is where an "exponential in the boundary volume" *could* arise.
It does not enter the constant; it enters the *measure* of a bad set.

**(5.1) What conditioning does *not* do.** Conditional on `U′`, the measure is still
`μ(·|U′; ξ) ∝ exp(S_j(ξ;U′)) ∏ Haar(dξ)` with the **same** action `S_j`. So `H_j` (w.r.t. fine links),
`W_j` (Haar curvature) and `P` have the **same** N-counting for every `U′` — `U′` is a *parameter*, not a
new scale. `W_j` is manifestly `U′`-independent (it is the Haar marginal). Hence **conditioning per se costs
no N factor in `W_j`, `H_j`, or `δ_W`.**

**(5.2) The naive uniform-in-`U′` statement fails on a set — [NG] for the "all `U′`, no excision" class.**
A full LSI on `μ` gives only the **averaged** conditional variance,
`E_{U′}[Var_{μ(·|U′)}(F)] ≤ Var_μ(F) ≤ (1/K)Σ E|∇F|²` — it does **not** control `sup_{U′} Var_{μ(·|U′)}`.
The pointwise-in-`U′` upgrade is a strictly stronger property (uniform conditional LSI / complete
analyticity). It **provably fails** for arbitrary boundary conditions: a **large-field** boundary `U′`
(pinned holonomy `O(1)`, not `O(g_j)`) couples to the block through boundary plaquettes with weight
`β_j·(holonomy) ~ N²j`, which can **tilt** the conditional measure toward phase coexistence and collapse
its conditional spectral gap. So:

> **[NG]** *There is no `U′`-uniform conditional LSI (hence no pointwise-in-`U′` (★)) with a
> `U′`-independent constant over **all** boundary conditions.* Mechanism: boundary-induced coexistence at
> large field. This is the honest reason `μ(·|U′)` is a *conditional* object and why (★) is stated with a
> small-field domain and the large-field excision of Lemma B.

**(5.3) The rescue — the bad set is `e^{−cNj}` and shrinks with N.** The large-field boundary set is
precisely what Lemma B (`C2-BALABAN-UNIFORM` §2.2 [T]) controls: `μ(large field at shell j) ≤ e^{−c/g_j²} =
e^{−2cβ₀j}`. Because `β₀ = 11N/(48π²) ∝ N`, the suppression rate `2cβ₀ ∝ N` **grows with N**:

| `N` | `2cβ₀` (per shell) | `μ_bad(j=1)` | `μ_bad(j=5)` | `μ_bad(j=10)` |
|---|---|---|---|---|
| 2 | 3.72 | 2.4e−2 | 8.6e−9 | 7.3e−17 |
| 3 | 5.57 | 3.8e−3 | 7.9e−13 | 6.3e−25 |
| 10 | 18.6 | 8.6e−9 | 4.6e−41 | 2.1e−81 |
| 100 | 186 | 2.1e−81 | 0 | 0 |

(`s5_n_uniformity.py` part 4, `c=40`, `L=2`, `d=4`.) So the **only** exponential in the problem is the
*measure of the exceptional boundary set*, and it is **exponentially small in N** (and in `j`). Off that
set — i.e. on the small-field domain where the tower actually lives (`g_j→0`) — the conditional measure is a
small perturbation of Haar, the LSI is `U′`-uniform, and `δ_W ≥ cN − O(1)`. **The conditioning does not cost
a factor in the constant; it confines a coexistence obstruction to a measure-`e^{−cNj}` set that is
negligible for large N.**

> **Conditioning status.** **[T]** Lemma B measure bound, `rate ∝ N`. **[NG]** all-`U′` uniform
> conditional LSI (boundary coexistence). **[C]** `U′`-uniform (★) on the small-field domain, since the
> large-field complement is `e^{−cNj}`-small. **No exponential in `δ_W`.**

---

## 6. Comparison to SZZ `K_S = N(1/2−24|β|)` — and the regime caveat

| | SZZ (`C1-SZZ-NSCALING`) | block-conditional (this note) |
|---|---|---|
| object | spectral gap / LSI of the **full, unconditional** Langevin | conditional coarse **effective Hessian** floor `δ_W` |
| coupling regime | **strong** `\|β\|<1/48` (`β_{SZZ}`) | **weak** tower `g_j→0`, `β_{SZZ,j}=11j/(12π²)≫1/48` |
| N-source of the `O(N)` | Ricci gain `N/2` − Hessian `8(d−1)N\|β\|` (both `O(N)`) | Haar curvature `W_geo ≥ N/12` (`O(N)`) − `O(1)` defect |
| N-scaling | **linear**, radius `1/48` N-free | **linear**, per-color `c_W ∈ [1/12, ≈0.13]` N-free |
| per-color constant | `1/2 − 24\|β\|` | `c_W = W_geo/N` (measured ≈0.125 at `N=2`); reservoir `E[Φ'']/N≈0.17` |
| `a`-scaling | fixed `a`, strong coupling | **`a⁰`**, `a`-uniform at fixed physical block |

**Answer to Q2.** The block decomposition's constant is `≥ c_W·N` with `c_W` **N-independent** (lower bound
`1/12` from `E[Φ'']`, measured `W_geo/N≈0.125`). It is **not** obtained *from* SZZ (SZZ is off-regime in the tower: witness part 1
shows `β_{SZZ,j} ≈ 0.093j > 1/48` at every shell), but **independently** from `W_geo`. **Conditioning on
`U′` does not cost a factor** (`1/N`, exponential, or boundary-volume) **in the constant**: the only
exponential is the `e^{−cNj}` measure of the large-field boundary set, which is *favorable* in N. The block
per-color coefficient is *smaller* than SZZ's (`≈0.13` vs `1/2`) only because `W_geo`'s proven bound `N/12`
and the measured `W_geo/N≈0.125` are weaker (Weyl/Vandermonde) estimates than SZZ's exact Ricci computation
— not because of any N-loss.

---

## 7. Status ledger and verdict

| item | statement | label |
|---|---|---|
| Hessian cost `O(N)` per mode, adjoint dim = multiplicity | SZZ Lemma 4.1 + trace structure | **[T]** |
| Haar curvature reservoir `E[Φ''] ≥ (N/12)‖v‖²` per mode, coupling-free | ch05 Lemma B | **[T]** |
| `W_geo = E[Φ'']−Var(Φ') ≥ 0`; `W_geo/N ≈ 0.125` measured | ch05 Thm / numerics | **[T]** (`≥0`) / **[R]** (`O(N)`) |
| projection `P`: `Lip C_blk = O(1)`, N- and a-independent | Lemma SC | **[T]** |
| boundary plaquettes `O(N)` per plaquette, support `O(a)` | `L3-band-decay` §1 | **[T]** (support) / **[O]** (rate) |
| `δ_W ≥ c·N`, N-free radius — *preserves C1's linear-N* | this note, conditional on (★) | **[C]** |
| `δ_W = Θ(a⁰)`, `a`-uniform — survives `a→0` at fixed block | this note, conditional on (★dec) | **[C]** |
| no `1/N`, no exponential in the constant | §3–§5 | **[T]** (mechanism) / **[C]** (synthesis) |
| all-`U′` uniform conditional LSI (no excision) | boundary coexistence | **[NG]** |
| bad-boundary set `μ ≤ e^{−cNj}`, shrinks in N | Lemma B | **[T]** |
| equipartition-defect coefficient (= the one-sided variance) | the open content | **[O]** |
| the block-conditional statement (★) itself | open, N/a-benign | **[O]** |

**Verdict.** The block-conditional version of (★) **does not lose N-uniformity and does not lose
a-uniformity.** N: `δ_W ≥ c_W·N` with `c_W` N-independent (the SZZ linear-N cancellation survives, because
the adjoint dimension `N²−1` is multiplicity and the bound is `O(N)` not `O(N²)`). a: `δ_W = Θ(a⁰)` at fixed
physical block (because `W_geo` is coupling-free and the bare-Hessian growth `~Nj` is cancelled by the field
amplitude `g_j`). Conditioning on `U′` adds **no** factor to the constant; the only exponential is the
`e^{−cNj}` measure of the large-field boundary set, which is negligible and shrinks with N. **The two genuine
obstructions are the `O(1)` equipartition-defect coefficient (the one-sidedness `Var≥0` subtracts) and the
`O(a)`-band clustering rate `(★dec)` — both analytic, both N/a-benign in their bookkeeping.** They are the
same `(★)∧(★dec)` as C2/G1-L3. The obstruction to L3 is one-sided hard analysis, **not** an N- or
a-combinatorial degradation. Clay 0%.

---

## 8. Verification ledger (this note)

| item | check | status |
|---|---|---|
| SZZ `K_S = N(1/2−24\|β\|)`, Lemma 4.1 `\|Hess\|≤8(d−1)N\|β\|`, radius `1/48` | read from `C1-SZZ-NSCALING.md` §1–§3 (itself verified against arXiv:2204.12737) | **[V]** (carried) |
| `W_geo≥0`, `E[Φ'']≥(N/12)‖v‖²`, Lemma SC `C_blk=O(1)` | `ch05-theorems.tex` §§the W_geo theorem, Lemma SC | **[V]** (in-repo) |
| tower flow `g_j²=1/(2β₀j)`, `β_j=4Nβ₀j`, `β₀=11N/(48π²)`; Lemma B `e^{−c/g_j²}` | `C2-BALABAN-UNIFORM.md` §2.1–§2.2 | **[V]** (in-repo) |
| `β_{SZZ,j}=11j/(12π²)>1/48 ∀j≥1`; reservoir`/N=1/12`; `E[H]·g_j²` j-independent; `μ_bad` vs N | witness `research/v3/w7/s5_n_uniformity.py` (`py -I`), parts 0–4 | **[V]** numerically |
| conditioning: averaged vs uniform-in-`U′` conditional LSI gap; boundary-coexistence `[NG]` | standard disintegration argument (§5.2); no external source downloaded | **[C]** analysis |

*Refs.* `C2-BALABAN-UNIFORM.md` §2.4, §2.1–2.2; `C1-SZZ-NSCALING.md` §1–§4; `docs/paper-v2/ch03-existence.tex`
§L3, §uv-face, eq. `dim`; `ch05-theorems.tex` (W_geo theorem, Lemma SC, Feshbach lemma); `G1-L3.md` §2;
`G1-L3-UV-FORCE-BOUND.md` §1–§6; `G1-L3-TOWER-POINCARE.md` §2–§5; `G1-UNIFORM-ZA.md` §1–§3;
`PROOF-COARSE-POINCARE.md` (a)–(c); `L3-band-decay.md` §1–§2. Shen–Zhu–Zhu CMP 400 (2023) 805,
arXiv:2204.12737. No third-party sources downloaded; the witness is self-contained (numpy only).
