# E2-FLIP-MINORANT — attack the E₂ flip as a *minorant of `m₁`* head-on: construct a computable uniformly-positive lower bound `Φ`, or prove none exists in the lower-bound class

> 2026-11-24 · Direct attack on the **minorant structure** of `E2-FLIP-NOGO §0` / Theorem 0: a flip is *exactly* an explicit functional `Φ` with `0 < c ≤ Φ ≤ m₁`, computable from accessible data, uniform in `(L,a)`, acid-clean, non-circular. This note does **not** re-run the four-horn data-mode dichotomy (`E2-FLIP-NOGO`); it takes the **five concrete LOWER-bound methods** named in the task — (i) spectral/eigenvalue enclosure, (ii) resolvent/Fredholm determinant, (iii) cluster/KP, (iv) comparison/domination, (v) transfer-matrix/positivity — and asks of each, rigorously: *is the output a genuine `c ≤ m₁`, computable, uniform, acid-clean, and free of `E₂`/gap input?*
> **Labels (charter §2):** **[T]** 真定理 · **[C]** 条件定理 · **[NG]** no-go for a *named class* · **[O]** 精确缺失步 · **[J]** 判断. External facts carried from the repo **[R]**; the fresh checks **[V-here]** are `research/unconditional/w7/e2_minorant_witness.py`. **No QED on anything unproven; this does NOT prove Clay and does NOT prove no proof can exist** (`∀class ≠ ∀possible`, `KEY-EXISTS §5`).

---

## Verdict up front

> ### **SPLIT VERDICT — one genus of minorant is genuinely CONSTRUCTED (on the strong-coupling family); the full `(L,a)`-uniform minorant across the crossover is IRREDUCIBLE [NG-for-class + [O]].**
>
> **(A) CONSTRUCTED [T] — the only genuine minorant:**
> On the strong-coupling family `F_sc = {(β,a): β<β₁(N)=O(1/N)}` there is an **explicit constant** `m₀(β)>0` — computed from the graded KP/character expansion, with **leading value** `ln(1/u(β))/a`, `u(β)=c_f(β)` the positive fundamental character coefficient — such that the minorant
> `Φ_sc(β,a) := m₀(β)/a`  satisfies  **`0 < Φ_sc(β,a) ≤ m₁(β,a)`**, **uniform in the spatial volume `L`** (`PROOF-P3-UNCONDITIONAL §5` residue `f₁²≥r₀`; the gap `m₁>0` is Osterwalder–Seiler '78 / SZZ '23). `u∈(0,1)` for `β<β₁`, so the rate is positive and explicit. This is a real, computable, `L`-uniform minorant — **the strongest `Φ` that exists**. It is `[T]` on `F_sc`.
>
> **(B) IRREDUCIBLE [NG] for the class `C_lb` of certified lower-bound methods, uniform across the crossover / into `a→0`:** every one of the five methods, extended off `F_sc`, satisfies exactly one of: **(α) it is an UPPER bound, not a lower bound** (Rayleigh/min–max, moment ratio, finite-matrix Gershgorin on a truncation); **(β) it needs a spectral-separation oracle `ρ ≤ E₂`** — and by Single-Scale Covariance `E₂=ρ_r·E₁`, so `ρ` is a fixed multiple of the very `m₁` being bounded (Theorem **ENCL**; covers Gershgorin/Brauer/Kato–Temple/Lehmann–Goerisch); **(γ) it is vacuous or non-computable** (trace/partition minorant is `≤0` until `t≳1/E₂`, and needs the full `Tr e^{-tH}`; Fredholm-determinant coefficients *are* the moments — Theorem **DET**); **(δ) it reduces to proving the Poincaré inequality = `H1`** (transfer-matrix/positivity/Cheeger/LSI — Theorem **POS**); **(ε) it needs an order-theoretic domination whose only universal member is the free field, gap `0`** (comparison/domination; the "stable quadratic minorant" engine of the cluster expansion fails at Savvidy — Theorem **CMP**).
>
> **Exact stuck step.** A uniform minorant needs a *uniform positive spectral-separation at `k→0`* on `(vac)^⊥` — i.e. either a certified `E₂` (which SSC makes `=ρ_r·m₁`), or a stable (sign-definite) quadratic minorant of the Savvidy Hessian (which is indefinite), or a Poincaré inequality `Var≤(1/c)ℰ` (which is `H1`). **One of these three — all the same open object** — is the gap.

---

## §0 The minorant object, and the direction bar (why "lower" is the whole content)

**Setup** (`E2-FLIP-NOGO §0`, `BRIEF §1`). `H ≥ 0` self-adjoint transfer Hamiltonian, unique vacuum `Ω` (`PROOF-GAUSS-INVARIANT-PF`), `m₁ = E₁ − E₀` with `E₀=0` after shift. (On `(vac)^⊥`, `H` has bottom `E₁=m₁`, second level `E₂`.) `ν` = KL spectral measure, `ν≥0`, `supp ν⊂[m₁²,∞)`, atom `f₁²=ν({m₁²})`.

**Definition (minorant / flip functional).** `Φ` is a *minorant on a family `F`* iff **(LB)** `Φ(μ) ≤ m₁(μ)` for all `μ∈F`; it is a *flip functional* iff in addition `inf_F Φ ≥ Δ>0`. Theorem 0 (`E2-FLIP-NOGO`) then forces `Φ(μ_free)≤0` at the free-massless endpoint — the acid test is *implied* by (LB), not an extra filter. So:

> **"Construct a flip functional" ≡ "produce a computable uniformly-positive lower bound `0<c ≤ Φ ≤ m₁`."** The task is exactly a *lower-bound* problem, and the only question is the *method*.

**The direction bar (rigorous; the trap the task flags).** For self-adjoint `H`, *every* variational instrument that minimizes/maximizes a Rayleigh quotient or a moment ratio bounds eigenvalues **from above**:

- **Courant–Fischer / Rayleigh–Ritz [T].** `E_k = min_{dim V=k} max_{0≠ψ∈V} R(ψ)` (and `= max_{codim V=k-1} min_{0≠ψ∈V}R(ψ)`), and `R(ψ)≥E₁` for all `ψ∈(vac)^⊥`. Hence for a trial subspace `V` (or a single trial) the computed Ritz value ≥ the true eigenvalue, and `E₁ ≤ inf_ψ R(ψ)` with equality only at the exact eigenvector. **Every single RQ value is `≥ m₁`; the inf is `m₁` but is not certifiable without the eigenvector** (`[V-here]` A: 200k random trials, `min R = 1.000000 = E₁`, all `R≥E₁`). A "minorant built from a Rayleigh quotient" would need `Φ=R≤m₁`, i.e. `R≤E₁`, which fails identically. **RQ is an UPPER bound.**
- **Moment ratio `M_{k+1}/M_k ≥ m₁²` [T]** (`GAP-G2-edge`): positivity ⇒ `m₁² ≤ (M_a/M_b)^{1/(a−b)}` — always upper.

So the minorant must come from a genuine *lower*-bound method. There are exactly five in the corpus's reach; §1 does each. The `[V-here]` direction table:

| method | output direction | source of the sign |
|---|---|---|
| Rayleigh / min–max / Ritz | **UPPER** on `E_k` | `E_k ≤ Ritz_k` (min–max realizes the eigenvalue from above) |
| moment ratio `M_{k+1}/M_k` | **UPPER** on `m₁²` | `ν≥0` |
| Gershgorin/Brauer on a *truncation* | UPPER on the operator | compression raises eigenvalues; tail unbounded |
| Temple/Weinstein/Kato | **LOWER** *given* `ρ≤E₂` | `Var≥(ρ−R)(R−E₀)` under `E₀≤ρ≤E₁` |
| Fredholm-det coefficients | moment data ⇒ UPPER/edge-blind | elementary symmetrics = moments |
| trace/partition `−(1/t)ln Z_⊥` | LOWER but `≤0` until `t≳1/E₂` | Tauberian (`Z_⊥≥e^{−tE₁}`) |
| cluster/KP | **LOWER**, `L`-uniform, only `β<β₁` | stable quadratic minorant `H″≥ε(−∂²)−C` |
| comparison/domination | LOWER iff model gap `≤m₁` **and** order holds | free model has gap `0` |
| transfer-positivity / Cheeger / LSI | the gap itself (`H1`) | gap `=` Poincaré constant |

---

## §1 The five lower-bound methods, one theorem each

### (i) Spectral / eigenvalue ENCLOSURE — needs a separation oracle `ρ≤E₂`  **[NG for the class]**

This is the class of *a posteriori* eigenvalue bounds: Gershgorin / Brauer / Weyl discs, Kato–Temple, Weinstein, Lehmann–Goerisch, and the modern "certified approximate eigenpair + resolvent bound" enclosures. All are functions of *(an approximate eigenpair `(ψ,R)`, its residual `Var=‖(H−R)ψ‖²`, and a lower bound `ρ` on the **next** level)*.

> **Theorem ENCL (minorant ⇒ separation ⇒ `E₂`). [T, + SSC-conditional]**
> Let `H` be self-adjoint on `(vac)^⊥`, bottom `E₁=m₁`, second `E₂`. For any normalized `ψ` with `R=⟨ψ,Hψ⟩`, `Var=‖(H−R)ψ‖²`, the **Temple bound**
> `m₁ ≥ R − Var/(ρ − R)`  is sound (`≤ m₁` for every `H` with that `E₂`) **iff `ρ ≤ E₂`**; and it is strictly increasing in `ρ`. Hence: (1) the *sharpest sound* Temple minorant is `R − Var/(E₂ − R)`, which consumes `E₂`; (2) any certified `ρ` satisfies `ρ ≤ E₂`, so a sound minorant is a function of `(R,Var,ρ≤E₂)`. **No function of `(R,Var)` alone lower-bounds `m₁`** (impossibility lemma, `TEMPLE-WEINSTEIN §3`: fix `R,Var`, move an atom to `−N`).
> *Proof.* Temple's inequality requires `ρ≤E₂` (proof: `Var≥(ρ−R)(R−E₀)` under `E₀≤ρ≤E₁`; here bottom `E₀→m₁`, second `E₁→E₂`). `∂_ρ[R−Var/(ρ−R)]=Var/(ρ−R)²>0`. The impossibility lemma is `LEVER §0B`. ∎ (`[V-here]` B: with `{0,1,1.5}`, `ρ=E₂` gives `b=1.0000≤E₁` sound; `ρ=2E₁` or `3.0` gives `b=1.0006…1.10>E₁` — an **overclaim**, a false *too-high* "lower" bound, growing with admixture.)
> **Circularity (SSC).** By Single-Scale Covariance (`E2-REGRESS §0.4`, `G2D-ISOLATION-I2 §2`), `E₂ = ρ_r·E₁` with `ρ_r=m_{0⁺⁺*}/m₁` a fixed `β,L,a`-independent number. Thus `ρ ≤ E₂ = ρ_r·m₁`: the separation oracle is a **fixed multiple of `m₁`**. A non-vacuous threshold (`Φ≥c`) pins `ρ`, `R ≈ m₁` and hence `m₁` itself to within the gap scale — the minorant *is* the gap, one rung up. **`[NG]` for the enclosure class; the exact input `ρ≤E₂` is `I2=H*`.**

*Gershgorin/Brauer specifically:* they bound eigenvalues of a **known finite matrix** from its entries. The YM transfer matrix is infinite-dimensional (continuous link variables); any truncation `P_N H P_N` has eigenvalues `μ_k ≥ E_k` (**upper**, Rayleigh–Ritz), so Gershgorin on the truncation bounds the *truncation* above the operator; converting to a lower bound on `H` needs the compression tail — the same `E₂`/UV oracle. There is no finite-matrix Gershgorin lower bound on the infinite operator.

### (ii) Resolvent / Fredholm DETERMINANT / trace — moments or Tauberian  **[NG]**

**DET (the determinant's analytic data = moments). [T via Theorem α]** With `T=e^{−aH}`, eigenvalues `λ_n=e^{−aE_n}∈(0,1]`, the Fredholm determinant `det(1−zT)=Σ_{n}(−1)^n z^n t_n`, `t_n=Tr(∧^nT)` = the **elementary symmetric functions** `= e_n(λ_1,λ_2,…)`. By Newton's identities these are fixed polynomial combinations of the power sums `p_k=Σλ_i^k=Tr(T^k)=Σ_n e^{−akE_n}` — **exactly the moment data** of `ν`. By the linear-data lemma (`E2-FLIP-NOGO Theorem α`: `δ₀+δ₂` vs `2δ₁`; `[V-here]` order-3 extension), finitely many `t_n` cannot pin the edge zero `z₁=1/λ₁=e^{aE₁}`; Maclaurin's inequalities (elementary symmetric means are log-concave) give **upper** bounds on the top eigenvalues. A *lower* bound on the first determinant zero in modulus needs the **full spectrum** (the edge is not a finite-data functional). So the determinant route `⊂` Horn α — **edge-blind**.
**Trace / partition minorant.** The one genuine *lower*-bound inequality in this family is `Z_⊥(t):=Tr_{⊥}e^{−tH}=Σ_{n≥1}e^{−tE_n} ≥ e^{−tE₁}`, giving `m₁ ≥ −(1/t)·ln Z_⊥(t)`. It is a real lower bound (`[V-here]` C). But: (a) `Z_⊥(t) → e^{−tE₁}` only for `t ≳ 1/(E₂−E₁)` — the bound is positive *only* once `t` **resolves `E₂`**; at every smaller, accessible `t` it is `≤0` (`[V-here]` C: edge `1`+continuum gives `−9.01,−3.52,−1.01` at `t=0.5,1,2`, turning `+0.32,+0.71` only at `t=5,10`; massless toy `≤0` throughout). (b) `Tr e^{−tH}` for the YM transfer operator is UV-divergent / not numerically evaluable in a way that certifies the excited trace. This is the trace-space twin of Horn γ: *finite-`t` = vacuous, `t→∞` = the answer, resolution needs `E₂`.* **`[NG]`.**

### (iii) Constructive CLUSTER / KP — a real minorant, only `β<β₁`; transport = the wall  **[T on F_sc; NG off it]**

This is the **one place a genuine, computable, `L`-uniform lower bound exists**. The high-temperature (character) expansion `e^{(β/N)Re tr U_p}=Σ_r c_r(β)χ_r(U_p)` has `c_r(β)≥0`; with `u(β):=c_f(β)∈(0,1)`, the graded Kotecký–Preiss criterion `u(β)K<1` (`K≤CN`) holds for `β<β₁=min(1/48,c/N)`, and the connected `trF²` correlator is dominated by a single flux tube `a₀u^{t/a}=a₀e^{−m₁t}`, `m₁=ln(1/u)/a`. Hence (OS '78 / SZZ '23 for the gap; `PROOF-P3-UNCONDITIONAL §5` for the residue)

> **`0 < Φ_sc(β,a) := m₀(β)/a ≤ m₁(β,a)`, uniform in `L`, for `β<β₁(N)`** (leading rate `ln(1/u(β))/a`).  `[T]`

`[V-here]` D computes `u(β)=c_f(β)` for SU(2) on `β∈[0.005,0.2]`: `u∈(0,1)` (e.g. `u=1.0×10^{−2}` at `β=0.02`), so `ln(1/u)∈[2.30,5.99]>0`. **This is the strongest `Φ` that exists.** It is a genuine, explicit, computable, `L`-uniform lower bound.
**Why it does not close the task.** `F_sc` excludes the crossover `β≈β_c≈2.47` and the continuum `a→0`; the physical statement needs `β→∞` along the AF curve. **Extending `Φ_sc` = the "transport across the Savvidy crossover," which is `transport≡confinement`** (`BRIEF §4`, `SMALLCIRCLE-ADIA`, `CHATGPT-METHOD-RGENTRY`). Two impossibilities, both `[T]`-level:
- *Analyticity radius.* The cluster expansion converges only on the strong-coupling disk `|β|<β₁`; the Savvidy crossover lies outside it, and no zero-free region across it is known (`MS-BLIND §1.3`).
- *Engine failure (the sharp form).* The cluster expansion's engine is a **stable quadratic minorant** `H″(A) ≥ ε(−∂²) − C` with `ε>0` uniform (`MS-BLIND §1.5`). At Savvidy the transverse Hessian `−D²+2adF_B` is **indefinite** (`λ_{0,−}=k_∥²−gB`, bottom `−gB<0`, extensive negative spectrum `≈(gB)²` per volume), so `ε>0` is impossible and **no sign-definite polymer weight exists** — the expansion cannot even begin at the crossover (`D-BBS-LOWRANK-KATO`: the low-rank+Kato transfer empties the UV cone but is soft at `k=0`, `PᵀHP|_{\rm coarse}≈−0.3`). **`[NG]` for extending the minorant; the exact stuck step is the stable quadratic minorant at `k→0`.**

### (iv) COMPARISON / domination by a reference model — order direction + acid test ⇒ `c=0` or import the mass  **[NG]**

A comparison minorant needs a reference model `M` with (a) a **computable** gap `gap(M) ≤ m₁` and (b) an **order relation** making `M` provably "floppier" than YM: `μ_M ⪯ μ_YM` in a sense that implies `DirichletForm_M ≤ DirichletForm_YM` on matched functions (then `Poincaré(M) ≤ Poincaré(YM) = m₁`, so `Φ=gap(M)`).

> **Theorem CMP. [T for the order direction; NG for the class]**
> (1) *Direction.* Every off-the-shelf domination inequality (Brascamp–Lieb, Prékopa, stochastic/convex order, log-Sobolev comparison) propagates a *lower* Poincaré constant from the *more log-concave* measure to the less. To get `gap(M)≤m₁(YM)` one needs `M` provably *at least as log-concave* as `μ_YM` **and** one needs to *know* that order — which requires the `Γ₂`/curvature of `μ_YM`, i.e. the Savvidy Hessian. That Hessian is **indefinite** (`Hess −D²+2adF_B`; `CLAUDE.md §5`, `W4` 7.1% negative after gauge-quotient), so YM is **not** log-concave and no positive-curvature reference survives with `c>0`.
> (2) *Acid test.* The **only** universally valid domination of any OS measure is by the free field, whose gap is `0`. By Theorem 0, `Φ(μ_free)≤0`; so a comparison minorant that is true for *all* admissible `μ` evaluates `≤0` — **`c=0`, vacuous**. To escape, the reference must import YM structure (the gluon condensate `⟨¼trF²⟩`), and the needed inequality is exactly **(SGL)** `m⁴ ≥ κ⟨¼trF²⟩, κ>0` (`MS-BLIND §4`) — which is **the gap** (it makes the condensate a mass floor). Nothing known implies (SGL); positivity gives only the reverse sign (`m² ≤` moment ratio). **`[NG]`; the exact object is (SGL).**

This is the same object as the Savvidy "stable quadratic minorant": a positive-quadratic *lower* model `A≤H` with `spec(A)∩(0,c)=∅` would dominate `H` and give the gap; the natural `A` (the fluctuation Hessian) is indefinite, and adding a positive mass term to fix it **is** the gluon condensate.

### (v) Transfer-matrix POSITIVITY — the gap = Poincaré constant = `H1`  **[O]**

The transfer matrix is positive and RP; the gap of `H` **is** the Poincaré constant of the measure `μ`, `m₁ = λ₁ = inf_{f} ℰ(f)/Var_μ(f)`, `ℰ(f)=⟨f,Hf⟩=‖∇f‖²_{L²(μ)}` (standard Dirichlet-form identity). Lower-bounding it by positivity means proving the inequality `Var_μ(f) ≤ (1/c)ℰ(f)` **for all `f`** — i.e. a Poincaré/coarse-mode anti-concentration inequality.

> **Theorem POS. [T reduction] The transfer-positivity minorant *is* `H1`.** Every positivity tool (Cheeger `gap≥h²/2`, LSI `gap≥2α_{LSI}`, Doeblin/coupling contraction, Brascamp–Lieb when applicable) lower-bounds `m₁` precisely by proving a Poincaré/anti-concentration estimate `Var(f)≤(1/c)ℰ(f)`. That estimate — uniform over `f`, at `k→0`, gauge-invariant — **is `H1`** (`PROOF-COARSE-POINCARÉ`; residual-dictionary `D.1`). Moreover, RP/positivity **alone** is scale-covariant (`GAP-G2-edge`: RP theories with arbitrarily small gaps exist) and forces only `m=0` atoms (`GOLDS-POS-MASS`); it cannot force `m₁>0` without the anti-concentration content. So the positivity route neither adds a 6th species nor shortcuts `H1` — **it is `H1`.** `[O]` (the open target), not a solution.

---

## §2 The constructed minorant, stated fully (`Φ_sc`)

> **Theorem (strong-coupling minorant). [T]** Fix `a>0`, `N`. There is `β₁=min(1/48,c/N)>0` and an explicit constant `m₀(β)>0` (equal to the KP flux-tube rate `ln(1/u)/a` up to the `1−O(u)` decoration factor, `PROOF-P3-UNCONDITIONAL §3–§4`) such that for all `β<β₁` and every spatial volume `L` (and strip length `t`),
> ```
>  0 <  Φ_sc(β,a)  :=  m₀(β)/a   ≤   m₁(β,a;L) ,
> ```
> with `u(β)=c_f(β)` the fundamental character coefficient, `u∈(0,1)` and `u≍β`. `Φ_sc` is explicit, computable, and `L`-uniform.
> *Hypotheses / imports.* Osterwalder–Seiler connected-composite + Kotecký–Preiss expansion (standard); the SU(N) multiplicity bound `C_g≤CN` (proven, numeric factor only shifts `β₁`); `PROOF-P3-UNCONDITIONAL §2–5`. The **gap** `m₁>0` is the classical strong-coupling gap (OS '78, SZZ '23); `PROOF-P3` independently supplies the residue `f₁²≥r₀>0`. (The bare `ln(1/u)/a` is the leading value; the dressed mass differs by the `O(u)` self-energy, which the KP radius controls, so the rigorous bound is the explicit `m₀(β)` with that rate — the witness D checks `u∈(0,1)`.)
> *Acid status.* `F_sc` does **not** contain the free-massless endpoint (that is `β→∞`). So the acid test is *not applicable on `F_sc`* — and by Theorem 0 it cannot be: any extension of `Φ_sc` to a family containing the free endpoint must fall to `≤0` there. `Φ_sc` is a genuine `[T]` minorant **on its restricted family**, not a flip functional on the full family.
> *What it is NOT.* It is not uniform in `β` across the crossover, does not reach `a→0`, and is not a Clay gap proof.

So the honest answer to "construct such a `Φ`" is: **`Φ_sc` exists and is as strong as the proven strong-coupling region allows; it cannot be made the full-family minorant by any of (i)–(v) without the gap.**

---

## §3 The no-go for `C_lb` (assembled) + minorant rigidity

**Definition.** `C_lb` := minorants produced by a *certified lower-bound method*: `C_lb = (i) ∪ (ii) ∪ (iii) ∪ (iv) ∪ (v)` — enclosure (needs `(R,Var,ρ)`, `ρ≤E₂`), determinant/trace (moment or Tauberian data), cluster/KP (needs a stable quadratic minorant), comparison (needs an order + reference gap), positivity (Poincaré/`H1`).

> **Theorem (five-method dichotomy). [NG for `C_lb`; uniform statement]** For every `Φ∈C_lb` with (LB) and uniform-positivity on the full `(L,a)` family `F`:
> `Φ` is **(α)** an upper bound / direction-wrong, or **(β)** circular (needs `ρ≤E₂=ρ_r m₁` — Theorem ENCL — or the free-field `c=0` — Theorem CMP — or the full spectrum / `t≳1/E₂` — Theorem DET), or **(γ)** reduces to `H1` (Theorem POS). There is **no fourth**: a uniform, non-circular, acid-clean, computable positive minorant.
> *Assembly.* (i)=β, (ii)=β/α (Theorem DET), (iii)=the restricted `[T]` (no escape off `F_sc`), (iv)=β/α (Theorem CMP), (v)=γ (Theorem POS). Each uses only the *method's own* structure; no cross-talk assumed. ∎

> **Theorem 0′ (minorant rigidity — no "cheap" minorant). [T]** Any `Φ` with `0<c≤Φ≤m₁` uniformly on `F` satisfies `sup_F Φ/m₁ ≤ 1` and, by the scale dichotomy (`E2-FLIP-NOGO §3`), a *dimensionful* uniform `Φ` is forced onto `Λ_QCD` with `Φ=c'·Λ_QCD ≤ m₁=c·Λ_QCD`, i.e. `c'≤c=m₁/Λ_QCD` — **the minorant's coupling-dependence is the gap's.** There is no functional that stays uniformly positive and *strictly below* `m₁` by a coupling-independent factor: "find a minorant" has no content beyond "find a lower bound," which is the gap. (This sharpens Theorem 0: the difficulty is not *weakness* of the minorant — weakness is allowed — but *existence*.)

---

## §4 Fresh witness `[V-here]` (`research/unconditional/w7/e2_minorant_witness.py`)

Run `py -3.12 e2_minorant_witness.py`. It re-verifies the load-bearing **directions** (tiny toy spectra; no heavy lattice job):
- **(A)** Rayleigh quotient: 200k random trials on `(vac)^⊥`, `min R = 1.000000 = E₁`, **all `R≥E₁`** ⇒ RQ is an **upper** bound; a "minorant from RQ" violates `Φ≤m₁`.
- **(B)** Temple: `ρ=E₂` gives `b=1.0000≤E₁` (sound); `ρ=2E₁`/`3.0` give `b=1.0006…1.10>E₁` (**overclaims** = false too-high "lower" bounds), growing with `φ₂` admixture. Only `ρ≤E₂` is sound.
- **(C)** trace/partition minorant: with edge `1`+continuum, `−(1/t)ln Z_⊥ = −9.01,−3.52,−1.01` at `t=0.5,1,2` (**vacuous ≤0**), becoming `+0.32,+0.71` only at `t=5,10` (= resolving `E₂`); massless toy `≤0` throughout. Confirms Tauberian: finite-`t` vacuous, `t→∞` = the answer.
- **(D)** strong-coupling: SU(2) `u(β)=c_f(β)∈(0,1)` on `β∈[0.005,0.2]` (`u≈1.0×10^{−2}` at `β=0.02`), `ln(1/u)∈[2.30,5.99]>0` ⇒ rate `>0` and `Φ_sc>0` — the genuine `[T]` minorant on `F_sc`.

These are checks of already-banked lemmas (min–max, Temple validity condition, Tauberian direction, character-expansion positivity), not new numerics; the substantive partials (`Temple 0/3474`, strong-coupling gap OS/SZZ, Savvidy `PᵀHP` min `−0.3`) are inherited `[R-internal]`. No fabricated bound.

---

## §5 Honest boundary — what this does and does NOT establish; the flip condition

- **Does NOT prove Clay; does NOT prove no proof exists.** `[NG]` is for the named method-class `C_lb` and the uniform statement (`∀class ≠ ∀possible`, `KEY-EXISTS §5`, `CONVERGENCE-redteam §d`). A not-yet-invented 6th-species lower-bound method could live outside `C_lb`; by Theorem POS/(L-2) of `E2-FLIP-NOGO` it would *be* `H1` or a genuinely non-positivity, non-circular, sign-controlling inequality — the charter §1 target.
- **Does NOT undersell `Φ_sc`.** On `F_sc` the minorant is explicit, computable, `L`-uniform, and `[T]`. This is real and it is the strongest `Φ` that exists; it simply cannot cross the Savvidy crossover by (i)–(v).
- **The single exact stuck step**, restated: a uniform positive spectral separation on `(vac)^⊥` — equivalently a certified `E₂` (SSC: `=ρ_r m₁`), or a stable (sign-definite) quadratic minorant of the indefinite Savvidy Hessian at `k→0`, or the Poincaré inequality `Var≤(1/c)ℰ` (`H1`). **One object, three names** (`GAP-G2-edge`, `SUPPORT-NOT-ATOM`: support/isolation/`H1` are one wall).
- **Flip condition (unchanged, sharpened for `C_lb`).** This verdict flips iff someone exhibits a *lower-bound method* `Φ∉C_lb` producing `0<c≤Φ≤m₁`, uniform in `(L,a)`, acid-clean, without an `E₂`-type separation oracle, without a free-field-collapsing comparison, and not equal to `H1`. By §1 that object must be a **stable, sign-definite, gauge-invariant quadratic minorant of the YM Dirichlet form at `k→0`** — the "6th species." None is known; finding it is the missing machine.

---

## §6 References

**Internal (absolute paths; the banked lemmas this note assembles):**
- `D:\Dev\math\CLAUDE.md` (§1 build-don't-stop; §2 labels; §5 wall = Savvidy concavity / key = `G″(0)` / `H1`); `D:\Dev\math\docs\BRIEF-FOR-MODEL.md` (§2 five conditions; §4 transport≡confinement; §5 subgoals).
- `E2-FLIP-NOGO.md` (Theorem 0, the four-horn `C0` no-go, scale dichotomy §3, loopholes L-1/L-2/L-3); `E2-REGRESS-CIRCULARITY.md` (§0.2 L1; §0.4 SSC `E₂=ρE₁`; §0.6 Temple two horns); `KEY-EXISTS-OR-IMPOSSIBLE.md` (§5 `∀class≠∀possible`, §7 calibration).
- `TEMPLE-WEINSTEIN.md` (§1 Temple/Weinstein, §3 impossibility lemma — Theorem ENCL); `TEMPLE-RIGOROUS-H.md` (§1 gap-agnostic `H`, §2 kinematic-edge killed); `PROOF-EDGE-E2.md` (Temple needs `E₂`; every bottom functional blind to `λ₁−λ₂`); `BS-VARIATIONAL-LOWER.md` (Temple needs `ρ≤E₂`; kinematic `ρ=4m_g²` fails); `LEVER-nonconvex-edge.md` (§0A free-field adversary, §0B linear-data lemma, §1–4 direction/hypothesis/insensitivity).
- `GAP-G2-edge.md` (moment ratio upper; residue–level-gap lemma); `SUPPORT-NOT-ATOM.md` (support vs atom); `H1-IS-IT-THE-GAP.md` (squeeze `f₁²/χ≤m₁²≤M₀/χ`, `G″=1/χ`); `PROOF-COARSE-POINCARÉ.md` (H1 = coarse Poincaré = `m_eff²`, Theorem POS); `GOLDS-POS-MASS.md` (positivity pins only `m=0`).
- `PROOF-P3-UNCONDITIONAL.md` (§2–5 KP criterion, `m₁=ln(1/u)/a`, `f₁²≥r₀` — the constructed `Φ_sc`); `MS-BLIND.md` (§1.5 stable quadratic minorant fails at Savvidy — Theorem CMP/§1(iii); §4 (SGL)); `D-BBS-LOWRANK-KATO.md` (low-rank+Kato, soft at `k=0`, `PᵀHP≈−0.3`); `G2D-ISOLATION-I2.md` (§2 SSC, §3 vNW pointwise).
- `INEQ-T2-variational-residue.md` (residue dual is an infimum = upper); `INEQ-R4-optimal-transport.md`, `INEQ-R5-donsker-varadhan.md` (the atom is usc; sup-duals degenerate to 0); `W3-adjacent-4D.md` (scale-free ⇒ `Λ_QCD`).
- `D:\Dev\math\research\unconditional\w7\e2_minorant_witness.py` (`[V-here]` A–D; run `py -3.12 e2_minorant_witness.py`).

**External — [standard]/[R] (textbook; not re-verified this session):** Temple *Proc. London Math. Soc.* (1928); Weinstein *Phys. Rev.* **40** (1932) 797; Kato *J. Phys. Soc. Japan* **4** (1949) 334; Kato *Perturbation Theory for Linear Operators* (1966); Lehmann–Goerisch (1970s); Courant–Fischer; Reed–Simon IV Thm X.12, XIII.5; Weyl (1912); Gershgorin (1931), Brauer (1946); Kotecký–Preiss *CMP* **103** (1986) 19; Osterwalder–Seiler *Ann. Phys.* **110** (1978) 440; Seiler LNP **159** (1982); Shen–Zhu–Zhu *CMP* **400** (2023) 805, arXiv:2204.12737; Savvidy (1977); Brascamp–Lieb *Adv. Math.* **22** (1976) 366; Gross log-Sobolev (1975); Bakry–Émery (1985); Bernstein/Widder (completely monotone). No fabricated bound; no QED.

---

*Status: **SPLIT** — **(A) `[T]` MINORANT CONSTRUCTED on `F_sc`** (`Φ_sc=m₀(β)/a` with leading rate `ln(1/u(β))/a`, explicit, computable, `L`-uniform, `β<β₁`, from OS'78/SZZ'23 + PROOF-P3); **(B) `[NG]` for the certified lower-bound class `C_lb` (enclosure / determinant–trace / cluster-KP / comparison / positivity) uniform across the crossover, assembled method-by-method** (Theorem ENCL: minorant⇒`ρ≤E₂`⇒SSC gives the gap; Theorem DET: determinant coefficients=moments, trace minorant Tauberian-vacuous; §1(iii): cluster minorant real only `β<β₁`, stable-quadratic-minorant engine fails at Savvidy; Theorem CMP: comparison ⇒ free-field `c=0` or (SGL); Theorem POS: positivity ⇔ Poincaré ⇔ `H1`). **Exact stuck step:** a uniform positive spectral separation on `(vac)^⊥` = certified `E₂` = stable quadratic minorant at `k→0` = Poincaré / `H1` — one object. NOT a proof of Clay; NOT a global impossibility theorem. Clay stays open.*
