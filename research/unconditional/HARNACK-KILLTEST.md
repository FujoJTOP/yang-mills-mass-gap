# HARNACK-KILLTEST — is a non-circular, L-uniform vacuum Harnack provable at β=2.5?

> 2026-11-24 · SU(2) Wilson lattice gauge theory, 4D. Target: the **local vacuum Harnack bound** on the Perron–Frobenius ground vector `h₀` that the Wasserstein/Ollivier route (`WASSERSTEIN-OLLIVIER.md`) flagged as its one surviving gap — the `h₀`-circularity of the Doob kernel `P = T(h₀·)/λ₀`.
> Companion witness: `research/unconditional/w7/harnack_killtest.py` (`py -3.12`, single job). All numbers below are `[V-here]`.
> Labels (charter §2): **[T]** true theorem · **[C]** conditional · **[O]** exact missing step · **[NG]** named-class no-go · **[V-here]** computed this round · **[R]** external, cited. **No QED.**

---

## 0. VERDICT — one line

**The absolute local Harnack bound on `h₀` IS provable non-circularly and IS `L`-uniform at β=2.5 — `[T]`, with the explicit constant `C = (8β_s+2β_t)\,sin(ψ/2) ≤ 10β ≈ 25` — so the Wasserstein doc's `§4` "circular / needs H1 to prove" diagnosis is *wrong on the mechanism*. But this is NOT a genuine opening, for a different and decisive reason: the bound is (i) **non-small** — it is `Θ(β)`, and the only way to turn it into a contraction (Dobrushin row sum `6(1−e^{−C})<1`) needs `C<0.182`, i.e. `β≲1/55`, which is *exactly the strong-coupling cell `β<1/48` already owned by P3*, while at β=2.5 the row sum is `≈6 ≫ 1` (no contraction); and (ii) **vacuous** — it is finite for the gapless *compact* free field (compact U(1) Coulomb phase) and, for a free mode, the Harnack constant is *linear in the gap* (`C_H = ½M²·Δ`), so it is **anti-diagnostic** for a gap lower bound. An *absolute* bound `|Δ\,log h₀|≤C` does not determine `h₀`, so the Doob kernel stays `h₀`-dependent. **The circularity was a red herring; the wall is the same local→global contraction = H1.** `[T]/[V-here]/[O]`. Clay(G2) 0%.

**一句话**：非循环、`L`-一致的 `h₀`-Harnack 在 β=2.5 **是可证的**（`C≤10β`，证明只用 `Th₀=λ₀h₀` 与 `h₀>0`，不预设隙）——所以 Wasserstein 册 §4 说它"循环"是**误诊**。但它**不是开口**：常数是 `Θ(β)`（β=2.5 行和≈6，无收缩）、且**对自由场空洞**（紧致无质量场也满足；自由模的 Harnack 常数 ∝ 隙，方向相反）。真正的墙仍是局部→全局收缩 = H1。

---

## 1. Setup — what `h₀` is, and what "local Harnack" means

Let `Λ = Λ_s × Λ_t` be a finite 4D lattice, spatial volume `V = L³` sites, spacing `a`. Let `E_s` be the spatial oriented links (`|E_s| = 3L³`), and `𝒜_s = G^{E_s}` the spatial gauge-field space with `G = SU(2)`.

**Transfer matrix.** The isotropic Wilson action
`S = β Σ_p (1 − ½ Re Tr U_p)`, `β = 4/g²` (SU(2)), factorizes over time slices as `e^{-S} = ∏_t T(U_t, U_{t+1})` with the (finite-`a`, Lagrangian) transfer kernel
`T(U,W) = e^{(β/2) M(U)}\, e^{(β/2) Σ_{ℓ} Re Tr(U_ℓ W_ℓ†)}\, e^{(β/2) M(W)}`,  `M(U) = Σ_{spatial p} Re Tr U_p`,  `U,W ∈ 𝒜_s`.  (Convention: spatial plaquettes split half–half between the two adjacent slices; the precise split only changes the `O(β)` constant, not its structure.)

**Perron–Frobenius ground vector `h₀` — `[T]` (restated from `PROOF-GAUSS-INVARIANT-PF`).** `T(U,W)>0` for all `U,W` (electric heat-kernel positivity × `e^{-a·magnetic}>0`), and `T` is compact and positivity-improving on `L²(𝒜_s, Haar)`. Hence `λ₀>0` is simple and the unique positive eigenvector `h₀>0` solves
> `(T h₀)(U) := ∫_{𝒜_s} T(U,W)\,h₀(W)\,dW = λ₀\,h₀(U)`.   `[T]`

**Local Harnack bound.** For `U,V ∈ 𝒜_s` differing on **one** spatial link `ℓ₀` (`U_{ℓ}=V_{ℓ}` for all `ℓ≠ℓ₀`), the bound is
> `| log h₀(U) − log h₀(V) | ≤ C(β)`, with `C` **independent of L** and explicit in `β`.   (★)

---

## 2. The non-circular mechanism — a kernel-ratio Harnack lemma

> **Lemma (`[T]`, exact).** Let `T(U,W)>0`, `h₀>0`, `T h₀ = λ₀ h₀`. Then for any `U,V`,
> `| log h₀(U) − log h₀(V) | ≤ sup_{W} | log T(U,W) − log T(V,W) |`.  (♥)

**Proof.** Since `T(V,W)h₀(W)>0`, define the probability `ν(dW) = T(V,W)h₀(W)\,dW / \big(λ₀ h₀(V)\big)` (the denominator is `>0`). Then
`h₀(U)/h₀(V) = ∫ T(U,W)h₀\,dW / (λ₀h₀(V)) = ∫ [T(U,W)/T(V,W)]\,ν(dW) = E_ν[ T(U,W)/T(V,W) ]`.
A positive random variable `X` obeys `inf log X ≤ log E[X] ≤ sup log X`; applied to `X = T(U,W)/T(V,W)` this gives `log(h₀(U)/h₀(V)) ∈ [inf_W log ratio, sup_W log ratio]`. Swapping `U↔V` gives the other side. ∎

**This proof uses only `T h₀ = λ₀ h₀` and `h₀>0` — no spectral gap, no correlation decay, no volume.** Hence **[T] and non-circular**: the Wasserstein doc's claim (`§4`) that proving (★) "requires controlling how a one-link change propagates into `h₀`, i.e. that correlations decay (massive)" is **incorrect**. The bound is a *compactness + positivity* fact, not a clustering fact.

**Exact SU(2) evaluation of the right-hand side — `[T]`/`[V-here]`.** Write `U = q_U`, `V = q_V` as unit quaternions; the identity `½ Tr(U W†) = ⟨q_U, q_W⟩` (4D dot) is verified to machine precision `[V-here]` (max error `1.1×10⁻¹⁶` over 2000 pairs), with distance `cos ψ = ⟨q_U,q_V⟩`, `|q_U−q_V| = 2 sin(ψ/2)`, diameter `ψ=π`. A one-link change affects only:
- the **4 spatial plaquettes** of the `U`-slice through `ℓ₀`: `½ReTr U_p = ⟨q_{U_{ℓ₀}}, q_{S_p}⟩` (product of unit quaternions), so `|Δ(½ReTr U_p)| = |⟨q_U−q_V, q_{S_p}⟩| ≤ |q_U−q_V| = 2sin(ψ/2)`; summed over 4 plaquettes and doubled into `M = Σ ReTr U_p = 2Σ½ReTr`, `|Δ((β/2)M)| ≤ (β/2)·4·2·|q_U−q_V| = 8β\,sin(ψ/2)`;
- the **1 temporal plaquette** `ℓ₀` to the `W`-slice: `Δ((β/2)ReTr(U_{ℓ₀}W†)) = β\,⟨q_U−q_V, q_W⟩`, so `sup_W |Δ| = β\,|q_U−q_V| = 2β\,sin(ψ/2)` (sup attained at `q_W ∥ q_U−q_V`; verified `[V-here]`, sampling ratios 0.996–0.9997).

Bounding `|log T(U,W)−log T(V,W)| ≤ |spatial| + sup_W|temporal|` yields

> `C(ψ) = (8β_s + 2β_t)\,sin(ψ/2)`,  isotropic `C(ψ)=10β\,sin(ψ/2) ≤ 10β`.   **(★′) `[T]`**

Every factor is a **local count** (4 spatial + 1 temporal plaquette) — **no `L`, no `a_t` in the constant** beyond `β`. Hence:

| property of (★′) | at β=2.5 | status |
|---|---|---|
| finite | worst-case `C_max=10β=25`; quarter-turn `17.7` | `[V-here]` |
| `L`-uniform | no `L` in `C` | `[T]` |
| non-circular | uses only `Th₀=λ₀h₀`, `h₀>0` | `[T]` |
| `β`-dependence | `≤10β` (at `β<1/48`: `<0.21`) | `[V-here]` |

**β-table `[V-here]`:** `C_max(β)/β = 10`; `C_max = 0.208` at `β=1/48`, `=25` at `β=2.5`, `=30` at `β=3.0`. (The mean-field *typical* per-link change is smaller — `≈0.815`, `O(1)`, reproduced in `[V-here]` §D and matching `WASSERSTEIN-OLLIVIER §3b`'s `0.82` — but the *worst-case sup* used by Dobrushin is `10β`.)

**Remark on the Hamiltonian limit (honesty).** For the finite-`a` Lagrangian kernel used above, (★′) is a clean `[T]`. In the `a_t→0` Hamiltonian kernel `T = e^{-a_t H}`, PF gives `T(U,W)/ (h₀(U)h₀(W)e^{-a_tE₀}) = 1 + O(e^{-a_tΔ})`, so the *ratio* converges to the `h₀`-ratio plus a **gap-controlled** correction — sharpening (♥) to the exact Harnack there would need `Δ`. As an *upper* bound the RHS is still finite by compactness, but the explicit `O(β)`, `L`-uniform form (★′) is stated for the finite-`a` kernel, which is the object the mass gap and `P3` are defined on. This nuance does not affect the kill below (it only helps the bound at finite `a`).

---

## 3. The kill — a finite Harnack does NOT give a contraction at β=2.5

To be *useful* to the Ollivier/W₁ route, a Harnack bound must convert to a **Dobrushin contraction**. Two kernels with log-density ratio in `[−C,C]` have `TV ≤ 1−e^{−C}`; a site with `6` neighbours (3 spatial directions × 2, in `d_s=3`) then has Dobrushin row sum `R = 6(1−e^{−C})`, and contraction `⟺ R<1 ⟺ C < −ln(5/6) = 0.1823`.

> **Contraction threshold — `[V-here]`.** With `C=10β`: contraction iff `β < 0.0182 < 1/48`. **At β=2.5: `R ≈ 6(1−e^{−25}) ≈ 6 ≫ 1` — NO contraction** (Dobrushin condition violated by `≈6×`).

| β | `C=10β` | `TV≤1−e^{−C}` | row sum `6(1−e^{−C})` | contraction? |
|---|---|---|---|---|
| 1/60 | 0.167 | 0.154 | 0.921 | **YES** |
| 1/48 | 0.208 | 0.188 | 1.128 | NO |
| 0.10 | 1.00 | 0.632 | 3.79 | NO |
| 2.50 | 25.0 | 1.000 | 6.00 | NO |

**Reading.** The kernel-ratio Harnack *does* yield a genuinely non-circular Dobrushin contraction — but **only in the strong-coupling cell `β ≲ 1/55 < 1/48`**, i.e. *exactly where the gap is already a theorem* (`P3`, `PROOF-P3-UNCONDITIONAL`; independently `GIRSANOV.md`, `LOC-BLIND.md` reach the same "only non-circular Dobrushin cell = `β<1/48`"). At the crossover `β=2.5`, `α_s~O(1)`, `C=Θ(β)` is **O(1), no small parameter** (`WASSERSTEIN-OLLIVIER §3b` said exactly this), the row sum is `≥1`, and the global contraction fails — the same local→global wall as every prior route. **So the bound is useless precisely at the crossover where the gap is not known.**

---

## 4. Free-field acid test — the mechanism is VACUOUS `[V-here]`

A tool must **fail** on the free massless field, else it is vacuous (charter §1). Test a single free mode: transfer `T(x,y)=exp(−(x−cy)²/(2s²))`, `c=e^{−aω}` (`ω` = mode frequency = **the gap scale**), PF vector `h₀(x)=exp(−ωx²/2)`.

1. **Kernel-ratio, non-compact field (`x,y ∈ ℝ`).** `log(T(x,y)/T(x',y)) = −(x−x')(x+x'−2cy)/(2s²)` grows **linearly in `y`**; `sup_y` **diverges**. `[V-here]`: for `x=0,x'=1`, the sup on `|y|≤Y` is `4.47, 14.4, 34.4, 74.4, 154, 314` for `Y=1,2,4,8,16,32` — unbounded. The Harnack bound **fails on the non-compact massless field** — but because the field space is **non-compact**, *not* because it is gapless.

2. **Harnack of `h₀` on a window `|x|≤M`.** `C_H = sup| \tfrac{ω}{2}(x²−x'²)| = \tfrac{M²}{2}\,ω = \tfrac{M²}{2}\,Δ`. **`C_H` is linear in the gap.** `[V-here]`: `ω=1.0→C_H=0.5`, `ω=0.1→0.05`, `ω→0→C_H→0`. The **massless** limit has the **smallest** Harnack constant (flat `h₀`) yet **no gap**.

3. **Compact free field (compact U(1) heat kernel, Coulomb phase).** The kernel-ratio sup is **finite** for every `β` (`[V-here]`: `≤ 2β`), even small — yet compact U(1) is **gapless** (photon).

**Acid test FAILED.** The kernel-ratio/Harnack mechanism is **vacuous** as a gap certificate:
- it is finite for **gapless compact** theories and infinite only for **non-compact** ones — it measures **compactness / field-range**, not the gap;
- for free modes `C_H ∝ Δ`, so *small* `C_H` ⇔ *small* gap — the **wrong sign** for a gap **lower** bound (and (★′) is an upper bound `C≤10β`, satisfied by `Δ→0`).

So a finite Harnack constant coexists with a gapless theory, and its size does not track the gap. `[V-here]`

---

## 5. What this corrects — circularity vs vacuity

`WASSERSTEIN-OLLIVIER.md` contains a tension. `§3b`/`§3a`/`§6` correctly conclude that the constant is `O(1)` at β=2.5, that the **local** row sum is `<1` (`0.096`) but the **worst-case/global** row sum is `≥1` (`1.09–3.27`), and that the route "dies on the local→global upgrade" (= H1). But `§4` attributes the obstruction to **circularity**: it says proving the local Harnack "requires controlling propagation … = H1." **This round corrects `§4`:** the absolute local Harnack is a *one-line* consequence of `T h₀=λ₀h₀` and `h₀>0` (Lemma (♥)), with explicit `L`-independent constant `≤10β` — it is **[T] and non-circular**, and it does **not** need correlation decay. The obstruction is **not** circularity of the *Harnack*; it is:
- **(a) the Harnack→contraction upgrade** (Dobrushin `q<1`, equivalently the W₁ global contraction) — which fails at β=2.5 (`R≈6`), = H1; and
- **(b) vacuity** — a finite Harnack holds for the gapless compact free field (`§4`).

Equivalently: even if `h₀` were *handed to us*, the W₁ contraction still fails globally at β=2.5 (row sum `≥1`), so the "`h₀`-circularity" was never the true gate — the gate is the local→global contraction. (This is consistent with `§6`'s own admission that points (i) and (iii) "are the same crossover control," and with `GIRSANOV.md`'s "only non-circular Dobrushin cell = `β<1/48`.")

---

## 6. Exact missing step `[O]`

For the absolute Harnack to be a *genuine opening* one would need, instead of (★),

> **[O-Harnack²]** either (i) a **relative** Harnack `|log h₀ − log h̃| ≤ C` against a **computable** proxy `h̃` (e.g. the mean-field `∏p_λ`, `λ=6β`) that does **not** invoke the gap — this would de-circularize the Doob kernel by controlling `P = T(h₀·)/λ₀` up to `e^{2C}` against the *known* proxy kernel; **and** (ii) the **same** proxy kernel must satisfy a **global** W₁/Dobrushin contraction `q<1` at β=2.5 (row sum `<1`), i.e. the local→global upgrade.

Neither (i) nor (ii) is a consequence of the absolute bound; both reduce to the crossover contraction **H1** (the very object). The mean-field proxy has row sum `≥1` too (`§3`; `WASSERSTEIN-OLLIVIER §3a`), so (ii) is the non-perturbative crossover control. **[O]**

---

## 7. Scorecard (charter's 5 hard conditions, `BRIEF-FOR-MODEL.md §2`)

| # | condition | result |
|---|---|---|
| 1 | non-perturbative | **PASS** — pathwise kernel ratio, no loop expansion |
| 2 | sign-controlling lower bound | **FAIL** — Harnack is an *upper* bound on `log`-variation; gives no gap lower bound |
| 3 | non-circular | **PASS (for the Harnack itself)** — `[T]`, only `Th₀=λ₀h₀`; but the *upgrade* to a contraction is circular (=H1) |
| 4 | non-SUSY + **fails on free massless** | **FAIL (vacuous)** — finite & even `→0` on the gapless compact free field; free-mode `C_H∝Δ` (anti-diagnostic) |
| 5 | `(L,a)`-uniform bound | **PASS** — `C≤10β`, no `L`; but useless at β=2.5 (row sum `≈6`) |

**Sharpest honest line.** The local vacuum Harnack bound on `h₀` is **not** the circular object the Wasserstein route feared: it is a compactness+positivity fact, provable in one line, with an explicit `L`-independent constant `C=(8β_s+2β_t)sin(ψ/2)≤10β` (`[T]`). But a finite/`L`-uniform Harnack bound is **vacuous** for the gap — it is `Θ(β)` (no small parameter at `α_s~O(1)`; Dobrushin row sum `≈6≫1` at β=2.5), finite for the *gapless* compact free field, and, for free modes, *proportional* to the gap (wrong direction for a lower bound). The true gate is the local→global W₁ contraction = **H1**, which the Harnack bound does not touch. **Verdict: (b) — no genuine non-circular L-uniform Harnack opening at β=2.5 — but the exact reason is *vacuous + Θ(β) non-smallness + row-sum≥1 at the crossover*, NOT L-blowup (it is L-uniform) and NOT literal circularity (it is provable).** Clay(G2) 0%; no QED.

---

## 8. Standard results cited `[R]`

- **Perron–Frobenius / Krein–Rutman** for positivity-improving compact kernels: Jentzsch; Krein–Rutman (in `PROOF-GAUSS-INVARIANT-PF`).
- **Dobrushin contraction coefficient** and the row-sum criterion: R. L. Dobrushin, *Theory Probab. Appl.* **13** (1968) 197; **15** (1970) 458. TV vs log-density-ratio bound `TV≤1−e^{−‖log(dP/dQ)‖_∞}` is elementary (Pinsker/Csiszár).
- **Harnack inequalities for Markov kernels / semigroups:** a kernel-ratio (sup-norm) Harnack is the discrete analog of the parabolic Harnack inequality (Moser, *Comm. Pure Appl. Math.* **17** (1964) 101; Fabes–Stroock); the "log of the ground state varies by ≤ sup log kernel-ratio" step is the Doob-`h`/ground-state-transform Harnack (Doob, *Classical Potential Theory*; for the elliptic ground-state pinching see also Bakry–Gentil–Ledoux, *Analysis and Geometry of Markov Diffusion Operators*, Springer 2014, §Lindqvist/ground-state Harnack).
- **Bakry–Émery / Li–Yau alternative routes** (checked, not used): Bakry–Émery `Γ₂`, Li–Yau gradient estimates require `Hess⪰−K`; the project's Savvidy obstruction (`IDEA-wasserstein.md`, `GIRSANOV.md`) makes `Γ₂` indefinite at β=2.5, so the *curvature* Harnack is dead where the *kernel-ratio* Harnack (this note) survives — but the latter is vacuous, as shown.
- **Strong-coupling / Dobrushin-uniqueness window** `β<1/48` (P3): Osterwalder–Seiler, *Ann. Phys.* **110** (1978) 440; SZZ (Schor), and the project's `PROOF-P3-UNCONDITIONAL.md`. Cross-check: `GIRSANOV.md`, `LOC-BLIND.md`, `P4-RP-uniform-gap.md` (quantitative PF constants `~e^{−cV}` for the *spectral* gap; the *local* Harnack is `L`-independent — consistent, as they are different objects).

*Internal refs.* `WASSERSTEIN-OLLIVIER.md` (§3a/§3b/§4/§6), `PROOF-GAUSS-INVARIANT-PF.md`, `GIRSANOV.md`, `LOC-BLIND.md`, `P4-RP-uniform-gap.md`, `PROOF-P3-UNCONDITIONAL.md`, `docs/BRIEF-FOR-MODEL.md §2`, `docs/PROOF-SKELETON.md §2`. Witness: `research/unconditional/w7/harnack_killtest.py`.
