# SLOT-PF-FIXPOINT — can the gauge-invariant Perron–Frobenius structure *itself* be the `k→0` sign-definite self-consistency? The single positive cone is **sign-blind to `λ₁`**; strict positivity forces `λ₁<λ₀` only *per finite `(V,a)`*, and its one sign-definite lower face (Cheeger/Poincaré) **is `H1`**

> 2026-11-24 · Route: **gauge-invariant transfer-matrix PF / Collatz–Wielandt / Doob / Hilbert-metric**. Target: the ONE open slot of `PROFILE.md` §3 — a **`k→0` sign-definite self-consistency, no gapless fixed point, uniform at the crossover, non-positivity, nonabelian, acid-clean**. Witness: `research/unconditional/w7/slot_pf_fixpoint.py` (`py -3.12`, single job, pure numpy/scipy).
>
> Labels (charter §2): **[T]** 真定理 · **[C]** 条件定理 · **[O]** 精确缺失步 · **[NG]** no-go for a *named class* · **[V-here]** computed this round · **[R]** inherited/banked.

---

## 判词 (Verdict, up front)

> **We make the PF structure its own self-consistency on `ρ:=λ₁/λ₀` in its three native forms — Collatz–Wielandt, Doob ratio-map, Hilbert-metric (Birkhoff) — and pin exactly what it can and cannot do.** **`[T]`** At each *finite* `(V,a)` strict positivity **does** force separation: `λ₀` simple, `Ω>0`, `λ₁<λ₀`, `ρ<1`, `Δ_V>0` (Jentzsch/Krein–Rutman; `PROOF-GAUSS-INVARIANT-PF`). **`[T]`** But the single PF cone is **sign-blind to `λ₁`**: every eigenvector `ψ` with `λ<λ₀` is `Ω`-orthogonal with `Ω>0`, hence **sign-changing**, hence **outside** the positive cone — so Collatz–Wielandt, Birkhoff/Hilbert contraction and Dobrushin-in-the-oscillation-cone (all of which act on `int C`) **cannot address `λ₁` at all**; the deflated operator `T_⊥=T−λ₀P_Ω` (whose norm **is** `λ₁`) is **sign-indefinite** and CW on it collapses. **`[V-here]`** `ρ=1` **is** admissible as a *uniform* limit with a strictly positive kernel: the free massless heat kernel is strictly positive (PF applies exactly) with `λ₀=1` simple and `ρ_L=e^{−a μ₁(L)}→1`; a two-block kernel with cross-coupling `δ>0` is strictly positive at every `δ` yet `1−ρ→0`. **`[NG]`** So strict positivity forces a gap **only per finite `(V,a)`, with no uniform modulus**, and it is **abelian-blind** (U(1) has the same strictly positive kernel — IL-ii fails). **`[O]`** The one positivity-native **sign-definite lower** bound is Cheeger `1−ρ ≥ Φ²/2` (conductance `Φ`, a ratio of positive measures); on the free field `Φ→0` (`[V-here]`), and proving `Φ≥δ` uniformly at `β≈2.5` **is exactly the coarse Poincaré `=H1`**. **Net: the PF self-consistency cannot be written as a nontrivial fixpoint `R(ρ)` at all (it is either the definitional gap — Schema D — or an upper-bound contraction — `WASSERSTEIN-OLLIVIER`); the exact failing step is the sign-blindness of the single cone ⇔ the absent second cone, and the uniform statement reduces to `H1`. Clay 0%. No QED.**

---

## §0 · Setting, labels, iron laws, acid test

**Object.** `T=e^{−aH}` on `H_phys=ker Gauss=L²(𝒜/G, μ̄)`, finite spatial volume `V`, spacing `a`. Kernel `T̄(ω,ω')>0` **for all** `ω,ω'` (`PROOF-GAUSS-INVARIANT-PF §1`, from the electric heat kernel times `e^{−a·magnetic}>0` — **full positive Haar, no Gribov/FP sign**). `T̄` self-adjoint, Hilbert–Schmidt at finite `V`. Eigenvalues `λ₀>λ₁≥λ₂≥…` (so `λ₀=e^{−aE₀}`). Set
```
   ρ := λ₁/λ₀ ∈ (0,1) ,      Δ := −(1/a) ln ρ = E₁−E₀   (the gap; the 0⁺⁺ overall gap in pure YM).
```
Target (`PROFILE §3`, `PROFILE-FACET-ANCHOR`): a relation `R(Δ,H)` that (P4) is **intrinsic/controller-free**, **nonlinear in the measure**, **IR-support-sensitive**, and **excludes `Δ=0`** by a **`k→0` sign-definite mechanism**, uniformly in `(V,a)` at `β≈2.5`.

**The three iron laws (`KEY-DIRECT-DERIVATION §0`).** Any lower bound must be **(IL-i) non-positivity** — the free massless field must *fail* it (`Φ(free)>0` ⇒ acid-dirty / vacuous); **(IL-ii) nonabelian** — else the U(1) massless no-go; **(IL-iii) no small parameter** — must hold at `α_s~O(1)`, `β≈2.5`. **Acid test:** clean iff `Φ(free massless)≤0`, *discriminating* iff it also gives `≤0` on the gapped-but-atom-free `:φ²:`.

---

## §1 · The PF self-consistency: three native forms, and why none is a genuine fixpoint in `ρ`

### 1.1 The exact PF facts — `[T]` (Jentzsch / Krein–Rutman, finite `V`)

Strictly positive, self-adjoint, compact `T̄` on the compact `𝒜/G` gives: **(i)** `λ₀=‖T̄‖` is **simple** with eigenvector `Ω>0` everywhere; **(ii)** the rest of the spectrum is a sequence `λ₁≥λ₂≥…→0`; **(iii)** hence `ρ<1` and `Δ_V>0`. This is the finite-volume gap and it is a genuine theorem. (`PROOF-GAUSS-INVARIANT-PF §2`; `INEQ-R1`; witness §A.) It yields `Δ_V` as an **upper** bound on the continuum gap (`Δ_∞=inf_V Δ_V`), never a lower one.

### 1.2 Form 1 — Collatz–Wielandt, and the deflation that kills positivity — `[T]` + `[V-here]`

For a **nonnegative** operator, the sign-definite self-consistency for the TOP eigenvalue is exact:
```
   λ₀ = inf_{f>0} sup_x (Tf)(x)/f(x) = sup_{f>0} inf_x (Tf)(x)/f(x),     achieved at f=Ω.        (CW)
```
`[V-here]` §B: on the Gaussian toy, `sup_{f>0} inf_x (Tf)_x/f_x = 1.704307 = λ₀` exactly (the Perron vector is a trial). **The naive hope** is to apply CW to the **deflated** operator on `(vac)^⊥`,
```
   T_⊥ := T − λ₀ P_Ω ,   P_Ω = |Ω⟩⟨Ω| .
```
The spectral fact is clean (`[T]`): `T_⊥Ω=0` and `T_⊥=T` on `Ω^⊥`, so `‖T_⊥‖=λ₁` (**exactly** the gap eigenvalue). But `T_⊥` is **sign-indefinite**: its kernel `T_⊥(ω,ω')=T̄(ω,ω')−λ₀Ω(ω)Ω(ω')` changes sign wherever `λ₀ΩΩ` exceeds `T̄`. `[V-here]` §B: `min(T_⊥)=−0.304`, **44 of 64 entries negative**. Consequently CW **does not apply** to `T_⊥` (CW needs `≥0`), and numerically the CW functional **collapses**: `sup_{f>0} inf_x (T_⊥f)_x/f_x = 0.000 ≪ λ₁=1.561`. **The deflation that exposes `λ₁` destroys the sign-definiteness that CW needs.** This is the precise structural failure of the task's proposed mechanism.

### 1.3 Form 2 — Doob ratio-map: triangular, not a feedback fixpoint — `[T]`

Normalize by the Perron vector, `P := T(Ω·)/λ₀` — a **stochastic** kernel on `𝒜/G` (`P1=1`, `spec(P)={λ_k/λ₀}`), the classical Doob `h`-transform. The nonlinear normalized power map
```
   Φ(h) := Th / ‖Th‖ ,      Φ(Ω)=Ω ,
```
has `ρ` as the spectral radius of its linearization `DΦ(Ω)=(I−ΩΩ*)T` on `Ω^⊥`; i.e. `ρ` is the **rate of the power iteration** toward the vacuum. **But note the triangularity:** `Ω` is fixed by the **linear** eigenproblem `TΩ=λ₀Ω` and depends on `λ₀` alone — **not** on `λ₁`. So `ρ` is *read off* the next eigenvalue and **never feeds back** into `Ω`. There is **no intrinsic relation `ρ=F(ρ)`** with `ρ` on both sides; the only "self-consistency" for `ρ` is the deflated **Rayleigh min–max**
```
   ρ = sup_{E_π g=0, g≠0} ⟨g,Pg⟩_π / ⟨g,g⟩_π ,          π = Ω²μ̄/Z  (stationary) ,
```
which is a variational characterization (Schema D), **not** a nonlinear feedback equation. **Consequence:** the slot's desired shape — an relation `R` whose solution set *intrinsically excludes* `ρ=1` — **cannot even be posed** in the PF structure: either it is the definitional gap (`1−ρ>0 ⇔` gap, Schema D tautology), or it is an upper-bound contraction (§1.4).

### 1.4 Form 3 — Hilbert-metric / Birkhoff contraction: also blind to `λ₁` — `[T]`

`T` is a contraction in the **Hilbert projective metric** `d_H` on `int C`: `d_H(Tf,Tg) ≤ κ(T)d_H(f,g)`, `κ(T)=tanh(D(T)/4)<1` (`D`=projective diameter of `T(C)`), giving geometric equidistribution of normalized iterates toward `Ω` (Birkhoff; Hopf). **But `d_H` is defined only on `int C`, and the `λ₁`-eigenvector `ψ₁` is sign-changing (§1.5), so `ψ₁ ∉ int C`.** Hence Birkhoff's contraction carries **no information about `λ₁`**. This is the same wall as `WASSERSTEIN-OLLIVIER` (whose `W₁` coefficient `q≥λ₁/λ₀`, sharp only on the free field) seen from the kernel side: the positive-cone contraction governs convergence to the **vacuum**, not the separation of the vacuum from the next state.

### 1.5 Central Lemma — sign-blindness of the single positive cone — `[T]`

> **Lemma SB.** Let `T` be strictly positive with Perron vector `Ω>0` (μ̄-a.e.). If `ψ` is any other eigenvector (`λ≠λ₀`, real, `T` self-adjoint), then `⟨ψ,Ω⟩_{L²(μ̄)}=0`; since `Ω>0`, `ψ` takes **both signs** (were `ψ≥0` a.e., `ψΩ≥0` and `∫ψΩ=0` force `ψ=0`). Hence `ψ∉int C`.
> **Corollary (blindness).** Every positive-cone functional — Collatz–Wielandt, Birkhoff/Hilbert, Dobrushin's oscillation-coefficient, PF vacuum-dominance — acts on `int C` and therefore **cannot see `λ₁`**; the deflated operator that exposes it (§1.2) is sign-indefinite. **The single PF cone pins the top of a channel (the vacuum); it is structurally incapable of measuring the separation of two states within one channel.**
> ∎ (`INEQ-R1` "second cone absent"; `PROOF-GAUSS-INVARIANT-PF §3b` "the glueball sits in `Ω^⊥`, the one face where the RP cone degenerates" — here sharpened from the *atom* to the *gap* `1−ρ` itself.)

This is the central reason the PF route cannot close: the gap `λ₁` needs a **second invariant cone** / total positivity (Karlin sign-regularity, `ψ₁` interior to a flag cone) on `𝒜/G`, and `𝒜/G=G^E/G` **carries no order** (`PROOF-F1SQ iii`), so such a cone does not exist; alternatively it needs a **sign-sensitive** (oscillation / `L²` / Cheeger) functional — which is §5.

---

## §2 · Gapless-fixed-point test — is `ρ=1` admissible? — `[T]` / `[V-here]`

### 2.1 At finite `(V,a)`: NO — `[T]`
PF gives `λ₀` **simple**, so `λ₁<λ₀`, `ρ<1`, `Δ_V>0` at **every** finite `(V,a)`. `ρ=1` (a degenerate ground level, two positive ground vectors) is **forbidden** by Jentzsch. **This is the one genuine positive result: strict positivity forces a finite-volume gap.**

### 2.2 Uniformly (over `V`, or `a→0`): YES — `[V-here]`
Strict positivity is **compatible** with `ρ→1`, for two explicit families:
- **Free massless field.** `T=e^{−aM}`, `M=−Δ_lat` on the `L³` torus, is a **manifestly strictly positive heat kernel** (`min T>0`, witness §C; `Ω`=constant mode). `λ₀=1` simple, `λ₁=e^{−a μ₁}`, `μ₁=4 sin²(π/L)`. Then
  `ρ_L=e^{−a μ₁} → 1` as `L→∞` (`[V-here]`: `L=2:ρ=0.018`, `4:0.135`, `8:0.557`, `→1`; energy gap `2 sin(π/L)=2.00,1.732,1.414,…→0`, matching `FINITE-VOL-CERTIFY`). **Strictly positive PF kernel, gapless in the continuum.**
- **Two-block kernel** (cross-coupling `δ>0`, strictly positive at every `δ`): `ρ=(1−δ)/(1+δ)`, `1−ρ=2δ/(1+δ)→0` while `min K=δ>0` (`[V-here]` §A table: `δ=10^{−1}:1−ρ=1.8×10^{−1}` … `δ=10^{−6}:1−ρ=2×10^{−6}`).

So `ρ=1` **is** a solution of "strictly positive PF" **in the uniform limit**. The PF self-consistency does **not** exclude the gapless fixed point; it is satisfied by the free field at `ρ→1`.

---

## §3 · Does strict positivity force separation? — split answer

- **Per finite `(V,a)`: YES — `[T]`.** `λ₁<λ₀` strictly (Jentzsch). A finite-volume gap **exists**.
- **Uniformly: NO — `[NG]`.** There is **no modulus** `1−ρ ≥ δ` derivable from positivity: the two-block family and the free field are explicit sequences with `1−ρ→0` and `min`-kernel `>0` (§2.2). Positivity supplies separation *without a rate*; the rate is a different (sign-sensitive) object (§5).

---

## §4 · Uniform attainment at `β≈2.5` — not a PF property — `[V-here]` / `[O]`

The only `[V-here]` positivity datum at `β=2.5` is **local** (single-link/plaquette), not the coarse mode:
- single-plaquette (2D SU(2) heat kernel) `λ_{1/2}/λ₀ = I₂(β)/I₁(β) = 0.507` at `β=2.5` (local gap `0.493`);
- single-link heat bath, aligned 6-staple background `β_eff=6β=15`, `I₂/I₁(15)=0.902` (local gap `0.098`) — **slow**, and `→1` as `β_eff` grows.

These are **local** relaxation ratios; the **coarse-mode** (the actual `ρ`) is **not** controlled by them. The block gap is a **fast** relaxation rate, silent at `k=0` (`PROOF-COARSE-POINCARE §b`: adiabatic elimination controls the memory kernel, never the coarse stiffness; `H2`'s variance bound runs the wrong way by `~e^{5N_p}`). No certified finite-volume interacting gap exists at `L≥2, β=2.5` (`FINITE-VOL-CERTIFY`: `dim_GI(L=2)=1.3×10⁷ ≫` dense ceiling). **Uniform attainment is not a PF property — it is the transport/crossover certificate `=H1`.** `[O]`

---

## §5 · The one sign-definite LOWER face — Cheeger — and its reduction to `H1`

The PF structure's *only* positivity-native **lower** bound on the gap is the conductance (Cheeger) inequality. For the reversible Doob chain `P` (stationary `π`), with `Φ := inf_A Q(A,A^c)/min(π(A),π(A^c))` the **conductance** (`Q` = Dirichlet form; `Φ` a **ratio of positive measures**, maximally sign-definite):
```
   1 − ρ  ≥  (1/2) Φ²            (Cheeger lower bound; [T])
   1 − ρ  ≤  2 Φ                 (Cheeger upper bound)  ⇒  g/2 ≤ Φ ≤ √(2g),  g:=1−ρ.
```
- **Acid (`[V-here]` §D):** free field `g=μ₁=4 sin²(π/L)`, `Φ∈[g/2,√(2g)]` → `0` as `L→∞` (`L=2:[2,2.83]` … `L=16:[0.076,0.552]`). The sign-definite lower mechanism manufactures **no uniform mass** on the free field (acid-clean but vacuous).
- **YM:** proving `Φ≥δ` uniformly at `β=2.5` is the **coarse Poincaré** inequality `Var_ν(φ)≤(1/λ)∫|∇φ|²` (Cheeger ⟺ Poincaré up to constants), which is **verbatim `H1`** (`PROOF-COARSE-POINCARE`: "`G″≥λ` uniformly ⟺ `m_eff²(k→0)≥λ` ⟺ `Var_ν(φ)≤(1/λ)∫|∇φ|²`"). So the PF self-consistency, pushed to its one sign-definite **lower**-bound face, **reduces exactly to `H1`.** `[O]`/`[NG]` for the uniform closure.

---

## §6 · Acid test + the three iron laws

Let `Φ_PF(V,a) := 1−ρ_V = 1−e^{−aΔ_V}` be the PF "gap functional."

- **IL-i (non-positivity): FAIL.** On the free massless field `Φ_PF^{free}=1−e^{−aE₁(L)}>0` at **every finite `L`** (the *box* gap — exactly `FINITE-VOL-CERTIFY`'s "a finite-volume gap is not the Clay gap"). So the PF functional certifies a **positive** gap on a continuum-gapless field: **acid-dirty at fixed `(V,a)`**; it only reaches `0` in the `L→∞` limit, where it is **vacuous**. The mechanism does **not distinguish** gapped YM from gapless free field at any finite `(V,a)`. `[T]`/`[V-here]`
- **IL-ii (nonabelian): FAIL.** The strict positivity holds for **U(1)** as well (free Maxwell has the same strictly positive kernel and a finite-V box gap). The PF mechanism is **abelian-blind** — and U(1) is massless. Positivity cannot be the nonabelian discriminator. `[T]`
- **IL-iii (no small parameter): fails as a bound.** Exact at all `β` in form, but the **uniform** coarse-mode gap at `β=2.5` is unproven (§4) — the very statement is the wall. `[O]`

So PF fails **two** of the three iron laws structurally (IL-i acid-dirty at finite `V`; IL-ii abelian-blind), and its IL-iii uniform form is `H1`. This is a sharper and more specific diagnosis than "the second cone is absent": **positivity is an abelian-blind, finite-`V`-box statement.**

---

## §7 · Verdict — constructed or the exact failing step

**`[T]` constructed positive core.** At finite `(V,a)`, strict positivity **does** force a gap (`λ₀` simple ⇒ `λ₁<λ₀` ⇒ `Δ_V>0`); this is real (`PROOF-GAUSS-INVARIANT-PF`, witness §A,C).
**`[NG]` for the uniform PF closure.** No sign-definite relation with `ρ` on both sides exists (triangular, §1.3); the single cone is **sign-blind to `λ₁`** (Lemma SB, §1.5); `ρ=1` is admissible as a uniform limit with strict positivity (free field, two-block; §2.2); PF is abelian-blind and acid-dirty at finite `V` (§6).
**`[O]` the exact missing step.** The one sign-definite **lower** face is Cheeger `1−ρ≥Φ²/2`; uniform `Φ≥δ` at `β=2.5` is the **coarse Poincaré `=H1`** (§5). Equivalently, what is missing is a **second invariant cone / total positivity** on the orderless `𝒜/G` — `ψ₁` interior to a flag cone with `λ₁` separated — which `PROOF-F1SQ iii` proves does not exist.

**Is `ρ=1` a solution?** finite `(V,a)`: **NO** (Jentzsch). uniformly: **YES** (free field / two-block). **Does strict positivity force separation?** per-`(V,a)`: **YES `[T]`**; uniformly: **NO `[NG]`**.

**Relation to prior no-gos (non-redundant steps).**
- vs **`WASSERSTEIN-OLLIVIER`**: that route uses the `W₁` contraction coefficient `q≥λ₁/λ₀` (an **upper** bound on `ρ`, sharp only on the free field — a *sharpness deficit*). Here the object is the **single PF cone** and the obstruction is **sign-blindness** (Collatz–Wielandt/Birkhoff undefined on the sign-changing side; deflation sign-indefinite) — a different, complementary *why*. Both reach `H1`, via distinct precise steps.
- vs **`INEQ-R1`/`PROOF-GAUSS-INVARIANT-PF §3b`**: those address the **atom `f₁²`** (PF cannot reach the edge *atom*). Here we address the **weaker gap `1−ρ`** and show even *that* is invisible to the single cone *uniformly*. A gap-level refinement.
- vs **`PROFILE-FACET-ANCHOR`/`META-GAPFORCING`**: consistent — PF is a self-anchoring **READ** (not FORCE) object; §5 makes the "READ⇒Poincaré=H1" reduction explicit through Cheeger.

**Sharpest line.** *Strict positivity is strong enough to make the vacuum unique and the finite box gapped, and nothing more: the very eigenvector that carries the gap is the one the positive cone cannot contain, the deflated operator that exposes it is sign-indefinite, and the free massless field — strictly positive, PF, `ρ→1` — passes every positivity test at every box size. The gap is a sign-sensitive isoperimetric (Cheeger/Poincaré) statement, and that is `H1`.*

**Honest boundary.** `[NG]` is for the *named class* (uniform PF / positivity-cone self-consistency), not a global impossibility (`∀class ≠ ∀possible`). `[T]` the finite-`V` gap + Lemma SB; `[V-here]` the numeric witness; `[O]` the uniform step = coarse Poincaré = `H1`. **Not a proof of Clay. Clay 0%. No QED.**

---

## §8 · Reproduce + references

```
py -3.12 research/unconditional/w7/slot_pf_fixpoint.py       # single job, ~2 s
```
Sections: A PF/Jentzsch + two-block `1−ρ→0`; B deflation sign-indefinite + CW collapse; C free massless heat kernel `ρ→1`; D Cheeger `Φ→0`; E SU(2) `β=2.5` local ratios; F verdict.

**Internal (absolute paths; banked labels carried).**
- `D:\Dev\math\research\unconditional\PROOF-GAUSS-INVARIANT-PF.md` (strictly-positive kernel; Jentzsch; `Ω>0`; glueball in `Ω^⊥`, the degenerate face).
- `D:\Dev\math\research\unconditional\INEQ-R1-perron-positive.md` (PF anchored at the vacuum; second cone absent; free field `f₁²=0`).
- `D:\Dev\math\research\unconditional\WASSERSTEIN-OLLIVIER.md` (`W₁` contraction `q≥λ₁/λ₀`, sharp only free; sharpness deficit).
- `D:\Dev\math\research\unconditional\PROOF-COARSE-POINCARE.md` (`G″≥λ` uniform ⟺ `m_eff²(k→0)≥λ` ⟺ `Var_ν(φ)≤(1/λ)∫|∇φ|²` = `H1`; block gap ≠ coarse stiffness).
- `D:\Dev\math\research\unconditional\FINITE-VOL-CERTIFY.md` (free box gap `2 sin(π/L)→0`; finite-V gap is not Clay; `dim_GI(L=2)=1.3×10⁷` barrier).
- `D:\Dev\math\research\unconditional\KEY-DIRECT-DERIVATION.md` (§0 iron laws / acid test/ Schur-convex).
- `D:\Dev\math\research\unconditional\PROOF-F1SQ.md` (`𝒜/G` orderless — no total positivity / Gantmacher–Krein cone).
- `D:\Dev\math\research\unconditional\PROFILE.md`, `PROFILE-FACET-ANCHOR.md` (the slot; Schema A–E table; Theorem SA; `READ≠FORCE`).
- `D:\Dev\math\research\unconditional\META-GAPFORCING-IDENTITY.md` (§3.1 identity-forcing no-go; §3.5 `READ≠FORCE`).

**External `[R]`.**
- Jentzsch (1912); Krein–Rutman (1948); Krasnosel'skii, *Positive Solutions of Operator Equations* (1964) — PF, CW.
- G. Birkhoff, *Lattice Theory* (1967); Hopf (1963) — Hilbert projective metric, contraction coefficient `tanh(D/4)`.
- R. L. Dobrushin, *Theory Probab. Appl.* **15** (1970) 458 — ergodicity coefficient.
- J. Cheeger, *Problems in Analysis* (1970) 195; P. Diaconis, D. Stroock, *Ann. Appl. Probab.* **1** (1991) 36 — Cheeger/conductance gap bounds.
- Osterwalder–Seiler, *Ann. Phys.* **110** (1978) 440; Lüscher, *CMP* **54** (1977) 283 — reflection positivity / transfer matrix.
- S. Karlin, *Total Positivity* (1968); Gantmacher–Krein — second cone / sign-regularity (the absent structure).

---

*Status: **PF-route facet — `[T]` finite-`V` gap (Jentzsch) + Lemma SB (sign-blindness of the single cone, sharpened from atom to gap) + `[V-here]` numeric witness + `[NG]` for the uniform PF / positivity-cone self-consistency class + `[O]` the uniform step = coarse Poincaré = `H1`.** `ρ=1` excluded per finite `(V,a)`, admissible as a uniform limit (free field / two-block). Strict positivity forces separation only per finite `(V,a)`, with no modulus, and is abelian-blind (IL-ii) and acid-dirty at finite `V` (IL-i). The deflated operator exposing `λ₁` is sign-indefinite and CW collapses; the only sign-definite lower face (Cheeger) is `H1`. **Not a proof of Clay; Clay 0%; no QED.***
