# WASSERSTEIN-OLLIVIER — Ollivier coarse-Ricci / W₁ contraction of the physical-time transfer: sound engine, but H1 in the strictly-weaker W₁ language

> 2026-11-24 · Test of a specific proposal: does the **Ollivier coarse-Ricci / 1-Wasserstein contraction** of the reflection-positive physical-time transfer (Doob-transformed) kernel `P = T(h₀·)/λ₀` give a `(L,a)`-uniform handle on the physical gap `Δ`?
> Companion witness: `research/unconditional/w7/wasserstein_ollivier.py` (`py -3.12`, single job). All numbers below are `[V-here]`.
> Labels (charter §2): **[T]** true theorem · **[C]** conditional · **[O]** exact missing step · **[NG]** named-class no-go · **[V-here]** computed this round · **[R]** external, verified. **No QED on anything unproven.**

---

## 0. VERDICT — 一句话

**The logical engine is sound and genuinely diameter-free `[T]`. The object it acts on (W₁ / Ollivier coarse-Ricci) is *strictly weaker* than the already-refuted W₂/CD/Bakry–Émery route, so the earlier no-go does *not* automatically transfer, and at β=2.5 the **local** coarse-Ricci curvature is indeed positive (`[V-here]`, local row sum 0.096 ≪ 1) where the W₂/BE curvature is Savvidy-negative. That is the genuine novelty. But the route is *still not a handle to Clay*, for three reasons that the numbers make explicit: (A) the Doob kernel needs `h₀` — circular `[O]`; (B) the Lipschitz bound gives `q ≥ λ₁/λ₀ = e^{−aₜΔ}`, so a proven W₁ contraction certifies only `Δ_cert = −(1/aₜ)log q ≤ Δ_true`, with equality **iff the contraction is exactly sharp** — false for an interacting SU(2) kernel, so it certifies a strict **sub-gap**, not `Δ`; (C) only the *infinitesimal* W₁ row sum is `<1`; the **global worst-case (TV/Dobrushin) row sum is 1.09–3.27 ≥ 1** (`[V-here]`), so a global all-pairs `(L,a)`-uniform W₁ contraction is unproven and is again the crossover control. **Net: H1 recollected in a strictly-weaker, non-redundant W₁ dialect, plus a sharpness deficit; not a new handle. Acid test passes (free massless gives `q→1`, no faked mass).** `[T]/[O]/[V-here]`

**为什么这不是"W₂ 已死 ⇒ W₁ 也死"的重复：** von Renesse–Sturm 的 **sharp equivalence**（`W₂-contraction(λ) ⟺ Hess V⪰λ`）是 **W₂** 定理；W₁ 没有这样的 sharp 等价（W₁ ⊂ W₂ 严格更弱），所以在 W₂/BE 被 Savvidy 杀的同一 β=2.5 点，W₁ 的**局部 coarse-Ricci 仍可为正**——这正是本轮实测（0.096）。此前 `IDEA-wasserstein.md`（W₂/CD 全测度被 Savvidy 证伪）与 `D2-COUPLING-DOEBLIN.md`（W₂ 耦合 ≡ `G″` ≡ H1）都是 W₂ 面；W₁ 面确实**不在**那两条 no-go 的射程内。

---

## 1. Soundness of the `W₁ ⇒ eigenvalue` step  — `[T]`, no diameter factor

**Claim.** If a Markov operator `P` on a metric space `(X,d)` satisfies `W₁(P(x,·),P(y,·)) ≤ q·d(x,y)` for all `x,y`, then every eigenvalue `λ` of `P` with a non-constant **Lipschitz** eigenfunction obeys `|λ| ≤ q`.

**Proof (Lipschitz-eigenfunction argument).** If `f` is 1-Lipschitz then `Pf` is `q`-Lipschitz: for `x,y`, `|Pf(x)−Pf(y)| ≤ W₁(P(x,·),P(y,·)) ≤ q d(x,y)` (Kantorovich–Rubinstein duality). Induct: `Lip(Pⁿf) ≤ qⁿ Lip(f)`. For an eigenfunction `Pf = λf` (complex `λ` allowed), `Lip(Pⁿf) = |λ|ⁿ Lip(f)`. With `0 < Lip(f) < ∞`, divide: `|λ|ⁿ ≤ qⁿ` ⇒ `|λ| ≤ q`. ∎

**Diameter / volume: NONE.** The chain uses only the pointwise pair inequality `|Pf(x)−Pf(y)| ≤ q d(x,y)` and `Lip`; **no `diam(X)`, no `vol`, no `L²` norm appears.** Hence if `q` is `(L,a)`-uniform, so is the eigenvalue bound. This is precisely where the direct Lipschitz route beats the standard Ollivier→**L²-Poincaré** route, which must compare `L²` to `W₁` and *thereby* imports a `diam²` factor (the `L²` variance is bounded by `diam² × Lip²`). The direct pair-Lipschitz inequality sidesteps that comparison. **Confirmed precisely: no diameter factor.**

**Caveats (each benign or flagged).**
1. **Lipschitz eigenfunctions.** The bound constrains eigenvalues carried by eigenfunctions with finite `Lip`. The transfer kernel on the compact gauge group is smooth (heat-bath × `e^{−a·magnetic}`), so its low-lying eigenfunctions are smooth, hence Lipschitz. The **constant** eigenfunction `φ₀=1` has `Lip=0` and is correctly excluded — consistent with `λ₀/λ₀=1 ≰ q<1`. So the bound applies to `j≥1`.
2. **Reversibility NOT required** for the Lipschitz step (works for non-reversible `P`, complex `λ`, via `|λ|ⁿ`). Self-adjointness/reflection positivity (the project's `PROOF-GAUSS-INVARIANT-PF`) is used only to have a real orthonormal eigenbasis and `E_j−E₀ = −(1/aₜ)log(λ_j/λ₀)`; PF gives `λ₀>0` and the transfer is self-adjoint in `L²(h₀² dU)`.
3. **Normalization.** `P = T(h₀·)/λ₀` is stochastic (`P·1=1`); `spec(P)={λ_j/λ₀}`. The Doob eigenfunction is `φ_j = ψ_j/h₀`. The bound is on `|λ_j/λ₀|`.
4. **Sign / ordering.** Physical gap `Δ = −(1/aₜ)log(|λ₁|/λ₀)` where `λ₁` is the second-largest in magnitude; RP/PF makes the relevant `λ₁>0`. Fine.

**Consequence (the crux, `[T]`).** The lemma gives `λ₁/λ₀ ≤ q`, i.e. `e^{−aₜΔ} ≤ q`. To *certify* `gap ≥ Δ` one needs the converse `q ≤ e^{−aₜΔ}`. Both hold together **iff `q = e^{−aₜΔ}` exactly** — i.e. the W₁ contraction must be **exactly sharp**. Then `Δ_cert = −(1/aₜ)log q ≤ Δ_true`, with equality only in the sharp case. **This sharpness requirement is the structural obstruction** (see §3d). The free field (§2) is the sharp case; the interacting SU(2) kernel is not.

---

## 2. Acid test — free field (must pass) — `[V-here]`

Free `U(1)`/Gaussian transfer of a mode of frequency `ω`: the Doob-transformed one-step map is the OU semigroup, and `W₁(P(x,·),P(y,·)) = e^{−aₜω}|x−y|` **exactly** (same variance, means separated by `e^{−aₜω}|x−y|`). Hence

> **`q^{free} = e^{−aₜω}`, `ω = √(m² + (2π/L)²)`.**  `[V-here]`, exact to 6 digits (`q_num ≡ e^{−aₜω}` in the table).

| L | m | aₜ | ω | q | (1−q)/aₜ |
|---|---|---|---|---|---|
| 8 | 0 | 0.1 | 0.7854 | 0.924465 | 0.755 |
| 64 | 0 | 0.1 | 0.0982 | 0.990231 | 0.0977 |
| 10⁶ | 0 | 0.1 | 6.3e−6 | 0.999999 | ≈ω→0 |
| 10⁶ | 0.3 | 0.1 | 0.3000 | 0.970446 | 0.296→m |

- **`m=0`:** `(1−q^{free})/aₜ → 2π/L → 0` as `L→∞`. The finite-volume box mass `2π/L` is **not** faked into a continuum mass — exactly right. **Acid test PASSED (non-vacuous): the tool manufactures no gap on the free massless field.**
- **`m>0`:** `→ m`. Recovers the true continuum mass.
- The `(1−q)/aₜ` entries carry an `O(aₜ)` correction (`(1−e^{−aₜω})/aₜ = ω − aₜω²/2+…`), so they approach `ω` only as `aₜ→0`; the **identity** `q^{free}=e^{−aₜω}` is exact at every `aₜ`.

The free field is the **sharp** case (`q^{free} = λ₁/λ₀`): it is the one kernel for which the W₁ bound closes exactly and certifies the full gap.

---

## 3. The SU(2) kill-test at β=2.5 — `[V-here]`

Setup. Wilson SU(2), `S=βΣ_p(1−½Re Tr U_p)`, physical coupling **β=2.5**. Single link with its `2(d−1)=6` staples at the (aligned, identity) background: `S_link(U)=exp((β/2)ReTr(U·S))`, `S=6I`, giving exponent `6β·u₀ = 15 u₀`, `u₀=½Re Tr U ∈[−1,1]`, Haar `∝ √(1−u₀²)du₀`. One staple flipped `I→−I` ⇒ `S=4I`, exponent `10 u₀`. Full-group diameter `d_G=π`, antipodes `±I`.

### 3a. Single-link contraction + influence row sum

| quantity | value |
|---|---|
| single-staple TV influence `a₁=TV(p₁₅,p₁₀)` | **0.1814** |
| single-staple W₁ influence (metric `|u₀−u₀'|`) | 0.0476 |
| all-6-flip range `TV(p₁₅,p₋₁₅)` | 1.0000 |
| **worst-case TV/Dobrushin row sum**, 6 staples `6a₁` | **1.088** |
| same, 18 link–staple incidences `18a₁` | 3.265 |
| local linear-response: `β·Var(u₀)` per staple | 0.0161 |
| **local (infinitesimal) W₁ row sum**, `6β·Var(u₀)` | **0.0963** |
| same, 18 incidence `18β·Var(u₀)` | 0.2889 |

**Reading.** At β=2.5 the links are strongly **frozen** (`⟨u₀⟩=0.902`, `std=0.080`), so the **local** response of one link to a neighbor's field is tiny: the infinitesimal W₁ row sum is `0.096 ≪ 1` — **local contraction holds; the coarse-Ricci curvature is positive here**, where the W₂/BE curvature is Savvidy-negative (`IDEA-wasserstein.md`). **This is the genuine, non-redundant observation: W₁ survives locally where W₂ dies.**
**But** the **worst-case TV/Dobrushin row sum is ≥ 1** (1.09 for the 6-staple count, 3.27 over 18 incidences). So a *global, all-neighbour* contraction — the only thing Ollivier's spectral bound can use — **fails at β=2.5**. Local positivity does not upgrade to a global `q<1`: far-pair / worst-case contraction is uncontrolled. There is **no small parameter** turning the worst-case row sum provably `<1`.

### 3b. Local vacuum Harnack constant of `log h₀`

Independent-link proxy `h₀(U)~∏_ℓ p_λ(u₀ℓ)`, `λ=6β=15`:

| quantity | value |
|---|---|
| bulk gradient `∂_{u₀} log p` at `u₀=0` | 15.00 `(=6β)` |
| gradient at the measure mean `u₀=0.902` | 10.17 `(O(β))` |
| **typical per-link `|Δ log h₀|`, 1σ link move** | **0.815 `(O(1))`** |
| range of `log p` on fixed window `u₀∈[−0.5,0.9]` | 20.3 `(O(β))` |

`log p` diverges only at the measure-zero poles `u₀=±1` (never sampled). **`L`-uniform? YES** — these are *local* per-link constants, independent of `L` (charter condition 5 is met). **The route does NOT die on `L`-uniformity.** It dies because the constant is **O(1)** (typical `Δlog h₀≈0.82`, gradient `O(β)`): **no small parameter at β=2.5** — the precise meaning of "the crossover is the wall."

### 3c. The minimal RP transfer with a genuine gap (2-link / single plaquette)

The proposal's "2-link / single-plaquette patch with the RP transfer kernel" is, on one spatial site, exactly **2D SU(2) YM** (after spatial gauge fixing the plaquette is `V_t V_{t+1}^{-1}`, `T(V,V')=exp((β/2)ReTr(VV'^{−1}))`), i.e. the SU(2) heat kernel. Its Peter–Weyl eigenvalues are **exact**: `λ_j/λ₀ = I_{2j+1}(β)/I₁(β)`, so

> `λ_{1/2}/λ₀ = I₂(β)/I₁(β)`: β=0.5 → 0.1237; β=1.0 → 0.2402; **β=2.5 → 0.5072** (`aₜΔ=0.679`); β=5.0 → 0.7193. `[V-here]`, no approximation.

Cross-check: `⟨u₀⟩_{von Mises, λ=β} = I₂(β)/I₁(β)` (recurrence `I₂ = I₀ − 2I₁/β`), and the script reproduces `0.5072` — confirms the heat-kernel normalization. **The single link itself (0+0D) has NO transfer spectrum** (heat-bath is rank-1 / full resampling); the 2D heat kernel is the minimal exactly-solvable RP transfer that actually carries a gap. This is the yardstick a W₁ contraction must be compared with: by §1 it must satisfy `q ≥ 0.507`, and certify `Δ` only if `q = 0.507` **exactly**.

### 3d. Sharpness deficit (why even a *proven* W₁ contraction cannot certify `Δ`)

`λ₁/λ₀ ≤ q` (Lipschitz lemma) **always**; `Δ_cert = −(1/aₜ)log q ≤ Δ_true`. The free field attains equality (§2); the interacting SU(2) kernel is generically **strict** (`q > λ₁/λ₀`). So a proven W₁ contraction yields a **strictly weaker** gap than the physical one. Certifying the physical `Δ` would require the Ollivier bound to be **exactly sharp** — which is itself knowing `Δ` (`[O]`). The 2D exact ratio `0.507` shows the target magnitude; no argument here gets a W₁ bound down *to* it.

---

## 4. Circularity — does the Doob transform smuggle in `h₀`? — `[O]`

**Yes.** `P = T(h₀·)/λ₀` is built **from the ground state `h₀`**, the very object we lack. The W₁ contraction `q` in §1–3 is a property of *this* `h₀`-dependent kernel. Two escape attempts fail:

- **Use `T` directly?** `spec(T)` gives `λ_j/λ₀` (ratio is basis-free), but `T` is not stochastic; its `W₁` contraction must be measured w.r.t. the metric adapted to the **physical** measure `h₀²dU` — still `h₀`.
- **Use a different reference metric (e.g. product Haar)?** Then `T`'s `W₁` contraction w.r.t. that metric does **not** bound the physical spectrum, because the Lipschitz argument in §1 requires the metric in which the eigenfunctions are Lipschitz *and* the reference measure is the invariant one. Mismatched metric ⇒ no spectral conclusion.

The local Harnack lemma (§3b) — localizing `log h₀` to one link — is exactly a **local regularity statement about `h₀`**. Proving it at β=2.5 without presupposing the gap requires controlling how a one-link change propagates into `h₀`, i.e. that correlations decay (massive) **and** the total (summed over the correlation volume `~1/Δ^d`) stays `O(1)` and `L`-uniform. **That propagation control is H1** (the same object). So the local-Harnack hypothesis is `[C]` **conditional on H1** — `[O]`, circular. (It is *not* independent new information; the numbers satisfy it in the independent-link proxy precisely because that proxy *omits* the correlations, which are the whole difficulty.)

---

## 5. Literature verification (charter §2 — no unverified propagation) — `[R]`

**(a) Nair 2026, 3D Yang–Mills.** Verified: **V.P. Nair, "Towards a Proof of Mass Gap in 3d Yang-Mills Theory," arXiv:2608.10133 [hep-th]** (v1 10 Aug 2026, v2 22 Aug 2026, 11 pp; CCNY). It is an **outline**, giving `λ ≥ m = e²c_A/(2π)` via an **operator inequality `T² − mT ≥ 0`** on the **Laplacian of the gauge-orbit space** `C=A/G*`, using the Karabali–Kim–Nair Hamiltonian (Schrödinger) representation; the gap is generated by the **WZW-action weight `e^{2c_A S(H)}`** in the inner product.
- **3D, not 4D: CORRECT.** **Orbit-space Laplacian spectral inequality: CORRECT.**
- **"Bakry–Émery": INCORRECT — and the opposite is true.** Nair explicitly *reviews and avoids* Singer–Ricci/Lichnerowicz and Bakry–Émery (obstructed by regularization) and proposes the Laplacian/WZW inequality *precisely to sidestep* the curvature route. So "same as the rejected orbit-curvature route" is **wrong on this citation**: the BE orbit-curvature route is Mondal, *JHEP* 12 (2023) 191 (and the working-paper bound below), **not** Nair. Nair's route is itself an attempt to avoid the Savvidy/BE obstruction — directly relevant, but it is an **outline** (regularization of Green functions, positivity of the renormalized potential, and the saturating state are all open) and **3D only**. Do **not** cite Nair as a BE result.
- Aside: **S. M. N. Sakib, "Towards the Yang-Mills Mass Gap via Gauge-Orbit Curvature: A Strong-Coupling Bound and Data-Driven Lattice Evidence," Cambridge Open Engage working paper, 21 Nov 2026** (single author; **not peer-reviewed**). It *does* use Bakry–Émery, improving the SZZ bound by `(N/2)Q_d`, widening the strong-coupling window from `β<1/48` to `β<195/7168` in 4D (factor 1.31) — still `O(10¹×)` below the physical `β≈0.95`, i.e. the same 40-order wall. A widened **strong-coupling** bound, **not** a gap. `[R]`

**(b) Bringmann–Cao–Hairer–Zhao, 2D nonabelian YM–Higgs.** Verified: **arXiv:2611.03640** [math.AP] (2 Nov 2026, ~350 pp), **"Construction of the 2D non-Abelian Yang–Mills–Higgs measure via stochastic quantisation"** (also INSPIRE 3222658). Content: **first construction** of the 2D nonabelian **YM–Higgs** measure via Parisi–Wu stochastic quantisation (covariant paracontrolled calculus); the gauge-covariant approximation converges to the Langevin dynamics (modulo gauge), which is **topologically irreducible and strong Feller ⇒ unique invariant measure**; non-perturbative, any fixed coupling.
- **G1 construction (globally defined process + unique invariant measure): CORRECT.**
- **2D, not 4D: CORRECT.** It gives **no 4D gap** and does **not** identify stochastic-time mixing with the physical Hamiltonian gap (it is 2D, and with a Higgs field). Matches `TRANSPORT-RIGOROUS-NUMERICS.md` §2.1 (G1, 2D). `[R]`

**(c) Harold D. Kirk III, 2026, truncated SU(2) projective contraction.** Verified: **"Quantitative Mass Gap Bounds from Projective Contraction for Truncated Lattice SU(2) Yang–Mills Theory," Zenodo record 18773860, 25 Feb 2026 (v1.0)** (ORCID 0009-0000-6014-8451); companion "Quantitative PF₂ Stability and Explicit Zero-Free Neighborhoods…", same date.
- **Truncated theory only: CORRECT** ("truncated lattice SU(2)"; minimally truncated). Uses Stroock–Zegarlinski LSI + Dobrushin–Shlosman analyticity + PF₂ cone stability; claims `m_lat(β) ≥ c β^{−p}` along the Gaussian trajectory, with extension to the full theory asserted via a "Hurwitz approximation argument."
- **High-rank extension left open: CORRECT** ("claimed pathways to SU(N)"; framework asserted group-agnostic). **Flag:** independent Zenodo preprint, **not** arXiv, **not** peer-reviewed; the full-theory extension is asserted, not rigorously closed. Do not propagate as a proven lattice gap. `[R]`

**Also noted, NOT to be propagated** (charter §2): self-published `hexagonmath.org` "explicit lattice YM mass gap for every compact simple group" (arXiv:2610.00341v5) — not peer-reviewed, unverified; excluded.

---

## 6. Exact missing step `[O]`

To promote this route from "sound engine on a strictly-weaker object" to a Clay handle would require **one lemma**, here in the W₁ dialect:

> **[O-W1]  `(L,a)`-uniform, all-pairs, *globally worst-case* Ollivier W₁ contraction of the physical `h₀`-Doob transfer kernel at rate exactly `q = e^{−aₜΔ}`, with `Δ>0` the true gap — at the crossover `β≈2.5`.**

Three simultaneous obstructions, each `[V-here]`-quantified:
1. **`h₀`-circularity** (§4): `P` is defined by the unknown ground state.
2. **Non-sharpness** (§3d): the Lipschitz lemma only gives `q ≥ e^{−aₜΔ}`; certifying `Δ` needs `q = e^{−aₜΔ}` *exactly*, generically false for interacting SU(2) — so one certifies `Δ_cert < Δ`.
3. **Global vs local** (§3a): only the **local** (infinitesimal) W₁ row sum is `<1` (`0.096`); the **worst-case** row sum is `≥1` (`1.09–3.27`). A global uniform contraction is unproven and is exactly the non-perturbative crossover control.

**[O-W1] ≡ H1** in the W₁/Markov-coupling dialect. The value of this round is the precise *localization* of the wall: the obstruction is **not** the W₂/BE/Savvidy curvature (W₁ locally escapes it, `[V-here]` 0.096), but **(i) the `h₀`-Doob circularity, (ii) the exact-sharpness requirement, and (iii) the local→global upgrade at the crossover.** Points (i)+(iii) are the same crossover control as every prior route; (ii) is a new, sharper statement that even a *proven* W₁ contraction certifies only a strict sub-gap.

---

## 7. Scorecard (charter's 5 hard conditions from `BRIEF-FOR-MODEL.md` §2)

| # | condition | result |
|---|---|---|
| 1 | non-perturbative | **PASS** — pathwise W₁, no loop expansion (same as D2) |
| 2 | sign-controlling lower bound `m_eff²(k→0)≥λ` | **FAIL@2.5** — certifies `Δ_cert ≤ Δ` only, and needs sharpness; local row sum positive but global ≥1 |
| 3 | non-circular (no presupposed gap) | **FAIL** — Doob needs `h₀`; local Harnack = H1 (§4) |
| 4 | non-SUSY + **fails on free massless** | **PASS (cleanest)** — `q^{free}→1` as `L→∞`, no faked mass (§2) |
| 5 | `(L,a)`-uniform bound | **PARTIAL/FAIL** — Harnack `L`-uniform (`[V-here]`), but global worst-case W₁ row sum `≥1`, not `L`-uniform |

**Sharpest honest line.** The engine is right and the object is genuinely new (W₁ is *strictly weaker* than the refuted W₂/CD/BE route, so von Renesse–Sturm does **not** apply; the β=2.5 local coarse-Ricci is **positive**, `[V-here]` 0.096, where W₂ is Savvidy-negative). But it is not a key: the Doob transform needs `h₀` (circular), the Lipschitz bound `q ≥ e^{−aₜΔ}` certifies only a **strict sub-gap** unless the contraction is *exactly sharp*, and the *global* worst-case contraction fails at β=2.5 (row sum `1.09–3.27 ≥ 1`). It relocates the wall precisely — from *curvature* (W₂/BE, dead at Savvidy) to *`h₀`-circularity + exact-sharpness + local→global upgrade* — and that relocation is a real, verifiable advance in coordinates. **Clay(G2) still 0%; no QED.**

*Refs.* Internal: `IDEA-wasserstein.md` (W₂/CD refuted at Savvidy), `D2-COUPLING-DOEBLIN.md` (W₂ coupling ≡ `G″` ≡ H1, sharp von Renesse–Sturm), `PROOF-GAUSS-INVARIANT-PF`, `TRANSPORT-RIGOROUS-NUMERICS.md` (G1/G2 split), `BRIEF-FOR-MODEL.md` §2. External `[R]`: V.P. Nair, arXiv:2608.10133 (2026, 3d YM, orbit-space Laplacian, **avoids** BE, outline); S.M.N. Sakib, Cambridge Open Engage working paper 21 Nov 2026 (BE strong-coupling widening `β<195/7168`, **not** a gap); Bringmann–Cao–Hairer–Zhao, arXiv:2611.03640 (2D nonabelian YM–Higgs, G1 + unique invariant measure); H.D. Kirk III, Zenodo 18773860 (2026, truncated SU(2), not peer-reviewed). Mathematical: Ollivier, *J. Funct. Anal.* **257** (2009) 689; von Renesse–Sturm, *CPAM* **58** (2005) 923 (W₂ sharp — does not transfer to W₁); Joulin–Ollivier, *J. Funct. Anal.* **262** (2012) 2050 (W₁ vs spectral, diameter-free direction).
