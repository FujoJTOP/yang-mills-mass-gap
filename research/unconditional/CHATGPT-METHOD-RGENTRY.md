# CHATGPT-METHOD-RGENTRY — Audit of the RG-ENTRY / Kotecký–Preiss bridge proposal

> 2026-11-23 · Assessment, not a proof. This note audits the "nonperturbative RG bridge into a strong-coupling KP basin" corridor (RG-ENTRY 3.x) and the Faizal–Shabir 2026 realization of it.
>
> Labels: **[T] 真定理** (unconditional within its stated domain) · **[C] 条件定理** (conditional theorem) · **[NG] no-go** (for a named class of approaches, not for the problem) · **[O] 缺失步 / open**.
> FS audit labels: **[V]** verified by me from the full PDF (arXiv:2606.19362v1, text extracted with `pdftotext`); **[R]** as-reported / not independently verified.
>
> Iron rules followed: no "QED" on an unproven step; every claim typed; the free-massless acid test is applied explicitly.

---

## Verdict up front

**Same wall, different corridor.** The RG-ENTRY framework is a legitimate *conditional theorem* and a cleaner way to organize the problem than the project's positive-geometry routes, but its single open hypothesis — O2-Bridge — is the constructive form of confinement, not a weaker problem. Fixed `J` controls only the number of RG steps, not the per-step estimates through the `O(1)`-coupling crossover. I found **no independent handle**: no monotone/Lyapunov functional, no compactness argument, and no contraction over the `J` fixed steps that avoids the per-scale estimate. The known candidates fail for structurally independent reasons (Schur order, local→global, edge→support, indefinite Savvidy Hessian, free zero mode).

**The single exact open obstacle** is:

> `sup_{a,L} 𝔄_κ(R_b^J S_{a,L}) < ε_KP` at `ℓ_* = 1/(cΛ)`, `J = ⌈log_b(M/c)⌉` —
> equivalently, a **common contracting/KP polydisc** for the UV (AF) and IR (strong-coupling) RG trajectories, proved without assuming the gap.

In the Faizal–Shabir paper this is hidden in Theorem 10.8's premise ("at a common intermediate scale … both trajectories lie in the contracting polydisc") and in Lemma F.10's explicit requirement that the two regulators be **based on the same local bare action at the same bare coupling**. The paper's Eq (10.72) also contains a genuine harmonic-divergence error: it asserts `Σ|ε_k| < ∞` while `ε_k` contains `O(g_k²) ∼ 1/k`.

**FS 2026 is real, published, and not a proof.** arXiv:2606.19362 [hep-lat], Fortschr. Phys. **74** (2026) e70097, DOI 10.1002/prop.70097. It is a well-organized conditional construction; the crossover is assumed, not crossed. The Zenodo record 22840696 cited in the prior internal note is **not** an audit of Faizal–Shabir; it is a different G-field proof-architecture paper (V. Liashkevych), and its public metadata/passport do not mention Faizal or Shabir. The FS-specific criticisms attributed to it there are **[R]** and should be re-sourced or dropped. This is a correction to `RG-ENTRY-VERIFICATION.md` §4 under CLAUDE.md §2.

---

## §1 The object on the table

The RG-ENTRY 3.x conditional theorem has the following named inputs:

| Input | Statement | Status |
|---|---|---|
| **H-RG** | exact, gauge-invariant **spatial** block-spin `R_b^sp` preserving time reflection positivity and observable correlators | **[O]** not constructed in general |
| **H-coeff** | canonical gauge-invariant spin-network polymer activity `𝔄_{κ,τ}` (rep-dim weight `τ`, time-diameter `ζ`) | **[O]** definitional/constructive input |
| **O2-Bridge** | `∃ ℓ_*=1/(cΛ)`, `κ,τ,ε_KP`: `sup_{a,L} 𝔄_κ(R_b^J S_{a,L}) < ε_KP`, `J=⌈log_b(M/c)⌉` finite | **[O]** the crux |
| **Lemma-KP** | activity `< ε` (Kotecký–Preiss) ⇒ uniform exponential clustering ⇒ (per local operator + H-RS(vac)) `spec(H) ⊆ [m_*,∞)` | **[T]** given the activity bound; standard KP/OS machinery |
| **H-uniform** | rates/limits uniform in `a↓0`, `L_phys↑∞` | **[O]** |
| **H-RS(vac)** | vacuum-representation cyclicity of local gauge-invariant operators | **[T] / legitimate scope** for the standard Clay formulation |
| **H-RS(all)** | gap in all electric/magnetic/topological/superselection sectors | **[O]**, stronger, **not required** by Clay as standardly stated |

Two of these are themselves open constructive problems: **H-RG** (an exact gauge-invariant block-spin with all the stated properties) and **H-coeff** (a canonical activity with the required positivity). The FS paper's "admissible class" is an assumption of H-RG, not a construction; its horizon projector/Landau representative is a gauge-fixed device whose exact gauge invariance and reflection positivity are part of the hypothesis.

The correct logical reading is:

- **Lemma-KP ⇒ gap** (with H-RS(vac), H-uniform): a **[C] 条件定理**, not a proof.
- **O2-Bridge ⇒ the hypothesis of Lemma-KP**: so O2-Bridge is a **sufficient** condition for the mass gap.
- **Gap ⇒ O2-Bridge** is **not known**. O2-Bridge is strictly stronger than the gap as a mathematical statement: a gapped but strongly interacting theory need not have a KP-small effective action at the gap scale. So the framework is not a literal equivalence; it is a sufficient reformulation.
- **O2-Bridge ⇒ "no nontrivial 4D CFT in the IR"**: if the effective action at a finite physical scale is KP-small, the short-range expansion gives exponential clustering and no massless CFT. This is the constructive negation of the Clay alternative (iii).
- **"No nontrivial 4D CFT" ⇒ O2-Bridge** is also not known, but proving either requires the same nonperturbative control of the RG orbit. The project's slogan "transport ≡ confinement" is a **difficulty** statement, and it is the correct one.

The free-massless acid test is passed **in principle**: the free massless field has a zero mode, its effective action at any finite scale is not KP-small (the massless mode is relevant; the activity is `O(1)`), so the hypothesis fails and the tool is not vacuous. The free massive field has a gap and is KP-small at scales `≫ 1/m`; the tool distinguishes the two. The problem is that the same discriminator that is correct for the free field is exactly what is unproven for YM.

**Scope note.** Lemma-KP gives the spectral gap `spec(H) ⊆ {0} ∪ [m_*,∞)`; it does **not** by itself give the KL residue atom `f₁²>0` or its isolation `I2`. The latter is the project's separate strong-coupling transport target (`TRANSPORT-MACHINE` Stage 3), a related but distinct open problem. The Clay mass gap is the spectral statement, so the RG-ENTRY route is aimed at the right object; do not read it as automatically giving the residue formulation.

---

## §2 Task 1 — Is O2-Bridge equivalent to confinement, or genuinely more tractable?

### 2.1 The logical skeleton

Let `C` denote the brief's statement "the RG orbit reaches the trivial/short-range IR, not a nontrivial 4D CFT," understood constructively as: there is a finite scale at which the effective action is in the trivial/KP basin. Then:

```
O2-Bridge  ⇒  C          [C]  (KP-smallness is the constructive form of trivial IR)
C          ⇒  gap        [C]  (Lemma-KP + H-RS(vac) + H-uniform)
gap        ⇒  O2-Bridge  [O]  NOT KNOWN
O2-Bridge  ⇒  gap        [C]  (Lemma-KP)
```

So O2-Bridge is **not logically equivalent** to confinement: it is a sufficient, constructive strengthening. But the only known routes to O2-Bridge are:

1. prove the gap and then run a cluster expansion at large distances (circular); or
2. prove the strong-coupling basin is reached across the crossover (confinement).

Neither is a weaker problem.

### 2.2 Why fixed `J` does not simplify

`J = ⌈log_b(M/c)⌉` is finite and independent of `a,L`, so the *number* of RG steps is controlled. It does **not** control the *norm* of the full effective action at each step. The RG map is:

- a contraction in the polymer norm only inside the weak-coupling polydisc `|g|≤ε, K≤ρ` (FS Thm 10.4);
- a convergent KP/character expansion only inside the strong-coupling domain `β<β_0` (FS Thm 3.1 / Osterwalder–Seiler);
- **not contractive at the crossover**. The linearized map has a marginal coupling direction and relevant directions; the Savvidy Hessian is indefinite (`−D²+2ad F_B`, negative transverse magnetic modes), so the BBS large-field trapping argument fails. The number of steps is finite, but the per-step norm can grow, and the constants in any finite-`J` Gronwall estimate depend on the `O(1)` coupling at the crossover. This is exactly the point already recorded in `RG-ENTRY-VERIFICATION.md` §2 and `BRIEF-FOR-MODEL.md` §4.

**Label: [NG] for "fixed `J` is a simplification" as a strategy.** It is not a theorem that no finite-step argument can ever work; it is a no-go for the known estimates.

### 2.3 Hunt for an independent handle

I considered the natural candidates. None exists in the current framework.

| Candidate handle | Why it does not close O2-Bridge | Status |
|---|---|---|
| **Monotone/Lyapunov RG functional** (c/a-theorem analogue) | 4D a-theorem is a statement about CFT fixed points, not a Wilsonian Lyapunov function along a massive flow; it does not force an atom. In the project's Schur-order analysis, monotone/total objects are Schur-concave, while atom-forcing is Schur-convex; they lie on opposite sides. | **[NG]** for known species |
| **Relative entropy / entanglement monotonicity** | Monotone under RG, but blind to the support edge; it bounds susceptibility/`H1`, not the atom `H*`. The capacity/Poincaré route is explicitly necessary-not-sufficient (`C-BBS-SUFFICIENCY.md`). | **[NG]** |
| **Exact Ward identity `G″·Var(a)=1`** | Exact and sign-controlling for `H1`, but `H1` is the inverse susceptibility, not the gap (`H1-IS-THE-GAP` correction). It does not see the second-level isolation `I2`. | **[NG]** for the atom |
| **Prékopa `W_geo ≥ 0`** | A genuine floor (`PROOF-JACOBIAN-SAVES`), but not a monotone along RG and does not force the atom; it is compatible with a massless CFT. | **[NG]** |
| **Capacity / Dirichlet / Cheeger** | Gives `λ_1 ≥ Φ²/4`, i.e. `χ ≤ 4/Φ²`: a uniform `H1`, not `H*`. The capacity dictionary dies on the full gauge-quotiented measure at the Savvidy cone. | **[NG]** |
| **Compactness / tightness over `J` steps** | Tightness gives subsequential limits of measures, not that the limit is KP-small. KP-smallness is an **open** condition; the crossover can approach its boundary, and the limit could be a nontrivial CFT. | **[NG]** |
| **Contraction over `J` fixed steps** | The Gaussian polydisc contraction and the strong-coupling KP basin are disjoint in the known rigorous regimes; no map is known that sends one into the other in finite steps. The linearization at the crossover is not contractive. | **[NG]** |
| **BBS trapping region** | The large-field half needs a positivity (`+12λφ²` in φ⁴) that YM does not have. The Savvidy Hessian is indefinite; WSTC substitutes only on the value tail, not the Hessian sign. | **[NG]** (`C-BBS-MARGINALLOG` §3) |
| **Large-N / master-field / Chatterjee** | Rigorous at fixed small `β` (strong coupling); `N=∞ → 3` is an asymptotic series with no uniform control; does not cross the crossover. | **[NG]** |
| **Topological/degree argument** | The space of effective actions is infinite-dimensional and not compact in the polymer norm; a nontrivial CFT fixed point is a possible limit. No topological necessity forces the KP basin. | **[NG]** |
| **Numerical/lattice evidence** | Supports confinement, but is not a proof; finite-volume gaps do not by themselves give an `(L,a)`-uniform lower bound, and the continuum extrapolation is the open step. | not a handle |

**Honest conclusion for Task 1:** no independent handle exists in the current framework. This is a **no-go for the listed species** (monotone, compactness, contraction, capacity, Ward, Prékopa), not a mathematical theorem that no future tool can exist. The project's independent no-gos — Schur order, local→global, edge→support, and `C-BBS-SUFFICIENCY`'s necessary-not-sufficient — all point to the same place.

### 2.4 Verdict on Task 1

**Not a literal equivalence, but the same wall in difficulty.** O2-Bridge is a sufficient reformulation of "the RG orbit reaches the trivial IR," and it is strictly stronger than the gap (no known implication from the gap/confinement to O2-Bridge). The only known way to prove it is to prove the crossover/confinement. It is not more tractable than the original gap; it is a cleaner statement of it.

---

## §3 Task 2 — The KP domain: what is `β_KP`, and can `ℓ_*` be KP-small without the gap?

### 3.1 `1/48` is not `β_KP` (corrected reading, not re-litigated)

The prior correction stands (`RG-ENTRY-VERIFICATION.md` §0–§1):

- SZZ's `1/48` is the **Bakry–Émery / LSI** sufficient region for their Langevin semigroup, **not** the convergence radius of a Kotecký–Preiss polymer expansion.
- The parameter conversion is `β_W = N² β_SZZ`; for SU(2), `β_SZZ < 1/48 ⇒ β_W < 1/12`.
- The project's P3 theorem uses `β < min(1/48, c/N)` for its **graded-strip KP lemma**; that is a method-specific sufficient bound, not the true KP radius.
- `β ≈ 2.2` is a **crossover**, not a `T=0` bulk transition. This note does not revive the wrong version.

### 3.2 The actual KP convergence domain

There is **no universal `β_KP`**. It depends on the polymer norm, the activity, and the expansion. What is known:

- For the standard character/polymer expansion of the Wilson action, a sufficient domain is `β_W < β_0` with `β_0` explicit but small. FS's §3 gives `β_0 = 1/(2e C_0 c_1)` (their Eq (3.12)); Osterwalder–Seiler give a group-dependent `β < c_d/N²` for SU(N). The true radius is unknown.
- The actual KP domain is in the **deep strong-coupling region** (`β` small, `g` large), far below the physical crossover `β_W ≈ 2.2`.
- The KP domain is **not** the same as the BE/LSI domain; neither implies the other in general without additional hypotheses. A measure can satisfy an LSI/BE inequality without being in a KP basin, and a KP expansion needs its own activity bound.

**Label: [T] that `1/48` is BE/LSI, not KP; [O] the exact `β_KP`.**

### 3.3 Can the effective action at `ℓ_*=1/(cΛ)` be proven KP-small without the gap?

This is the crux. The answer is **no**, for a structural reason:

- If `c = O(1)`, then `ℓ_* ∼ 1/Λ` is the **gap/crossover scale**. The effective coupling there is `O(1)`; the theory is neither in the weak-coupling Gaussian basin (`g ≪ 1`) nor in the strong-coupling KP basin (`β ≪ β_0`, `g ≫ 1`). The KP-smallness of the full gauge-invariant effective action at this scale is exactly the unproven crossover estimate.
- If `c ≫ 1`, then `ℓ_*` is in the **UV**. The effective action is near Gaussian by asymptotic freedom, but this is the UV basin; it says nothing about the IR gap. Proving KP-smallness there would be a UV small-field constructive-RG statement (Balaban-type), not the mass gap.
- If `c ≪ 1`, then `ℓ_*` is in the **deep IR** (`ℓ_* ≫ 1/Λ`). A gapped theory should have a trivial effective action there, but proving that is precisely the gap: the cluster expansion at large distances is what one is trying to establish. The implication "gap ⇒ trivial effective action at `ℓ_*`" is not a known theorem; it is the same cluster-expansion problem.

FS's Lemma 10.5 does **not** solve this. It enters at **small `g_bare`**, i.e. the **weak-coupling/UV** side, after `k0` coarse-graining steps; its proof claims `K0 ≤ A_b g_bare²` and `g0 ≈ g_bare`. The phrase "by construction in the strong-coupling sector" is a **parameter-identity error**: small `g_bare` is weak coupling, not strong coupling. Section 10 is titled "Weak-Coupling Extension and Asymptotic Freedom"; its polydisc is the Gaussian one. The strong-coupling base of Section 3 is at `β < β_0` (large `g`, fixed `a`), which is a **different, disjoint** region. Lemma 10.5 therefore does not put the effective action at `ℓ_*` into the strong-coupling KP basin; it puts a UV trajectory into the Gaussian polydisc.

**Is this confinement restated?** In the sense of difficulty, yes. The statement "the effective action at `ℓ_*` is KP-small" is a sufficient condition for the gap; the only known ways to prove it are (i) prove the gap and run the cluster expansion, or (ii) prove the strong-coupling basin is reached across the crossover. Both are the same nonperturbative dynamics. It is not a definitional tautology, but it is not a weaker problem.

**Label: [O] the KP-basin statement at `ℓ_*`; [NG] for the claim that it is a shortcut around confinement.**

---

## §4 Task 3 — `H-RS(all)`: which sectors, what gap, and does Clay require them?

### 4.1 The sectors and their gaps

| Sector | Operator / label | Gap in the confining phase | Threat to the vacuum gap? |
|---|---|---|---|
| **Vacuum** | local gauge-invariant operators (`tr F²`, `T_{μν}`, …) | `Δ_vac = inf spec(H|_{H_vac⊖CΩ}) > 0` | This is the Clay gap |
| **Electric 1-form** | Wilson loop/line in rep `R` | static charge energy grows as `σR`; isolated charge has infinite energy; no normalizable one-particle state | No finite lower gap; sector is not part of the vacuum representation |
| **Magnetic 't Hooft** | 't Hooft loop / twisted boundary condition | 't Hooft loop has perimeter law (magnetic charges not confined); the twisted ground state is a flux/torelon at `E ∼ σL`, above the vacuum and growing with `L` | No lower gap; twisted sectors fix a quantum number, not a Hessian sign (`FORGE5-TWISTED`) |
| **Topological / θ** | instanton number, θ vacuum | θ is a Bloch label; RP pins `θ ∈ {0,π}`; topological charge is blind to spacetime translations | No new lower gap; θ is an orthogonal axis (`G1-UNIQUENESS-VACUUM` §2(iii)) |
| **Center `Z_N`** | center-twisted boundary conditions | at `T=0` the center is not a bulk order parameter; on finite `T⁴` the twisted sectors are boundary conditions; `Z_N` "vacua" are finite-`T` artifacts | No lower gap; uniformity across center sectors is a boundary-condition statement (`G1-OS4-GIVEN-GAP` §4) |

**Quantification.** In the vacuum sector the gap is `Δ_vac`. In the electric/magnetic sectors the relevant energy is the string energy `σR` (or `σL` for a torelon), which is positive and grows with separation/volume; there is no finite-volume-independent lower gap. The topological sectors share the bulk Hamiltonian and do not lower the gap. The center-twisted sectors have a flux/torelon gap `∼σL → ∞`; they are not a threat. Therefore the minimal gap over all sectors is the vacuum gap, **assuming confinement**; but proving the uniform gap across all sectors is a stronger, independent statement.

**Free Maxwell acid test.** On `T³`, free Maxwell has a massless photon in the **vacuum** sector (`m_phys = 0`), while a twist lifts the zero mode (`m_tw > 0`). Thus a sector-wise gap (even if true in the twisted sector) does **not** imply the vacuum gap. This is the exact counterexample to any attempt to prove the mass gap by proving gaps in flux sectors (`FORGE5-TWISTED` §3–§4). `H-RS(all)` is not a substitute for `H-RS(vac)`.

### 4.2 Does Clay require `H-RS(all)`?

The standard Clay/Wightman formulation asks for a quantum field theory satisfying the Wightman (or OS) axioms with a unique vacuum and a mass gap. The Wightman axioms are formulated on a **single** Hilbert space carrying a unitary representation of the Poincaré group, with a cyclic vacuum; the GNS construction from the vacuum state gives the **vacuum representation**. Superselection sectors are additional structure, not part of the axioms. The mass gap is a statement about the spectrum of the Hamiltonian on the vacuum representation. Therefore:

- **`H-RS(vac)` is legitimate and sufficient** for the mass gap as standardly formulated. The spectral measure proof (Lemma-KP / Theorem 8.6) needs only the vacuum-representation Reeh–Schlieder property. Reeh–Schlieder is a structural theorem of the Wightman/OS framework **once the theory exists**; it is not an extra dynamical input. The dynamical input is the gap/clustering.
- **`H-RS(all)` is stronger and not required** by Clay. It would be required only if one defines the "physical Hilbert space" as a direct sum over all electric/magnetic/topological sectors and demands a gap in each. That is a different, stronger problem.
- **Restricting to the vacuum representation is honest if explicitly stated.** It is not honest to use the phrase "Reeh–Schlieder" to pretend to cover charged sectors; the sectors are inequivalent representations, and a gap in one does not imply a gap in another.

FS's Theorem 8.6 is explicit: it uses "the Reeh–Schlieder property in the OS framework" for the local algebra acting on the vacuum, i.e. `H-RS(vac)`. It does **not** claim `H-RS(all)`. That is the correct scope for the mass gap; the area-law statement (Theorem 6.7) is a separate expectation-value statement about Wilson loops and does not require `H-RS(all)` either. FS passes this honesty check.

**Label: `H-RS(vac)` = [T]/legitimate scope; `H-RS(all)` = [O], stronger, not required. The vacuum-representation restriction is honest for Clay; it is a scope choice, not a cheat.**

---

## §5 Task 4 — Audit of Faizal–Shabir 2026

### 5.1 Identity and status

- arXiv:**2606.19362** [hep-lat], v1, 9 June 2026; cross-list hep-th.
- Fortschr. Phys. **74** (2026) 4, **e70097**; DOI **10.1002/prop.70097** (OA).
- Authors: Mir Faizal, Arshid Shabir.
- Companion papers: Int. J. Geom. Meth. Mod. Phys. (2026) 2650111–2650114, attached at the end of the arXiv PDF.
- **Status: published but contested; Clay unsolved.** The 2-year scrutiny rule applies. The paper is a conditional construction, not an unconditional proof.

### 5.2 Theorem 10.4 — local, not global entry

**Exact content [V].** Theorem 10.4 (One-step contraction near the Gaussian fixed point): there exist `ε, ρ > 0`, `q ∈ (0,1)`, `c₁,c₂,c₃ > 0` such that for all `(g_k,K_k)` in the polydisc `D(ε,ρ)`, the renormalization map satisfies

```
K_{k+1} ≤ q K_k + c₁ g_k²,
g_{k+1} = g_k − β₀ g_k³ + R_k,   |R_k| ≤ c₂|g_k|⁵ + c₃|g_k|K_k + c₃K_k²,
β₀ = 11N/(48π²).
```

The theorem is a **local** contraction near the Gaussian fixed point. The paper itself states [V]:

> "Theorem (10.4) does not by itself guarantee that an arbitrary admissible trajectory will enter the polydisc D(ε,ρ)."

So the local-vs-global criticism is **verified at the equation level**: Theorem 10.4 is not an entry theorem.

### 5.3 Lemma 10.5 — the parameter identity

**Exact content [V].** Lemma 10.5 (Entry in finitely many steps by one-parameter tuning): there exist `ḡ > 0` and `M = M(b)` such that if `0 < g_bare ≤ ḡ` and one chooses `k0 ≥ M` with the effective IR-scale condition `μ₀ b^{-k₀} ≤ a^{-1}`, then the `k0`-scale data `(g₀,K₀) ∈ D(ε/2,ρ/2)`, and after at most `m(b)` further steps the trajectory enters `D(ε,ρ)`. The proof claims [V]:

> "By construction in the strong-coupling sector, `(g₀,K₀)` is finite and reflection-positive; moreover, as `g_bare → 0` one has `g₀ → 0` and `K₀ → 0` …"

and

> "`K₀ = K_{k₀} ≤ A_b (g_bare)²`, `|g₀ − g_bare| ≤ B_b (g_bare)³`."

**The flag [V].** `g_bare → 0` is the **weak-coupling** limit; in Wilson normalization `β_W = 2N/g²`, so `g_bare → 0` is `β_W → ∞`. The "strong-coupling sector" phrase is therefore a **parameter-identity error**. Section 10 is explicitly the "Weak-Coupling Extension and Asymptotic Freedom"; its polydisc is the Gaussian/UV polydisc. The strong-coupling base of Section 3 is `β < β₀` (large `g`, fixed `a`), a **different and disjoint** region. Section 3 itself calls the small-`β` regime "small bare coupling"; in Wilson normalization small `β` is **strong** coupling, so the paper's "bare coupling" is used inconsistently across §3 and §10. Lemma 10.5 is an entry lemma at weak coupling/UV, not a bridge from strong coupling.

The proof also assumes the KP/BKAR expansion at the mesoscopic scale converges with activity controlled by powers of the covariance `C(k₀)` and by `g_bare`; that is the UV small-field regime. It is not a proof that the effective action at the IR/crossover scale is in the strong-coupling KP basin.

### 5.4 Theorem 10.8 — AF–IR identification and the same-bare-action premise

**Exact content [V].** Theorem 10.8 claims `S_n^{AF} = S_n^{IR}` for all `n`, hence a unique universal continuum limit with a gap, area law, and asymptotic freedom. The proof says [V]:

> "At a common intermediate scale `k` large enough that both trajectories lie in the contracting polydisc, we write the difference of their one-step kernels as `T_k` …"

and later [V]:

> "The initial difference at `k = 0` can be absorbed into the same summable error by adjusting the starting scale or, equivalently, by observing that admissibility provides a finite Lipschitz distance between the initial one-slice data of the two constructions."

**The verified gap [V].** The premise "both trajectories lie in the contracting polydisc" is **assumed**, not proved. More strongly, it is incompatible with the paper's own definitions: the strong-coupling base has `β < β₀` (large `g`), while the weak-coupling polydisc has `|g| ≤ ε` (small `g`); for `ε² < 2N/β₀` the two regions are disjoint. To put the strong-coupling trajectory into the weak-coupling polydisc one must run it to the UV; that is exactly the unproven transport. The "common intermediate scale" is the bridge, and the bridge is the wall.

**Lemma F.10 [V].** The UV-alignment lemma is stated as:

> "Lemma F.10 (UV alignment). Let `Ω₀, Ω₁` be admissible regulators **based on the same local bare action** (for instance the Wilson plaquette action at the same bare coupling), differing only in their admissible slice projector and blocking/finite-range-decomposition choices. Then … `d_adm(Ω₀,Ω₁)` can be made arbitrarily small …"

Theorem F.6 gives a bound proportional to `d_adm(Ω₀,Ω₁)`; to conclude equality one needs `d_adm → 0`, not merely `d_adm < ∞`. Lemma F.10 supplies the `d_adm → 0` only for **the same local bare action at the same bare coupling**. The AF trajectory (small `g_bare`, weak coupling) and the IR construction (strong coupling, `β < β₀`) are **not** based on the same bare action. The theorem's use of F.10 is therefore unjustified. This is exactly the "same-bare-action premise" criticism: **verified [V]**. If the paper intended the strong-coupling base to be a *blocked* version of the same weak-coupling bare action, then proving that the blocked effective action is in the strong-coupling KP basin is exactly O2-Bridge; the premise is still the missing step.

**Label: Theorem 10.8 = [C] conditional on the common-polydisc/same-bare-action premise; the premise is [O].**

### 5.5 Equation (10.72) — the harmonic divergence

**Exact content [V].** With `h_k := g_k^{-2}`, the paper writes

```
h_{k+1} − h_k = 2β₀ log b + ε_k,          (10.72)
ε_k = O(g_k²) + O(K_k/g_k) + O(K_k²/g_k²).
```

It then asserts:

> "The previously established decay `K_k → 0` together with `g_k → 0` implies `Σ_{k≥k₀}|ε_k| < ∞`."

**This is false as written [V].** `g_k ∼ (2β₀ log b · k)^{-1/2}`, so `g_k² ∼ 1/(2β₀ log b · k)` and

```
Σ_k g_k² ∼ (2β₀ log b)^{-1} Σ_k 1/k = ∞.
```

The `O(g_k²)` term is not summable; it is harmonic. The expansion of `g^{-2}` gives `ε_k = +2β₀ log b · g_k² + …`, with a nonzero positive coefficient, so there is no cancellation. Consequently `h_k = h_{k₀} + 2β₀ log b (k−k₀) + O(log k)`, not `+ O(1)`; the statement "the `O(1)` remainder has a finite `k → ∞` limit" is wrong.

**What survives.** The leading asymptotic `g_k ∼ (2β₀ log b · k)^{-1/2}` is unchanged; the error is in the subleading constant and in the summability bookkeeping. The true difference between two trajectories with different initial data is `O(k^{-3/2})` (the linearized flow gives `δg ∼ k^{-3/2}`), which is summable; the paper's stated bound `|g_j − g_j'| ≤ C(g_j² + g_j'²) ∼ C'/j` is too weak and not summable. So Eq (10.72) is a **repairable but real** error: the claimed `Σ|ε_k| < ∞` and the `O(1)` remainder are not valid as written, and the telescoping in Thm 10.8 needs a sharper estimate.

### 5.6 Additional verified issues

**(a) The `ε_k` name collision [V].** In Eq (10.72), `ε_k` is the running-coupling remainder `O(g_k²)+O(K_k/g_k)+O(K_k²/g_k²)`. In Eq (10.78), `ε_k` is the interlacing defect `‖E_k‖ ≤ ε_k` with `Σ ε_k < ∞`. They are **different sequences** with the same symbol in nearby sections. This is an expositional trap, not a mathematical error by itself, but it obscures which summability is used where.

**(b) Scale-direction inconsistency [V].** Section 10.1 and the interpolation paragraph in §10.2 use the effective IR scale `μ_k = μ₀ b^{-k}` (coarse-graining, IR direction), while Theorem 10.7 uses `μ_k = μ₀ b^k` (UV direction). The recursion `g_{k+1} = g_k − β₀ g_k³` decreases `g` with `k`, which is AF only if `k` is the UV direction (or if `μ` increases with `k`). In the IR direction, AF says `g` should **increase**. The paper mixes the two conventions. This is the "RG scale direction self-consistency" criticism: **verified [V]** as a textual inconsistency; the leading asymptotic statement can be read in the UV convention, but the "from the ultraviolet/from the infrared meet" language then does not match the displayed flow.

**(c) The companion Part 1 admission [V].** The attached Part 1 paper states explicitly:

> "In this Part 1, our analysis is confined to the strong-coupling domain `0 < β < β(N)` at fixed lattice spacing `a`, and the existence of a spectral gap is established only in this regime. The critical question of whether this nonzero gap survives as `β` increases (approaching the continuum limit) is deliberately left unanswered here. It will be addressed in Part II and Part III via a multi-scale renormalization group (RG) analysis."

This is a direct admission that the strong-coupling gap is **fixed-`a`** and the transport to the continuum is the open question. The main paper's Sections 5–10 claim to supply that transport; the audit above shows the bridge is assumed, not proved.

**(d) Section 6 "choose K such that a_K lies in the strong-coupling domain" [V].** The continuum area-law proof (Prop. 6.5) coarse-grains until `a_K` is in the strong-coupling domain, then invokes the strong-coupling area law at the coarse mesh. This is the same O2-Bridge in the area-law channel: it assumes the effective action at the coarse physical scale is in the strong-coupling basin. The phrase "choose `K` such that `a_K` lies in the strong-coupling domain" is an existence assumption, not a proved consequence of the UV flow.

### 5.7 The Zenodo record 22840696 — correction

The prior internal note (`RG-ENTRY-VERIFICATION.md` §4) attributes specific FS criticisms to "independent audit Zenodo 22840696." I fetched the record. **The attribution is not corroborated:**

- Zenodo record **22840696**, DOI **10.5281/zenodo.22840696**, title: **"II. Yang–Mills Mass Gap in G-Field Theory: Proof Architecture and Pathways, Constructive Partial Results, and the Remaining Model-Specific Frontier"**, author **V. Liashkevych** (metadata transliteration: Volodymyr Lyashkevuch), published **19 September 2026** (v2), CC BY 4.0, ~38.4 MB in four parts.
- The record is a **G-field proof program**, not an audit of Faizal–Shabir. Its public landing page and the English "Scientific Passport" PDF contain **zero occurrences of "Faizal" or "Shabir"** (I downloaded and searched the passport). It does not contain the specific criticisms (tightness, covariance, clustering, Wilson-carrier hierarchy, Part (3) normalization reversal) attributed to it in the internal note.
- The record's public description states that the 4D pure-Yang–Mills identification "has not yet been completed" and that no mass-gap result for standard Yang–Mills is claimed. Its "no-bypass" conclusion is that "finite-regulator positivity cannot automatically be promoted to continuum noncollapse" — a caution that is consistent with our verdict, but it is **not** an FS-specific audit.
- **Status of the attributed FS criticisms: [R], not [V].** They may come from a different source or be unsupported; they should not be cited as record 22840696 without a page reference. This is a citation correction under CLAUDE.md §2.

The FS audit above stands on its own **[V]** evidence; the Zenodo correction does not rescue or damage the paper, it only fixes the provenance of the prior internal note.

### 5.8 FS audit scorecard

| Audit point | Status | What is actually true |
|---|---|---|
| Thm 10.4 local-vs-global entry | **[V]** | local one-step contraction in `D(ε,ρ)`; paper admits it does not guarantee entry |
| Lemma 10.5 `g_bare` "strong-coupling" | **[V]** | small `g_bare` is **weak** coupling/UV; "strong-coupling sector" is a parameter-identity error |
| Thm 10.8 AF–IR initial alignment | **[V]** | common-polydisc premise assumed; Lemma F.10 requires same local bare action at same bare coupling; AF/IR are different bare actions |
| Eq (10.72) `Σg_k²` | **[V]** | `ε_k` contains `O(g_k²) ∼ 1/k`; `Σ|ε_k| < ∞` is false; `h_k` has `O(log k)`, not `O(1)` |
| `ε_k` name collision | **[V]** | Eq (10.72) running remainder vs Eq (10.78) interlacing defect |
| RG scale direction | **[V]** | §10.1 and the §10.2 interpolation use `μ₀ b^{-k}` (IR); Thm 10.7 uses `μ₀ b^k` (UV); `g` decreases with `k` is AF only in the UV reading |
| Strong-coupling base fixed-`a` | **[V]** | companion Part 1 explicitly leaves the `β`-increase/continuum transport unanswered |
| Section 6 `a_K` strong-coupling domain | **[V]** | same O2-Bridge in the area-law channel |
| Zenodo 22840696 as FS audit | **[R]/corrected** | record is a different G-field paper; no FS mention in landing page/passport |
| "Clay solved" | **NO** | published but contested; conditional construction; Clay unsolved |

---

## §6 Task 5 — Final verdict

### 6.1 Different corridor or same wall?

**Different proof architecture, same wall.** The RG-ENTRY corridor is not a reheated version of the project's positivity/curvature/spectral routes: it is a constructive-RG/KP route, and it has a clean conditional-theorem structure. It passes the free-massless acid test in principle, and it correctly isolates the crossover as the single bottleneck. But:

- Its open hypothesis (O2-Bridge) is the constructive form of "the RG orbit reaches the trivial IR," i.e. confinement in the brief's sense.
- Fixed `J` does not simplify the per-step estimate; there is no monotone, compactness, or contraction handle.
- The Faizal–Shabir realization assumes the bridge in Theorem 10.8 and in Lemma F.10, and contains a genuine Eq (10.72) summability error.
- The strong-coupling and weak-coupling basins are disjoint; the "common intermediate scale" is the missing bridge, not a proved overlap.

So the corridor is **legitimate as a conditional framework**, but it is **not more tractable**. It is the same wall with a better-organized proof skeleton.

### 6.2 The single exact open obstacle

> **Prove, uniformly in `a` and `L`, that the exact gauge-invariant effective action after `J` block-spin steps at the physical crossover scale `ℓ_* = 1/(cΛ)` is Kotecký–Preiss-small — equivalently, that the UV (AF) and IR (strong-coupling) RG trajectories share a common contracting/KP polydisc — without assuming the gap.**

In the RG-ENTRY notation: `sup_{a,L} 𝔄_κ(R_b^J S_{a,L}) < ε_KP`. In FS's paper: the missing premise of Theorem 10.8 (`both trajectories lie in the contracting polydisc`) plus the same-bare-action requirement of Lemma F.10, with the Eq (10.72) bookkeeping repaired. In the project's language: this is the machine that transports the strong-coupling gap/residue across the Savvidy crossover and takes the continuum limit — the machine the charter says does not yet exist.

### 6.3 What would change the verdict

- A proof of O2-Bridge (or a proof that it is false at `ℓ_*` on the lattice, e.g. the KP norm of the blocked effective action at `β ≈ 2.5` is `O(1)`), or
- a genuinely independent handle: a monotone/Lyapunov functional that forces the atom, a compactness argument that certifies KP-smallness, or a contraction over `J` steps that is uniform at the crossover, or
- a proof that the strong-coupling and UV trajectories are based on the same bare action with a summable defect (the FS Thm 10.8 premise), which would make the identification a theorem.

Absent one of these, the corridor remains conditional and the wall remains.

---

## §7 Scorecard

### 7.1 The five hard conditions (`BRIEF-FOR-MODEL.md` §2)

| # | Condition | RG-ENTRY status | Reason |
|---|---|---|---|
| 1 | Nonperturbative | **PASS (form)** | exact RG + polymer/KP, no perturbative expansion assumed |
| 2 | Sign-controlling lower bound | **CONDITIONAL** | Lemma-KP gives a gap lower bound **if** the activity bound holds; the activity bound is the missing step |
| 3 | Non-circular | **FAIL (in practice)** | O2-Bridge is the constructive form of confinement; the only known proofs use the gap/confinement |
| 4 | Non-SUSY + free-field acid test | **PASS** | free massless field is not KP-small (zero mode); free massive is; the tool discriminates |
| 5 | `(L,a)`-uniform lower bound | **CONDITIONAL/OPEN** | H-uniform is a named hypothesis; the uniform activity bound is exactly O2-Bridge |

### 7.2 Framework hypotheses

| Hypothesis | Status |
|---|---|
| H-RG (exact gauge-invariant spatial block-spin, RP-preserving) | **[O]** not constructed in general; FS's "admissible class" assumes it |
| H-coeff (canonical gauge-invariant spin-network activity) | **[O]** definitional/constructive input |
| O2-Bridge | **[O]** the single exact obstacle |
| Lemma-KP | **[T]** given the activity bound (standard KP/OS) |
| H-RS(vac) | **[T]/legitimate** for the standard Clay formulation |
| H-RS(all) | **[O]**, stronger, **not required** by Clay |
| H-uniform | **[O]** |

### 7.3 What is genuinely established

- **Lemma-KP** is a true theorem given the KP activity bound; the KP/OS machinery is standard.
- **FS Theorem 3.1** (strong-coupling character/polymer convergence at `β < β₀`, fixed `a`) is a true theorem in its domain, essentially Osterwalder–Seiler.
- **FS Theorem 8.6** (Tauberian gap from exponential clustering) is a true conditional theorem; it correctly uses `H-RS(vac)`.
- **FS Section 5 interlacing inequality** (`T_{k+1} = V_k T_k^b V_k − D_k + E_k`, `D_k ≥ 0`, `‖E_k‖ ≤ ε_k`) is a genuine conditional theorem: if `Σ ε_k < ∞`, the gap persists. The open part is the uniform summability of `E_k` across the crossover.
- **The H-RS split** is correct and important: `H-RS(vac)` is what the spectral proof needs; `H-RS(all)` is an independent, stronger statement.
- **The free-massless acid test** is passed in principle by the KP-smallness hypothesis.

---

## §8 Sharpest line and exact stuck step

**Sharpest line.** O2-Bridge is not a bridge; it is the far bank renamed. The RG-ENTRY theorem is a valid conditional theorem, but its one open hypothesis is the constructive form of confinement, and the Faizal–Shabir "identification" assumes the two banks are the same polydisc. The UV bank (`g` small) and the strong-coupling bank (`β` small, `g` large) are disjoint in the known rigorous regimes; the "common intermediate scale" is the unproven crossover. Fixed `J` controls the number of steps, not their norms; no monotone, compactness, or contraction handle is known; the free-massless field fails the hypothesis exactly where YM is unproven.

**Exact stuck step.**

```
sup_{a,L} 𝔄_κ( R_b^J S_{a,L} )  <  ε_KP      at  ℓ_* = 1/(cΛ),  J = ⌈log_b(M/c)⌉,
```

i.e. the (`a,L`)-uniform KP-smallness of the full gauge-invariant effective action at the crossover scale, equivalently a common contracting/KP polydisc for the UV (AF) and IR (strong-coupling) RG trajectories. In FS's paper: the common-polydisc premise of Theorem 10.8 plus the same-bare-action requirement of Lemma F.10 (and the Eq (10.72) `Σg_k²` bookkeeping). No independent handle exists in the current framework; this is the single exact missing step.

---

## §9 Refs

**Internal (absolute paths):**

- `D:\Dev\math\CLAUDE.md` — iron rules, current frontier.
- `D:\Dev\math\docs\BRIEF-FOR-MODEL.md` — the wall, 5 hard conditions, 11-species death map, "transport ≡ confinement."
- `D:\Dev\math\docs\PROOF-SKELETON.md` §0 — G1/G2 ledger, the single species (variance upper bound).
- `D:\Dev\math\research\unconditional\RG-ENTRY-VERIFICATION.md` — prior audit; corrections `1/48` vs KP, H-RS split, FS status.
- `D:\Dev\math\research\unconditional\TRANSPORT-MACHINE.md` — four-stage machine; Stage 3 = H1/I2; no weakening.
- `D:\Dev\math\research\unconditional\C-BBS-SUFFICIENCY.md` — capacity/Poincaré = H1, necessary-not-sufficient.
- `D:\Dev\math\research\unconditional\C-BBS-UVTRAP.md` — UV trap, thermal shell = H1.
- `D:\Dev\math\research\unconditional\C-BBS-ARRHENIUS.md` — thermal-shell Arrhenius, exact residual.
- `D:\Dev\math\research\unconditional\C-BBS-MARGINALLOG.md` — BBS large-field positivity breaks; Savvidy indefinite Hessian.
- `D:\Dev\math\research\unconditional\G1-OS4-GIVEN-GAP.md` — OS4/unique vacuum, center sectors.
- `D:\Dev\math\research\unconditional\G1-UNIQUENESS-VACUUM.md` — θ/topological sectors, OS4 vs gap.
- `D:\Dev\math\research\unconditional\FORGE5-TWISTED.md` — 't Hooft/twisted sectors; free Maxwell counterexample.
- `D:\Dev\math\docs\YANG-MILLS-CONDITIONAL-REDUCTION.md` §4.1 — P3 strong-coupling residue theorem, scope, `1/48` is the Ursell radius.

**External:**

- Osterwalder–Seiler, *Ann. Phys.* **110** (1978) 440 — strong-coupling OS axioms, area law, mass gap.
- Kotecký–Preiss, *Commun. Math. Phys.* **103** (1986) 19 — abstract polymer expansion.
- Shen–Zhu–Zhu, *Commun. Math. Phys.* **400** (2023) 805; arXiv:2204.12737 — BE/LSI strong-coupling region (`1/48` is theirs, not KP).
- Bauerschmidt–Brydges–Slade, *LNM* **2242** (2019) — constructive RG / BBS.
- Balaban, *Commun. Math. Phys.* **109** (1987) 249; **116** (1988) 1; **122** (1989) 175/355 — multiscale gauge RG; stops at the crossover.
- Osterwalder–Schrader, *Commun. Math. Phys.* **31** (1973) 83; **42** (1975) 281 — OS reconstruction.
- Faizal–Shabir, arXiv:**2606.19362**; *Fortschr. Phys.* **74** (2026) e70097; DOI **10.1002/prop.70097** — the paper audited in §5.
- V. Liashkevych, Zenodo record **22840696**, DOI **10.5281/zenodo.22840696** (2026-09-19) — the G-field proof-architecture paper; **not** an FS audit (see §5.7).
- 't Hooft, *Nucl. Phys.* **B138** (1978) 1; **B153** (1979) 141 — monopoles, flux sectors.
- Komargodski–Schwimmer, *JHEP* **12** (2011) 099 — 4D a-theorem (CFT fixed points, not a massive-flow Lyapunov function).
- Zamolodchikov, *JETP Lett.* **43** (1986) 730 — 2D c-theorem (contrast: no 4D analogue used here).
- Jaffe–Witten, *Quantum Yang–Mills Theory* (Clay problem description) — Wightman/OS axioms, mass gap; vacuum representation.

---

*Status of this note: **[C]/[NG]/[O]**, not a proof of Clay. It records the exact missing step and the verified audit points; it does not claim the wall has fallen.*
