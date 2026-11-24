# E2-FLIP-L1-VNW — attacking the UNIFORMIZATION of the one genuine non-H1 escape: same-channel von Neumann–Wigner *pointwise* isolation of the lightest `0⁺⁺`. Can the level ordering be made `(L,a,β)`-uniform without presupposing the gap?

> 2026-11-24 · Head-on attack on **escape L-1** of `E2-FLIP-NOGO` / `E2-REGRESS-CIRCULARITY`: the same-channel vNW isolation reads **LEVEL ORDERING** (not the spectral measure `ν`), so it is genuinely *outside* the class `C0` (the four-horn dichotomy) and is the **only** surviving non-`H1` loophole. That loophole is **[T]-shaped but SUB-UNIFORM**. This note attacks the one thing that would turn it into a proof of Clay: making the isolation **uniform** in `(L,a,β)` — an **absolute** margin `E₂−E₁ ≥ κ`, not a relative one. Four concrete uniformization mechanisms are built and tested: (i) symmetry / quantum-number / radial-excitation rigidity; (ii) a 't Hooft / non-invertible anomaly forbidding the crossing; (iii) a different `J^{PC}` channel (`2⁺⁺`) with provable uniform ordering; (iv) the two-particle (multi-glueball) threshold keeping the atom below the cut uniformly in `L`.
>
> **Labels (charter §2):** **[T]** 真定理 / **[C]** 条件定理 / **[NG]** no-go for a *named class* / **[O]** 精确缺失步 / **[J]** 判断 / **[V-here]** fresh computation this session (`research/unconditional/w7/e2_flip_l1_vnw_witness.py`, run this session, output in §A). External/library facts are **[R]** (standard/as-reported; re-verified no external source this session). **No QED on anything unproven. This document does NOT prove Clay and does NOT prove that no tool can exist (`∀class ≠ ∀possible`, `KEY-EXISTS §5`).**

---

## Verdict up front

> ### **IRREDUCIBLE — uniformizing the vNW isolation is exactly `ρ = H*` PLUS the gap itself, and no symmetry / anomaly / channel / threshold supplies the missing absolute scale.**

**The precise reason, one paragraph.** vNW delivers a **relative** ordering: `E₂ > E₁`, i.e. `E₂ ≥ (1+ε)m₁` with `ε>0` *pointwise*. Temple's `ρ`-slot needs the opposite kind of object: a **constant** `C` (independent of the unknown `m₁`) with `C ≤ E₂` — an **absolute** floor. By Single-Scale Covariance (`E₂ = ρ_true·m₁`, `ρ_true` a fixed number) these are not the same statement: **`E₂ ≥ C` uniformly  ⟺  `m₁ ≥ C/ρ_true`  = the gap.** So "uniform vNW isolation" in the Temple-sound (absolute) sense **is** the gap, up to the ratio `ρ_true` — it cannot be an *input* to proving the gap; it presupposes it. And a certified **relative** ratio `E₂ ≥ (1+c)E₁` fed to Temple is **an exact tautology**: with the true ratio, Temple returns `m₁ ≥ m₁` (zero information); with a weaker certified ratio it returns a fixed **multiple** `K·m₁`, still homogeneous, still **no absolute floor** ([T], §1 + [V-here]). The four uniformization mechanisms each fail to break this: **(i)** symmetry labels *sectors*, and `0⁺⁺`/`0⁺⁺*` are the *same* sector ⇒ any symmetry acts identically on both, fixing at best a protected *ratio*, never an absolute spacing; **(ii)** an anomaly is *scale-invariant* by definition (a dimensionless 't Hooft bracket) ⇒ can never manufacture a dimensionful gap, and a same-sector pair carries no anomaly charge to distinguish it; **(iii)** a `2⁺⁺`-channel ordering is the **wrong sector** (a `2⁺⁺` bottom does not lower-bound the overall `m₁ ≤ m_{2⁺⁺}`) **and** the `0⁺⁺`/`2⁺⁺` pair is *exactly degenerate* at leading strong coupling with **no** cross-channel inequality existing (`PROOF-IN11`); **(iv)** the margin below the two-glueball cut is `2m₁−m₁ = m₁` = the answer, and `L`-uniformity needs `m₁` on both ends (small `L`: gap `∝1/L` ≠ `Λ`; large `L`: Lüscher corrections `∝ e^{−m₁L}`). **The exact stuck step:** a `β,a,L`-uniform **ratio** `ρ ≥ 1+c` (= the relative, acid-clean half) **AND** the absolute scale `m₁ ≥ μ` (= the gap). vNW supplies the pointwise `ρ>1` (genuine, credited), **none** of the uniformity. This is the same `O(u)`-decoration sign wall (`SUB2`, `Cu` vs `Be^{−δt}`) as `H*`, in the level-ordering coordinate.

**What is *not* claimed.** This is **[NG] for the four named uniformization classes**, **[T-conditional on SSC]** for the structural reason, **[O]** for the exact missing step, **[J]** for "no second scale exists." It is **not** a global impossibility theorem, and it does **not** under-sell the pointwise vNW break, which is real (`[T]`-shaped) and acid-clean — it is simply **sub-uniform**.

---

## §0 What "uniform vNW isolation" must mean, and the direction requirement

### 0.1 The two objects (and a direction check on the task's phrasing)  [T]
vNW same-sector no-crossing: two eigenvalues of one Hermitian `H(β)` in the same irrep cannot cross as `β` varies (von Neumann–Wigner 1929; Reed–Simon IV §XIII.1). Given strong-coupling cost-grading `m₁ < m_{0⁺⁺*}`, the ordering is sign-constant on the connected `β`-line (confinement is a crossover), so `m₁ < m_{0⁺⁺*}` **pointwise wherever both are bound** (`G2D §3`, `PROOF-NO-DISSOLUTION-KINEMATIC §3`). This is **relative isolation**: `E₂ ≥ (1+ε)E₁`, `ε>0`.

Temple (`TEMPLE-WEINSTEIN §1`) on `H|_{(vac)^⊥}`: **if `ρ ≤ E₂` and `m₁ ≤ R < ρ`, then `m₁ ≥ R − Var/(ρ−R)`.** To output an `E₁`-lower bound (the gap) the `ρ`-slot must hold a **constant lower bound on the second eigenvalue** `ρ ≤ E₂`, `ρ` independent of the unknown `m₁`. That is **absolute isolation**: `E₂ ≥ C`, `C` a fixed number.

> **Direction lemma [T].** An **upper** bound on `E₂` (always variationally available: eigenvalues are min-max = inf-of-max, hence *upper*-computable) **cannot lower-bound `E₁`.** Only a **lower** bound `ρ ≤ E₂` (Temple/Weinstein) can. So vNW's sole possible route into a proof is to *certify `ρ ≤ E₂` uniformly* (= READ-A below); "an `E₂` upper bound from vNW" is vacuous for the gap. *(The task's parenthetical "E₂ upper from vNW + E₁ lower" is therefore read, correctly, as "vNW fills Temple's `ρ ≤ E₂` slot ⟹ Temple outputs an `E₁` lower bound.")*

### 0.2 Two readings of "uniform isolation"
- **(READ-R) Relative:** `E₂ ≥ (1+c)E₁` for a certified constant `c>0`, uniform in `(L,a,β)`. This is the **ratio** `ρ = m_{0⁺⁺*}/m₁ ≥ 1+c`. **Acid-clean** (free field: no `0⁺⁺` tower, §4). It is real progress but **not** a Temple anchor for the absolute gap (§1).
- **(READ-A) Absolute:** `E₂ ≥ C` for a fixed number `C>0` (equivalently `E₂−E₁ ≥ κ`, `κ>0`, uniform). This *is* what Temple needs. §1 shows READ-A **⟺ the gap**.

### 0.3 Lemma L1 (recap)  [T]
*(`TEMPLE-WEINSTEIN §3`, `E2-REGRESS §0.2`.)* No function of `(R, Var)` alone lower-bounds `inf supp ν`; a valid trial bound must smuggle a **hard edge** `ρ ≤ E₂`. This is why READ-A (a hard absolute edge) is indispensable and READ-R (a ratio) is not enough.

---

## §1 The logical spine — absolute isolation is the gap; relative isolation is Temple-inert

### 1.1 Theorem A (READ-A ⟺ gap)  [T-conditional on SSC]
Suppose Single-Scale Covariance: `E₂ = ρ_true·m₁` with `ρ_true` a fixed `(β,L,a)`-independent number (`[J]`, empirical `ρ≈1.5–1.6`, `G2D §2`, `RATIO-RHO §5`). Then the uniform absolute-isolation statement
`(READ-A)  ∃C>0: E₂ ≥ C  ∀(L,a,β)`
is **equivalent** to the uniform gap statement
`m₁ ≥ C/ρ_true  ∀(L,a,β).`
*Proof.* `E₂ = ρ_true m₁` gives `E₂ ≥ C ⟺ m₁ ≥ C/ρ_true`. Both directions. ∎

> **Consequence.** READ-A is **not an input** that could prove the gap — it **is** the gap, written with `ρ_true` in the exponent. Equivalently (`G2D`'s correction): uniform READ-A ⟺ `gap ∧ ratio` (it needs *both* an absolute scale and a bounded ratio). It can no more supply the gap than "the gap" can. This is the sharpest single sentence of the L-1 attack: **the sound Temple anchor *requires* READ-A, and READ-A *is* the gap.**

### 1.2 Theorem R (READ-R is Temple-inert)  [T; witness §A]
Suppose only **relative** isolation is certified: `E₂ ≥ (1+c)E₁` (`c>0` constant), and feed Temple the best available slot `ρ = (1+c)E₁`. Use a trial `ψ = √(1−s)ψ₁ + √s ψ₂` in the eigenbasis (true energies `E₁=m₁`, `E₂=ρ_true m₁`, `ρ_true ≥ 1+c`). Then `R = m₁[1+s(ρ_true−1)]`, `Var = s(1−s)(ρ_true−1)²m₁²`, and
```
m₁ ≥ R − Var/(ρ − R)  =  m₁ · K,   K := 1 + s(ρ_true−1) − s(1−s)(ρ_true−1)² / [ (1+c) − 1 − s(ρ_true−1) ],
```
a **fixed multiple** `K·m₁` of the unknown mass — **independent of `m₁`'s value**. In the extremal case `c = ρ_true−1` (certified ratio equals true ratio), `K = 1` exactly: **Temple returns `m₁ ≥ m₁`, an exact tautology.** In no case is there an **absolute** floor: `m₁ → 0` remains admissible.
*Proof.* Direct algebra (above); the homogeneity is manifest — every term is `∝ m₁`. **[V-here]** `research/unconditional/w7/e2_flip_l1_vnw_witness.py` (output §A) confirms `out/m₁ = 0.99744` identically for `m₁ ∈ {0.1,1,3,10,100}` (certified `1.4`, true `1.5`), and `out = m₁` exactly when certified `=` true `1.5`. ∎

> **Consequence.** Even if one *proved* the ratio `ρ ≥ 1+c` (the hardest open sub-target `SUB2`/`H*`), it would **not** yield the gap through Temple: READ-R is a scale-covariant (homogeneous) statement and admits `m₁→0`. The ratio is genuinely useful only **bolted onto an already-known absolute scale** — e.g. the strong-coupling SZZ gap at `β<1/48` — which is why `G2D §5` says ratio + SZZ closes uniform I2 *at strong coupling*; it does **nothing** for the continuum `a→0` and nothing at general coupling.

### 1.3 The exact identity of the wall
Combining §1.1–§1.2:
```
sound uniform Temple anchor  =  READ-A  =  gap ∧ ratio           (§1.1)
                              =  (m₁ ≥ μ)  ∧  (ρ ≥ 1+c).
vNW supplies:                    ρ > 1  pointwise only,  NO uniform c,  NO μ.      (§3)
∴ uniformizing vNW  =  (the open ratio 1+c)  ∧  (the gap μ)   =  H*  in level-ordering coordinates.   ⟹  [O]
```
**Nothing is gained by "uniformizing" that was not already both halves of `H*`.**

---

## §2 The four uniformization mechanisms, each attacked

### Route (i) — symmetry / quantum-number barrier / radial-excitation rigidity  → [NG]
- **Symmetry labels *sectors*, not the splitting.** Ground `0⁺⁺` and excited `0⁺⁺*` are the **same** sector: same `J^{PC}=0⁺⁺`, both `A1g`, both N-ality 0, both trivial under the 1-form center `Z_N^{(1)}` (which charges *lines*; every glueball is `N`-ality 0 ⇒ blind, `PROOF-IN11 §(ii)`), both carry no generalized (higher-form / non-invertible) charge, both built from the same single-plaquette operator. Any symmetry — invertible or not — acts **identically** on the two, so it can constrain *neither* their splitting *nor* an absolute floor. The **most** a symmetry can do is **order** them (= vNW no-crossing, the pointwise fact), not give a margin.
- **"Rigidity" fixes ratios, not spacings.** A protected level ratio (hydrogen `E_n=−Ry/n²`; Landau ladder `m_n²=(2n+1)m₁²`; integrable-Toda S-matrix) is `E_n = (pure number)·Λ`: it fixes `ρ_true` (READ-R) but the *absolute spacing* `E₂−E₁ = (ρ_true−1)Λ` is still `∝Λ` — a fixed multiple of the single scale. **The Landau-ladder `ρ=√3` is refuted** as a glueball ratio: it is the **single-gluon** fluctuation spectrum in a constant background, **tachyonic** at the bottom, and realised verbatim by a **free** charged particle (acid-test-vacuous) — `RATIO-RHO §2`. Pure 4D YM has **no** exact rigidity: not integrable, no Coulomb-branch modulus, not SUSY.
- **Verdict.** At best a symmetry supplies a protected **ratio** `ρ_true` (READ-R). It cannot supply READ-A's absolute `C`; and even a proven `ρ_true` fails §1.2 (Temple-inert). **[NG]** for "symmetry uniformity of the *absolute* isolation."

### Route (ii) — 't Hooft anomaly / non-invertible symmetry forbidding the crossing  → [NG]
- **An anomaly is scale-invariant.** An 't Hooft anomaly is the failure of a symmetry to be gauged, measured by a **dimensionless** descent bracket `𝒜` (a cocycle on the symmetry category). Being dimensionless, it **cannot** produce a dimensionful gap `C` — the same reason a symmetric phase can be gapped or gapless: anomalies fix *content*, not *energies* (`E2-REGRESS §4.3`, `PROOF-1FORM-ANOMALY`).
- **It constrains sector content, not intra-sector spacing.** Anomaly matching says *which representations / which generalized charges* must appear in the IR. `0⁺⁺` vs `0⁺⁺*` are the **same** representation and the **same** (trivial) charge under every generalized symmetry ⇒ the anomaly is reproduced by the single `0⁺⁺` ground state and constrains **nothing** about the excited `0⁺⁺*`. At `θ=0` the electric 1-form anomaly is trivially matched by the gapped confining vacuum; at `θ=π` a mixed anomaly forces **at most** a gapped TQFT (`PROOF-1FORM-ANOMALY`) — a *gap* statement that is still not a *splitting* statement, and trivially scale-free. Non-invertible symmetries (`Rep Z_N` at `θ=π`, etc.) **forbid operators / couplings**, not energy splittings between two same-symmetry states.
- **The only no-crossing content is vNW itself**, and vNW gives relative `ρ>1`, not absolute `C`. **[NG]** — an anomaly can forbid the *crossing* (already vNW), never *produce* a dimensionful uniform floor; it is sector-and-scale-free.

### Route (iii) — a different `J^{PC}` channel (`2⁺⁺`) with provable uniform ordering  → [NG] (double failure)
- **(a) Wrong sector, wrong direction.** Clay needs `spec(H) ∩ (0,Δ) = ∅` over **all** channels, i.e. a lower bound on `m₁ = min_{channels} m_c`. A lower bound in the `2⁺⁺` channel bounds `m_{2⁺⁺}`, and `m₁ ≤ m_{2⁺⁺}` — the **wrong direction** to lower-bound `m₁`. A cross-channel bottom never lower-bounds the global minimum.
- **(b) No cross-channel inequality exists.** Positivity / Källén–Lehmann is **intra-sector**; at `θ=0`, `P` and `C` make every **mixed** two-point kernel vanish identically, so there is **no** cross-channel inequality (`PROOF-IN11 §(ii)`). Worse, `0⁺⁺` and `2⁺⁺` are **exactly degenerate** at leading strong coupling (`a_{0⁺⁺}=a_{2⁺⁺}=1`, both single-plaquette operators); which is lighter is a subleading `O(1)` **dynamical** splitting that **no symmetry fixes** (`PROOF-IN11 §(i)`). Even the "which channel is lightest" fact is not a theorem.
- **(c) Even within `2⁺⁺`, uniform ordering = SSC.** `E^{(2⁺⁺)}₂−E^{(2⁺⁺)}₁ ≥ κ` ⟺ `m_{2⁺⁺} ≥ κ/(ρ₂−1)` by SSC ⇒ bounds `m_{2⁺⁺}`, not `m₁`, and is itself "the gap *in the `2⁺⁺` channel*," the same wall. **[NG]** on all three counts.

### Route (iv) — two-particle threshold keeping the atom below the cut uniformly in `L`  → [NG]
- **The margin IS the answer.** The lightest `0⁺⁺` sits below the two-glueball cut `θ=2m₁` by `θ−m₁ = m₁` (`PROOF-NO-DISSOLUTION-KINEMATIC §1`; `PROOF-DISSOLUTION-SPECTRAL-DENSITY §1`). "Uniform margin below the cut `≥ κ`" is **verbatim** "`m₁ ≥ κ`" — the gap (`E2-FLIP-KINEMATIC §1.3`, reason γ). *(The true nearest competitor is the same-channel excited atom `≈1.5m₁`, also `∝m₁` — `PROOF-ISOLATION §4.4`.)*
- **`L`-uniformity is circular on both ends.** Small `L ≪ 1/Λ`: the finite-volume gap is the **perturbative / Casimir scale** `∝ 1/L` (`[R]` standard, Lüscher femto-universe) — a **different scale** from `Λ`; there is no `L`-independent constant margin spanning small→large `L`. Large `L`: Lüscher finite-volume corrections to the atom and the threshold are `O(e^{−m₁ L})`, controlled uniformly only by knowing the decay rate `m₁` — circular.
- **Pinning is conditional, not proven.** The "atom below every cut" lemma assumes the gap is open (`m_{0⁺⁺} = m_g > 0`); it **inherits** the gap, does not prove it (`PROOF-NO-DISSOLUTION-KINEMATIC`, "exact condition," conditional (K-A)). **[NG]** — margin `=m₁`; `L`-uniform control needs `m₁`; the no-dissolution antecedent is the gap.

---

## §3 What vNW *does* achieve — and its precise ceiling (do not under-sell; do not over-sell)

- **Genuine pointwise break [T].** `m₁ < m_{0⁺⁺*}` at every `β` where both are bound, **without** presupposing `m₁>0` (it presupposes only the strong-coupling *ordering*). It is **outside class `C0`** — it reads the **ordering** of two levels, not a functional of `ν` — so it survives all four horns (α finite-data, β rigid, γ smooth, δ extremal) of `E2-FLIP-NOGO`. This is a real, non-circular, acid-clean partial. **Credited.**
- **Ceiling = relative, not absolute.** It gives `ρ>1` pointwise, **never** a uniform `c>0`, and **never** an absolute `C`. By §1.2 even a *proven* `c` is Temple-inert; by §1.1 a *proven absolute* `C` **is** the gap. So the pointwise break, however genuine, sits strictly **below** both halves of READ-A.
- **Dissolution is the OPE-blind disappearance of the *weight*, which vNW cannot see.** The excited `0⁺⁺*` "dissolves" not by a crossing (vNW forbids that) but by its residue `f₁²→0` at a **fixed below-threshold** pole — the free-scalar mechanism, invisible to OPE (`s→∞`), allowed by upper-semicontinuity of `f₁²` (`PROOF-DISSOLUTION-SPECTRAL-DENSITY`). vNW orders levels **while they are Kato-isolated**; it predicts **nothing** about their persistence. So even the pointwise statement lives only on the stretch where both states remain discrete, and stops where the atom weight vanishes — precisely the sub-uniformity `G2D §3` flags.

---

## §4 Free-massless acid test — precise (and what it does *not* discriminate)  [T; §A]

Free massless (`g=0`, `β→∞`): the `0⁺⁺` channel (the `¼tr F²` two-**photon** cut) is a continuum from `s=0` with **no isolated `0⁺⁺` bound state** — there is no tower to order, `E₁=E₂=0`.
- **READ-A tool:** needs `C ≤ E₂ = 0`; any `C>0` **fails** ⇒ returns nothing. **Acid-clean** (cannot manufacture a gap on the free field). ✔
- **READ-R tool:** reads `E₂ ≥ (1+c)E₁` as `0 ≥ 0` (trivial); but no `0⁺⁺` state exists to which it applies ⇒ **vacuous**, certifies no gap. ✔
- **Critical honesty point.** The free field gives **no isolation** — correct. **But the acid test does *not* separate the real pointwise vNW (genuine) from a would-be uniform one**: **both** are acid-clean, because **both** fail identically on the free field (no tower). The discriminator is **not** the acid test; it is **Single-Scale Covariance** (whether the margin can be absolute). Passing the acid test is **necessary** but does **not** imply closing Clay. The L-1 attack's entire content is the SSC failure of the *uniform* margin, **not** any acid-test shortfall.

---

## §5 The verdict and the exact stuck step

### 5.1 Verdict
> **IRREDUCIBLE.** Uniformizing same-channel vNW isolation is **provably equivalent, by SSC, to `gap ∧ ratio = H*`** (§1.1). No symmetry, anomaly, `J^{PC}` channel, or threshold supplies the missing **absolute scale** (§2). The pointwise vNW break is real and acid-clean but is **relative isolation**, which Temple cannot convert into an absolute floor (§1.2). **[NG] for the four named uniformization classes; [T-conditional on SSC] for the equivalence; [O] for the exact missing step.**

### 5.2 The exact stuck step (one box)
```
need (sound, uniform Temple anchor):   READ-A:  E2 >= C ,  C a NUMBER indep of m1 ,  (L,a,beta)-uniform
have (vNW, genuine):                   READ-R:  E2 > E1 ,  i.e. E2 >= (1+eps) m1 ,  eps>0 POINTWISE
bridge (single step):                  (1+eps) m1 >= C uniformly  <=>  m1 >= C / (uniform eps)
                                       requires BOTH  a uniform eps (= ratio rho >= 1+c)  AND  m1 >= mu.
equivalent (by SSC, E2 = rho_true m1):  READ-A  <=>  m1 >= C/rho_true  = the GAP.
∴ uniformizing vNW presupposes the gap; it cannot prove it.   EXACT MISSING STEP = H*  (ratio 1+c  +  gap mu).
```
The **relative** half (`ρ ≥ 1+c`) is acid-clean, `β,L`-invariant, and is the sharpest genuine open sub-target — but proving it is the `O(u)`-decoration sign estimate `SUB2` (`Cu` vs `Be^{−δt}`) = the scale-invariant face of `H*` (`RATIO-RHO §4`). The **absolute** half (`m₁ ≥ μ`) is the gap (continuum transport, untouched). vNW supplies the **pointwise** `ρ>1`; nothing else.

### 5.3 What would flip this to FLIP  [precise]
Exhibit **one uniform absolute isolation** `E₂ ≥ C` with `C` independent of `m₁`, **provable without assuming** `m₁ ≥ C/ρ_true` — i.e. a **second dimensionful scale** `Λ₀ ≢ multiple-of-Λ_QCD` with a proven `m₁ ≥ f(Λ₀)` — **or** a symmetry / anomaly that produces a **dimensionful** intra-sector gap (impossible for a scale-invariant anomaly, §2-ii). None exists; pure YM at `T=0`, infinite volume, `θ=0`, no fermions has **one** scale (`G2D §2`, `E2-FLIP-KINEMATIC §6`). The vNW pointwise break is **not** such a `Λ₀`.

---

## §6 Honest boundary — what this document does NOT say

- It does **not** prove Clay and does **not** prove no proof exists. It proves a **[NG] for the four named uniformization mechanisms**, a **[T-conditional on SSC]** equivalence (READ-A ⟺ gap), and records the **[O]** exact missing step. `∀class ≠ ∀possible` (`KEY-EXISTS §5`, `CONVERGENCE-redteam §d`).
- The **SSC equivalence (§1.1)** is **[T]-conditional**: `E₂ = ρ_true m₁` with `ρ_true` fixed (`[J]`, empirical `ρ≈1.5–1.6`, `RATIO-RHO §5`). If `ρ_true` were unbounded one could have `E₂ ≥ C` with `m₁→0`; pure YM gives no mechanism for a mass hierarchy at the gapless point (QCD-like, one scale). The *dimensional* content (`ρ_true` bounded ⇒ READ-A ⟺ gap) is `[T]`; the boundedness is `[J]`.
- It does **not** under-sell the pointwise **vNW** break (§3): genuine, non-circular, outside `C0`, acid-clean. It simply is **READ-R**, not **READ-A**.
- The small-`L` "`∝1/L` Casimir scale" (§2-iv) is `[R]`/standard (Lüscher finite-volume), not re-derived here.
- The **[V-here]** witness (§A) is a two-level check of *already-banked* lemmas (Temple soundness, `L1`, SSC); it is not evidence about the world.
- **Do not conflate:** "uniform vNW isolation" in the **relative** (READ-R) sense being achievable-but-sub-uniform is *not* a gap proof (§1.2); only **absolute** (READ-A) isolation is, and that **is** the gap (§1.1).

---

## Appendix A — [V-here] numerical witness (`w7/e2_flip_l1_vnw_witness.py`, run this session)

Pure-Python, no deps. Reproduced output (essence):
```
(1) RELATIVE isolation rho_slot=(1+c)*m1, certified 1.4, true 1.5:
    m=0.1..100 : Temple_out/m = 0.99744  IDENTICAL for all m  -> homogeneous, NO absolute floor.
(1b) certified ratio == true ratio (1.5):  Temple_out = m  exactly ( == 'm >= m' tautology ).
(2) ABSOLUTE isolation rho=C=1.5 (fixed):  floor appears only where C<=E2=1.5*m i.e. m>=1.0,
    which was ASSUMED to write C<=E2  -> CIRCULAR.  (m=0.5: C>E2 -> the 0.5038 row is an
    UNSOUND overclaim, reproducing E2-FLIP-KINEMATIC Model A.)
(3) FREE massless E1=E2=0: absolute tool needs C<=0 -> any C>0 FAILS (acid-clean, inert);
    relative tool reads 0>=0 trivial, no tower -> vacuous. Free field gives NO isolation.
(4) SCALE COVARIANCE: Temple_out/lambda = 0.99744 for lambda=0.01,1,100 (identical).
```
These witness `L1` + Temple soundness + SSC (the exact absolute-vs-relative distinction), **not** new theorems.

---

## §7 References

**Internal (absolute paths; the banked lemmas this note attacks/assembles):**
- `D:\Dev\math\CLAUDE.md` (§1 build-don't-stop; §2 labels; §5 wall/machine); `D:\Dev\math\docs\BRIEF-FOR-MODEL.md` (§1 corrected `H*`; §2 five constraints; §3 death map).
- `D:\Dev\math\research\unconditional\E2-FLIP-NOGO.md` (§5 (L-1) the vNW pointwise loophole this note attacks; §1 class `C0`); `E2-REGRESS-CIRCULARITY.md` (§0.4 SSC; §0.6 Temple two horns; §2.4 vNW corner; §5 exact missing steps); `E2-FLIP-KINEMATIC.md` (§1.3 margin `=m₁`; §2 cross-channel table; §4 vNW not-this-escape; §6 flip condition).
- `D:\Dev\math\research\unconditional\G2D-ISOLATION-I2.md` (§2 SSC; §3 vNW pointwise; §5 ratio `=H*`); `RATIO-RHO.md` (§2 LLR refuted — free/tachyonic/wrong object; §3 no symmetry protection; §4 ratio `=H*`; §5 `ρ_lat≈1.5`); `KEY-EXISTS-OR-IMPOSSIBLE.md` (§5 `∀class≠∀possible`; §7 calibration).
- `D:\Dev\math\research\unconditional\PROOF-NO-DISSOLUTION-KINEMATIC.md` (§1 pinning margin `=m_g`; §3 same-channel vNW vs cross-channel; exact condition K-A); `PROOF-DISSOLUTION-SPECTRAL-DENSITY.md` (§1 operative `Λ`=excited atom; §atom `f₁²→0` OPE-blind); `PROOF-ISOLATION.md` (§4 (P3) residue; operative `Λ≈1.5m_g`); `PROOF-IN11.md` (§(i) `0⁺⁺`/`2⁺⁺` exact leading-order degeneracy, `a_c=1`, no symmetry orders; §(ii) no cross-channel inequality, positivity intra-sector); `GAP-G2-edge.md` (L2 residue–level); `TEMPLE-WEINSTEIN.md` (§1 Temple, §3 L1); `LEVER-nonconvex-edge.md` (§0B linear-data lemma); `PROOF-1FORM-ANOMALY.md` (1-form anomaly blind to local `trF²`; `θ=π` gapped TQFT); `SCHUR-CONVEX.md` (Family-I upper); `CONVERGENCE-redteam.md` (§d).
- `D:\Dev\math\research\unconditional\w7\e2_flip_l1_vnw_witness.py` (`[V-here]` §A; run `py -3.12 e2_flip_l1_vnw_witness.py`).

**External — [standard]/[R] (textbook; not re-verified this session):**
- von Neumann–Wigner (1929) — same-irrep level repulsion / no crossing; Kato, *Perturbation Theory for Linear Operators* (1976) VII §3; Reed–Simon IV §XIII.1, Thm XIII.5 (Temple/Weinstein); Temple, *Proc. London Math. Soc.* (1928); Weinstein, *Phys. Rev.* **40** (1932) 797; Courant–Fischer (min-max = upper).
- Osterwalder–Seiler, *Ann. Phys.* **110** (1978) 440; Källén (1952); Lehmann (1954); Lüscher finite-volume / femto-universe (small-`L` finite-volume gap); Shen–Zhu–Zhu, *CMP* **400** (2023) 805, arXiv:2204.12737 (strong-coupling gap); Morningstar–Peardon, *PRD* **60** (1999) 034509 (quenched spectrum `0⁺⁺`1710<`2⁺⁺`2390<`0⁻⁺`2560 MeV).
- Gaiotto–Kapustin–Seiberg–Willett, arXiv:1412.5148; Cordova–Ohmori, arXiv:1910.04962; Cordova–Dumitrescu–Intriligator, arXiv:1802.04790 (1-form / non-invertible anomalies, dimensionless 't Hooft brackets).

---

*Status: **[NG] for the four named uniformization mechanisms** (symmetry / anomaly / different-`J^{PC}` / two-particle-threshold, each failing to supply an absolute scale); **[T-conditional on SSC]** that uniform absolute vNW isolation ⟺ `gap ∧ ratio = H*` (so it presupposes the gap) and that certified *relative* isolation is Temple-inert (`m₁ ≥ m₁`); **[O]** exact missing step = uniform ratio `ρ≥1+c` (relative, acid-clean, the `O(u)`/`SUB2` wall) **plus** the absolute `m₁≥μ` (continuum transport). The pointwise vNW break is real, outside `C0`, acid-clean, but **sub-uniform**; the free-massless acid test passes yet does **not** discriminate it from a would-be uniform one (SSC does). **Verdict: IRREDUCIBLE — uniformizing L-1 is `ρ=H*` plus the gap. NOT a proof of Clay; NOT a global impossibility theorem (`∀class ≠ ∀possible`). Clay stays open.***
