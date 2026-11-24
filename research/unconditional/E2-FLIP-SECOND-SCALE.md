# E2-FLIP-SECOND-SCALE — is there a *genuine independent second dimensionful scale* in 4D pure SU(N) YM that beats the Scale Dichotomy (Thm S), or is the theory dimensionally rigid (single-scale)?

> 2026-11-24 · Adversarial hunt, on the **one named escape** from the E₂ no-go: `E2-FLIP-NOGO §3` (Theorem S: a uniform dimensionful flip functional is forced onto `Λ_QCD`, and `Φ≤m₁=c·Λ` forces knowing `c=m₁/Λ` = the gap) and `E2-FLIP-KINEMATIC §2/§6` (the kinematic-E₂ escape dies "because pure YM lacks an independent second dimensionful scale"). The task: **find such a scale, or confirm rigidity.** Six candidates examined — topological/θ, rank/center/level (N), deformation (double-trace / center-stabilized / 't Hooft-twist), finite volume L, defect/'t Hooft loop/domain wall, anomaly/trace-anomaly.
>
> **Labels (`CLAUDE.md §2`):** **[T]** 真定理 / **[C]** 条件定理 / **[NG]** no-go for a *named class* / **[O]** 精确缺失步 / **[J]** 判断 / **[V-here]** fresh computation this session. External/banked facts are **[R]** (as-reported from the repo; not re-verified this session). **No QED on anything unproven. This does NOT prove Clay, and it is NOT a global impossibility theorem** (`∀class ≠ ∀possible`, `KEY-EXISTS §5`).

---

## Verdict up front

> ### **RIGIDITY-CONFIRMED [NG + T-dimensional].** No genuine independent second dimensionful scale exists in pure 4D SU(N) YM at `T=0`, `θ=0`, infinite volume, zero fermions. **Every candidate either (A) is `∝Λ` — and then it is *not independent* and cannot anchor a floor (Lemma SA below); or (B) is a genuine external input — and then it fails the Clay requirement** (L-uniformity, monotonicity direction, or connection-to-pure-YM). There is no third horn.

**One paragraph.** The single-scale fact is dimensional, not empirical: pure YM has **one dimensionless coupling** `g` (equivalently one `Λ_QCD` by dimensional transmutation) plus **discrete data** (`N`, `θ/2π`) plus **regulators** (`a`, `L`). Every dynamical mass — `√σ`, `m_{0⁺⁺}`, `T_c`, `χ_t^{1/4}`, `⟨T^μ_μ⟩^{1/4}` — is a **fixed pure-number multiple of that one scale**, `Q=q·Λ`. The sharp point the six-candidate hunt turns on: *an anchor `m₁ ≥ α·Q` with `Q∝Λ` is, by homogeneity, the pure-number inequality `α ≤ m₁/Q` — i.e. it presupposes the very ratio `m₁/Λ` it was meant to bound* (**Lemma SA**, §1). So being `∝Λ` is not a technicality; it is fatal **regardless of whether `Q` lies inside or outside the E₂ class `C0`** — which is exactly what extends Theorem S beyond the spectral class. The only `Q` that escapes Lemma SA are **external scales we set by hand** (`L`, `T`, a deformation coupling `κ`). Those are checked in §2 and each fails for a *structural* reason, not a missing estimate: `L`/`T` are monotone in the **wrong direction** (`Δ(finite volume) ≥ Δ(ℝ⁴)`, so a finite-region bound is an *upper* bound on the target); deformations define a **different theory** whose pure-YM limit is the gap itself. Hence the dimensional rigidity is confirmed for the full family, and the kinematic-E₂ escape cannot be re-run: there is no `Λ₀ ≢ multiple-of-Λ` to feed it.

**What is and is not claimed.** `[T]`: the dimensional inventory (one coupling ⇒ one scale up to pure numbers) and the homogeneity reduction of Lemma SA. `[NG]` for the named class: "a flip functional / Temple-`E₂` anchor built from any dimensionful `Q` that is `∝Λ` *or* is an external scale whose pure-YM limit is taken." **Not** claimed: that no future mathematics can prove the gap — only that the *specific* "second-scale anchor" shape is closed, six-for-six, with a structural reason each.

---

## §0 The escape's exact logical form — and why it is *not* the same object as Theorem S

### 0.1 What Theorem S already closed (and what it left)
`E2-FLIP-NOGO §3` (Theorem S) says: a uniform dimensionful flip functional `Φ` is, on an `(L,a)`-uniform gauge-invariant **local spectral** functional, a function `Ψ(Λ_QCD)` of the only scale; `(LB)` then requires `Ψ(Λ)≤m₁=c·Λ`, so proving `m₁≥Δ` needs `c=m₁/Λ` — circular. **Its scope is `Φ ∈ C0`** (functionals of the KL spectral measure / local 2-pt data; `E2-FLIP-NOGO §1`).

### 0.2 What the second-scale escape actually proposes — and the gap in Theorem S's scope
The escape leaves `C0`: take a dimensionful quantity `Q` that is **not** a functional of the `trF²` spectral measure — the **string tension** `√σ` (1-form sector), the **topological susceptibility** `χ_t^{1/4}`, the **deconfining temperature** `T_c`, a **defect/lattice** scale, an **anomaly** scale — and anchor a *kinematic* threshold to `Q`: "`m₁ ≥ f(Q)`" with `f` explicit, `Q` proven `>0`, and `Q` **not** `∝Λ`. If such `Q` exists, Theorem S's "only `Λ`" inventory is bypassed and the kinematic-E₂ (`E2-FLIP-KINEMATIC`) can be re-run with `Q` in the Temple slot.

So the escape needs, precisely, a quantity `Q` with **three** properties:
- **(P1) genuine:** `Q>0` provable in pure YM **without** assuming the gap;
- **(P2) not `∝Λ`:** `Q/m₁` is *not* an unknown pure number that must itself be bounded (else circular);
- **(P3) anchor:** an explicit `f` with `m₁ ≥ f(Q)` **uniform in `(L,a)`**, `f(Q)>0`, and **acid-clean** (`f(Q)=0` on the free massless field).

`P1+P3` is the E₂-flip condition; `P2` is the new demand that would beat Theorem S. The six candidates are exactly the attempts at `P2`.

### 0.3 The spine (why the hunt is decidable)
Pure 4D SU(N) YM is defined by **one dimensionless coupling** `g₀(a)` on the lattice (continuum: one `Λ`), **discrete data** `N∈ℤ₊`, `θ/2π∈[0,1)`, and **regulators** `a`, `L`. There is **no other continuous dimensionful parameter** and **no other continuous dimensionless coupling**. Hence:
> **[T] (dimensional inventory).** Every dynamical dimensionful observable of pure YM is `Q = q(g,N,θ)·Λ^{d_Q}` with `q` a pure number. The `g`-dependence of `q` is only through `Λ` (the single scale) and the running; there is no second independent dimensionful input.

Two immediate consequences, which the six candidates exhaust:
- `q` may be **large or small** (a numerical hierarchy is allowed — e.g. `m_{0⁺⁺}/√σ≈3.9`), but it is still **one scale up to pure numbers**. "Numerical separation" is not "independence."
- The only quantities **not** of the form `q·Λ^{d}` are **external inputs** `L`, `T`, deformation couplings `κ`, which are *set by hand* and are **absent from the `T=0`, `∞`-volume, undeformed theory** that Clay concerns.

The rest of this document proves that horn (A) (`∝Λ`) cannot anchor (**Lemma SA**) and audits horn (B) (external) candidate-by-candidate.

---

## §1 Lemma SA — the scale-anchor dichotomy  [T + J-dictionary]

This is the structural heart: it upgrades "Theorem S is only about `C0`" to "**no `∝Λ` quantity, in or out of `C0`, can anchor**."

**Lemma SA (scale-anchor insufficiency).** Let `Q = q·Λ^{d}` be a dimensionful observable of pure YM (`d>0`, `q` a pure number), and let `m₁ = c·Λ` be the gap. Consider an anchor of *homogeneous* form
```
        (Anchor)        m₁  ≥  α · Q^{1/d}          (α a "provable" dimensionless constant).
```
Then:
- **(SA-1) [T].** `(Anchor)` ⇔ `c ≥ α·q^{1/d}`, a **purely dimensionless** inequality. Consequently: (i) the *maximal* `α` for which `(Anchor)` is true is `α_max = c/q^{1/d} = m₁/Q^{1/d}` — the unknown **ratio**; (ii) proving `(Anchor)` for any fixed `α>0` is proving `c ≥ α q^{1/d}`, i.e. a **lower bound on the gap ratio `c=m₁/Λ`** (up to the measured `q`); and (iii) `(Anchor)` is **scale-invariant** — under `Λ→λΛ` both sides `∝λ` and the inequality is unchanged, so a `∝Λ` anchor **cannot detect** the absolute scale, only restate a ratio.
- **(SA-2) [J/O, the dictionary].** Even securing `Q>0` (i.e. `q>0`) is **not** a scale measurement — it is a **dimensionless condensate/density lower bound** (`Q ≥ q₀ Λ^{d}`), which is the charter's `⟨¼trF²⟩`-type object, one of the equivalent faces of the key (`CLAUDE.md §5`: "`G″(0)>0` … / 胶子凝聚 `⟨¼trF²⟩` — 同一个对象的多个名字"). Bounding it below is a nonperturbative order-parameter estimate of the same difficulty class as the gap.

*Proof of SA-1.* Substitute `m₁=cΛ`, `Q^{1/d}=q^{1/d}Λ`; `Λ>0` cancels. Scale-invariance is the same cancellation. ∎ *(The only escape from SA-1 is `Q` **not** `∝Λ` — i.e. `Q` an external scale, or a quantity whose ratio to `m₁` is fixed by a theorem not reducible to dimensional analysis. §2 shows none exists.)*

**Why this is strictly stronger than Theorem S.** Theorem S's `(D)` considered `Φ=Ψ(Λ)` with `Φ∈C0`. Lemma SA makes **no membership assumption** on `Q`: the string tension, `χ_t`, `T_c`, `⟨T^μ_μ⟩` are all *outside* `C0` (they are not functionals of the `trF²` KL measure), yet SA-1 still kills them, because it uses only **dimensionality**, not the spectral class. This is the precise sense in which the second-scale hunt cannot rescue the escape by "leaving `C0`."

**The only loophole Lemma SA names.** A `Q` that is **not** `∝Λ`: (a) an **external input** (`L`,`T`,`κ`), or (b) a `∝Λ` quantity whose ratio to `m₁` is pinned by a **symmetry/kinematic theorem** (not by dimensional analysis) — e.g. same-channel von Neumann–Wigner. §2/§4 show (a) fails *uniformity/direction/connection* and (b) is the known **sub-uniform** vNW corner (`G2D-ISOLATION-I2 §3`), which is a *level-ordering* statement, not a *second scale*, and never yields a uniform `κ`.

---

## §2 The six candidates — audited, each with its acid test and its named failure

For each: **is it `∝Λ`?** (→ SA-1 kills it) · **if not, does it survive (P3) uniformity/direction?** · **acid test** (`Q` or the bound must *fail* on the free massless field, per charter §1) · **verdict + label** · **banked anchor**.

### (i) Topological / θ-angle scale — "Λ_θ"  [NG]
- **What it is.** `θ` is an **angle**: dimensionless, `θ/2π∈[0,1)`. The only `θ`-dependent *dimensionful* object is the vacuum-energy curvature `E''(θ)|_{θ=0}=χ_t` (topological susceptibility), dimension 4, and its scale is `χ_t^{1/4}≈qΛ` with `q≈0.4–0.7` (`[R]` lattice).
- **∝Λ?** **Yes.** `χ_t = q_χ Λ⁴`, `q_χ` a pure number. By **SA-1** it cannot anchor (`χ_t^{1/4}` in the Temple slot reduces to `c ≥ α q_χ^{1/4}` = the gap). The would-be "`Λ_θ`" is `χ_t^{1/4} ∝ Λ` — *same* single scale.
- **At `θ=0` (Clay's theory) `θ` is not even a knob** — it is a fixed phase, generating no scale. `THETA-PI` (banked): θ=π offers no leverage to θ=0; its "scale" content is `χ_t`/wall-tension, `∝Λ`.
- **Acid test.** Free massless Maxwell has **no finite-action instantons**; its "topological susceptibility" is a **UV contact term** (`∝δ⁴(0)`), not a dynamical `Λ`-scale — so a `χ_t`-anchor is either **UV/regulator**-valued or `∝Λ`; it cannot be the YM-specific, free-failing lever of the charter.
- **Verdict: [NG].** `∝Λ` (SA-1). Also needs `χ_t>0` (instanton density), itself a `∝Λ⁴` condensate lower bound (SA-2).

### (ii) Rank / center / level structure — "a scale `∝N^a Λ`"  [NG]
- **∝Λ?** **Yes.** `N` is a dimensionless integer; `N^aΛ = (N^a)·Λ` — the factor `N^a` is a **pure number**, so by **SA-1** `m₁ ≥ α N^a Λ` reduces to `c ≥ α N^a`, still a pure-number gap inequality. The `[V-here]` witness (F3, Appendix) checks this explicitly: `m₁/(N^a√σ) = (m₁/√σ)/N^a` is unchanged — the `N`-dressing buys nothing.
- **Center `Z_N` / level:** `Z_N` is **discrete**; there is **no 4D Chern–Simons level** (CS is 3D; θ is the 4D analogue and is dimensionless). Large `N` is a **dimensionless control parameter**: masses are `c(N)Λ`, `c(N)=c_∞+O(1/N²)`. Large-`N` could, in principle, let one *compute* the ratio `c_∞` — but that is **proving the planar gap** (open; `PROOF-LARGE-N-RIGOROUS`: rigorous large-`N` gives only the *strong-coupling* gap, no continuum `'t Hooft` construction, and `N=∞→3` has no theorem — the `1/N` series is asymptotic, `O(c/9)` uncontrolled). And even success would yield the **ratio** `c_∞`, not a second scale.
- **Acid test.** A pure-number factor cannot be a scale at all: it fails the acid test *trivially* by producing the same ratio on the free field (`c_free=0`). Vacuous as a lever.
- **Verdict: [NG].** `N` is dimensionless; `N^aΛ` is `∝Λ` (SA-1). The real large-`N` content is "compute the ratio," which is the planar gap, not a scale.

### (iii) Deformation scale — double-trace / center-stabilized / 't Hooft-twist  [NG]
- **What it is.** A double-trace term `δS= (κ/2)∫O²` (`O=trF²`) has `κ` of dimension `[mass]⁴`, so `Λ_dt:=κ^{1/4}` is a **genuine independent dimensionful input** (it escapes SA-1 — it is *not* `∝Λ`). Same for a center-stabilizing deformation coupling or a heavy adjoint mass. This is the **strongest candidate** and deserves the careful kill.
- **Does it survive "connect to pure YM"?** **No.** The deformation defines a **different theory**. To recover pure YM one takes `Λ_dt→0` (deformation IR-irrelevant) or `→∞` (decouples). The anchor lives in the *deformed* theory: `m₁^{def} ≥ f(Λ_dt)`. Then:
  - As `Λ_dt→0`, uniformity demands `m₁^{def} ≥ Δ>0` **uniformly as the deformation is removed** — but the deformation is an **infimum**: every finite-deformation reference gap is an **upper bound** on the undeformed one (`MONOTONE-COUPLING §4`, banked in `SMALLCIRCLE-ADIA §5`). The `Λ_dt→0` limit **is** the gap; the bound degenerates (`f(Λ_dt)→0` unless you already control the undeformed gap).
  - The center-symmetric small-volume phase **is not pure YM**: for pure YM with periodic BC the one-loop GPY potential favors the **center-broken (deconfined)** configuration (`SMALLCIRCLE-ADIA §2`, `[T]` one-loop). Center stabilization (double-trace, adjoint fermions, twist) is a *different theory* — `QCD(adj)`, deformed YM — and "(CP)" = adiabatic continuity is **false** for pure YM (GPY transition) and **fails even in the closest analogue** (four phases in `QCD(adj)`, Cossu–D'Elia). So the deformation→pure-YM bridge is the transport wall.
- **'t Hooft twist:** a **discrete** choice of bundle/flux — **dimensionless**, no scale. **Acid test on the deformation:** the deformed *free* field still has no bound state (`m₁=0`); the mechanism is silent on the free field (fine) but that only means it does not fabricate a gap — it does not supply the anchor in pure YM.
- **Verdict: [NG].** The independent scale is **real but belongs to a different theory**; its removal (or the twist, or the stabilization) cannot be connected to pure YM without the gap. Uniformity/direction (`MONOTONE-COUPLING`, `SMALLCIRCLE-ADIA`) kill it.

### (iv) Finite physical volume L  [NG]
- **What it is.** `L` is a genuine independent dimensionful scale — an **IR regulator we control**.
- **∝Λ?** No (it is external) — so SA-1 does not apply. It must be killed by **(P3) uniformity/direction**, and it is — *structurally*:
  - **Uniformity forbids `1/L`** (`E2-FLIP-NOGO §3-C` explicitly: uniformity kills `1/a,1/L`). Any threshold at `c/L →0` as `L→∞`; a threshold that is a floor after `L→∞` is `∝Λ`.
  - **Direction is backward.** By min–max, `Δ(L)` is **non-increasing** in `L` and `Δ(ℝ⁴)=inf_L Δ(L)` (`SMALLCIRCLE-ADIA §3` Lemma). So `Δ(L) ≥ Δ(ℝ⁴)`: the finite-volume gap sits at the **more-gapped** end (a box raises the lightest mode) and any *computed* finite-`L` value is an **upper** bound on the target. To floor `Δ(ℝ⁴)` you need `Δ(L) ≥ δ` for **all** `L` including `L→∞` — which *is* the gap. A lower bound at a *single* finite `L` is not a bound on the **infimum** `inf_L Δ(L)`; it points the wrong way.
  - **The "minimal-gap volume" `L\*` is itself `∝1/Λ`**: the deconfinement/crossover circumference is `L_c=1/T_c`, `T_c∝Λ` — so locating the minimum gap reintroduces `c=m₁/Λ`.
- **Acid test.** The **KK gap `1/L` is a finite-volume artifact present for the free massless field too** — reading it as a dynamical mass is the classic vacuity trap (`SMALLCIRCLE-ADIA §4`, "spurious-`1/r` trap"). Applied correctly (zero-KK-mode mass), the free field is gapless.
- **Verdict: [NG].** External but **backward** (`Δ(ℝ⁴)=inf_L Δ(L)`) and **non-uniform**; unavailable as a floor. This is the same wall on the volume axis.

### (v) Defect / interface / 't Hooft loop / domain wall  [NG]
- **What it is.** 't Hooft loops (magnetic `Z_N^{(1)}` defects), center vortices, domain walls, interfaces.
- **∝Λ?** **Yes.** The 't Hooft-loop (vortex) tension `σ̃ = q·Λ²`, the domain-wall tension (at `θ=π`) `∝Λ³`, interfaces likewise — all `∝` powers of `Λ` (SA-1). At `θ=0` (Clay) there is a **unique vacuum** → **no domain wall** (walls exist only at the `θ=π` level crossing, which Clay does not concern; `THETA-PI`).
- **Sector mismatch (independent second kill).** The 't Hooft loop / center vortex is **1-form** (`Z_N^{(1)}`); the residue `f₁²=|⟨0|¼trF²|G⟩|²` and the `0⁺⁺` mass are **0-form** local data. `trF²` is N-ality-0 and **provably blind** to the 1-form symmetry (`PROOF-1FORM-ANOMALY`): a point cannot carry a 1-form charge. So even a proven `σ̃>0` lands in the wrong channel; `σ̃>0 ⇒ m_{0⁺⁺}>0` only through the **measured** `m_{0⁺⁺}≈3.9√σ` (`BS-THRESHOLD §4`, `PROOF-CENTER-VORTEX`), not as a theorem.
- **P1 also fails non-circularly:** at `θ=0` the electric `Z_N^{(1)}` carries **no 't Hooft anomaly**, so the anomaly **does not force** area law (`PROOF-1FORM-ANOMALY`, `PROOF-CENTER-VORTEX §4`); vortex condensation at `β=2.5` is the confinement problem re-expressed (proven only `β<1/48`/`β<0.706`).
- **Acid test.** A defect tension is a confinement object; the free field has no `σ̃` (`σ̃=0`) — but the mechanism is silent, and it is `∝Λ` and cross-sector anyway.
- **Verdict: [NG]**, doubly: `∝Λ` (SA-1) **and** 1-form vs 0-form sector mismatch.

### (vi) Anomaly-induced / trace-anomaly scale  [NG]
- **What it is.** The trace anomaly `T^μ_μ=(β(g)/2g)·¼trF²` (`θ`-type dilaton `f_G`); 't Hooft anomalies; the conformal anomaly.
- **∝Λ?** **Yes — it *is* `Λ`.** `⟨T^μ_μ⟩=(β/2g)⟨¼trF²⟩ ∝ Λ⁴`; the trace anomaly is the **mechanism of dimensional transmutation**, i.e. the *definition* of the single scale, **not** a second one. The "dilaton decay constant" `f_G` is, by the exact trace identity, **the residue rescaled**: `f₁²=f_G² m_G⁴/(β/2g)²`, so `f_G≥c ⇔ f₁²≥r₀` — verbatim the gap/residue (`DILATON-FLOOR`: "the dilaton decay constant is the unprotected order parameter in dilaton clothing"). `'t Hooft` anomalies are **dimensionless/discrete**, and the 1-form anomaly is **line-sector** (blind to `trF²`; `PROOF-1FORM-ANOMALY`).
- **Acid test.** The free massless field is **conformal** (`β=0`) ⇒ `⟨T^μ_μ⟩=0`. So a `⟨T⟩^{1/4}`-anchor gives `m₁≥0` on the free field (correct number, **certifies nothing** — vacuous), and on YM gives `m₁ ≥ α⟨T⟩^{1/4}` = `c ≥ α·(const)` — the gap ratio (SA-1, circular).
- **Verdict: [NG].** The anomaly scale **is** `Λ` (SA-1); the anomaly residue **is** the gap residue (`DILATON-FLOOR`).

---

## §3 Assembly — the dichotomy, six-for-six

**Theorem (second-scale rigidity, for the named family). [NG + T]**
Let `Q` range over the dimensionful quantities of pure YM at `T=0, θ=0, ∞`-volume, `N` fixed. Then `Q` is one of:
- **(A) `Q∝Λ`** (candidate classes i, ii, v, vi, and any internal observable such as `√σ, m_{G}, T_c, χ_t^{1/4}, ⟨T^μ_μ⟩^{1/4}, σ̃`). By **Lemma SA-1** any homogeneous anchor `m₁≥αQ^{1/d}` is the pure-number inequality `c ≥ αq^{1/d}`, i.e. it presupposes `c=m₁/Λ` — **circular**; by **SA-2** even `Q>0` is a `∝Λ^d` condensate bound. **No anchor.**
- **(B) `Q` external** (candidate classes iii, iv; also `T>0` if one leaves `T=0`). `Q` escapes SA but is **absent from the target theory** and fails **(P3)**: `L`,`T` are **backward-monotone** (`Δ(reg) ≥ Δ(ℝ⁴)`, `Δ(ℝ⁴)=inf Δ`) and **non-uniform** (`1/L`,`1/T→0`); deformations define a **different theory** whose pure-YM limit (`κ→0`, `r→∞`) **is** the transport/gap (`MONOTONE-COUPLING`, `SMALLCIRCLE-ADIA`); twists are **dimensionless**.
There is **no third horn** — the inventory §0.3 is exhaustive (`[T]`: one coupling + discrete data + regulators). Hence **no genuine independent second scale furnishes an `(L,a)`-uniform, direction-correct, acid-clean anchor**: the kinematic-E₂ escape cannot be re-run.

*Why "six accidents" is really "one dichotomy":* the split is exactly **(A) `∝Λ` vs (B) non-`∝Λ`**, and each horn has one structural obstruction (`SA`-circularity vs `uniformity/direction/connection`). The six candidates are *instances*, not six separate failures — which is what makes the verdict structural (`[T]`-dimensional + `[NG]`-family) rather than a list.

**Witness `[V-here]` (Appendix).** `e2_second_scale_witness.py` verifies on the standard `[R]` ratio table: (F1) every `∝Λ` anchor reduces to the pure ratio `m₁/Q` (= "best `α`") — no absolute floor; (F2) the anchor inequality is invariant under `Λ→λΛ` (it cannot see the scale); (F3) `N^aΛ` is still `∝Λ`; and the converse — a **known external** `Q_ext` (not `∝Λ`) **does** give a floor, confirming the dichotomy's (B) branch is the only door, and it is the door §3(B) closes. All arithmetic, **not** new QCD evidence; values `[R]`.

---

## §4 What would flip this — precise flip condition, and the one honest live corner

**Flip condition (stated so it can be refuted).** Exhibit a quantity `Q` with **(P1)** `Q>0` provable non-circularly; **(P2)** `Q` **not** `∝Λ` (so SA-1 does not apply); **(P3)** an explicit `f` with `m₁≥f(Q)` **uniform in `(L,a)`**, and **acid-clean** (`f(Q)=0` on the free massless field). Equivalently: a **`Λ₀` provably `≢ multiple-of-Λ`** with a **kinematic threshold fixed by `Λ₀`** *and* a **theorem `m₁≥f(Λ₀)`**. Per §0.3/§3, pure YM at `T=0,θ=0,∞` has **no such `Λ₀`**; the only non-`∝Λ` scales are external and fail (P3). **If such `Q` existed, Theorem S is beaten and `E2-FLIP-KINEMATIC` can be re-run.**

**The one genuine *non-scale* corner (do not under-sell, but do not mis-file).** Same-channel **von Neumann–Wigner** isolation (`G2D-ISOLATION-I2 §3`, `E2-FLIP-NOGO` L-1) genuinely breaks the *pointwise* isolation circularity **without** a second scale and **without** assuming the gap: two `0⁺⁺` levels cannot cross, so `m₁<m_{0⁺⁺*}` pointwise wherever both are bound. It is **not** a "second scale" (it is a **level-ordering**, not a dimensionful `Q`), it is **sub-uniform** (strict inequality, never a uniform `κ`), and it **dissolves at large `β`** (Kato isolation of the excited `0⁺⁺*` fails). Uniformizing it is the ratio `ρ=m_{0⁺⁺*}/m₁≥1+c` — verbatim `H*`/`SSC`, the scale-invariant face of the same wall. So it is real but does **not** rescue the second-scale escape.

---

## §5 Honest boundary

- **[T] what is rigorous:** the dimensional inventory (§0.3: one dimensionless coupling + discrete data + regulators ⇒ one scale up to pure numbers) and Lemma SA-1 (homogeneity ⇒ a `∝Λ` anchor is a pure-number/gap inequality; scale-invariant). These are dimensional algebra, not dynamics.
- **[NG] for a named family, not `∀possible`:** "anchor built from `Q∝Λ` **or** from an external scale whose pure-YM limit is taken." The vNW pointwise corner is **not** in this family and is a real (sub-uniform) partial; it is not a scale.
- **[J]/[O] the dictionary step (SA-2):** "`Q>0` is a `∝Λ`-condensate lower bound = one face of the key" is the charter's dictionary (`CLAUDE.md §5`), a judgment/known-equivalence, **not** a fresh theorem.
- **Empirical `[R]`:** the lattice ratios (`m_{0⁺⁺}/√σ≈3.9`, `χ_t^{1/4}≈0.18 GeV`, `T_c≈0.27 GeV`, `T_c/√σ≈0.61`) are standard/as-reported, **not re-verified** this session; they illustrate single-scale constancy, they are not evidence *for* rigidity (rigidity is `[T]`-dimensional). `√σ` is the scheme-independent "one scale"; `Λ_MSbar` is a scheme-dependent label of the same scale (ratios to it carry a scheme conversion) — this does not affect any conclusion.
- **This does NOT prove Clay** and does NOT prove no proof exists (`KEY-EXISTS §5`, `∀class≠∀possible`). It closes the **"second-scale anchor"** shape for the named family, six-for-six with structural reasons, and thereby **confirms the dimensional rigidity** the `E2-FLIP-KINEMATIC` escape needed to beat. Clay stays open.

---

## Appendix — `[V-here]` witness (`research/unconditional/e2_second_scale_witness.py`)

Pure-Python, `numpy`-free, single process (memory rule). Output (essence):

```
(F1) anchor Q, standard [R] ratios (GeV): m1=m_0++=1.710, sqrt(sigma)=0.440
     anchor Q        Q      Q/sqrt(sig)   m1/Q (= best alpha)
     sqrt(sigma)   0.440       1.000            3.886
     m_2++         2.400       5.455            0.713
     chi_t^(1/4)   0.180       0.409            9.500
     T_c           0.270       0.614            6.333
     Lambda_MSbar  0.250       0.568            6.840
     -> best alpha is the unknown ratio; proving m1>=alpha Q is proving
        alpha <= m1/Q = the gap.  No absolute floor.
(F2) Lam -> lam*Lam (lam=0.1,1,10): 'm1 >= 0.1 Q' truth INVARIANT;
     ratio m1/Q unchanged (e.g. 3.886).  A scale-covariant anchor cannot
     see the absolute scale -> it restates a ratio, cannot produce a floor.
(F3) Q = N^a*sqrt(sigma):  m1/Q = (m1/sqrt(sigma))/N^a  (pure)  -> no anchor.
CONVERSE: a KNOWN external Q_ext (NOT prop. to Lam) DOES give m1 >= alpha Q_ext
          (a genuine floor) -> confirms the ONLY door is an external scale,
          which the document closes (backward / different-theory / non-uniform).
ACID: free field m1=0 -> every anchor gives 'm1>=0' (vacuous, charter-consistent).
```

These are **witnesses of Lemma SA** (dimensional algebra), **not** new evidence about QCD. The structural failures of external scales are argued in §2–§3 (monotonicity direction; different-theory/degenerate-limit; dimensionless twist). Run: `py -3.12 -I e2_second_scale_witness.py`.

---

## §6 References

**Internal (absolute paths; the banked verdicts audited here):**
- `D:\Dev\math\CLAUDE.md` (§1 build-don't-stop, acid test; §2 labels; §5 single-scale / condensate = key); `D:\Dev\math\docs\BRIEF-FOR-MODEL.md` (§0 task types; §2 five constraints; §3 death map).
- `research\unconditional\E2-FLIP-NOGO.md` (§3 Theorem S — the barrier being attacked; §5 L-1 vNW, L-2 `H1`); `E2-FLIP-KINEMATIC.md` (§0.3 SSC; §2 "no second scale"; §6 flip condition); `E2-FLIP-DYNAMICAL.md` (sibling: Savvidy lever irreducible; §4 exhaustion of families); `E2-REGRESS-CIRCULARITY.md` (§0.4 SSC, §6 flip condition).
- `G2D-ISOLATION-I2.md` (§2 SSC — "to escape, one needs a threshold set by a different scale … absent in pure YM"; §3 vNW pointwise; §5 ratio = `H*`).
- `THETA-PI.md` (θ=π: no leverage to θ=0; θ scale = `χ_t`/wall tension `∝Λ`); `DILATON-FLOOR.md` (`f_G` = residue in disguise; trace anomaly = `Λ`; `f₁²=f_G²m_G⁴/(β/2g)²`).
- `PROOF-1FORM-ANOMALY.md` (1-form anomaly line-sector, `trF²` blind; θ=0 no anomaly); `PROOF-CENTER-VORTEX.md` ('t Hooft loop/vortex `∝Λ²`, 1-form, `m_{0⁺⁺}∝√σ` measured); `BS-THRESHOLD.md` (§4 `√σ` 1-form mismatch; `m_{0⁺⁺}≈3.9√σ`).
- `SMALLCIRCLE-ADIA.md` (§2 small-r pure YM center-broken/GPY; §3 backward monotonicity `Δ(ℝ⁴)=inf_r Δ(r)`; §4 spurious-`1/r`; §5 deformation/twist ≠ pure YM, adiabatic continuity = (CP) false); `PROOF-LARGE-N-RIGOROUS.md` (large-`N` = strong-coupling gap only; `N→∞→3` no theorem); `MONOTONE-COUPLING.md` (deformation infimum = upper bounds).
- `KEY-EXISTS-OR-IMPOSSIBLE.md` (§5 `∀class≠∀possible`); `SCHUR-CONVEX.md` (positivity Schur-concave). `e2_second_scale_witness.py` (`[V-here]`).

**External — `[R]` (standard/as-reported; not re-verified this session):**
- Morningstar–Peardon, *Phys. Rev. D* **60** (1999) 034509 (glueball/`√σ` ratios); Shen–Zhu–Zhu, *CMP* **400** (2023) 805, arXiv:2204.12737 (strong-coupling lattice gap, `β<1/48`); Gross–Pisarski–Yaffe, *Rev. Mod. Phys.* **53** (1981) 43 (GPY small-circle center breaking).
- 't Hooft, *Nucl. Phys. B* **138** (1978) 1; *Nucl. Phys. B* **153** (1979) 141 (twisted sectors / 't Hooft loop); Tomboulis–Yaffe, *CMP* **100** (1985) 313 (vortex free energy ⇒ `σ>0`); Gaiotto–Kapustin–Seiberg–Willett, arXiv:1412.5148 (generalized symmetries).
- Shifman–Ünsal, *Phys. Rev. D* **78** (2008) 065004, arXiv:0802.1232 (double-trace); Ünsal–Yaffe, *Phys. Rev. D* **78** (2008) 065035, arXiv:0803.0344 (center-stabilized / volume independence); Ünsal, *Phys. Rev. D* **80** (2009) 065001, arXiv:0709.3269 (magnetic bions); Cossu–D'Elia, *JHEP* **0907** (2009) 048 / arXiv:0904.1353 (four phases, against smooth continuity).
- Standard lattice value `χ_t^{1/4}≈0.18 GeV` `[R]` (as banked; no specific arXiv asserted here); Jaffe–Witten, *Quantum Yang–Mills Theory*, Clay (2000).

---

*Status: **RIGIDITY-CONFIRMED** — **[NG] for the named family** "second-scale anchor" + **[T] dimensional inventory + Lemma SA**. No genuine independent second dimensionful scale exists in pure `T=0,θ=0,∞`-volume YM: every candidate is either `∝Λ` (Lemma SA: anchor = pure-number gap inequality, circular) or external and failing L-uniformity / monotonicity-direction / connection-to-pure-YM. The kinematic-E₂ escape cannot be re-run; the one real non-scale corner (vNW pointwise) is sub-uniform and `≡ H*` when uniformized. Witness `[V-here]` verifies the dimensional algebra on `[R]` ratios. NOT a proof of Clay; NOT a global impossibility theorem (`∀class ≠ ∀possible`). **Clay stays open.***
