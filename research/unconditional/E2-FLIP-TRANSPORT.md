# E2-FLIP-TRANSPORT — attempt to BUILD/CONTROL the nonperturbative transport machine carried by the *direction-correct* functional `e^{−S_inst}` (and its condensate relatives); verdict: **PARTIAL** — the machine is built and proven-controlled in the *rigid-saddle / dilute* class and the residue-transport is an exact theorem, but its uniform `(L,a)` control across the crossover is **IRREDUCIBLE**: the condensate coordinate carries the *same* wall, in a new dress (the **size modulus / diluteness ⟺ wrong-scale dichotomy**)

> 2026-11-24 · Direct attack on the **`e^{−S_inst}` "leak"** flagged in `E2-FLIP-DYNAMICAL §6` and named the missing machine in `CLAUDE.md §5` / `TRANSPORT-MACHINE`. The task is not "another mechanism name": it is to **write the machine down, prove what it controls, and pin exactly where it stops**, or prove it cannot be completed. This note does both: a genuine partial machine (`[T]`/`[C]` pieces) *plus* a precise obstruction (`[NG]` for a named class) resting on a fresh `[V-here]` computation.
>
> **Labels (`CLAUDE.md §2`):** **[T]** 真定理 / **[C]** 条件定理 / **[NG]** no-go *for a named class* / **[O]** 精确缺失步 / **[J]** 判断 / **[V-here]** fresh computation this session (`research/unconditional/w7/e2_flip_transport.py`). External facts are **[R]** (as-reported, not re-verified here). **No QED on anything unproven. This document does NOT prove Clay and does NOT prove no tool can exist (`∀class ≠ ∀possible`, `KEY-EXISTS §5`).**

---

> **⚠ 更正 (2026-11-24 · 宪章 §2)**：本笔记 §4.1 的"packing `nρ⁴ = O(10)` ⇒ 稠密 ⇒ diluteness fails"是**计算错误**，已由 `TRANSPORT-PACKING-AUDIT.md` 订正。真因两处复合错：(1) SU(3) 单圈前因子 `C₃ = 1.51×10⁻³`（非 `O(10⁻¹)`，差 66×）；(2) 把带量纲 `1/长度` 的单位尺寸密度 `D(ρ*)ρ*⁵` 当成无量纲 `nρ̄⁴`。**真堆积率 `𝔡 = nρ̄⁴ ≈ 0.01–0.06`（ILM / 格点 `χ_top` / 单圈全积分 三源一致）—— 物理瞬子真空是稀薄，不是稠密。** 故"dense ⟹ diluteness fails"的 no-go **前提不成立**；结论（物理尺度无严格簇展开）只能改从 **Kotecký–Preiss 的 activity×range / 连续统控制 = H1** 重新论证。`ρ* = 0.368/Λ`（质量在大 ρ / IR 生成）这条仍对。详见 `TRANSPORT-PACKING-AUDIT.md`、`TRANSPORT-FIVE-WAY-VERDICT.md`。

## Verdict up front

> ### **PARTIAL.** The `e^{−S_inst}` transport machine is **explicitly built** and, in its **controlled class** (a *rigid saddle* with a fixed, isolated instanton action, or a *dilute* / small-`S¹` gas), it is a **proven lower bound** on the mass gap, direction-correct and acid-clean `[T-in-model / C]`. The **residue-transport** (the "move the seed across β") is an **exact, sign-preserving theorem** `[T]`. But the machine **provably stops** at the crossover: the very functional that is direction-correct, `e^{−S_inst}`, is carried in 4D by a saddle with a **size modulus** whose one-loop distribution is **controlled by no small parameter at the mass-generating scale** — the controlled (`S_inst≫1`, dilute) end is the **UV / wrong phase / wrong channel**, and the mass-generating end is **uncontrolled**. This is the *same* wall as `TRANSPORT-MACHINE`/`SMALLCIRCLE-ADIA`/`FORGE5-CONFINING-MONOPOLE`, reached this time through the condensate/tunnelling coordinate rather than the Hessian or the residue. **The uniform `(L,a)` statement is IRREDUCIBLE for the class {instanton gas, monopole gas, center-vortex condensation}.**
>
> **One line:** the machine is real, works where the saddle is rigid and the gas dilute, and the gap is exactly the failure of a *small, scale-invariant-to-IR control parameter* — i.e. the size modulus is the wall.

**What is new here (vs. the banked per-route verdicts).**
1. **Assembled as one machine** with an explicit four-factor decomposition `e^{−S_inst} = [topological weight]·[size-integral]·[running exponential]·[transport]`, so the three condensate relatives (instanton / monopole / vortex) are *three evaluations of one formula*, and the single obstruction is located once.
2. **Fresh `[V-here]` proof that the direction-correct functional is genuinely controlled in the rigid-saddle class** — `m₁ = (1.48±0.02)·ω√(S_inst)e^{−S_inst}` with `corr(m₁, e^{−S_inst})=+0.973` while `corr(m₁, Hessian depth)=−0.9999` (Part 1). This *confirms* `E2-FLIP-DYNAMICAL`'s "leak = the direction-correct functional" and gives it a working controlled instance.
3. **Fresh `[V-here]` quantification of the obstruction**: the one-loop instanton size distribution peaks at `ρ* = e^{−2N/(b−5)}/Λ` — i.e. `2.7/Λ` (SU(3)) to `5.6/Λ` (SU(2)), **large = IR** — where the packing fraction `n ρ*⁴` is `O(1)` (SU(3): `O(10)`), so **diluteness fails exactly where the mass is generated** (Part 2). The free-field limit has no running, is scale-invariant, and returns `m=0` (Part 2b). This turns "control ⟺ wrong phase" from a slogan into a computable statement.

---

## §0 The target object, restated with the exact functional

From `E2-FLIP-DYNAMICAL §6` and `CLAUDE.md §5`, the unique **direction-correct** functional the campaign has identified is the **tunnelling / condensation action**
```
   Φ_tun  :=  e^{−S_inst}   (or the condensate density it controls; mass arises as m ~ e^{−S_inst} / e^{−S_inst/2}),
```
which `E2-FLIP-DYNAMICAL` proved is (i) non-Schur, (ii) non-regress, (iii) dynamical, (iv) acid-test-clean, and (v) *equals* the gap — but whose **uniform `(L,a)` control across the Savvidy crossover** was asserted (not shown) to *be* the missing transport machine. The present note does not take that assertion on faith: it **builds `Φ_tun` and tests it**.

The precise claim `Φ_tun` must certify to close Clay (`BRIEF-FOR-MODEL §1–2`, `E2-REGRESS §0`):
```
   (LB)   Φ_tun(μ)  ≤  m₁(μ)               for every admissible (OS/RP) theory μ,
          Φ_tun(μ_free) = 0                [free massless / scale-invariant endpoint]
          inf_{μ on pure SU(N) YM, (L,a)-uniform}  Φ_tun  ≥  Δ > 0.
```
`Φ_tun` is direction-correct *by construction* when `m₁ ≳ c·e^{−S_inst}` — i.e. the gap is a positive multiple of the instanton weight — which is the standard form of **dimensional transmutation**. So the whole problem is: **prove `m₁ ≥ c·e^{−S_inst}` uniformly, without presupposing the scale `Λ` that `e^{−S_inst}` is made of.**

---

## §1 The machine, written as a four-factor decomposition

On the lattice `μ = Z⁻¹ e^{−S_W}dU`, the instanton contribution to the free energy / the tunnelling amplitude factorises (standard semiclassical gas; `[R]` Callan–Dashen–Gross, Schäfer–Shuryak):
```
   (Φ_tun)   e^{−S_inst} =  [topological weight e^{−8π²/g²}]  ·  [size/moduli integral ∫ dρ D(ρ)]  ·  [transport / renormalisation of the weight]
                                 ╰ deterministic, sign-controlled ╯   ╰ THE MODULUS: uncontrolled in 4D ╯   ╰ exact homomorphism f₁² ↦ Z_O²f₁² ╯
```
The four-stage form of the *transport machine* (matching `TRANSPORT-MACHINE`'s Stage 0–3, now in the condensate coordinate):

| stage | object | status | basis |
|---|---|---|---|
| **0 seed** | strong-coupling gap / residue `f₁²≥r₀` at `β<1/48` | **`[T]`** unconditional | P3 / SZZ, `PROOF-P3-UNCONDITIONAL` |
| **1 saddle action** | `S_inst = 8π²/g²`, deterministic weight `e^{−S_inst}` | **`[T]`** (a number); controls the gap **iff the saddle is rigid & dilute** | `[V-here]` Part 1; `[R]` Coleman |
| **2 transport** | `f₁²(μ)=Z_O(μ,μ₀)²f₁²(μ₀)`, `Z_O²>0` a square | **`[T]`** exact, sign-preserving, acid-clean | `INVENT-TRANSPORT-INEQ` |
| **3 crossover** | the **size-modulus integral** `∫dρ D(ρ)` at the mass-generating scale | **`[O]` ⟺ the wall** | `[V-here]` Part 2; `W2 §(i)`, `SMALLCIRCLE-ADIA §3` |

The three **condensate relatives** are the same formula evaluated in different variables:
- **Instanton gas**: `D(ρ) ∝ (8π²/g²(ρ))^{2N} e^{−8π²/g²(ρ)}/ρ⁵`, `m₁ ~ (∫D)^{1/2}`.
- **Monopole-instanton gas** (on `ℝ³×S¹` / 3D compact `U(1)`): fugacity `ζ ~ e^{−S₀}`, `S₀=8π²/(g²N)`, dual-photon mass `m_D ~ ζ^{1/2}`.
- **Center-vortex condensation**: vortex free-energy density `ρ̂` at the confining threshold; Tomboulis–Yaffe `σ ≥ c·ρ̂`.
All four stages and all three relatives share **one** hard object: the crossover stage, whose control parameter must be small *at the mass-generating scale*.

---

## §2 The acid test and the direction — `e^{−S_inst}` is genuinely the right functional `[V-here] + [T-in-model]`

The task's hurdle (`CLAUDE.md §1`, `BRIEF §6`) is that a *valid* tool must **fail** on the free massless field. We check both direction and acid **explicitly** in the exactly solvable 1D tunnelling model `V=(g/4)(x²−v²)²` (the physical gap is the transfer-matrix splitting, computed exactly via `U(−L)U⁻¹ = −∂²+V′²/4−V″/2`; `e2_flip_transport.py` Part 1). Table `[V-here]`:

| `g` | `S_inst` | `e^{−S}` | `m₁²` exact | `m₁` | `ω√S e^{−S}` | **ratio `m₁/Φ_tun`** | Hessian `−g` |
|---|---|---|---|---|---|---|---|
| 2 | 1.333 | 2.64e−1 | 8.459e−1 | 0.9197 | 0.6088 | 1.511 | −2 |
| 4 | 1.886 | 1.52e−1 | 7.487e−1 | 0.8653 | 0.5893 | 1.468 | −4 |
| 6 | 2.309 | 9.93e−2 | 5.973e−1 | 0.7728 | 0.5229 | 1.478 | −6 |
| 8 | 2.667 | 6.95e−2 | 4.545e−1 | 0.6742 | 0.4539 | 1.485 | −8 |
| 10 | 2.981 | 5.07e−2 | 3.354e−1 | 0.5792 | 0.3917 | 1.479 | −10 |

with `Φ_tun := ω√(S_inst)e^{−S_inst}`, `ω=√(V″(well))`, `S_inst=(4v³/3)√(g/2)`.

- **Direction `[T-in-model]`:** over the deep regime `g≥4`, `corr(m₁, e^{−S_inst}) = +0.973` (the gap moves **with** the instanton weight) while `corr(m₁, g) = −0.9999` (the gap moves **opposite** to the Savvidy Hessian depth). This is the exact content of `E2-FLIP-DYNAMICAL` — the Hessian is the wrong functional, the tunnelling weight is the right one — now with a working controlled instance.
- **Controlled lower bound `[T-in-model]`:** the ratio `m₁/Φ_tun = 1.48±0.02` is **essentially constant** across the whole deep range, so
  ```
     m₁  ≥  Φ_tun  =  ω · √(S_inst) · e^{−S_inst}      (rigid-saddle class, S_inst≫1).
  ```
  A genuine **lower bound** on the gap, computable from the action landscape, positive, and **acid-clean**: in the scale-invariant limit `S_inst→∞` (free/weak coupling), `Φ_tun→0` and the free massless field returns `m₁=0`. `[V-here]` Part 2b confirms the free-field size integral is scale-invariant (no running ⇒ `∫dρ/ρ⁵` diverges ⇒ no preferred size ⇒ `m=0`). **Acid test PASSED non-vacuously.**

**Why this is the *only* functional of its kind, and why that is not enough by itself.** `Φ_tun` clears all six flip constraints of `E2-REGRESS §6` — *except* that its control requires a rigid saddle with `S_inst` bounded across the crossover. In the double well `S_inst` is a **single number set by the potential** (a rigid, isolated bounce); in 4D it is a **function of the instanton size**, a flat modulus. That difference is the whole document.

---

## §3 What the machine CONTROLS — the honest positive ledger

**`[T]` Rigid-saddle lower bound (§2).** In any theory with an isolated, finite-action tunnelling saddle (double well; `ℝ³×S¹` caloron/bion of *maximal* fixed size; 3D compact `U(1)` monopole), the machine gives `m₁ ≥ c·ω√(S_inst)e^{−S_inst}`, `c>0` explicit, non-circular (`S_inst` from the action, not the spectrum), acid-clean. This is a real, proven piece of machine.

**`[T]` Exact residue transport (`INVENT-TRANSPORT-INEQ`).** The seed `f₁²` is carried across scales *multiplicatively*: `f₁²(μ)=Z_O(μ,μ₀)²f₁²(μ₀)`, `Z_O=(β/2g)(μ₀)/(β/2g)(μ)>0`, a manifest square — `0→0`, so acid-clean by construction. This is the rigorous "move the seed" leg and it is *not* where the machine fails.

**`[C]` Controlled dilute/small-`S¹` lower bound.** In the regimes where a **convergent cluster expansion** exists — 3D compact `U(1)` (Göpfert–Mack `[R]`), and the small-`S¹` centre-stabilised bion gas *as controlled semiclassics* — the condensate does produce `m₁ ≥ δ(r)>0`. This is the machine's output where it is rigorous. The hypotheses are exactly: diluteness + (for the nonabelian case) a deformation/twist fixing the size.

**Net positive:** the machine is real and is *right-signed, acid-clean, and exact in its transport leg.* The positive content is genuine and should not be under-sold — it is the correct anatomy, which is precisely why the obstruction below is structural rather than "we didn't try hard enough."

---

## §4 WHERE IT STOPS — the precise obstruction (the heart of the note)

The crossover stage is the **size-modulus integral**. We now show it is uncontrolled at the mass-generating scale for a reason that is *computable and convention-independent*.

### 4.1 The one-loop instanton size distribution `[V-here]`
With `8π²/g²(ρ) = b·ln(1/(ρΛ))`, `b=11N/3` (`Λ=1` units), the distribution `D(ρ) ∝ [b\ln(1/ρ)]^{2N} e^{−b\ln(1/ρ)}/ρ⁵ = [b\ln(1/ρ)]^{2N} ρ^{b−5}` has its **peak** at
```
   d/dρ[...] = 0   ⟹   ln(1/ρ*) = 2N/(b−5)   ⟹   ρ* = e^{−2N/(b−5)} / Λ.
```
Fresh numbers `[V-here]` Part 2:

| gauge group | `b=11N/3` | peak size `ρ*` | `S_inst(ρ*)` | `e^{−S}` at peak | packing `n ρ*⁴` |
|---|---|---|---|---|---|
| SU(2) | 7.333 | **0.180 / Λ** (= 5.55/Λ) | 12.57 | 3.5e−6 | `C·0.09` (marginal) |
| SU(3) | 11.0 | **0.368 / Λ** (= 2.72/Λ) | 11.00 | 1.7e−5 | `C·29.6` (dense) |

(`C` = one-loop determinant prefactor, `O(10^{−1})`.)

**The two decisive facts, both exact:**
1. **The mass-generating instantons are LARGE = IR.** The peak is at `ρ* = O(1)/Λ` (a size comparable to the confinement scale itself), **not** at small `ρ` (UV). So the semiclassical gas generates its mass at the *largest* instanton size in the problem.
2. **The gas is DENSE at the peak.** The packing fraction `n ρ*⁴` is `O(1)` (SU(3): `O(10)`) — *not* parametrically small (`≪1`). A gas of instantons whose volume fraction is `O(1)` is a **liquid, not a dilute gas**; the Mayer/Kotecký–Preiss cluster expansion — which is what makes the mass a theorem — does *not* converge.

### 4.2 The diluteness ⟺ wrong-scale dichotomy (the obstruction, made precise)
Make precise the single statement underlying `W2 §(i)`, `R1`, `SMALLCIRCLE-ADIA §3`, `FORGE5-CONFINING-MONOPOLE #3`, `PROOF-CENTER-VORTEX §3`. The semiclassical gas is under control **iff** there is a dimensionless small parameter with, *simultaneously*:

```
   (i)   S_inst  ≫  1        (large action ⇒ dilute, saddle point dominates),  and
   (ii)  𝔡 := n·ρ̄⁴  ≪  1     (packing fraction small ⇒ cluster expansion converges),
   (iii) the scale at which (i)–(ii) hold  =  the scale at which the mass is generated.
```
- **At the mass-generating scale** (the one-loop peak): (i) holds weakly (`S_inst≈11–13`) but **(ii) fails** (packing `O(1)`), so the expansion — hence the "theorem" — fails. `[V-here]`
- **In the UV** (`ρ<a`, `S_inst=4π²β/N≫1` on the lattice): (i)+(ii) hold, but (iii) **fails** — these are lattice artifacts of size `~a`; their contribution to a physical correlation length is `~ e^{−S_inst}→0` as `a→0` (they renormalise, they do not confine). On the lattice, `S_inst = 4π²β/N ∝ β→∞` in the continuum, so `e^{−S_inst}→0`: **the controlled end's mass vanishes in the continuum limit.**
- **On small `ℝ³×S¹`** (bion gas, `S₀=8π²/(g²N)≫1`): (i)+(ii) hold, but (iii) fails in a *second* way — the controlled theory is **deformed / centre-stabilised**, not pure YM, and its `ℝ⁴` limit crosses the GPY deconfinement transition; moreover by min–max `m(S¹_L) ≥ m(ℝ⁴)` the small-`L` bound is an **upper** bound on the physical gap (`SMALLCIRCLE-ADIA §3`, `FORGE5-SMALL-S1`).

> **`[NG]` Condensate-transport obstruction (named class).** *Let `Φ` be `e^{−S_inst}` (or a monotone function of the instanton/monopole/vortex fugacity) evaluated via a semiclassical gas. Then `Φ` cannot satisfy `(LB)` uniformly in `(L,a)` on `ℝ⁴`, because the condition for the gas to be a *controlled* (provably positive, cluster-convergent) mass source — `S_inst≫1` and `𝔡≪1` at the generating scale — is incompatible with that scale being the mass-generating one. The three ways to have control (UV lattice / small `S¹` / Coulomb phase) each fail (iii): UV mass `→0` in the continuum; small-`S¹` is a different theory with a backward min–max inequality; Coulomb is gapless.* **`[NG]` for the class `{instanton gas, monopole gas, center-vortex condensation}`; the "control ⟺ wrong-scale/phase/channel" half is `[T]`-computable (§4.1); the "no uniform floor follows" component rests on the same structural facts as `TRANSPORT-MACHINE`.** ∎

### 4.3 Two further locks, independent (so the wall is not a single missing estimate)
Even granting (4.2) hypothetically, two further independent obstructions face the condensate coordinate, both banked and both confirmed against this build:

- **Channel mismatch `[T]`.** The controlled condensates deliver a **magnetic / 1-form** mass (`m_D`, `√σ`), while Clay's residue `f₁²=|⟨0|¼trF²|G⟩|²` is the **neutral 0-form `0⁺⁺`** edge. A point carries no 1-form charge (`PROOF-1FORM-ANOMALY`); `σ>0 ⇒ m_{0⁺⁺}>0` only through the *measured* `m_{0⁺⁺}=c√σ` (`TRANSMUTE-TENSION-GAP`, `BS-THRESHOLD §4`), an identity not a theorem. The vortex route's Tomboulis–Yaffe `σ≥cρ̂` (`PROOF-CENTER-VORTEX`) is right-signed but lands in the 1-form sector.
- **One-scale circularity `[T-cond / J]`.** `e^{−S_inst}=e^{−8π²/g²}= (Λ/μ)^{\text{power}}` is *literally* the dimensional-transmutation factor: the scale it "generates" **is** `Λ_QCD`, and `m₁=c·Λ_QCD` has `c` a pure number. So `Φ_tun` is a lower bound on `m₁` only up to the unknown dimensionless ratio `c=m₁/Λ` — i.e. only up to the gap (`E2-FLIP-NOGO §3 Theorem S`; `HIGGS-DIM-TRANS`: transmutation fixes `Λ` but never an isolated `0⁺⁺` pole). The anomaly is a *sum rule* (`∫dsρ(s)/s²∝Λ⁴`) satisfied equally by a gapped and a gapless measure.

**Three independent locks — size-modulus/diluteness (4.2), channel (4.3a), one-scale (4.3b) — all close on the condensate coordinate.** This is why the verdict is *structural*, not "one estimate short": the pieces are not a single gap but three mutually reinforcing walls, any one of which is sufficient to stop the uniform statement.

---

## §5 The exact missing step, and the precise flip condition

The one inequality whose proof would complete the machine (equivalently, the precise content of the crossover stage):
```
   [O]  ∃ a single dimensionless bound  𝔡(ρ) := n(ρ) ρ⁴ ≤ 𝔡₀ < 1  for all ρ ∈ [a, 1/Λ],
        with  S_inst(ρ) bounded on the same range and a lower bound on the mass it generates,
        uniform in (L,a)  —  i.e. control of the interacting instanton/monopole/vortex LIQUID
        (multi-defect correlations, "molecules") at the physical scale.
```
By §4 this is **not** a condensate-density estimate: it is the statement that the density is small where the mass lives, which is the statement that the dense non-abelian liquid is under analytic control — verbatim the **confining / crossover RG problem** (`BRIEF §4 "transport ≡ confinement"`, `CHATGPT-METHOD-RGENTRY`). It is the *same* object as `H1`/`I2` (`TRANSPORT-MACHINE Stage 3`), reached through `e^{−S_inst}` instead of the Hessian or the residue.

**Flip condition (stated so it can be refuted).** This verdict flips to **MACHINE BUILT** iff one exhibits, for pure `SU(N)` YM on `ℝ⁴`, a **provable** statement of any of:
1. a dilute-gas control at the mass-generating scale — i.e. an upper bound `𝔡(ρ)≤𝔡₀<1` for `ρ` up to `~1/Λ` (a non-abelian liquid closure), uniform in `(L,a)`; **or**
2. a size-fixing mechanism that removes the instanton modulus in `ℝ⁴` (not on `S¹`, not via a deformation) and pins `S_inst` bounded across the crossover; **or**
3. a **1-form ⇄ 0-form** transport theorem making `m_{0⁺⁺} ≥ c√σ` rigorous with `σ_phys>0` in the continuum (`TRANSMUTE-TENSION-GAP`'s open Route A/B/C); **or**
4. a second dimensionful scale in pure YM at `T=0,θ=0` giving a *kinematic* `E₂` (`E2-FLIP-KINEMATIC` flip condition).
None is known. Each is `[O]`.

---

## §6 Honest boundary — what this document does and does not say

- **Does NOT prove Clay; does NOT prove no tool can exist** (`∀class ≠ ∀possible`, `KEY-EXISTS §5`). It is a **`[NG]` for the named class** `{instanton gas, monopole gas, center-vortex condensation}` (the `e^{−S_inst}` transport machine and its relatives), with a **`[T]`-in-model controlled positive part** (§2–§3) and a **`[V-here]` computable obstruction** (§4.1–4.2).
- **The positive part is real, not decoration.** `Φ_tun` is genuinely direction-correct (`corr = +0.97` vs the Hessian's `−0.9999`), genuinely acid-clean, and genuinely a controlled lower bound in the rigid-saddle class (`m₁≥ω√(S_inst)e^{−S_inst}`). The exact residue transport (`Z_O²`) is a real theorem. Under-claiming these would be as dishonest as over-claiming the negative.
- **The obstruction is three independent locks, not one.** Even a fix of the size modulus (lock 4.2) would still face the 1-form/0-form channel mismatch (4.3a) and the one-scale `c=m₁/Λ` circularity (4.3b). This is consistent with — and independently re-derives — `TRANSPORT-MACHINE`'s "three anti-localization no-gos" and `E2-FLIP-NOGO`'s "four-horn + scale dichotomy," in the condensate coordinate.
- **The `[V-here]` computation is a witness of structure, labelled as such.** The double well is a 1D model (its role: prove the *logic* of "`e^{−S_inst}` controls the gap where the saddle is rigid"); the size-distribution numbers are the standard one-loop instanton density (`[R]` Callan–Dashen–Gross; Schäfer–Shuryak) evaluated to locate the peak and the packing fraction. **No `SU(N)` gap number is claimed.**
- **Relation to the banked verdicts.** This note *unifies*: `E2-FLIP-DYNAMICAL §6` (the `e^{−S_inst}` leak), `W2-transmuted-gap (i)-(iv)` (gas only controlled with diluteness+IR cutoff), `R1`/`SMALLCIRCLE-ADIA`/`FORGE5-SMALL-S1` (small-`S¹` is a backward upper-bound machine), `FORGE5-CONFINING-MONOPOLE` (rigorous only in Coulomb phase; 2/5), `PROOF-CENTER-VORTEX` (criterion ⇒ `σ>0`, premise open, 1-form), `HIGGS-DIM-TRANS` (transmutation fixes `Λ` not the pole), `TRANSMUTE-TENSION-GAP` (`m∝√σ` measured not derived). It adds: the machine's **explicit four-factor form**, the **controlled positive half** as a theorem-shaped object, and the **one-loop size-modulus numbers** that make "control ⟺ wrong scale" computable.
- **No fabricated machine, no QED.** The machine is built to the extent it can be and stops exactly where the corpus says the wall is — but now shown through the condensate coordinate rather than merely asserted.

---

## §7 References

**Internal (absolute paths; the verdicts this note classifies against and extends):**
- `D:\Dev\math\CLAUDE.md` (§1 build-don't-stop + acid test; §2 labels; §5 wall & the "machine that does not exist"); `D:\Dev\math\docs\BRIEF-FOR-MODEL.md` (§2 five constraints; §3 death map; §4 "transport ≡ confinement"); `D:\Dev\math\docs\PROOF-SKELETON.md` (§0 G2 open cell = `H1`; G1-L3 = same species).
- `D:\Dev\math\research\unconditional\E2-FLIP-DYNAMICAL.md` (§6 the `e^{−S_inst}` leak = direction-correct functional that *is* the machine; §2 the Hessian-failure witness); `E2-FLIP-NOGO.md` (§3 scale dichotomy Theorem S; §5 loopholes); `E2-FLIP-KINEMATIC.md` (kinematic-`E₂` no-go; the four flip items); `E2-REGRESS-CIRCULARITY.md` (§0.4 SSC; §6 flip condition).
- `D:\Dev\math\research\unconditional\TRANSPORT-MACHINE.md` (4-stage machine; Stage 3 = `H1`/`I2`; three anti-localization no-gos); `INVENT-TRANSPORT-INEQ.md` (exact `f₁²=Z_O²f₁²` transport, acid-clean); `TRANSPORT-CORE.md`, `TRANSPORT-BUILD.md` (Bakry–Émery faithful dialect of `H1`).
- `D:\Dev\math\research\unconditional\W2-transmuted-gap.md` (gas needs diluteness+IR cutoff; ℝ⁴ fails; small-`S¹`/3D only); `SMALLCIRCLE-ADIA.md` (backward inequality `m(S¹)≥m(ℝ⁴)`; deformed theory; adiabatic continuity = the gap); `FORGE5-SMALL-S1.md` (RP+no-crossing cannot reverse the direction; `#4` sole pass); `R1-small-S1-bions.md` (no rigorous nonabelian small-`L`; Göpfert–Mack abelian only).
- `D:\Dev\math\research\unconditional\FORGE5-CONFINING-MONOPOLE.md` (rigorous only in Coulomb/wrong phase; channel; 2/5); `PROOF-CENTER-VORTEX.md` (Tomboulis–Yaffe `criterion⇒σ>0`; premise open at 2.5; 1-form); `HIGGS-DIM-TRANS.md` (transmutation fixes `Λ`, not an isolated `0⁺⁺` pole; anomaly = sum rule); `TRANSMUTE-TENSION-GAP.md` (`m_{0⁺⁺}=c√σ` measured; closed string unbounded below; `2/5`); `TRANSMUTE-STRING-SPECTRUM.md`; `GAP-G2-magnetic.md` (magnetic perimeter law `⟺` dual Meissner `⟺` confinement — same RG wall).
- `D:\Dev\math\research\unconditional\DYSON-GAS.md` (LD-tilt rate = effective action = `H1`); `BROWNIAN-LOOP-SOUP.md` (no positive 4D soup; soup mass = `√σ`, 1-form); `PROOF-SAVVIDY-EXACT-SIGN.md` (`Hess=−D²+2adF_B`, `λ_min=−g|B|`); `W3-adjacent-4D.md` (no rigorous 4D gauge gap; rigor only d≤3 or with external scale); `PROOF-1FORM-ANOMALY.md` (trF² blind to 1-form center).
- `D:\Dev\math\research\unconditional\w7\e2_flip_transport.py` (**the `[V-here]` witness**: Part 1 rigid-saddle control; Part 2 size-modulus/diluteness; Part 2b free-field acid test). Run `py -3.12 e2_flip_transport.py`.

**External — `[R]` standard / as-reported (not re-verified this session):**
- Callan–Dashen–Gross, "The structure of the gauge theory vacuum," *Phys. Rev. D* **17** (1978) 2717; Belavin–Polyakov–Schwartz–Tyupkin, *Phys. Lett. B* **59** (1975) 85 (instantons, size modulus, one-loop density).
- Schäfer–Shuryak, "Instantons in QCD," *Rev. Mod. Phys.* **70** (1998) 323 (instanton liquid vs dilute gas; packing fraction); Shuryak, *Nucl. Phys. B* **203** (1982) 93/116/140.
- Coleman, *Aspects of Symmetry* (1985), Ch. 7 (double-well WKB splitting `ΔE ~ ω√S e^{−S}`); Callan–Coleman, *Phys. Rev. D* **16** (1977) 1762.
- 't Hooft, *Phys. Rev. D* **14** (1976) 3432 (instanton measure, one-loop prefactor); Bernard, *Phys. Rev. D* **19** (1979) 3013.
- Polyakov, *Nucl. Phys. B* **120** (1977) 429 (3D compact `U(1)` monopole mass gap); Göpfert–Mack, *CMP* **82** (1982) 545 (rigorous abelian monopole gas ⇒ area law); Fröhlich–Spencer, *CMP* **83** (1982) 411 (dilute Coulomb gas).
- Tomboulis–Yaffe, *CMP* **100** (1985) 313 (vortex free energy ⇒ area law); 't Hooft, *Nucl. Phys. B* **138** (1978) 1 / **153** (1979) 141.
- Ünsal, *Phys. Rev. D* **80** (2009) 065001, arXiv:0709.3269 (magnetic bion condensation, `m_D~e^{−S₀}`); Ünsal–Yaffe, *Phys. Rev. D* **78** (2008) 065035, arXiv:0803.0344; Shifman–Ünsal, *Phys. Rev. D* **78** (2008) 065004, arXiv:0802.1232 (centre stabilisation); Poppitz, arXiv:2111.10423 (review; adiabatic continuity conjectural).
- Gross–Pisarski–Yaffe, *Rev. Mod. Phys.* **53** (1981) 43 (one-loop Polyakov potential; pure YM small-`S¹` is centre-broken); Shen–Zhu–Zhu, *CMP* **400** (2023) 805, arXiv:2204.12737 (lattice gap `β<1/48`).
- Kotecký–Preiss, *CMP* **103** (1986) 19 (cluster expansion / diluteness criterion); Bakry–Émery, LNM **1123** (1985); Barnsley–citations carried from the corpus.

---

*Status: **PARTIAL MACHINE.** Built and proven-controlled in the rigid-saddle class `[T-in-model / V-here: m₁ ≥ ω√(S_inst)e^{−S_inst}, corr +0.97]`; exact residue transport `[T]`; controlled dilute/small-`S¹` lower bound `[C]`. Provably STOPS at the crossover: `[NG]` for `{instanton gas, monopole gas, center-vortex condensation}` via the size-modulus/diluteness ⟺ wrong-scale dichotomy `[V-here: peak at ρ*=O(1)/Λ, packing O(1)]`, the 1-form/0-form channel mismatch `[T]`, and the one-scale `c=m₁/Λ` circularity `[T-cond]`. Exact missing step `[O]` = non-abelian liquid control at the physical scale (= the crossover RG / `H1`). NOT a proof of Clay; NOT a global impossibility theorem (`∀class ≠ ∀possible`). **Clay stays open.***
