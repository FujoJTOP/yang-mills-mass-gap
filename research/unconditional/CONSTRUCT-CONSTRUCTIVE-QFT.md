# CONSTRUCT-CONSTRUCTIVE-QFT — the best constructed 4D gauge measure, what it proves, and the exact missing crossover estimate

> **VERDICT (no QED).** "Build the theory" is the honest long-line target. The state of the art is a
> **patchwork of constructed measures, one per regime**, none of which reaches the crossover with a gap:
> **(1) existence** of the 4D OS theory is, separately, essentially done (OS reconstruction + Balaban UV +
> the machine-checked Lean construction) — but **existence only, gap excluded**; **(2) a rigorous mass gap
> exists for `β<1/48`** (Shen–Zhu–Zhu, cluster/LSI) — a small-`β` (absolute) region; **(3) the rigorous
> continuum control at `β→∞`** (Balaban) is **UV, not IR**; **(4) 3D** stochastic quantization of YM-Higgs
> is constructed (CCHS), but **4D is the critical dimension** where the counterterms are the `u=2`
> renormalon. **At the crossover `β≈2.5` there is no constructed measure with a proven gap.** Every
> constructive-QFT method is an *absolute-convergence* method and needs a **small parameter**; at the
> crossover the only available small parameters — `1/N` (at `N=3`), a small circle, strong coupling — are
> all `O(1)` at the physical point. **The exact missing estimate is a uniform-in-`(L,a)` lower bound on the
> `0⁺⁺` spectral gap at `β≈2.5`, i.e. `χ_lattice ≲ a²` = `H1`.** Clay 0%.

---

## §0 · Labels, acid, iron laws (carried)

**[T]** 真定理 · **[C]** 条件定理 · **[O]** 精确缺失步 · **[NG]** no-go for the named class · **[V-num]** new numeric.
**Acid**: the free massless field (gap `→0`) must be a solution; the gapped-atomless `:φ²:` must not be certified.
**IL-i** non-positivity · **IL-ii** nonabelian · **IL-iii** no small parameter at `α_s~O(1)`.

---

## §1 · The constructed-measure landscape (what EXISTS, what it PROVES)

| regime | constructed object | what it proves | gap? |
|---|---|---|---|
| **Existence (G1)** | 4D OS theory via OS reconstruction of the lattice measure + Balaban UV control (+ machine-checked Lean, OpenAI 2026 `[T]`) | the Wightman/OS theory **exists**, all OS axioms | **NO** (explicitly excludes the gap) |
| **Strong coupling `β<1/48`** | lattice measure; cluster/LSI expansion (Osterwalder–Seiler; Shen–Zhu–Zhu, CMP 400 (2023) 805) `[T]` | reflection positivity + LSI + Poincaré ⇒ **mass gap, area law, clustering**, `L`-uniform | **YES**, but `λ≡g²N>96` |
| **Continuum `β→∞`** | Balaban RG (CMP 109/116/122) `[T]` | UV stability, renormalization, controlled effective action | **NO** (UV, not IR) |
| **Infrared crossover `β≈2.5`** | — | — | **none rigorous** |
| **3D stochastic quantization** | CCHS (2022) `[T]` | rigorous 3D YM-Higgs measure + FP analysis | 3D only; **4D critical** |
| **4D stochastic quantization** | Hairer regularity structures `[NG]` | needs subcriticality + finitely many counterterms | **4D marginal** ⇒ counterterms = `u=2` renormalon |
| **Large-`N` master field** | Singer (1994) `[T]`, **2D only** | exact 2D master field `⟨W⟩=exp(−λA/2)` | **4D unconstructed** |
| **Area law (large `β`)** | Cao–Nissim–Sheffield (arXiv:2505.16585) `[T]` | Wilson **area law**, N-independent, into Hooft regime | area law, **not the transfer-matrix gap** |

**The honest picture:** the theory is constructed *at both ends* (existence + strong-coupling gap on one
side, UV continuum control on the other) with a **notch in the IR middle** exactly at `β≈2.5`. No single
constructed object covers the crossover with a gap.

---

## §2 · The most complete constructed object, and what it proves

The **most complete** constructed object is the **OS theory itself** (G1): by OS reconstruction of the
lattice measure plus Balaban's UV renormalization, one has a genuine 4D Wightman field theory (and this is,
separately, machine-checked). **What it proves:** existence, all OS axioms, reflection positivity, the
transfer matrix. **What it does not prove:** that the transfer matrix has a spectral gap above the vacuum
in the `0⁺⁺` channel. The massless free scalar is a perfectly constructed OS theory with `Δ=0`; existence
does not force a gap.

The most complete **gap** statement is the strong-coupling theorem (SZZ): a rigorous, `L`-uniform mass gap,
but confined to `λ>96` (`β<1/48`). Between `1/48` and `2.5` the measure is constructed (via G1) but its
`0⁺⁺` gap is **not bounded below by any rigorous estimate**.

---

## §3 · The three "small parameters" — and why none transfers to the physical crossover

Constructive QFT is **absolute-convergence**: it needs a small parameter to make a cluster/expansion
convergent. At the crossover `α_s~O(1)`, the only candidate small parameters are deformations, and each is
`O(1)` at the physical point:

| small parameter | where it exists | value at the physical point | obstruction |
|---|---|---|---|
| **strong coupling** `λ` (cluster/LSI) | `λ>96` (SZZ) | `λ_cross≈7.3` | needs `λ` **large**, physical `λ` is `O(1)`; ratio `96/7.3≈13` |
| **`1/N`** (large-`N` master field) | `N→∞` | `N=3 ⇒ 1/N≈0.33` | `O(1)` correction; and the 4D master field itself is unconstructed (Singer 2D) |
| **small circle** `Λ L` (deformed / Unsal–Yaffe) | `L≪1/Λ` | decompactification `L→∞ ⇒ ΛL≫1` | the limit back to `R⁴` is the transport wall (`I2`) |
| **small `β`** (stochastic quant / Balaban) | `β→0` or `β→∞` | `β=2.5` is **interior** | no absolute convergence in the interior |

**None of the four is small at `β≈2.5, N=3, L=∞`.** This is the precise sense in which "the crossover is
beyond constructive QFT": it is not that the constructions are weak, it is that **their convergence
parameter is `O(1)` exactly where the creature lives**. (`BASE-PHASE-RIGOROUS` gives the same conclusion
machinery by machinery: cluster `≈13×` short; Dobrushin an order of magnitude; Pirogov–Sinai wrong object;
Peierls circular.)

---

## §4 · The exact missing estimate

What a "construction at the crossover" would have to deliver, stated with no slack:

> **A uniform-in-`(L,a)` lower bound on the `0⁺⁺` transfer-matrix spectral gap `m_{0⁺⁺}(β,a)` at
> `β≈2.5`**, equivalently (Lemma G, `PNT-CHI-CERTIFICATE`) a uniform upper bound `χ_phys ≤ 1/λ` with
> `χ_phys = χ_lattice/a_ℓ²`, equivalently `χ_lattice ≲ a_ℓ²` — dimensional transmutation.

Every constructive route reduces to this one estimate:
- cluster/LSI: the expansion radius would have to be uniform to `λ=7.3` — `[NG]` (structural, `≈13×`).
- stochastic quant: the 4D SPDE would have to be made well-posed with all-orders controlled counterterms
  — `[O]`, the `u=2` renormalon in the SPDE coordinate.
- large-`N`: the 4D master field would have to be constructed and the `0⁺⁺` mode linearized — `[O]`.
- Balaban: would have to be extended from UV to an IR gap — `[NG]` (it does not produce `m₁`).
- area law (CNS): would have to cross the area-law→gap bridge in the `0`-form channel — `[O]`.

All five are the **same estimate** in five coordinates. This is `H1`.

---

## §5 · Honest bottom line

- **Nothing is "built" at the crossover.** The constructed-measure landscape covers existence (G1), strong
  coupling (`β<1/48`), and the UV (`β→∞`); the IR crossover is a **notch**, and the exact missing estimate
  is a uniform `(L,a)` `0⁺⁺` gap bound at `β≈2.5` = `H1`.
- **"Build the theory" is well-posed but not a shortcut.** Constructing the full measure does **not** by
  itself give the gap (the free massless scalar is a constructed gapless OS theory). The gap is a property
  of the *specific* nonabelian measure, and every constructive method's convergence parameter is `O(1)`
  exactly at the crossover.
- **What is genuinely available:** existence (`[T]`); the strong-coupling gap (`[T]`); Balaban UV (`[T]`);
  3D stochastic quantization (`[T]`); the area law at large `β` (`[T]`). Together: the base with a notch.
- **The long-line target, stated honestly:** not "construct the measure" (done, for G1), but **"construct a
  small parameter at `α_s~O(1)`"** — a deformation whose convergence survives decompactification/large-`N`
  back to the physical point. None exists; inventing it is the long-line construction.
- **Clay = 0%.** Base stands with a notch; the notch is `H1`.

*Refs:* `BASE-PHASE-RIGOROUS.md`, `ROOF-STOCHASTIC.md`, `PNT-CHI-CERTIFICATE.md`, `ROOF-LARGEN-MASTER.md`,
`INVENT-SYNTHESIS.md`, `PROFILE.md`; witness `research/unconditional/w7/construct_smallparams.py`.
Citations: Shen–Zhu–Zhu arXiv:2204.12737 / CMP 400 (2023) 805; Balaban CMP 109/116/122;
Cao–Nissim–Sheffield arXiv:2505.16585; CCHS (2022); Singer (1994); Osterwalder–Seiler (1978).
