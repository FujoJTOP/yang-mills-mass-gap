# ICA-TEST — the "Irreversible Constraint Accumulation / ratchet" distance is the **conditional entropy** of the fine state given its coarse image: its monotonicity **is** the data-processing inequality, it is **larger for the gapless field than for the gapped one** (wrong sign), and all four kill tests land. **KILLED: known monotone, wrong direction, not necessary, not sufficient.**

> 2026-11-24 · Adversarial test of a proposed sixth-species object. **No QED. Clay 0%.**
> **Conjecture under test (ICA).** For a coarse-graining sequence `T_1, T_2, …`, let
> `D_n = dist(T_n⋯T_1(x), R_n)` where `R_n` is the set of states reachable by *allowed local recovery
> operations*. **ICA conjecture:** `D_{n+1} ≥ D_n − ε_n` (a ratchet: the "damage" cannot be undone
> locally). **Question:** is `D_n` a *new* monotone, or is it a known monotone (relative entropy,
> Wasserstein, RG monotone, free energy) renamed — and does it survive the four acid tests?
> **Reads:** `Q-WAVE-SYNTHESIS.md`, `T-WAVE-SYNTHESIS.md`, `docs/paper-v2/ch06-tools.tex`,
> `unconditional/{MONOTONE-COUPLING,TRANSPORT-BLIND,FALSIFY-ENTROPY}.md`, `v3/U5-ENTANGLEMENT-BINDING.md`,
> `v3/R5-EXACT-ANCHORS.md` §1.1, `v3/T5-CENTER-VORTEX.md`.
> **Witness:** `research/v3/w7/ica_witness.py` (`py -3.12 -I`, numpy only, deterministic, exit 0; output §6).
> **Labels.** **[T]** true theorem/identity · **[C]** conditional/model statement · **[NG]** no-go for a named
> class · **[O]** exact missing step · **[V]** numerical witness here · **[std]** standard tool used as such.

---

## 0 · Verdict (up front; no QED, Clay 0%)

**ICA is killed by both heads of its own criterion.** (i) `D_n` **reduces to a known monotone** — the
conditional entropy `H(x | y_n)` of the fine configuration given its `n`-fold coarse image (equivalently the
mutual-information / accumulated data-processing deficit); the ratchet `D_{n+1} ≥ D_n` is the
elementary **data-processing inequality** `H(x|y) ≤ H(x|f(y))` for `y_{n+1} = f(y_n)`, which is
**model-independent** and therefore cannot see a mass gap. (ii) `D_n` **fails all four acid tests**: it is
**maximal on the free massless field** (kill 1), tracks a **free energy** on exactly-solvable 2D YM which
has an area law but **no `0⁺⁺` atom** (kill 2), keeps growing with **no gap** (kill 3), and vanishes
exactly on 1D gauge theory which **does** have a transfer gap (kill 4). Worse, the sign is **backwards**: a
finite mass makes the field locally recoverable and *shrinks* `D` (`D_1: 599.9→435.6` nats from `m=0` to
`m=1`, §6), so "gap ⇒ large/persistent `D`" is false at the level of a single free field. The only
mass-dependent content is an `O(1)` subleading term on top of a `∝` number-of-modes UV background — the
same **gap-blind leading term** every entropy monotone has already been shown to carry (U5 §3, H1-T4 §3,
FALSIFY-ENTROPY). This is the **Schur-concave / wrong-direction** profile (MONOTONE-COUPLING §2,
TRANSPORT-BLIND): monotone objects are *totals* that bound the spread; the gap lives at the concentrated end.

**One line.** *`D_n` is the conditional entropy `H(x|y_n)`; its ratchet is the DPI, so it is a known monotone
and model-blind; it is largest on the free massless field (acid fail), equals a free energy on 2D YM with no
atom (kill 2), grows without a gap (kill 3), and is `0` on gapped 1D (kill 4); a finite mass shrinks it —
the sign is reversed. KILLED. No QED. Clay 0%.*

---

## 1 · Precise definition of `R_n` and `D_n` on a lattice gauge model

**The model.** Lattice torus `Λ = (ℤ/Lℤ)^d`, gauge group `G ∈ {U(1), SU(N), Z_N}`, link variables
`U_ℓ ∈ G`. Configuration space `Ω = G^{E}`, gauge group `𝒢 = G^{V}` acting by conjugation; physical space
`Ω/𝒢`. Gibbs measure `μ ≡ μ_β = Z^{-1} e^{-S_W} dU`, `S_W` the Wilson action. The `0⁺⁺` atom of the source
`O = ¼ tr F²` is `f₁² = ν({m₁²})` on the Källén–Lehmann measure `ν`.

**Coarse-graining (deflation) `T`.** A deterministic, gauge-covariant, local map `T: Ω_{n} → Ω_{n+1}`. Canonical
choice: **block average** in the flat chart `X_ℓ = log U_ℓ` (the deflation the corpus already uses,
`docs/paper-v2/app.tex`): level-`n` link ↔ block of `b^n` fine links, `T_n` the orthogonal projection onto
block-constant links. `T_n` is linear (flat chart) and preserves the zero-mean (constant-mode-fixed) subspace.
Write `y_n = T_n⋯T_1(x)`, and `μ_n = (T_n⋯T_1)_# μ_0`.

**Allowed local recovery operations.** Two equivalent formalizations; we test both.

- **(D) Distributional (information-theoretic) — the object we compute.** With `y_n = T_n⋯T_1(x)`, set
  ```
     D_n := H( x | y_n )        (conditional entropy;  = H(x) − I(x;y_n))                  [D-def]
  ```
  the **Bayes-optimal residual entropy** of the fine configuration given its `n`-fold coarse image — the
  irreducible fine detail the coarse image fails to determine. (In a classical Gibbs model the infimum of
  `D(μ_0 ‖ Λ_# μ_n)` over *all* recovery maps is trivially `0` — Bayes inversion `Λ(x|y) = p(x|y)` reproduces
  `μ_0` exactly — so the non-trivial object is precisely this residual entropy, i.e. the mutual-information
  deficit.) Restricting recovery to the ICA "allowed **local** operations" can only make the residual
  **larger**, `D_n^{loc} ≥ H(x|y_n)`, so every acid kill below holds a fortiori. For a Gaussian field
  `D_n = ½ log det Σ_{x|y_n} + const` (Schur complement), and `D_n = 0` iff `x` is a function of `y_n`.
- **(P) Point (geometric).** `R_n := { y ∈ Ω_n : T^{-n}(y) is a single gauge orbit }` — the states from which
  the fine configuration is *exactly* recoverable. Then `D_n = dist(y_n, R_n)` in the gauge-invariant
  (bi-invariant / flat-chart) metric. For Gaussian models this distance to the exactly-recoverable subspace
  `{x : x ∈ σ(y_n)}` equals (D) (Lemma 1).

**The ratchet.** `D_{n+1} ≥ D_n − ε_n`. In form (D) the residual is non-decreasing; in form (P) the distance
to the exactly-recoverable set is non-decreasing; `ε_n` is the one-step slack (taken `0` for exact recovery).

> **Lemma 1 (identification + ratchet) [T].** For every model (gauge or not), every `n`, writing
> `y_{n+1} = T_{n+1} y_n` (a deterministic function of `y_n`),
> ```
>   D_n = H( x | y_n ) = H(x) − I(x;y_n) = Σ_{k=1}^{n} δ_k ,    δ_k := I(x;y_{k-1}) − I(x;y_k) ≥ 0 ,
> ```
> the **conditional entropy** = **mutual-information deficit** = accumulated **data-processing deficit**.
> Hence `D_{n+1} ≥ D_n` **exactly** (`ε_n=0`): conditioning on a coarser image cannot lower entropy
> (`H(x|y) ≤ H(x|f(y))`, the data-processing inequality). For a Gaussian field `D_n = ½ log det Σ_{x|y_n} +
> const`, which is also the point-version distance to the exactly-recoverable subspace. ∎

**Identification with known monotones [T].** `D_n` is, term by term:
- the **conditional entropy** `H(x|y_n)` — a textbook monotone;
- equivalently the **mutual-information / accumulated data-processing deficit** (`H(x) − I(x;y_n)`; Lindblad/Petz);
- in the point/Gaussian form, the **free-energy difference per block** (`ΔF` of the RG step) when `R_n` is the
  exactly-invertible set among a product measure (2D, §2.2), or the **RG relevant-coupling distance**
  (distance to the critical surface) if `R_n` is instead taken to be the critical manifold (§3).

There is no fourth option: each is a **known** object, and none is gap-forcing.

---

## 2 · The four kill tests on exact models

### 2.1 · Kill (1) — free massless field: `D_n` is maximal, not `0` [NG]+[V]

Exact Gaussian free scalar, `C = (m²−Δ)^{-1}` on the zero-mean subspace (`L=512`, `b=16`, one step). The
one-step ICA distance is `D_1 = H(x|y_1)`:

| `m` | `D_1` (nats) | `D_1(0) − D_1(m)` |
|---|---|---|
| 0.00 | **599.92** | 0.00 |
| 0.05 | 598.83 | 1.09 |
| 0.10 | 595.72 | 4.20 |
| 0.30 | 569.97 | 29.95 |
| 1.00 | 435.64 | 164.27 |

The **free massless field (`m=0`) gives the largest `D_1`**, ~`600` nats — it is emphatically *not* `0`. In
2D (16×16, block 2×2) the same distance is `D_1 = H(x|y_1) = 130.5` nats at `m=0` (mutual information
`I = 1005.9`). A tool that returns a large positive value on the free massless field is **vacuous** by the
charter's acid rule (`ch06 §acid`): it cannot distinguish the free field from YM. **[NG]** (A gap-forcing tool
must return `0`/`≤0` here.)

### 2.2 · Kill (2) — 2D YM holonomy: `D_block` is a free-energy difference, blind to the atom [T]+[V]

2D pure YM is exactly solvable and is a **product measure over independent plaquettes** (Migdal
exact solution; for the finite-group instance this is exact — on a disk every plaquette assignment is
realized by links). For `Z_N` with single-plaquette law `q(s) ∝ e^{β cos(2πs/N)}`, a `b×b` block has `b²`
independent plaquettes; its coarse holonomy is their group product `q^{⊞ b²}`. The ICA distance is exact:
```
   D_block = b² H(q) − H(q^{⊞ b²})   =   entropy / free-energy difference of the merged plaquettes.
```
Witness (`N=2,3,4`, `β=0.5,1,2`): `D_block > 0` for every entry (e.g. `N=2, β=1`: `D_b=2,3,4 = 0.826, 2.599,
5.152`; `N=4, β=0.5`: `3.917, 10.545, 19.825`). Two structural facts:

- **`D_block` is a known thermodynamic quantity** (a free-energy/entropy difference of merging `b²`
  independent local variables), **not** a new order parameter.
- 2D pure YM **confines** (`σ_fund = g²C₂/2 > 0`, exact Migdal) **yet has `f₁²(trF²)=0` exactly**: in 1+1D
  `trF² = trE² = H` is the (conserved) Hamiltonian density, diagonal in the character basis, so
  `⟨0|trF²|G⟩ = C₂·δ_{0R} = 0` (R5 §1.1, T5). Thus `D_block > 0` **without** a `0⁺⁺` atom: ICA tracks the
  free energy / string tension, not the atom. **[T]** for the identities, **[V]** for the numbers.

### 2.3 · Kill (3) — ICA but no gap: `D_n` keeps growing on the gapless field [T]+[V]

The ratchet holds and `D_n` grows at `m=0` (`L=512, b=16`, successive halvings): `D_1=599.9, D_2=652.1,
D_3=683.8, D_4=702.4`, increments `52.2, 31.6, 18.6` — **all positive, no early saturation** (at fixed `L`
`D_n ↑ H(x)`; `H(x)` itself diverges as `L→∞` for the massless field, the log-divergent low-mode entropy).
The free massless field has **no gap**, yet `D_n` grows monotonically throughout — so "ICA holds with growing
`D`" does **not** imply a gap. `D_n` is the DPI deficit, `> 0` for any lossy channel on any model.

### 2.4 · Kill (4) — gap but ICA `=0`: 1D gauge theory, and the massive free field [T]+[V]

- **Exact 1D `Z_N` (or any `G`) gauge on a circle.** 1D has **no local gauge-invariant degree of freedom**:
  gauge-fixing leaves a single global holonomy, and blocking the links loses nothing (the coarse holonomy
  *is* the fine holonomy; the interior links are pure gauge). Hence `D_n = 0` **identically** for every `N,
  β`, yet the transfer matrix `∝ −Δ_G` on `L²(G)` has a **gap** (first Casimir). This is `gap ∧ ICA=0`.
  (Equivalently `f₁²(1D, trF²) = 0`, R5 §1.1.)
- **Massive free field.** `D_1` **shrinks** with `m` (`599.9 → 435.6`, §2.1): a finite correlation length
  `ξ=1/m` makes the field locally recoverable, so the *gapped* theory has the *smaller* ICA distance. At
  `ξ ≪ block` the discarded detail is exponentially small — the wrong side for a certificate.

So ICA is **not necessary** for a gap (§2.4) and **not sufficient** (§2.3).

---

## 3 · The strongest steelman — "distance to the RG critical surface" — and why it also dies [C]/[NG]/[O]

The most charitable non-degenerate reading: take `R_n` to be the **RG critical/stable manifold** (not the
exactly-recoverable set), so `D_n = dist(y_n, critical surface) = |renormalized mass|` (the magnitude of the
relevant coupling). This *does* have the right sign for a free field — gapless ⇒ `D=0`, massive ⇒ `D` grows —
and it passes the free-massless acid test. It dies on three counts:

1. **It is the known RG monotone** — the renormalized mass / relevant eigenvalue (Zamolodchikov-type
   `c`-function logic). Computing `D_n` *is* computing the beta function of the mass; in a nonperturbative
   gauge theory that is the hard part, not a new lever. **[T]/[std]**
2. **"Relevant coupling grows" ⇏ gap.** A relevant deformation can drive a crossover or a first-order
   transition with a gapless or metastable spectrum; conversely a gap can be generated by a
   marginally-irrelevant strong-coupling effect (as in 4D YM). So test (3) fires. **[C]/[NG]**
3. **4D YM's mass is emergent, not UV-relevant.** At the asymptotically-free Gaussian fixed point the pure
   gauge theory sits *on* the critical surface (`g→0`, marginal coupling); the mass is generated by
   **dimensional transmutation** at the strong-coupling scale `Λ`, not turned on as a UV-relevant operator.
   Hence the *computable-regime* distance to the critical surface is `≈ 0` while the true gap is `> 0` — test
   (4) fires exactly in the regime where the tool is usable. Getting a nonzero computable `D_n` requires
   already knowing `Λ`, i.e. the gap (circular). **[O]**

So the steelman is either a rename of the mass (circular read) or fails tests (3),(4); it is not a sixth
species.

---

## 4 · Why the sign is wrong: `D_n` is a Schur-concave total (the banked obstruction, in one more costume) [T]

The witness makes the sign explicit: `D_1` is **largest for the gapless (spread-out) field and shrinks as the
mass concentrates it**. This is the corpus's established monotone anatomy:

- **Monotone objects are totals; the atom is a coherence.** `D_n = H(x|y_n)` counts *how much detail is
  spread across discarded modes* — a basis-independent aggregate. The atom `f₁²` is a basis-dependent
  coherence in the operator-`O` basis (U5-4, `ch07 §schur-master`). An aggregate of the coarse/fine
  relative entropy cannot lower-bound a coherence; the sign runs the wrong way (MONOTONE-COUPLING §2,
  TRANSPORT-BLIND §1/§3).
- **The mass enters only subleadingly.** The `~600`-nats `D_1` is the **UV mode count** (`∝` number of
  discarded modes, gap-blind); the entire `m`-dependence is the `O(1)` `D_1(0)−D_1(m)` column. This is
  precisely the area/volume-law blindness U5 §3 and H1-T4 §3 flag: the leading entropy is identical for
  gapped and gapless theories; only an `O(1)` subleading piece carries IR information, and no theorem
  converts it into an atom floor (FALSIFY-ENTROPY: convex/density-requiring/Markov — all blind).
- **Free field is the counterexample to every floor.** Just as a universal floor `w₁ ≥ f(S)` fails on the
  free chain (large `S`, `w₁=0`, U5-3), a gap-forcing use of `D_n` fails on the free massless field (large
  `D`, no gap, §2.1/§2.3).

---

## 5 · Ledger, exact missing step, verdict

| item | statement | label |
|---|---|---|
| identification | `D_n = H(x|y_n) = H(x)−I(x;y_n) =` Σ data-processing deficits | **[T]** |
| ratchet | `D_{n+1} ≥ D_n` (`ε_n=0`) is DPI `H(x|y)≤H(x|f(y))`; model-independent | **[T]** |
| point form | `=` distance to exactly-recoverable subspace `=` Gaussian Schur entropy / free-energy diff | **[T]** |
| kill 1 (free massless) | `D_1 = 599.9` nats 1D, `130.5` nats 2D (`m=0`) → acid FAIL | **[NG]**+**[V]** |
| kill 2 (2D YM) | `D_block = b²H(q)−H(q^{⊞b²}) > 0`, a free-energy difference; 2D has `σ>0`, `f₁²(trF²)=0` | **[T]**+**[V]** |
| kill 3 (ICA, no gap) | `D_n` grows (599.9→702.4, `m=0`) with no gap → not sufficient | **[T]**+**[V]** |
| kill 4 (gap, ICA=0) | 1D exact `D_n=0` with transfer gap; massive field shrinks `D` → not necessary | **[T]**+**[V]** |
| steelman (RG critical surface) | `= |mass|`, a known RG monotone; relevant coupling ⇏ gap; 4D mass emergent | **[C]**/**[NG]**/**[O]** |
| sign | `D` largest at the gapless end; mass shrinks it → wrong direction (Schur-concave total) | **[T]**+**[V]** |
| **ICA verdict** | **KILLED** — reduces to a known monotone **and** fails all four acid tests | **[NG]** |
| Clay | untouched; no new lever | **[O] · 0%** |

**Exact missing step (unchanged by ICA).** The gap still needs the `0⁺⁺` atom `f₁² = ν({m₁²})` with
isolation `I₂` on the Källén–Lehmann measure — equivalently `w₁ = f₁²/M₀ ≥ w₀ > 0` (T-WAVE §5, U5 §8) —
an **operator-specific coherence**, which no conditional entropy / recovery / free-energy / RG monotone
supplies (`H*` of T-WAVE-SYNTHESIS). ICA re-expresses the wall (as "the RG increments are not summable") and
adds a sharp wrong-direction no-go; it does not lower it.

**Verdict.** *ICA is a rename of the conditional entropy `H(x|y_n)` (= mutual-information deficit =
accumulated DPI deficit = free-energy difference / RG relevant-coupling distance). Its ratchet is the
data-processing inequality — true for every channel and every model, hence gap-blind. On the exact free
massless field it is maximal (acid fail); on exact 2D YM it is a free energy with no atom; it grows without a
gap and vanishes with a gap (1D); and a finite mass *shrinks* it (wrong sign). All four kill tests land.
KILLED. No QED. Clay 0%.*

---

## 6 · Witness output (`research/v3/w7/ica_witness.py`, `py -3.12 -I`, numpy only, exit 0)

```
(A) Gaussian free scalar, 1D L=512, block b=16:  D_1 = H(x|y_1)
        m         D_1    D_1(0)-D_1(m)
     0.00    599.9189           0.0000
     0.05    598.8316           1.0872
     0.10    595.7180           4.2009
     0.30    569.9703          29.9485
     1.00    435.6445         164.2744
   leading ~600 nats = UV mode count (gap-blind); mass effect is the O(1) column.
   ratchet at m=0: D = 599.9, 652.1, 683.8, 702.4 (increments 52.2, 31.6, 18.6; no saturation)
(A') 2D Gaussian 16x16, block 2x2:  D_1 = H(x|y_1) = 130.5 (m=0), 128.4 (m=0.3), 109.8 (m=1)
(B) exact 1D Z_N gauge on a circle: D_n = 0 identically (transfer gap, no local DOF)
(C) exact 2D Z_N gauge, D_block = b^2 H(q) - H(q^{*b^2}) > 0 for all N,beta, b (free-energy diff)
      e.g. N=2,beta=1: D_b=2,3,4 = 0.826, 2.599, 5.152 ; N=4,beta=0.5: 3.917, 10.545, 19.825
```

Witness is deterministic and reproducible; the identities are §1–§4, the numbers are corroboration **[V]**.

---

*Refs (internal).* `Q-WAVE-SYNTHESIS.md` (§1 four faces; `χ≤1/λ` ceiling); `T-WAVE-SYNTHESIS.md` (§1
predicates `A,W,I,G`; §5 missing `w₁≥w₀`); `docs/paper-v2/ch06-tools.tex` (§acid three-adversary test;
§toolwave aggregate-upper-bound no-go); `unconditional/MONOTONE-COUPLING.md` (§2 wrong direction: gap ↓ in
β, ↑ in m); `unconditional/TRANSPORT-BLIND.md` (§1/§3 `Γ₂`=H1 lossless; Haar floor = compactness artifact);
`unconditional/FALSIFY-ENTROPY.md` (convex/density-requiring/Markov classes blind); `v3/U5-ENTANGLEMENT-BINDING.md`
(§3 no universal floor, free chain `S∝log l`, `w₁=0`; §5 operator-blindness); `v3/R5-EXACT-ANCHORS.md` §1.1
(2D `trF²=H` conserved, `f₁²=0`; 1D `f₁²=0`); `v3/T5-CENTER-VORTEX.md` (§2 `σ>0 ⇏ atom`, 2D exact);
`unconditional/INVENT-NEW-MONOTONE.md` (Perelman W: monotone on gapless side, sign-blocked on gapped side).

*Refs (external, standard; used as textbook [std], not re-derived).* Fawzi–Renner, *Quantum conditional mutual
information and approximate Markov chains*, CMP **340** (2015) 575, arXiv:1112.4452 (the recoverability family; class. Gibbs shadow = conditional entropy, §1);
Sutter–Fawzi–Renner, *Universal recovery map for approximate Markov chains*, Proc. R. Soc. A **472**
(2016), arXiv:1504.07297; Petz, *Sufficient subalgebras and the relative entropy of states of a von Neumann
algebra*, CMP **105** (1986) 123; Lindblad, *Completely positive maps and entropy inequalities*, CMP **40**
(1975) 147 (data processing); Zamolodchikov, JETP Lett. **43** (1986) 730 (`c`-function / RG monotone);
Migdal, *Recursion equations in gauge theories*, Sov. Phys. JETP **42** (1975) 413 (exact 2D YM product
measure). No third-party downloads; no `[C]` upgraded to `[T]`; no QED; Clay 0%.
