# SL-TEST — the spectral-identity cluster (Lineage / Event / Shadow): an exact **2×2** separates every one of them from the mass gap, and each collapses onto a **named banked object** — the cluster is the campaign's three coordinates **W / I / H1** in projective / spectral / counterfactual dress, not a sixth species

> **2026-11-24 · Wave SL.** Target: `m₁² = inf supp ν > 0` with edge atom `f₁² > 0` (the `0⁺⁺`
> Källén–Lehmann atom `W`); Clay is `W ∧ H1` (or isolation `I ∧ (ρ²≤C²)`), `CLAUDE.md` §5,
> `W-WAVE-CHARTER.md` §0. Three candidate *spectral-identity* engines are put on trial:
> **(A) Spectral Lineage** (projective system of finite-resolution algebras `Aₙ`, states `ωₙ`,
> intervals `Iₙ`, restrictions `r`; *exist iff evidence survives refinements*),
> **(B) Unfinished Spectral Event** (candidate projections `Pₙ` with allowed/forbidden refinements;
> *three fates*), **(C) Counterfactual Spectral Shadow** (*remove a candidate object's traces; ask if
> the remaining observables satisfy the same relations*). Each is given its smallest **exact 1D/2D
> model**, then run through a fixed **four-case battery**: free massless · 2D YM holonomy · a case
> with the property but no gap · a gapped case without the property.
>
> **Labels (charter §2).** **[T]** true theorem/identity · **[C]** conditional/reading · **[NG]** no-go
> for a named class · **[O]** exact missing step · **[V]** exact numeric here. **No QED. Clay 0%.**
> **Witness:** `research/v3/w7/sl_test.py` (`py -3.12 -I`, numpy 2.4.6 / scipy 1.18.1, one job, no
> downloads; output quoted §6).

---

## 0 · Verdict up front (no QED; Clay 0%)

1. **The four-case battery is exact and separates ATOM from GAP in both directions `[T]+[V]`.** The
   common control set (§2) is an exact 2×2: **free massless** (no gap, no atom); **2D YM SU(2)** (gap
   `3e²/8`, tr F²-atom `f₁²=0` because tr F²=H is conserved); **`wδ_m + (1−w)Unif[0,Λ]`** (atom weight
   `w>0` but `inf supp=0`, no gap); **gapped atomless** (2-particle threshold `c√(s−m²)`, `inf supp=m²>0`,
   no atom). A positive control (free massive scalar: atom at `m²` **and** gap) fixes the fourth corner.
   So **atom and gap are independent**, exactly the campaign's known dichotomy (R6 / U6 / W-WAVE §0),
   here re-derived with exact numbers.

2. **(A) Spectral Lineage `[NG]` — reduces to the inverse limit (compactness) + `W`.** The
   *existence* of a lineage is the finite-intersection property on the weak-\* compact set of consistent
   features: **automatic**, with no gap content (Banach–Alaoglu). The only quantitative content is that
   the lineage carry **uniform mass `≥w>0`** — which *is* the atom floor `W`. Free massless exhibits an
   `O(1)`-weight lineage whose **location `m_n=λ₁(L)→0`** (the `R1` face): lineage present, floor zero.
   2D YM has a gap but **no** tr F²-lineage (`f₁²=0`). ⇒ *neither necessary nor sufficient.*

3. **(B) Unfinished Spectral Event `[NG]` — the "three fates" are the spectral-flow trichotomy; it
   tracks isolation `I` or renames the gap, blind to the atom weight.** Fires = isolated eigenvalue
   (a spectral-flow crossing / a discrete lower edge); Forbidden = rank drop; **Unfinished = the
   eigenvalue approaches the threshold with no finite certificate = `H1` (weak-\* discontinuity)**. By
   `FALSIFY-FLOW.md`, the same crossing supports `f₁²=0` or `f₁²=1` — the event is **spectrum, not
   overlap**. Two readings, both `[NG]`: point-projection → **operator-algebra generation**; uniform
   spectral gap → **spectral-flow existence** (the gap restated). 2D YM is the gapped control; the thin
   atom is the "fires-without-gap" separator.

4. **(C) Counterfactual Spectral Shadow `[NG]` — reduces to operator-algebra generation /
   persistent homology + `H1`.** "Removing the candidate's traces" = the **atomless bump adversary**
   of `ch09` Thm 4 / `R3 [C3]`: every finite datum is invariant while the atom drops. So "the remaining
   observables satisfy the same relations" is exactly **weak-\* continuity / `H1`**; the algebraic reading
   ("is `O` determined by the relations") is **operator-algebra generation**; the filtration reading
   (essential vs boundary bar) is **persistent homology**, whose value-data no-go (`PERSISTENT-HOMOLOGY.md`)
   applies. "No shadow" ⇔ the relations *determine* `O` ⇔ generation — independent of any gap.

5. **The cluster is the banked separation in disguise `[C]`.** (A) → **`W`** (the edge weight / atom
   floor); (B) → **`I`** (isolation of the first excitation) *or* the gap by spectral flow; (C) →
   **`H1`** (finite-data continuity). These are precisely `W-WAVE-CHARTER §0`'s three coordinates, of
   which the gap needs a **conjunction**. So the cluster reconstructs the campaign's own decomposition;
   the still-missing machine is the conjunction — `[O]`, unmoved.

---

## 1 · The four-case battery (exact; the common control)

All three ideas are judged against the same exact models. "Atom" always means the edge atom
`f₁²=ν({inf supp ν})`; "gap" means `m₁²=inf supp ν>0`.

| model | exact construction | gap? | atom? |
|---|---|---|---|
| **free massless** | periodic `1D` Laplacian `λ_k=2−2cos(2πk/L)`, `λ₁→0` | **no** (`m₁²→0`) | **no** (atomless in continuum) |
| **2D YM SU(2)** | holonomy transfer, `H|R⟩=(e²/2)C₂(R)|R⟩`, `C₂(j)=j(j+1)` | **yes** (`Δ=3e²/8`) | **no** — tr F²=H diagonal ⇒ `f₁²=0` |
| **`wδ_m+(1−w)Unif[0,Λ]`** | `w=0.3, m=1, Λ=3` | **no** (`inf supp=0`) | **yes** (weight `0.3` at `m`) |
| **gapped atomless** | `ρ(s)=c√(s−m²)` on `[m²,m²+W]` (2-particle threshold) | **yes** (`inf supp=m²>0`) | **no** (a.c. `√`-edge) |
| *(positive control)* free massive | `δ`-atom at `m²` | **yes** | **yes** |

Exact witness numbers (§6): `m₁²(L)=5.86e-1,1.52e-1,3.84e-2,…,6.02e-4` for `L=8…256` (→0 like `4π²/L²`,
ratio `0.95→0.9999`); 2D YM `E=(0,0.375,1.0,1.875,3.0)`, `Δ=0.375`; `wδ_m+Unif`: mass in `[−d,d]` is
`(1−w)·d/Λ` (continuum reaches 0) while mass in `[m−d,m+d]→0.300` (the atom); gapped atomless
`M(d)=(2/3)d^{3/2}/Z`, vanishing like `d^{3/2}` at `d=0.2,0.05,0.01,0.003`.

**Reading.** The two mutable corners (rows 3 and the 2D-YM/gapped-atomless rows) put **property-without-gap**
and **gap-without-property** on the table simultaneously. Any "spectral identity" that is not the
literal pair `(f₁²>0 at the edge) ∧ (uniform χ)` will fall into one of the empty cells; the battery is
the acid test for *which*.

---

## 2 · (A) Spectral Lineage

**Formalisation.** A projective system `{(Aₙ, ωₙ, Iₙ, rₙ)}` with `rₙ:A_{n+1}→Aₙ` (coarse restrictions),
states `ωₙ∘rₙ=ω_{n+1}`, nested intervals `Iₙ`, and a **lineage** = a consistent family `χₙ` of candidate
edge-atoms (`rₙχ_{n+1}=χₙ`, `supp χₙ⊆Iₙ`). "Existence iff evidence survives refinements":
`ωₙ(χₙ) ≥ w` for all `n`.

**Claim A1 `[T]` (existence = compactness).** The set `E = {(χₙ) : rₙχ_{n+1}=χₙ, supp χₙ⊆Iₙ}` is a
projective limit of nonempty compact sets; by finite-intersection + Banach–Alaoglu it is **nonempty
whenever each `Eₙ≠∅`**. So "a lineage exists" asserts only that the *bottom of the spectrum* exists —
true of **every** theory with a nonempty spectrum, gap or not. It is the **inverse limit**, and carries
no quantitative content.

**Claim A2 `[T]` (the content is `W`).** The only non-vacuous version fixes `w>0`: a lineage of atoms
each of mass `≥w` at the edge. That is exactly the atom floor `f₁²≥w₀` (`W`). It is *location-sensitive*:
in free massless the finite-`L` lineage has `O(1)` mode weight but **location `m_n=λ₁(L)→0`** (`R1`) —
lineage present, floor zero.

**Battery.**

| model | lineage (uniform mass `w>0` at edge)? | gap? |
|---|---|---|
| free massless | **no** (mass in any `[0,d]`→0; location→0 with `O(1)` mode weight) | no |
| 2D YM | **no** (`f₁²=0` for tr F²; the theory is gapped for a *different* operator) | yes |
| `wδ_m+Unif` | **yes** (atom `w=0.3` at `m`) | **no** (`inf supp=0`) |
| gapped atomless | **no** (no atom) | yes |

**Verdict `[NG]`.** Lineage *existence* = inverse limit (automatic). Lineage *positivity* = `W`. The
gap additionally needs `m₁²>0` (location, `R1`) and isolation — i.e. `W∧H1`. So (A) is `W` in
projective dress, not the gap; free massless (lineage-with-zero-location) and `wδ_m+Unif`
(atom-without-floor) are the two separating adversaries. **Kill clause triggered: inverse limit.**

---

## 3 · (B) Unfinished Spectral Event

**Formalisation.** A tower of candidate projections `Pₙ` (rank-1 projectors on the lowest excitation in
`ℋₙ`) with a refinement relation: `P_{n+1}` *allowed* if it is compatible with `Pₙ` under the embedding.
**Three fates:** (i) **fires** — a consistent family with limit `P_∞≠0` (a discrete spectral edge
survives); (ii) **forbidden** — a finite `n` with no allowed refinement (rank drop, a finite certificate
of absence); (iii) **unfinished** — every finite `Pₙ≠0` and each single step is allowed, but no
consistent family / no uniform isolation (the eigenvalue approaches the threshold without crossing).

**Claim B1 `[T]` (the trichotomy is the spectral-flow trichotomy).** On the minimal tower
`H(λ)=diag(λ,1)`, candidate `P=1{eig<½}` (witness §6): `λ=0.2` → rank 1, isolation `0.8` = **fires**;
`λ=0.6` → rank 0 = **forbidden**; `λ_n=½−1/n → ½⁻` → rank 1 for every `n` but isolation
`1.00,0.75,0.625,…,5.0e-1 → ½` — **unfinished**, no finite certificate. This is exactly the flow of an
eigenvalue through a threshold: crossing (+1) / bracket / tangential touch (0-at-threshold). So (B) is a
**spectral-flow / eigenvalue-isolation** object.

**Claim B2 `[NG]` (blind to the atom weight).** `FALSIFY-FLOW.md`, re-run (§6): on the *same* fired
`H(λ=0.2)`, the edge weight is `f₁²=|<0|O|1>|² = 1.000` for `O=σₓ`, but `0.000` for `O=σ_z` and
`O=𝟙`. The event fires regardless of the overlap. **Spectrum ≠ overlap**: the event cannot certify the
atom `f₁²`.

**Battery — two honest readings.**

| model | (B1) point-projection *fires*? | (B2) uniform spectral gap? | gap? |
|---|---|---|---|
| free massless | no (rank grows / edge→threshold) | no (eigenvalue→0) | no |
| 2D YM | **yes** (isolated `R=½`, discrete) | **yes** (`Δ=3e²/8`) | yes |
| thin atom `ν_ε` | **yes** (isolated `ε`-atom vs `s₂`) | no (`ε→0`) | **no** |
| gapped atomless | **no** (`√`-edge, no point) | **yes** (`m²>0`) | yes |

**Verdict `[NG]`.** Reading B1 (point projection) fires for the thin atom yet the floor is `ε→0` —
law of a *discrete spectral point*, blind to its value = **operator-algebra generation** of a limit
projection. Reading B2 (uniform isolation) is the gap by spectral flow = **spectral-flow existence**;
by B2/`FALSIFY-FLOW` it still cannot produce the residue `f₁²`. The "new" third fate, **unfinished**, is
verbatim `H1` (no finite certificate = the weak-\* discontinuity of the atom). **Kill clauses triggered:
spectral-flow existence + operator-algebra generation (+ H1 in the third fate).**

---

## 4 · (C) Counterfactual Spectral Shadow

**Formalisation.** Algebra of observables `𝔄` with a relation set `R` (Ward, positivity, OS, moment
constraints). Candidate `O∈𝔄`. *Remove `O`'s traces*: project to `𝔄/⟨O⟩` (set `O`'s moments/correlators
to the shadow's). Property: the remaining observables **still satisfy `R`** ("a shadow exists").

**Claim C1 `[T]` (finite `R` ⇒ the shadow is the bump; = `H1`).** Exact (`R3 [C3]`, reproduced §6):
`ν_A = 0.45δ_{1.0}+0.55δ_{3.3}` (`f₁²=0.45`) vs its atom-smoothed `ν_B(h)`: `ΔM₀=ΔM₁=0` **exactly**;
`ΔM₂=3.8e-4, 1.5e-5, 1.5e-7` at `h=0.05,0.01,0.001`; all finite data agree while the edge atom
`0.45→0`. By `ch09` Thm 4 the atomless measures are weak-\* dense, so **every** finite relation `R` is
satisfied by the shadow that deletes the atom. "Shadow exists" is therefore **`H1`/weak-\* continuity** —
the relations cannot certify the atom they are asked to certify.

**Claim C2 `[T]` (algebraic / filtration reading).** "`O` is not determined by `R`" = `O`∉ the algebra
generated by `R`-expressions on the rest = **operator-algebra generation**. In the resolution filtration
(`R` = boundaries), "remove the feature, relations survive" = the feature is a **boundary** (inessential
bar) vs **essential** — i.e. **persistent homology**; `PERSISTENT-HOMOLOGY.md`'s value-data / sign-inverted
no-go applies. "No shadow" ⇔ `R` determines `O` (e.g. the exact `tr(O)=1` rank-1 control, §6) —
determinacy, not spectrum.

**Battery.**

| model | shadow exists (relations survive)? | gap? |
|---|---|---|
| free massless | yes (free relations linear; nothing to delete) | no |
| 2D YM | yes (tr F² conserved; holonomy area/Casimir relations intact) | yes |
| `wδ_m+Unif` | **yes** (delete the `w`-atom → `Unif` still satisfies finite `R`) | **no** |
| gapped `tr(O)=1` control | **no** (`O` determined; deletion breaks `R`) | yes |

**Verdict `[NG]`.** Shadow-existence = `H1`/weak-\* (the bump adversary); its negation = generation.
The property holds without a gap (`wδ_m+Unif`) and fails without a gap's absence (the determined `O`),
so it is **orthogonal to the gap** and reduces to **operator-algebra generation / persistent homology +
H1**. **Kill clauses triggered: operator-algebra generation + persistent homology (+ H1).**

---

## 5 · The cluster is the three banked coordinates, not a sixth species `[C]`/`[O]`

Put the three reductions side by side against `W-WAVE-CHARTER §0` (`gap = W ∧ H1`, or via isolation `I`):

| idea | what its property *is*, exactly | banked coordinate | kill object |
|---|---|---|---|
| **(A) Lineage** | consistent family of atoms of **uniform mass `≥w`** | **`W`** — the edge atom floor `f₁²≥w₀` | **inverse limit** |
| **(B) Event** | isolated discrete edge surviving refinements | **`I`** — isolation of first excitation (or the gap, by flow) | **spectral-flow existence / operator-algebra generation** |
| **(C) Shadow** | finite relations invariant under deleting the atom | **`H1`** — the weak-\* ceiling/continuity | **operator-algebra generation / persistent homology** |

The campaign already knows (`W-WAVE-CHARTER §0`, `Q-WAVE-SYNTHESIS §4`, `U6-BINDING-LOGIC`) that `W`,
`I`, `H1` are **independent outputs** and the gap requires a **conjunction** (`W∧H1` or `I∧(ρ²≤C²)`).
The spectral-identity cluster reconstructs exactly that decomposition: each idea isolates one factor and
names it in a new language (projective system / projection tower / counterfactual), and **none supplies
the other factor**. The battery's two separating corners (`wδ_m+Unif` for A/C; 2D YM and gapped-atomless
for A; thin atom for B) are the exact witnesses that the *property alone* is not the gap.

**Exact missing step `[O]`** — unchanged: a **uniform-in-`(L,a)` edge-atom floor with isolation**
(`f₁²≥r₀` **and** `m₁²`-isolation) on a gauge-invariant carrier, conjoined with the `H1` ceiling. That is
`H*` (`Q-WAVE-SYNTHESIS §4`). No member of the cluster is that machine; each is one of its inputs
restated. **Flip condition:** exhibit a refinement-stable, non-compactness functional that *forces* the
edge floor `≥r₀` at `β=2.5` (and `=0` on the free field) without presupposing the location/isolation —
none exists, and `R3`'s Stone–Weierstrass bound (`sup_ν|E−g|≥f₁²/2`) says a finite-data one cannot.

---

## 6 · Witness — `research/v3/w7/sl_test.py` (exact; quoted)

`py -3.12 -I sl_test.py` (exit 0; numpy 2.4.6 / scipy 1.18.1; one job, no downloads). Abridged output:

```
[CASE 1] free massless 1D:  L=8 m1^2=5.86e-1 ; L=256 m1^2=6.02e-4 ; ratio to 4pi^2/L^2 -> 0.9999
          edge atom read M(d)=0.2,0.05,0.01,0.003 for d=0.2,0.05,0.01,0.003  -> 0 like d  (no atom; continuum proxy rho=1)
[CASE 2] 2D YM SU(2): E=(0, 0.375, 1.0, 1.875, 3.0); gap Delta=0.375=3e^2/8
          trF^2=H diagonal => <0|trF^2|first-excited>=0 => f1^2=0 ; Var(H)_vac=0  (GAPPED, ATOMLESS)
[CASE 3a] thin atom nu_e: eps=0.5..0.001 ; edge_atom w=e/2lam: 0.25 -> 5e-4 ;
          M0=1, <s>->3, chi=1/<s> bounded  => finite data blind to the vanishing atom (H1)
[CASE 3b] w delta_m+Unif: mass near m ->0.300 (atom) ; mass near 0 -> 0 like density*d ; inf supp=0 (NO gap)
[CASE 4] gapped atomless sqrt-edge: M(d)=1.7e-2,2.2e-3,1.9e-4,3.2e-5 -> d^{3/2} ; inf supp=m^2>0 (GAP)
[A] lineage: free massless location m_n= 5.86e-1..6.02e-4 -> 0 with O(1) weight ; 2D YM f1^2=0
[B] three fates H(lam)=diag(lam,1): lam=0.2 FIRES(iso 0.8) ; lam=0.6 FORBIDDEN ;
          lam_n=0.5-1/n: rank 1 all n, isolation 1.00,0.75,...,5.0e-1 -> UNFINISHED (spectral-flow)
    FALSIFY-FLOW: same fired H(0.2): f1^2 = 1.000 (sigma_x), 0.000 (sigma_z), 0.000 (identity)
[C] bump adversary: dM0=dM1=0 exact ; dM2=3.8e-4,1.5e-5,1.5e-7 ; edge atom 0.45->0 (Thm 4)
    2x2 SUMMARY: atom and gap independent (row 3b: atom,no gap ; rows 2/4: gap,no atom)
```

`sl_test.py` computes everything by exact closed form / dense eigendecomposition of `2×2` and `L×L`
matrices; no Monte Carlo, no third-party data.

---

## 7 · Ledger

| item | statement | label |
|---|---|---|
| four-case battery | free massless / 2D YM / `wδ_m+Unif` / gapped-atomless, exact | **[T]+[V]** |
| atom ⟂ gap | separating corners: `wδ_m+Unif` (atom,no gap), 2D YM & gapped-atomless (gap,no atom) | **[T]+[V]** |
| 2D YM tr F² atom | `tr F²=H` conserved ⇒ `f₁²=0` while gapped `Δ=3e²/8` | **[T]** (`W`-charter acid #2)+**[V]** |
| **(A) existence = inverse limit** | nonempty projective limit = finite-intersection / Banach–Alaoglu, automatic | **[T]** |
| **(A) content = `W`** | uniform-mass lineage = edge atom floor; free massless location→0 (`R1`) | **[T] → [NG]** for gap |
| **(B) trichotomy = spectral flow** | fires / forbidden / unfinished = crossing / bracket / threshold-touch | **[T]** |
| **(B) blind to `f₁²`** | same fired `H(λ)`: `f₁²=1` (`σₓ`) vs `0` (`σ_z`,`𝟙`) | **[T]+[V]** (FALSIFY-FLOW) |
| **(B) readings = generation / flow** | point-projection = operator-algebra generation; uniform gap = spectral-flow existence; unfinished = H1 | **[T] → [NG]** for gap |
| **(C) shadow = bump = H1** | finite relations invariant under atom deletion (Thm 4 / R3 [C3]) | **[T]+[V]** |
| **(C) negation = generation** | "no shadow" ⇔ relations determine `O`; filtration reading = persistent homology | **[T] → [NG]** |
| **cluster = `W`/`I`/`H1`** | the three ideas are the three banked coordinates, each one factor | **[C]`/`[O]** |
| **gap = `W∧H1` (or `I∧…`)** | conjunction; the missing machine | **[O]** · **0%** |
| Clay | unbuilt | **[O]** · **0%** |

---

## 8 · Sharpest line

**All three spectral-identity engines pass every exact model and still cannot tell a mass gap from its
absence.** The four-case battery is a clean 2×2: `wδ_m+(1−w)Unif[0,Λ]` has an atom of weight `0.3` and
`inf supp=0` (property, no gap), while 2D YM (`tr F²=H` ⇒ `f₁²=0`) and the `√`-threshold measure
(`inf supp=m²>0`, atomless) have the gap and not the atom — so atom and gap are independent, the campaign's
central dichotomy, exactly. Each idea then collapses onto a *named* banked object: **(A)** reduces to the
inverse limit (existence = Banach–Alaoglu, automatic) plus `W`; **(B)**'s three fates are the spectral-flow
trichotomy, blind to the vertex `f₁²=⟨0|O|G⟩` (`f₁²=1` or `0` on the same crossing), i.e. operator-algebra
generation or the gap renamed; **(C)**'s "shadow" is the atomless bump of `ch09` Thm 4 — finite data is
invariant while the atom vanishes — i.e. `H1`/weak-\* continuity, or generation. The cluster is precisely
the campaign's three separated coordinates `W`, `I`, `H1` in projective / spectral / counterfactual dress;
the gap needs their **conjunction**, which is `H*` — the still-missing machine. **No sixth species; no QED;
Clay 0%.**

---

*Refs.* **Internal:** `CLAUDE.md` §1,§5 (acid test; the wall; missing machine); `research/v3/R3-ATOM-FUNCTIONAL.md`
(atom functional, stone–Weierstrass bound, [C1]–[C5]); `docs/paper-v2/ch09-weakstar.tex` (Thm 4, Lemmas 1–4,
weak-\* topology); `research/v3/X2-FACTORIZATION.md` (independent-axes pattern, sign audit); `W-WAVE-CHARTER.md`
§0 (`W`,`I`, gap = conjunction); `Q-WAVE-SYNTHESIS.md` §4 (H1 ceiling, `H*`); `research/unconditional/FALSIFY-FLOW.md`
(spectrum ≠ overlap); `FALSIFY-TOPOLOGICAL.md` (atom = weak-\* discontinuity); `PERSISTENT-HOMOLOGY.md`
(value-data, sign-inverted); `JONES-SUBFACTOR.md` (operator-algebra, mass-blind); `INVENT-NEW-INDEX.md`
(rigidity no-go); `D6-DECISION.md` (exact solvable-model playbook). **Witness:** `research/v3/w7/sl_test.py`.
No `[C]` upgraded to `[T]`; no fabricated theorem. **No QED; Clay 0%.**
