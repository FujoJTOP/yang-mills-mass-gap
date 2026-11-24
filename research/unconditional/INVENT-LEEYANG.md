# INVENT-LEEYANG — a complex-β analytic-structure functional: a genuinely NEW *type* (phase-sensitive, non-order-parameter, non-Markov, non-topological), but acid-dirty and atom-blind — the fix (IR-only zeros) is `H1`

> **VERDICT (no QED).** A Lee–Yang / Fisher functional — the distance `D(β_0)` from `β_0=2.5` to the
> nearest singularity of `Z(β)` in the complex `β` plane — is the **first candidate that satisfies the type
> spec** of the roof round: **phase-sensitive, non-order-parameter, non-Markov, non-topological.** That is a
> real positive: it proves the species can have this type. **But it fails the content tests:** `(i)` acid
> **FAIL** — the nearest singularity is the **kinematic `β=0` branch point**, present in every theory, so
> `D(free massless)=β_0>0` (it fires on the gapless field); `(ii)` **atom-blind** — `D` is a *gap-detector*
> (`D=β_0+m²/λ_max`, monotone in `m²`, independent of `f₁²`), so it fires on the gapped-atomless `:φ²:`;
> `(iii)` it gives **no** renormalon-free `χ_lattice ≲ a²` bound — the analyticity radius is
> **kinematic/UV-dominated**, and isolating the IR Fisher zeros from the `β=0` singularity *is* the
> UV/IR separation at the crossover `= H1`. **New type, not new content. Clay 0%.**

---

## §0 · Labels, iron rules, acid test

**[T]** 真定理 · **[C]** 条件定理 · **[O]** 精确缺失步 · **[NG]** no-go for the named class · **[V-num]** new numeric.
**IL-i** non-positivity (fail on free massless) · **IL-ii** nonabelian · **IL-iii** no small parameter.
**Acid** (`SRF-FREE-FIELD-ACID`): the functional must be `≤0` (`=0`) on the **free massless** field, **and
discriminatingly** `≤0` on the gapped-atomless `:φ²:` composite (`f₁²=0`).

---

## §1 · The functional

For `Z(β)` the lattice partition function (complex `β`), define

```
D(β_0) := dist(β_0, Sing Z)        (nearest singularity of log Z in Cβ)
ρ(β_0) := density of Fisher zeros on a disc of radius r around β_0,   r→0
```

`D` (resp. `ρ`) is **phase-sensitive**: at a real transition the zeros pinch the real axis (`D→0`,
`ρ→∞`); in a single phase (crossover) they stay at finite distance (`D>0`). This is exactly the
"phase-sensitive but **not** an order parameter" type the roof round asked for (`ROOF-TRANSFER` B2: no
order parameter because no transition — but the *analytic structure* still changes at a crossover).

**Claim to test:** the 0⁺⁺ gap is controlled by `D(2.5)` — a finite distance to the nearest zero ⇒ gap.

---

## §2 · Type classification — the round's genuine positive `[V-here]`

Against the four *dead* paradigms of `ROOF-SYNTHESIS`:

| paradigm | dead because | Lee–Yang/Fisher |
|---|---|---|
| order parameter | no transition to condense (B2) | **not** an order parameter (no condensate) ✓ |
| topological index | rigidity (constant on β) | **not** an index (can vary continuously) ✓ |
| Markov functional | Theorem BE (homogeneous/zero/sign-indefinite) | **not** a gradient-flow functional ✓ |
| 't Hooft dual | 0-form impossible (linking/center) | **not** a dual (no loop algebra needed) ✓ |

⇒ **Lee–Yang/Fisher is a genuinely NEW species in *type*.** It is the first candidate outside all four
dead paradigms, and it is phase-sensitive by construction. *This is the round's real yield:* the species
the roof needs **can** be an analytic-structure functional — the type exists.

---

## §3 · Acid test — FAIL `[V-num]`

Witness `w7/invent_leeyang.py` (free Gaussian, `Z(β)=∏_{λ>0}(m²+βλ_k)^{-1/2}×(zero-mode)`):

| model | nearest singularity | `D(β_0=2.5)` | fires? |
|---|---|---|---|
| **free massless** (`m²=0`, gap `=0`) | kinematic branch point `β=0` | **`2.5000`** | **YES (wrong)** |
| `:φ²:`-like massive (`m²=0.1`, `f₁²=0`) | `β=−m²/λ_max<0` | `2.5250` | YES (wrong) |
| `:φ²:`-like massive (`m²=1`) | `β=−0.25` | `2.7500` | YES (wrong) |
| 1D Ising (gapped) | zeros on `Im β` axis (`Re β=0`) | `2.5000` | YES |

**Root cause `[T]`:** the nearest singularity of any lattice `Z(β)` is the **kinematic `β=0` branch
point** (the free/UV structure), which is **present in every theory, gapped or not.** So `D(free
massless)=β_0>0` — the functional **fires on the gapless field** ⇒ **acid FAIL.** A finite analyticity
radius is a property of *positivity of the transfer matrix* (all theories), **not** of the IR gap.

---

## §4 · Gap-detector, not atom-certifier `[NG]`

`D = β_0 + m²/λ_max` for the massive Gaussian — **monotone in `m²`, independent of `f₁²`.** So `D` is a
**gap-correlate**, and it fires on the gapped-atomless `:φ²:` (`f₁²=0`) exactly as it fires on a gapped
theory with an atom. This is the known "gap-detector ≠ atom-certifier" distinction: `D` cannot tell a
KL atom from a two-particle continuum edge. **`[NG]`** for the atom-certifying class.

---

## §5 · No renormalon-free `χ_lattice ≲ a²` bound `[O]`

The roof's missing strut is a scheme-independent `χ_lattice ≲ a_ℓ²` (`ROOF-UNIFORM-CHI`) /
`σ_w(2.5)>c₀` (`ROOF-PEIERLS-0FORM`) / `R_vac` bound (`ROOF-TRANSPORT-SCALAR`). The Lee–Yang `D` does
**not** provide it:

- The analyticity radius is set by the **kinematic `β=0`/UV** singularity, not the physical condensate;
- The **physical (IR) Fisher zeros** in 4D gauge are **not isolated** — the zero-density `ρ` is
  **UV-divergent** and dominated by lattice artifacts, so no scheme-independent `a→0` limit exists;
- Anything that *would* isolate the IR zeros (subtract the kinematic part, or restrict to the IR
  effective action) already assumes the IR structure — **circular.**

**Exact failing step `[O]`:** separating the **IR Fisher zeros** from the **kinematic `β=0` singularity**.
That separation is *precisely* the UV/IR split at the crossover — **`H1` verbatim** in the analytic-
structure coordinate.

---

## §6 · Verdict

- **New *type*: YES `[V-here]`.** Lee–Yang/Fisher is the first candidate satisfying the roof round's type
  spec (phase-sensitive, non-order-parameter, non-Markov, non-topological). The species *can* be an
  analytic-structure functional — **not** an order parameter, index, or gradient-flow functional.
- **New *content*: NO.** Acid FAIL (kinematic `β=0` fires on free massless); atom-blind (fires on
  `:φ²:`); no renormalon-free `χ` bound (IR-zero isolation `=H1`).
- **The genuine handle:** the fix is an **IR-only** analytic-structure functional — one that sees only the
  zeros *created by the interaction*, not the kinematic ones. That would be acid-clean (the free field has
  no interaction zeros) and atom-sensitive (the interaction zeros know about `f₁²`). **Constructing the
  IR-only zero functional is equivalent to the UV/IR separation at the crossover = `H1`.**

**Scorecard: type 4/5 (new species type found), content 0/5 (acid-dirty, atom-blind, kinematic).**
**Clay 0%; no QED.**

*Refs:* `ROOF-SYNTHESIS.md` (the four dead paradigms; B2 single-phase), `ROOF-UNIFORM-CHI.md`,
`ROOF-CONDENSATE-FLOOR.md`, `PROFILE.md`; witness `w7/invent_leeyang.py`.
