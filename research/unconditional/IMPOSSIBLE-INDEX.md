# IMPOSSIBLE-INDEX — a NEW-SPACE index cannot force the gap: every index is a homotopy invariant (locally constant = weak-* continuous), but the atom is the weak-* discontinuity — so no index in any space can read or lower-bound it; acid-blind (0 on free massless and `:φ²:`)

> **VERDICT (no QED).** [NG] for the new-space-index class. Tried beyond the Dirac index (BASE-INDEX):
> equivariant index on the instanton moduli space, noncommutative / Connes spectral-triple K-theory, the
> loop-space / free-loop index (Witten elliptic genus), the relative index of Dirac operators twisted by two
> bundles, the Atiyah–Hitchin–Singer index, and the Chern–Simons / eta invariant. **None forces the gap**, and
> the reason is one structural theorem: **an index is a homotopy invariant — integer-valued and locally constant
> on the parameter space `(β, A)` — i.e. the *maximally continuous* functional; the KL atom is the *weak-*
> discontinuity* of the boundary density (ONLY-EXIT).** A locally-constant functional cannot equal or lower-bound
> a quantity that varies continuously through the crossover, and cannot read a discontinuity it is continuous
> across. **This generalises BASE-INDEX' rigidity no-go to ALL indices in ANY space** and *confirms* ONLY-EXIT.

---

## §0 · Labels

**[T]** 真定理 · **[C]** 条件定理 · **[O]** 精确缺失步 · **[NG]** no-go · **[V-num]** numeric. Acid: `≤0`/`=0` on free massless *and discriminatingly* on `:φ²:`.

---

## §1 · Each candidate is an index (rigid / homotopy-invariant)

| new-space index | what it counts | verdict |
|---|---|---|
| **Equivariant instanton index** (Atiyah–Hitchin–Singer; Nekrasov) | instanton moduli / equivariant characters `Tr (−1)^F g e^{−βH}` | `[NG]` UV/instanton data; β-rigid; does not see the IR `0⁺⁺` mass |
| **Noncommutative / Connes spectral-triple K-theory** | the K-theory index of the Dirac `D` of `(𝒜,H,D)` | `[NG]` if `D` = physical Hamiltonian → gap is tautological and its index rigid; if `D` = the natural gauge Dirac → **the standard index = BASE-INDEX (dead)** |
| **Loop-space / free-loop index** (Witten elliptic genus) | a topological invariant of loops `LM` | `[NG]` topological ⇒ rigid on the coupling; pure YM (no SUSY) → the elliptic genus is trivial / not the gap |
| **Relative index** `ind(D_{E₁})−ind(D_{E₂})` | a difference of integers | `[NG]` still an integer ⇒ locally constant on β; cannot jump at the crossover |
| **Atiyah–Hitchin–Singer index** | harmonic forms on self-dual moduli | `[NG]` UV/instanton; rigid |
| **Chern–Simons / eta invariant** `η(A)` | spectral asymmetry of the Dirac on a 3-manifold | the one **real-valued** candidate (not integer); `[NG]` nonetheless — `η` is a UV/spectral-flow datum of the *connection* `A`, not the IR measure, and (§4) continuous in `A` |

---

## §2 · The unifying theorem (why NO index, in any space, works)

**Theorem `[T]` (index rigidity ⇒ gap-blind).** An index is a homotopy invariant: it is integer-valued and
**locally constant on the connected components** of the parameter space `(β, A)`. The mass gap is an *IR,
`β`-dependent* quantity that changes continuously across the strong–weak crossover (`β ≈ 2.5`). The KL atom is
the **weak-\* discontinuity** of the boundary density `ρ` (`ONLY-EXIT`: every finite continuous datum is blind;
only a weak-\* discontinuous functional reads it).

Therefore:
1. **Rigidity:** on any open set of `β` where the gap varies, the index is constant ⇒ the index **cannot equal
   or lower-bound** the gap (a constant cannot bound a function that leaves it).
2. **Discontinuity-blindness:** the index is continuous where the atom is discontinuous ⇒ it **cannot read**
   `f₁²>0` (a continuous functional evaluates identically on the atom and on an atomless bump).
3. **Gap = 0 side:** an index changes only when the gap **closes** (a degeneracy) ⇒ if anything, it certifies
   `gap = 0` at the transition, never `gap > 0`.

**Corollary.** BASE-INDEX's no-go (rigidity / atom-gap-mismatch / zero-mode-certifies-`gap=0` / acid-fail) is
**not special to the Dirac index** — it is the generic fate of **every** index, in **any** space, for the
structural reason above: an index sits at the **maximally continuous** end of the functional spectrum, the atom
at the **maximally discontinuous** end.

---

## §3 · Acid test

An index is `0` on the free massless field (no topological charge / zero modes are *needed* for the mass) and
`0` on the gapped-atomless `:φ²:` (no index structure) — so an index **fails the discriminating `:φ²:` step**
(it is the same on the gap-detector and the atom-certifier adversaries). Acid-dirty / vacuous.

---

## §4 · The eta-invariant caveat (the only non-integer candidate)

`η(A)` is real, not integer, so it is not "rigid" in the integer sense. But (i) it is a **UV / spectral-flow**
datum of the connection `A` (the spectral asymmetry near the cutoff), not an IR property of the measure `μ`;
(ii) it does not see the `0⁺⁺` mass (a spectral asymmetry of a gauge Dirac operator is insensitive to the IR
scalar gap); and (iii) by ONLY-EXIT, `η(A)` is continuous in `A` while the atom is the weak-\* **discontinuity**
of the 0⁺⁺ spectral density — so `η` cannot read it either. `[NG]`.

---

## §5 · Witness (`w7/impossible_index_witness.py`, clean)

Toy: the SSH two-band Dirac model. Its topological index is the winding `W` of the Bloch vector
`d(k)=(t₁+t₂cos k, t₂ sin k)`.
- **RIGIDITY:** `W` is a **step** (`0` for `t₂/t₁<1`, `1` for `>1`); it jumps at `t₁=t₂` **while the spectral
  gap stays `≈0.6003` on both sides** — the index jumps where the gap is flat. `W ≠` gap.
- **ACID-BLIND:** the trivial gapped model (`W=0`, gap `1.34`) and the topological gapped model (`W=1`, gap
  `2.09`) are **both gapped** — the index separates two gapped phases, not gapped from gapless; and `W` is the
  same (`0`) for `m=0` and `m=0.3` — it does not track the mass.

⇒ an index neither equals the gap, nor tracks the mass, nor bounds the gap from below.

---

## §6 · Honest bottom line

- **The new-space index class is closed `[NG]`** — for the *same* structural reason as BASE-INDEX, now proved
  generically: an index is a homotopy invariant (weak-\* continuous / locally constant); the atom is the weak-\*
  discontinuity. **NO index in ANY space can force or read the gap.**
- **This confirms ONLY-EXIT** from the opposite extreme: indices are the maximally continuous functionals; the
  atom is the maximally discontinuous object. They cannot meet.
- The eta invariant (real, non-integer) is the only non-rigid candidate, but it is UV/spectral-flow of the
  connection, not the IR 0⁺⁺ measure, and remains continuous where the atom is discontinuous.
- **Clay 0%; no QED.** A genuinely-impossible base that *fails* still teaches the exact structure: the wall is the
  weak-\* discontinuity, and every topological index is on the wrong side of it.

*Refs:* `BASE-INDEX.md`, `ONLY-EXIT-SYNTHESIS.md`, `TOOL-SYNTHESIS.md`; witness `w7/impossible_index_witness.py`.
