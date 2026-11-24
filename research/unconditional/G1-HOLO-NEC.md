# G1-HOLO-NEC — bulk NEC cannot force the IR wall: it *is* the holographic a-theorem, which *admits* an IR AdS throat (gapless CFT). The wall is an input (=H1), not a bulk-positivity consequence

> **VERDICT [NG] (no QED, Clay 0%).** The question: does a bulk null-energy / c-theorem argument FORCE the
> pure-YM dual to have an IR wall (no AdS IR fixed point), thereby excluding a gapless conformal IR? **No.**
> The bulk NEC is, via Einstein, exactly the holographic **c-theorem** (Freedman–Gubser–Pilch–Warner 1999), and
> the c-theorem **admits an IR AdS throat** — a flow from UV-AdS to IR-AdS satisfies NEC with both endpoints AdS
> ⇒ a gapless conformal IR is NEC-compatible. A gapped wall is *also* NEC-compatible. So the NEC does **not**
> distinguish IR-CFT from wall. Forcing the wall needs an *additional* input — a dilaton potential that forces a
> singularity (the hard/soft wall, put in by hand to model confinement), or "the pure-YM dual has no IR CFT"
> (= G1 itself). Hence the NEC argument either **allows the CFT IR** (no exclusion) or **assumes the wall = H1**.

---

## §0 · Labels

**[T]** 真定理 · **[C]** 条件定理 · **[O]** 精确缺失步 · **[NG]** no-go · **[V-num]** new numeric. Acid: `≤0`/`=0` on free massless *and discriminatingly* on `:φ²:`.

---

## §1 · The NEC is the holographic a-theorem (and that is the problem)

Domain-wall ansatz `ds² = e^{2A(r)} η_{μν}dx^μdx^ν + dr²` (bulk dimension `d+1`). Einstein + the bulk
**null energy condition** `(T_{rr}+T_{tt}) ≥ 0` are *equivalent* to
`A''(r) ≤ 0` (A concave), i.e. the **holographic c-theorem**: the c-function
`c(r) ∝ 1/|A'(r)|^{d-1}` is **monotone decreasing toward the IR** (`[T]`, FGPW 1999; the NEC is the bulk
image of the boundary a-theorem).

**Key consequence:** the c-theorem is a *monotonicity*, **not** an exclusion of IR fixed points. It permits **both**
endpoints:
- **IR AdS throat (gapless CFT):** `A'` runs from `1/L_UV` (UV AdS) to `1/L_IR > 0` (IR AdS), *decreasing* ⇒
  `A'' ≤ 0` ⇒ NEC satisfied, with **both endpoints AdS**. So a gapless conformal IR is **NEC-compatible**.
- **Gapped wall:** `A' → 0` (a cap / singularity), also `A'' ≤ 0` ⇒ NEC satisfied.

Witness (`w7/g1_holo_nec.py`, exact) confirms both: a smooth UV-AdS→IR-AdS flow has `A'' ≤ 0`, `c` decreasing,
and both `A' > 0`; a wall flow has `A'' ≤ 0` with `A'(IR) → 0`. **Both NEC-satisfying — NEC does not decide.**

---

## §2 · Why NEC cannot force the wall

By the holographic dictionary, bulk NEC ⟺ boundary a-theorem ⟺ unitarity + monotonicity. The a-theorem
**allows an IR CFT** (`a_IR > 0`). Therefore **no bulk NEC argument can exclude an IR AdS throat** — it would be
excluding something the a-theorem permits. To force the wall you would need more than NEC:
- a **dilaton potential that forces a singularity** at finite `r` (the hard/soft wall) — but that potential is a
  **modeling choice** put in precisely to model confinement, not derived from "pure gauge, no matter";
- or the statement **"the pure-YM dual has no IR AdS throat"** — which **is (G1)** (the gap) itself.

Also: "pure gauge / no matter in the *bulk*" would give **pure gravity = pure AdS = a CFT** — the *opposite* of a
wall. So the wall is never a positivity consequence; it is either modelled in, or assumed.

**Holographic caveat (compounds it):** pure YM (non-SUSY, `N=3`, strong coupling) has **no controlled gravity
limit** (no large-N, no near-AdS, no SUSY). So even the bulk NEC setup is heuristic; the NEC argument inherits
this incompleteness.

---

## §3 · Verdict and the exact step

- **`[NG]`**: bulk NEC / c-theorem **does not force the wall** and **does not exclude a gapless conformal IR**.
- **`[T]`**: NEC ⟺ holographic a-theorem (monotonicity, admits IR AdS).
- **`[O]`**: forcing the wall requires an input beyond NEC — a symmetry-breaking dilaton potential (modeling
  choice) or "no IR CFT for pure YM" = **(G1) = H1** in the holographic/bulk coordinate.
- **Acid**: the NEC argument "succeeds" only by *assuming* the wall — exactly the circularity `STAR-COUNTEREXAMPLE`
  found for the β-function route. The holographic base re-derives `H1`; it does not go around it.

**Clay = 0%; no QED.** Files: `G1-HOLO-NEC.md`, `w7/g1_holo_nec.py` (+ `.log`).
