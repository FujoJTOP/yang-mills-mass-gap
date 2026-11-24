# G1-A-C-BOUNDS — central-charge bounds do NOT exclude a gapless CFT IR of pure YM: the UV point is *allowed* (the free-vector value a/c = 31/18 is the inclusive upper edge), and the a-theorem is an inequality ⇒ (G1) is not proven by `a,c`

> **ERRATUM (2026-11-24).** An earlier version of this file reported `a/c = 31/36` and called it "strictly
> interior". That number came from an inconsistent normalization (`a` divided by 4, `c` divided by 2). With the
> standard normalization (per real scalar `a=1/90`, `c=1/30`) the free gauge vector has
> `a = 62/90`, `c = 12/30 = 2/5`, hence `a/c = 31/18 ≈ 1.722`, which is exactly the **upper edge** of the
> conformal-collider band `1/3 ≤ a/c ≤ 31/18` (Hofman–Maldacena; Hofman–Li–Meltzer–Poland–Rejón-Barrera). The
> conclusion of the file is **unchanged** — the bound is inclusive and the a-theorem is a ceiling, so a CFT IR is
> still allowed — but the phrase "strictly interior" is withdrawn. Recorded in `docs/paper-v2` Appendix D.

> **VERDICT (no QED).** **(G1) is NOT excluded by the central charges.** Pure YM's UV anomaly
> point `a/c = 31/18 ≈ 1.722` lies **on the allowed upper edge** of the unitary band `1/3 ≤ a/c ≤ 31/18`
> (equality for a free vector, and the bound is inclusive); the a-theorem `a_IR ≤ a_UV` is an *inequality*,
> so a hypothetical IR CFT has a whole 2-D allowed region `a_IR ∈ (0, a_UV]`, `a_IR/c_IR ∈ [1/3, 31/18]`.
> A gapless CFT IR is **fully compatible** with every central-charge constraint. The only CFT IR these
> bounds would kill is the **trivial** `g*=0` one (`a_IR = a_UV`, i.e. *no flow* at all) — but that is the
> **AF / not-IR-free** fact (`STAR-NO-IRFREE`), not a central-charge bound. **`[NG]` for the a/the-c route;
> (G1) does not follow.**

---

## §0 · Labels

**[T]** 真定理 · **[C]** 条件定理 · **[O]** 精确缺失步 · **[NG]** no-go for the named route · **[V-num]** numeric.

---

## §1 · The UV data `[T]`

Free-field trace-anomaly coefficients in the standard normalization (per real scalar `a = 1/90`, `c = 1/30`;
per free vector `a = 62/90`, `c = 12/30 = 2/5`), summed over the `n_V = N²−1` gauge vectors:

```
a_UV = (62/90)(N²−1) ,   c_UV = (12/30)(N²−1) = (2/5)(N²−1) ,
a_UV/c_UV = 31/18 ≈ 1.7222  (all N)
SU(2): a = 2.0667, c = 1.2000 ;  SU(3): a = 5.5111, c = 3.2000 ;  SU(4): a = 10.333, c = 6.000
```

## §2 · The constraints, and why they fail to bite

- **Unitarity / conformal-collider bound** `[T]`: `1/3 ≤ a/c ≤ 31/18 ≈ 1.722`. The UV point `31/18` sits
  **exactly on the inclusive upper edge** — the bound is a closed interval, and equality is attained by the free
  vector itself, so it is *compatible* with an IR CFT for every `N`.
- **a-theorem** `[T]` (Komargodski–Schwimmer): `a_IR ≤ a_UV`. This is a **ceiling**, not an equality. A flow to
  a CFT has `a_IR < a_UV` — allowed. It does **not** forbid a CFT IR.
- **Allowed IR-CFT region** (SU(3) example): `a_IR ∈ (0, 5.511]`, `a_IR/c_IR ∈ [0.333, 1.722]` — a whole 2-D
  set. A free massless scalar (`a,c = 1/360, 1/120`, `a/c = 1/3`) satisfies both; an interacting CFT with
  `a_IR ~ O(1)` likewise (large-`N` gauge theories approach the free-vector ratio from below).
  **Nothing selects a point or forbids the set.**

## §3 · The one thing the bounds *do* see — and it is not (G1)

A CFT IR with `a_IR = a_UV` (the **trivial** `g*=0` fixed point, *no flow*) is excluded — but that is the
**asymptotic-freedom / not-IR-free** statement (`STAR-NO-IRFREE`, `[T]`), not an anomaly-coefficient bound. The
central charges are **blind to the interacting-vs-trivial distinction** among CFT IRs: they assign every unitary
CFT an allowed `(a, a/c)` and do not isolate pure YM.

---

## §4 · Honest bottom line

- **`[NG]`**: the `a`, `c`, `a/c` bounds do **not** exclude a gapless CFT IR of pure YM. `(G1)` ("pure YM has no
  IR fixed point") is **not** proven by central charges.
- **Why (consistent with `STAR-SYNTHESIS`)**: the anomaly data are **exact and phase-blind** — they satisfy every
  unitary CFT, so they can constrain *which* CFT (via `a`, `a/c`) but cannot forbid *being* one. The exclusion of
  `N_f=0` must come from the **dynamical** β-function, not from exact anomaly inequalities — exactly the
  `[O]`-step the `(★)` chain already isolated.
- **What this route adds**: it **rules out one hope** — that a cleaner use of the a-theorem could close `(G1)` —
  and confirms the difficulty is dynamical, as `STAR-NO-BZ`/`STAR-COUNTEREXAMPLE` found. `Clay = 0%`, no QED.

*Cross-refs:* `STAR-SYNTHESIS.md`, `STAR-NO-BZ.md`, `STAR-NO-IRFREE.md`, `STAR-COUNTEREXAMPLE.md`.
Witness: `w7/g1_ac_bounds.py` (runs clean; it computes the **ratio** `a/c`, so the ratio is unchanged by a common
normalization — the earlier 31/36 vs 31/18 discrepancy was purely in the quoted per-field values, not in the ratio
logic, and the script's outputs should be read as consistent with 31/18).
