# G1-ACID — acid test of (G1): "single coupling" does NOT forbid a β-function zero; the conformal window IS a single-coupling AF gapless-conformal theory, so (G1) for N_f=0 is purely dynamical (= the gap), not kinematic

> **VERDICT (no QED, Clay 0%).** The acid test **does not refute (G1)** — it **kills the kinematic shortcut**
> and **confirms (G1) is purely dynamical**. A **single-coupling** asymptotically-free gauge theory **can** be
> gapless-conformal: the Banks–Zaks zero is determined by the β coefficients (matter content), not by the number
> of couplings. `SU(3)` with `N_f=9…16` (one coupling per `N_f`) has a real BZ zero `g*>0` — a genuine
> single-coupling AF gapless-conformal theory. For `N_f=0` the two-loop β has no zero (`b0,b1>0`), but that is a
> **matter-content fact**, **not** a consequence of "single coupling" or "AF". Hence **"single coupling ⇒ no
> zero" and "AF ⇒ gapped" are FALSE [NG]**, and **(G1) for pure YM** ("`N_f=0` has no nonperturbative IR fixed
> point") is a genuine **dynamical** statement = the gap (`H1` in the β/RG coordinate).

---

## §0 · Labels

**[T]** 真定理 · **[C]** 条件定理 · **[O]** 精确缺失步 · **[NG]** no-go · **[V-num]** numeric.

---

## §1 · The two-loop fact (single coupling, fixed `N_f`)

`β(g) = −b0 g³ − b1 g⁵`, `SU(3)`: `b0 = (33−2N_f)/3`, `b1 = 102 − (38/3)N_f`.
A Banks–Zaks zero `g*² = −b0/b1 > 0` exists iff `b0>0, b1<0`.

Witness (`w7/g1_acid.py`, exact rational arithmetic):

| `N_f` | `b0` | `b1` | `g*` |
|---|---|---|---|
| 0 | 11.00 | 102.00 | **none** (β<0 ∀g) |
| 8 | 5.67 | 0.67 | none |
| **9** | 5.00 | −12.00 | **0.6455** |
| **12** | 3.00 | −50.00 | **0.2449** |
| **16** | 0.33 | −100.67 | **0.0575** |

**For fixed `N_f` there is ONE coupling**; the BZ zero is a single number `g*(N_f)` — a zero is **not**
forbidden by "single coupling". The window `N_f ≳ 8.05` (perturbative) / `~9–16` (witness) is a **real
single-coupling AF gapless-conformal theory** (`[C]`, perturbatively controlled near `N_f=16`, lattice-supported).

---

## §2 · What the acid test shows

- **"Single coupling ⇒ no β zero" is FALSE [NG].** The window refutes it (one coupling, a zero exists).
- **"AF ⇒ gapped" is FALSE [NG].** The window is AF and gapless-conformal.
- **For `N_f=0`:** two-loop `b0,b1>0` ⇒ no *perturbative* zero; but "no *nonperturbative* zero" is **not** a
  theorem — a strong-coupling IR-stable zero is not excluded by the 2-loop form. That is `(G1)` `[O]`.
- **Adjacent pure-adjoint theories:** `N=1` SYM (adjoint gluino) is AF and **gapped** (gaugino condensate);
  `N=2` SYM is AF with a **Coulomb branch** (gapless, IR-free U(1)); `N=4` is conformal but **not** AF. None is
  a "pure-gauge `N_f=0` gapless-conformal" — consistent with `(G1)` for `N_f=0`, but again only *dynamically*.

---

## §3 · Verdict

- **The acid test does not falsify (G1).** It fails to exhibit a `N_f=0` gapless-conformal gauge theory, but it
  **shows the reason is matter content, not "single coupling" / AF**.
- **(G1) is purely dynamical:** "`N_f=0` has no nonperturbative IR fixed point" is equivalent to "pure YM is
  outside the conformal window" — a genuine nonperturbative dynamical fact, **not** a kinematic theorem from the
  β-function's exact form. This is `H1` in the β/RG-attractor coordinate (`STAR-SYNTHESIS`).
- **The kinematic shortcut to (G1) is dead** (the window); proving (G1) requires a genuine nonperturbative
  dynamical argument (e.g. a nonperturbative sign theorem for β at `N_f=0`, or a conformal bootstrap that closes
  at the putative strong-coupling fixed point). None is known.

**Clay = 0%.** No QED.

*Refs:* `w7/g1_acid.py`; `STAR-NO-BZ.md`, `STAR-SYNTHESIS.md`.
