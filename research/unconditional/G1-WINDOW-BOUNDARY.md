# G1-WINDOW-BOUNDARY — N_f=0 is FAR outside the SU(3) conformal window, but "outside" is rigorous only *perturbatively*; the rigorous nonperturbative statement IS (G1) itself

> **VERDICT (no QED, Clay 0%).** The SU(3) Banks–Zaks conformal window, at 2-loop, is `N_f ∈ [9,16]`
> (`N_f^c = 306/38 ≈ 8.05`); at 4-loop `N_f^c ≈ 11`; 5-loop higher. **`N_f=0` sits far below — a distance of
> ~9–12 units of `N_f`, whereas the 2→4-loop shift of the boundary is only +3.** So `N_f=0` is not a borderline
> case: every perturbative order and the lattice step-scaling put it far outside. **But "N_f=0 is outside the
> window" is a theorem only in the *perturbative* sense** (for `N_f ≤ 8` the 2-loop β has no zero); the exact
> identities (trace anomaly, a-theorem, unitarity) do **not** forbid a nonperturbative IR fixed point, so the
> **rigorous nonperturbative** statement "`β_{N_f=0}` never reaches zero" **is (G1) itself, and remains `[O]`.**

---

## §0 · Labels

**[T]** 真定理 · **[C]** 条件定理 · **[O]** 精确缺失步 · **[NG]** no-go · **[V-num]** new numeric.

---

## §1 · The window, perturbatively `[T]`/`[V-num]`

SU(3) with `N_f` massless fundamental Dirac fermions. With `β(g) = −β₀g³ − β₁g⁵` (per `16π²` factors),

```
β₀(N_f) = 11 − (2/3)N_f            (>0  ⟺  N_f < 16.5, asymptotic freedom)
β₁(N_f) = 102 − (38/3)N_f          (<0  ⟺  N_f > 306/38 ≈ 8.05, Banks–Zaks zero g*² = −β₀/β₁ > 0)
```

- **2-loop window:** `N_f ∈ {9,…,16}` (AF upper edge 16, BZ lower edge `N_f^c = 8.05`).
- **4-loop** (van Ritbergen–Vermaseren–Larin, hep-ph/9701390; lower-edge analyses e.g. Ryttov–Sannino): the
  boundary is pushed **up** to `N_f^c ≈ 11` (some analyses 10–11). **5-loop** pushes higher (~12).
- **At `N_f=0`:** `β₀=11>0`, `β₁=102>0` ⇒ **no zero at any perturbative order**; the β-function is *most negative*
  (maximally asymptotically free). This is the farthest possible point from the window.

**Distance:** `N_f^c − 0 ≈ 9–12`, while the loop-by-loop shift of `N_f^c` is small (2→4 loop: `8→11`, i.e. **+3**).
Reaching `N_f=0` from the boundary would need a nonperturbative shift ~3–4× the known 4-loop correction.

---

## §2 · What is rigorous and what is not

- **`[T]` (algebraic):** for `N_f ≤ 8`, `β₀>0` **and** `β₁>0` ⇒ the 2-loop β has **no** zero. `N_f=0` is in this
  set. This is a genuine theorem — but it is *perturbative*.
- **`[C]` (semi-rigorous, distance):** the series corrections are monotone-increasing in the boundary
  (2-loop 8 → 4-loop 11 → 5-loop ~12) and far too small to close a gap of ~9–12 down to `N_f=0`. So `N_f=0` is
  robustly far outside *if the perturbative expansion is a guide at all*.
- **`[O]` (the honest wall):** the *nonperturbative* statement "`N_f=0` is outside the window" is **verbatim**
  "pure YM has no IR fixed point" — i.e. **(G1)**. The exact identities do **not** forbid an IR zero: at
  `β(g*)=0` the trace anomaly `(β/2g)⟨tr G²⟩=0`, the a-theorem allows `a_IR<a_UV`, and unitarity allows a
  conformal IR. So the exclusion **must** be dynamical, and no nonperturbative sign theorem for β exists.
- **Lattice `[C]` (not rigorous):** IR-conformal step-scaling (LSD, latKMI) places `N_f=12` conformal, `N_f=8`
  confining/near-conformal, `N_f=0` confining — unanimously `N_f=0` outside. Unanimous evidence ≠ proof.

---

## §3 · Honest bottom line

- **`N_f=0` is FAR outside the window** — a gap of ~9–12 units of `N_f`, and the higher-order boundary shifts
  (+3 per 2 loops) cannot bridge it. So the question "*is `N_f=0` conformal?*" is **not marginal**: every method
  (2-loop theorem, 4/5-loop, lattice) says no.
- **But this does not prove (G1).** "`N_f=0` outside the window" is rigorous only *perturbatively*; the
  nonperturbative version **is** (G1). The load-bearing `[O]` is exactly: *a nonperturbative proof that
  `β_{N_f=0}(g)` has no strong-coupling zero* (no IR fixed point). No such theorem exists.
- **Genuine value of the boundary location:** it shows `N_f=0` is a *far* boundary case, not a knife-edge — so the
  conformal-window framing does not reduce the difficulty; it relocates the wall to "prove the β-function of
  pure YM has no strong-coupling zero," which is (G1)/(H1) in the β coordinate. **Same wall, sharpest coordinate.**

*Labels:* `[T]` 2-loop no-zero for `N_f≤8`; `[V-num]` `N_f^c ≈ 8→11→12`; `[C]` distance robustness + lattice;
`[O]` nonperturbative no-zero = (G1); `[NG]` no rigorous nonperturbative proof. No QED.
