# G1-NONPERT-BETA — no nonperturbative sign theorem for `β_{N_f=0}`: the sign of β at strong coupling is not protected by any exact identity; G1 is the gap in the β-function coordinate

> **VERDICT (no QED).** Target (G1): prove `β_{N_f=0}(g) < 0` for **all** `g>0` (no IR fixed point) nonperturbatively.
> **No route yields a nonperturbative SIGN theorem.** Perturbatively `β<0` for small `g` is `[T]`; but for large
> `g` the sign is **dynamical**, and — decisively — **no exact identity forbids `β(g*)=0`**: at a putative IR
> fixed point the trace anomaly vanishes (`θ=(β/2g)trF²→0`), the a-theorem allows `a_IR>0`, and unitarity allows
> a conformal IR. So any sign theorem must be a **dynamical computation**, which is the gap itself. G1 is `H1` in
> the β-function/RG-attractor coordinate.

---

## §0 · Labels

**[T]** 真定理 · **[C]** 条件定理 · **[O]** 精确缺失步 · **[NG]** no-go. Acid: discriminates free massless / `:φ²:`.

---

## §1 · Perturbative sign `[T]`

`β(g) = −β₀ g³ − β₁ g⁵ − …` with `β₀ = 11N/(48π²) > 0`, `β₁ > 0` (for `N_f=0`). So `β(g) < 0` for small `g`
asymptotically. Witness: SU(2),SU(3),SU(4) — no positive zero on `g∈(0,1]`. This is only trusted at weak coupling.

---

## §2 · The four routes, and why each fails to give a nonperturbative sign theorem

**(i) a-theorem monotonicity `[NG]`.** `a_UV ≥ a_IR`, stationary only at fixed points. At a putative IR fixed
point `g*`, `a_IR = a_CFT(g*) > 0` is **allowed** (`a_UV > a_IR > 0`); a gapped IR has `a_IR = 0`, also allowed.
The a-theorem constrains the *flow*, **not the sign of β** — it excludes neither branch. Witness: both `a_IR=0`
(gapped) and `0<a_IR<a_UV` (CFT) satisfy the a-theorem.

**(ii) exact trace-anomaly β + OPE `[NG]`.** `θ^μ_μ=(β/2g)trF²`. The trace-anomaly coefficient `β/2g` is a
*dynamical* function of `g`; the OPE fixes its *operator structure*, not its sign. At `β(g*)=0` the anomaly
operator vanishes — exactly the CFT condition — so the anomaly identity is **satisfied** at a fixed point, not
violated. No sign.

**(iii) background-field / Schwinger–Dyson `[NG]`.** The SD/2PI equations give the *gap equation* (a self-consistency
for the mass/kernel), not a sign of `β(g)`. The exact β is defined by the renormalization of `g`, and no
nonperturbative SD closure pins its sign. (A reduced `F(κ)=λκ/(1+κ²)` has a positive fixed point only for
`λ>1` — the `k→0` tachyon — which is the gap, circular.)

**(iv) lattice nonperturbative running `[O]`.** The Schrödinger-functional step-scaling *measures* `g²(L)`; it
shows growth (AF → confinement) over the probed range but is **numerical, finite-range**, not a theorem. A
*rigorous* all-scale lower bound on the running is the crossover/transport problem = `H1`.

---

## §3 · The general obstruction `[T]`: no exact identity forbids `β(g*)=0`

At an IR fixed point `β(g*)=0`, **all** exact structures are satisfied:
- trace anomaly `θ=(β/2g)trF² = 0` (what a CFT wants);
- a-theorem `a_UV ≥ a_IR > 0` (allowed);
- unitarity / reflection positivity (a conformal IR is unitary);
- the conformal window is a **real** demonstration: `SU(3)` with `N_f=12,16` has a genuine Banks–Zaks zero — an
  asymptotically-free theory that is **gapless and conformal in the IR**.

⇒ **No sign theorem can be derived from exact structure.** The exclusion of a pure-YM IR fixed point is
**dynamical** (`N_f=0` outside the conformal window) = confinement = the gap.

---

## §4 · The one exact sign theorem is SUSY — and decoupling kills it `[T]/[NG]`

The NSVZ exact β for `N=1` SYM gives `β ∝ −(3N − N_f(1−γ)) g³/(1−g²N/8π²)` — for `N_f=0`, `β<0` **exactly**.
But this is **SUSY** (holomorphy/NSVZ); pure YM has no such protection, and the exact SUSY gap does not survive
decoupling to pure YM (`IMPOSSIBLE-SUSY-DECOUPLE`, `STAR-NO-BZ`). So the exact sign lives in a theory pure YM is not.

---

## §5 · Honest verdict

- **No nonperturbative sign theorem for `β_{N_f=0}` exists from exact structure.** The four routes are `[NG]`
  (a-theorem, anomaly, SD) or `[O]` (lattice all-scale bound).
- **The sign of `β` at strong coupling is dynamical**, and the witness (the artificial `β(g)=−β₀g³(1−g²/g*²)`,
  which has an IR-stable zero while its 2-loop part `−β₀g³` has none) shows 2-loop `β<0` **cannot** exclude it.
- **G1 ⟺ "pure YM is outside the conformal window" ⟺ confinement ⟺ `H1`**, now in the β-function/RG-attractor
  coordinate. It is the sharpest *re-statement* of the wall (no measure, no atom, no weak-*), and a concrete,
  falsifiable question — but **not a separate provable statement**: its proof would be the gap.
- **One genuinely new opening (untried technique pool):** the a-theorem + conformal bootstrap at the *putative*
  pure-YM fixed point — but the fixed point (if it existed) is strongly coupled, where the bootstrap is `[O]`.
- **Clay = 0%, no QED.**

*Refs:* `STAR-NO-BZ.md`, `STAR-SYNTHESIS.md`, `STAR-NO-DILATON.md`, `IMPOSSIBLE-SUSY-DECOUPLE.md`;
Komargodski–Schwimmer a-theorem (arXiv:1107.3987); Polchinski NPB 303 (1988) 226; NSVZ β (Novikov–Shifman–Vainshtein–Zakharov 1983).
