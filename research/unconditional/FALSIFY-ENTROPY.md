# FALSIFY-ENTROPY — nonlinear information functionals do NOT force the atom from finite data: they are convex (blind, TOOL-DUAL), density-requiring (diverge ⇒ presence-only ⇒ J→∞=H1), or Markov (blind, Theorem BE). The only-exit claim SURVIVES.

> **VERDICT (no QED).** Attempt to falsify "the only exit is YM-structure or `J→∞`=H1" by a **nonlinear
> information functional** (differential entropy `−∫ρlogρ`, relative entropy, Fisher information, Rényi,
> log-Sobolev) that lower-bounds the atom weight `f₁²` from finite data. **The falsification FAILS.** Every
> candidate falls into one of three dead classes: **(1) convex in ν** (relative entropy, Rényi `α>1`) ⇒
> blind (TOOL-DUAL: the atomless measure is feasible, exact duality); **(2) density-requiring** (differential
> entropy, Fisher, Shannon/Rényi limit) ⇒ an atom makes them **diverge** (−∞ entropy, ∞ Fisher), so they are
> **presence-only detectors needing the full density** = `ε→0` = `J→∞` = `H1`, not a finite-data forcing bound;
> **(3) Markov-semigroup** (log-Sobolev) ⇒ Theorem BE (homogeneous/zero/sign-indefinite) ⇒ blind. **No finite-data,
> non-convex, non-YM, non-`J→∞` functional in this family forces `f₁²`.** The only-exit claim is corroborated.

---

## §0 · Labels

**[T]** 真定理 · **[C]** 条件定理 · **[O]** 精确缺失步 · **[NG]** no-go · **[V-num]** new numeric. Acid: `≤0`/`=0` on free massless *and discriminatingly* on `:φ²:`.

---

## §1 · The three dead classes

### (1) Convex in ν ⇒ blind (TOOL-DUAL)
Relative entropy `D(ν‖ref)` and Rényi divergence `α>1` are **convex** in `ν`. The feasible set of measures
matching any fixed finite data (mass, `χ`, `m₁²`, finitely many moments) is **convex and contains the
atomless measure** (mass at `M/χ`). A convex functional therefore cannot certify `ν({m₁²})>0` uniformly over
the feasible set — the atom-weight infimum is `0` and (TOOL-DUAL) the LP duality gap is exactly `0`. So convex
information functionals inherit the integral-class blindness. **[NG]**

### (2) Density-requiring ⇒ diverge ⇒ presence-only ⇒ `J→∞`=`H1`
Differential entropy `H(ρ)=−∫ρlogρ`, Fisher `J(ρ)=∫ρ'²/ρ`, and the Shannon/Rényi limit require a **density**.
Native to an atom: `H=−∞` (a delta has no density; a width-`ε` bump has `H≈f₁²log(ε)` → `−∞`) and `J=∞`
(a width-`ε` bump has `J≈f₁²/ε²` → `∞`). So these functionals are **presence detectors by divergence** —
exactly the `TOOL-HARDY`/`L²` class: they see the atom only through the `ε→0` limit, i.e. they need the full
measure (`J→∞`), and even then give *presence*, not a *weight lower bound*. **[NG]** as an escape — this *is* `J→∞`=`H1`.

### (3) Markov-semigroup ⇒ Theorem BE ⇒ blind
The log-Sobolev / entropy-production functional of the gradient-flow semigroup is a Markov-semigroup functional;
`BASE-FIELDMEASURE`'s Theorem BE proves any such functional is homogeneous in the Hessian, zero at equilibrium, or
sign-indefinite — it cannot be the κ-independent source. **[NG]**

---

## §2 · Witness (`w7/falsify_entropy.py`, numpy, `py -3.12 -I`)

Atom `f²=0.5` at `s=1` + uniform continuum on `[1,4]`, smoothed to a Gaussian bump of width `ε`; `ε→0`
recovers the atom.

| ε | M0 | M1 | M2 | M3 | H(ρ) | Fisher(analytic `f²/ε²`) |
|---|---|---|---|---|---|---|
| 0.20 | 0.9970 | 1.7488 | 4.0199 | 11.1864 | **0.9561** | **12.5** |
| 0.10 | 1.0001 | 1.7501 | 4.0055 | 11.1417 | **0.6841** | **50.0** |
| 0.05 | 1.0001 | 1.7501 | 4.0017 | 11.1305 | **0.3805** | **200.0** |
| 0.02 | 1.0001 | 1.7501 | 4.0007 | 11.1273 | **−0.0426** | **1250.0** |

**Finite moments `M0..M3` change by only `5.3×10⁻³`** from `ε=0.20` to `0.02` ⇒ any finite-moment functional is
essentially **blind** to the atom. Meanwhile `H` decreases toward `−∞` and Fisher grows as `1/ε²` ⇒ these are
**divergence/presence** detectors that need `ε→0` (the full density = `J→∞`). So the witness confirms: finite
data cannot distinguish (moments flat), and the information functionals only "see" the atom by diverging at
infinite resolution — not a finite-data forcing bound.

---

## §3 · Honest bottom line

- **Falsification fails.** No nonlinear information functional gives a **finite-data, non-convex, non-YM,
  non-`J→∞`** lower bound on `f₁²`. Convex (KL, Rényi `α>1`) → blind; density-requiring (entropy, Fisher,
  Shannon) → divergence/presence at `ε→0` = `J→∞`; log-Sobolev → Theorem BE.
- **The only-exit claim is corroborated on this route:** to force (not merely read/detect) the atom, one needs
  a **nonlinear constraint that uses the YM structure** (not "ν is a positive measure") or the **infinite-data
  limit `J→∞`** = `H1`.
- **What it does yield:** entropy/Fisher are presence detectors in the same class as `TOOL-HARDY` (`L²`
  divergence) — a real but read/presence tool, not a forcing bound.
- **Clay 0%; no QED.** The invariant did not move.

*Refs:* `TOOL-SYNTHESIS.md`, `TOOL-DUAL.md`, `TOOL-HARDY.md`, `BASE-FIELDMEASURE`/Theorem BE, `PROFILE.md`;
witness `w7/falsify_entropy.py`. Labels [T]/[NG]/[O]/[V-num].
