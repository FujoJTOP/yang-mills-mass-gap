# TOOL-SYNTHESIS — the FORCING blind spot is structural: 7/8 complex-function-theory tools are blind (OPE-RHO) or read (A3); only the L²/singular-kernel route (HARDY) escapes OPE-RHO for atom-PRESENCE detection, but not for forcing the weight

> **VERDICT (no QED, Clay 0%).** Eight complex-function-theory tools were built to fill the identified
> TOOL BLIND SPOT — a *non-integral, forcing, continuum-robust* lower bound on the atom (or its isolation),
> the one thing neither the forcing tools (all blind) nor the reading tools (all read-not-force = A3) provide.
> **Seven of eight fail, and the failure is structural, not a tool-choice artifact.** The classification is
> exact: **values of `Π` → integrals → OPE-RHO (blind, upper bounds only); zeros/poles of `Π` → A3
> (read-not-force); positivity → upper-bounds an usc object; geometry → support-only (a point has zero
> capacity).** The **one genuine escape** is `TOOL-HARDY`: an `L²`/singular-kernel functional
> `N(ε)=∫_{ε<|x−m₁²|<δ}|G(x+i0)|²dx` that is **quadratic with a singular kernel** — outside the bounded-kernel
> OPE-RHO integral class — and is **acid-clean and continuum-robust** (unlike the WI complex-zero, destroyed
> by 0.01% continuum). It detects atom **presence** and recovers the weight from the divergence rate, but the
> weight needs the pole **location** (= the gap) and `ε→0` (= A3), and the `L→∞` survival = isolation = `H1`.
> **So: "detect whether the atom exists" is fillable (HARDY); "force a positive lower bound on the atom
> weight" is NOT — it needs a nonlinear YM-structure constraint or infinite data `J→∞` = `H1`.**

---

## §0 · Labels

**[T]** 真定理 · **[C]** 条件定理 · **[O]** 精确缺失步 · **[NG]** no-go · **[V-num]** new numeric. Acid: `≤0`/`=0` on free massless *and discriminatingly* on `:φ²:`.

---

## §1 · The eight tools

| tool (file) | verdict | mechanism |
|---|---|---|
| **Dual / LP** `TOOL-DUAL` | `[NG]` | strong duality gap = 0; the **atomless measure is feasible** (mass at `M/χ`), dual optimum 0 ⇒ no dual trick. **Meta: the blind spot is convex-structural.** Escape = nonlinear constraint or `J→∞` = H1 |
| **Potential** `TOOL-POTENTIAL` | `[NG]` | capacity / equilibrium measure are **support-only** (a point mass has zero capacity); reading the log-coefficient = A3 |
| **Conformal / extremal** `TOOL-CONFORMAL` | `[NG]` | every conformal invariant is support-only; the resolvent level-set radius `≈f₁²/R` = A3 |
| **Carathéodory** `TOOL-CARATHEODORY` | `[NG]` | Taylor coefficients **are integrals** (`∫dν/(s+z₀)^{k+1}`); positivity = **upper** bound only (atom usc) |
| **Nevanlinna–Pick** `TOOL-NEVANLINNA` | `[NG]` | point evals of `Π` are **integrals** (OPE-RHO); LP range `[0, f_max]` includes 0; **values → integral class, zeros → A3** |
| **Toeplitz / Szegő** `TOOL-TOEPLITZ` | `[NG]` | Christoffel is an **upper** bound (`λ_n ↓ f₁²`); no finite-moment lower bound exists (`Δ_min=0`) |
| **Hardy / Nehari (L² singular kernel)** `TOOL-HARDY` | **PARTIAL** | `N(ε)=∫|G|²` (quadratic, singular kernel) **ESCAPES OPE-RHO**: diverges iff atom, atomless finite; **continuum-robust** (coeff 0.49984/0.49984/0.49987 for continuum mass 0/0.25/1.0); distinguishes atom from equal-mass bump. **But weight needs location (=gap) + `ε→0` (=A3); `L→∞` survival = isolation = H1** |
| **Zero-counting / argument principle** `TOOL-ZEROCOUNT` | `[NG]` | count is not a faithful level counter; **positivity admits `s₂=s₁`** (continuum with no floor atom) ⇒ cannot exclude isolation 0; reads `s₁`, not the gap = `I2=H1` |

**7 / 8 blind or A3. 1 / 8 (HARDY) escapes OPE-RHO for presence only.**

---

## §2 · The structural classification (the round's real theorem)

The campaign's recurring dichotomy — **atom is usc, finite data is lsc** — is now *proved across the whole
complex-function-theory family*:

- **Values / moments / Taylor / point-evals of `Π`** — finite measurable **integrals** with bounded kernels ⇒
  **OPE-RHO applies verbatim** ⇒ blind to a thin atom; at best **upper** bounds (`Δ_max`, Christoffel, NP range).
- **Zeros / poles / residues of `Π`** (Padé, WI, argument principle) — the atom's *only* non-integral content,
  reachable only by analytic continuation ⇒ **A3 (read-not-force)**.
- **Positivity** (Carathéodory–Fejér, Christoffel) — **upper-bounds** an usc object; it cannot force a lower bound.
- **Geometry / potential / conformal / capacity** — **support-only**; a point mass has zero capacity, so a point
  atom is invisible to every geometric invariant.
- **Convex duality** — **exact** (gap 0) with the atomless measure feasible ⇒ inherits the blindness.

⇒ **No bounded-kernel/analytic-interpolation/geometric/positivity/convex tool can FORCE a lower bound on the
atom.** This is a structural theorem, not a failure of any one tool.

---

## §3 · The one genuine escape: the L²/singular-kernel route (HARDY)

`N(ε)=∫_{ε<|x−m₁²|<δ}|G(x+i0)|²dx` is **quadratic** in `ν` with a **singular (unbounded) kernel** — it is
**not** a bounded-kernel integral, so **OPE-RHO does not apply**. Consequences (all `[T]/[V-num]`):
- **Diverges iff atom** (`N·ε→2f₁⁴`), **finite iff atomless** ⇒ acid-clean, detects presence.
- **Continuum-robust** (recovery coefficient stable under continuum mass `0→1`; contrast the WI complex-zero destroyed by `0.01%`).
- **Distinguishes atom from an equal-mass smooth bump** ⇒ genuinely outside the integral no-go.

**But it does not fill the *forcing* blind spot:** the quantitative weight needs the pole **location** `m₁²`
(locating the pole = the gap), the `ε→0` limit (= A3), and the `L→∞` survival of the isolated pole =
isolation = `H1`. So HARDY gives a real, non-integral, continuum-robust **presence detector** — a genuine new
tool — but not a free forcing lower bound.

---

## §4 · Honest bottom line

- **No gap.** 8 tools, 7 `[NG]`, 1 partial escape. No QED. `Clay = 0%`.
- **The *forcing* blind spot is structural, not fixable by complex-function-theory tools.** The exhaustive
  classification (values=integral/OPE-RHO; zeros/poles=A3; positivity=upper-bound; geometry=support-only;
  duality=exact/zero) shows the forcing lower bound cannot come from this family. **Escape = a nonlinear
  constraint that uses the YM structure (not "ν is a positive measure"), or infinite data `J→∞` = `H1`.**
- **The one real yield:** `TOOL-HARDY` — a **non-integral, continuum-robust, acid-clean atom-presence detector**
  (`N(ε)` divergence), the first tool to genuinely escape the OPE-RHO integral no-go. It detects *presence*;
  it does not force the *weight* (A3) or establish *isolation* (`H1`).
- **Useful read-tools accumulated:** exact `s₁` zero ladder, Rouché zero-free strip, Re-slope functional
  (`WI-MODIFIED-FUNCTIONAL`), effective-mass plateau (P3), the SU(3) protocol.
- **What remains (unchanged):** `H1` = the atom weight + isolation, now blocked in yet more coordinates — and
  the tool search has *proved* that no complex-function-theory inequality forces it. **The invariant did not move.**

*Refs:* `TOOL-DUAL.md`, `TOOL-POTENTIAL.md`, `TOOL-CONFORMAL.md`, `TOOL-CARATHEODORY.md`, `TOOL-NEVANLINNA.md`,
`TOOL-TOEPLITZ.md`, `TOOL-HARDY.md`, `TOOL-ZEROCOUNT.md`; `SAVVIDY-WI-SYNTHESIS.md`, `CONSTRUCT-SYNTHESIS.md`,
`lean/YMStatement/H1Sharpening.lean`.
