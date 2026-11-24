# ONLY-EXIT-SYNTHESIS — 8 falsification routes all fail: the only forcing exit is `H1`; "YM structure" is not a loophole but `H1` in disguise, and the deep reason is the weak-* discontinuity of the boundary density

> **VERDICT (no QED, Clay 0%).** The claim to test: *forcing a lower bound on the KL atom (or its isolation)
> requires either (a) a YM-specific nonlinear structure or (b) infinite data `J→∞`.* Eight independent routes
> tried to **falsify** it by building a forcing tool that is **neither**. **All eight fail.** The claim is
> **confirmed** — and **sharpened**: "YM structure" is **not a separate exit** (proving the field→atom
> implication *is* `H1`), so there is **one exit, `H1`**. The deep unifying reason (FALSIFY-TOPOLOGICAL):
> **the atom is exactly the weak-* discontinuity of the boundary density; every finite datum is weak-*
> continuous and atomless measures are weak-* dense — so any functional that discriminates/forces the atom
> must be weak-* discontinuous, hence not finitely computable (= A3 read, or `J→∞` = `H1`).** Scope: confirmed
> within the 8-route falsification attempt; the weak-* argument is general but not claimed as a universal theorem.

---

## §0 · Labels

**[T]** 真定理 · **[C]** 条件定理 · **[O]** 精确缺失步 · **[NG]** no-go · **[V-num]** new numeric. Acid: `≤0`/`=0` on free massless *and discriminatingly* on `:φ²:`.

---

## §1 · The eight falsification routes (all fail)

| route (file) | attempted counterexample | why it fails |
|---|---|---|
| **HARDY upgrade** `FALSIFY-HARDY-UPGRADE` | push the `L²`/singular kernel from presence to **forcing the weight** without location/limit | `[T]` over any fixed window the kernel is **bounded** (OPE-RHO blind); the escape lives only at the unknown pole, and "singular kernel ⟺ location + `ε→0`" are the same thing — cannot decouple |
| **Higher singular kernel** `FALSIFY-HIGHER-SINGULAR` | `∫\|G^{(k)}\|²`, `∫\|G\|^p` — stronger singularity forces the weight | `[T]` order `p` only changes the divergence power `ε^{−(p−1)}`, not the need for the pole location + `ε→0` (= A3) |
| **Entropy / information** `FALSIFY-ENTROPY` | differential/relative entropy, Fisher, Rényi | convex ⇒ blind (TOOL-DUAL); density-requiring ⇒ diverges only at `ε→0` = infinite data, presence only; Markov ⇒ BE |
| **OPE coefficients** `FALSIFY-OPE` | force from perturbative-UV OPE coefficients | `[T]` `Π` is **linear in ν**, each OPE coefficient *is* the moment `M₋ₖ` ⇒ no nonlinear constraint; finite OPE = OPE-RHO blind, full = `J→∞` = H1 |
| **Nonlinear Ward/anomaly** `FALSIFY-WARD-NONLINEAR` | nonlinear combinations of exact Ward `κχ=1−η`, anomaly `A`, positivity, AF tail | `[T]` an atomless flat bump matches **every** moment with `f₁²=0` and satisfies all exact identities ⇒ any continuous function of them is blind; Ward is a tautology |
| **Topological / winding** `FALSIFY-TOPOLOGICAL` | winding / argument-principle of the Herglotz function (stronger than zero-counting) | `[T]` off-axis the winding is **identical** (blind); it discriminates only across ℝ, where the atom is a **weak-* discontinuity** ⇒ **non-finitely-computable** |
| **Spectral flow** `FALSIFY-FLOW` | spectral flow of a deformed transfer operator | `[T]` flow is a **spectrum** property, the atom is an **overlap** — independent (nonzero flow coexists with `f₁²=0`); flow is rigid on connected components |
| **Field-level** `FALSIFY-FIELD` | functionals of the gauge FIELD beyond ν (Wilson loops, `⟨tr Fⁿ⟩`, instantons, gradient flow) | `[T]` (A) correlator-functionals = ν = OPE-RHO blind; (B) beyond-ν data are 1-form/topological, **Schur-orthogonal** to `tr F²` = no kinematic info; the bridge is a dynamical theorem = `H1` |

**0 / 8 falsifications succeed. The claim stands.**

---

## §2 · The sharpening: one exit, not two

`FALSIFY-FIELD` closes the last hope — that "use the specific YM structure" is a separate loophole. It is not:
- **A field functional that is a function of the 0⁺⁺ correlator** is a function of ν ⇒ integral ⇒ blind.
- **A field functional using data beyond ν** (Wilson loops, topology) is 1-form/topological, Schur-orthogonal to
  the 0⁺⁺ scalar ⇒ carries **no kinematic** information on `f₁²`; linking it to `f₁²` requires a **dynamical
  theorem** about the specific YM measure — which **is** `H1`.

⇒ **"YM structure" ≡ `H1`.** There is **one** exit (`H1`), reachable by two *coordinates* — a YM-specific
nonlinear constraint (= `H1`) or infinite data `J→∞` (= `H1`) — not two exits.

---

## §3 · The deep reason (the round's theorem): weak-* discontinuity

`FALSIFY-TOPOLOGICAL` gives the unifying characterization of *why* every route fails identically:

> **The KL atom `f₁² = δ` is exactly the weak-* discontinuity of the boundary density of `Π`.** Every finite
> datum (moments, Taylor coefficients, point evals, bounded-kernel integrals, Wilson loops via their ν-image)
> is **weak-* continuous**; the atomless measures are **weak-* dense** (a bump of width `h` matches any finite
> data to `O(h)`, raising `f₁²` from `0` to its value). Hence **any functional that discriminates or forces the
> atom must be weak-* discontinuous** — and such a functional is **not finitely computable**: a finite sample is
> always matched by an atomless measure. Evaluating it requires the boundary value across ℝ (= A3, read-not-force)
> or the full limit `J→∞` (= `H1`).

This is why **all eight routes fail identically**, why the singular-kernel route (`TOOL-HARDY`) detects *presence*
but not the *weight*, and why the "atom blind / finite-data continuous" dichotomy (`PNT-OPE-RHO`, `PROFILE`) is
not an accident of any one tool — it is the **weak-* topology of measures**.

---

## §4 · The finite-L corollary (important for what the gap actually is)

`FALSIFY-FIELD` corollary: at **finite `L`**, the atom is **trivially** `f₁² > 0` (the transfer matrix is
compact, the spectrum discrete). So a finite-lattice "measurement" of `f₁²>0` is **vacuous** as regards Clay.
The Clay content is the **`a→0`, `L→∞` survival + isolation** of the atom = `H1`. This reconciles the earlier
"measured `G″=+0.28` / `f₁²>0`" observations: they are finite-volume facts, not the gap.

---

## §5 · Honest bottom line

- **The "only exit = YM structure or `H1`" claim is CONFIRMED, not falsified, by all 8 routes** — and **sharpened
  to a single exit `H1`** (YM structure that forces ≡ `H1`; infinite data ≡ `H1`).
- **The reason is a theorem, not a tool failure**: the atom is the **weak-* discontinuity** of the boundary
  density; finite data is weak-* continuous; atomless measures are weak-* dense ⇒ any forcing functional is
  **non-finitely-computable** ⇒ only A3 (read) or `J→∞` (= `H1`).
- **Scope (honest):** this is confirmed **within the 8 falsification routes**, and the weak-* argument is general;
  it is **not** claimed as a rigorous universal no-go theorem, but it is the strongest evidence to date and
  explains every prior result.
- **What remains:** `H1` — the `a→0`, `L→∞` survival + isolation of the `0⁺⁺` atom. **The only exit is the wall.**
- **Clay = 0%.** No QED.

*Refs:* `FALSIFY-HARDY-UPGRADE.md`, `FALSIFY-HIGHER-SINGULAR.md`, `FALSIFY-ENTROPY.md`, `FALSIFY-OPE.md`,
`FALSIFY-WARD-NONLINEAR.md`, `FALSIFY-TOPOLOGICAL.md`, `FALSIFY-FLOW.md`, `FALSIFY-FIELD.md`;
`TOOL-SYNTHESIS.md`, `PNT-OPE-RHO.md`, `SAVVIDY-WI-SYNTHESIS.md`, `PROFILE.md`, `lean/YMStatement/H1Sharpening.lean`.
