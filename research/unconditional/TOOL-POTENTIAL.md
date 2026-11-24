# TOOL-POTENTIAL — potential theory as a lower-bound tool for the atom `f₁²`: reading is A3, capacity is blind, precision is tied to the unknown atom ⇒ **[NG]** (no forcing bound); real yield = a clean *why*

> **VERDICT (no QED).** The logarithmic potential `U^ν(z)=∫log|z−s|dν(s)`, the equilibrium measure and the
> logarithmic capacity do **not** force a lower bound on a thin atom. Three exact facts: **(1)** reading the
> atom = the coefficient of the `log|z−s₁|` singularity = `f₁²` = the Stieltjes residue `−Res F` = **A3
> (read-not-force)**; **(2)** complex-point evaluations are **integrals** (continuous kernel `log|z−·|`) = the
> **OPE-RHO class**, and to resolve a thin atom `f₁²` the resolution must satisfy `ε ≪ f₁²` — data precision
> tied to the **unknown** atom ⇒ no uniform bound; **(3)** the **capacity / equilibrium measure are
> support-determined** (the log capacity of `[s₁,S]` is `(S−s₁)/4`, the equilibrium measure is arcsine),
> **independent of `f₁²`** — a thin atom at the floor does not change the support ⇒ **capacity is blind to the
> atom mass**. So the species is a *restatement* of A3 in potential-theory language, not a forcing tool.

---

## §0 · Labels
**[T]** 真定理 · **[C]** 条件定理 · **[O]** 精确缺失步 · **[NG]** no-go · **[V-num]** new numeric.

---

## §1 · The dictionary (potential theory ≡ Stieltjes transform)
`∇U^ν(z) = ½ F(z)` with `F(z)=∫dν(s)/(z−s)` the Stieltjes transform. So the log potential is the
Stieltjes transform integrated; an atom at `s₁` of mass `f₁²` is a `f₁² log|z−s₁|` singularity whose
coefficient is `f₁² = −Res_{s₁} F`. **This is exactly the A3 object.** The equilibrium measure `μ_K` of the
support `K` and the log capacity `cap(K)` are functionals of `K` **alone**.

## §2 · Three routes, all non-forcing
- **(A) Reading** `[T]`: `U(s₁+iε)=f₁² log ε + C_cont(s₁) + o(1)`, so `f₁² = lim U(ε₁)−U(ε₂) / log(ε₁/ε₂)`.
  This reads `f₁²` but is the **residue = A3 (read-not-force)** and needs **analytic continuation** to
  complex `z` (not physically measurable).
- **(B) Precision tied to the unknown** `[V-num]`: the finite-`ε` error is `O(ε)`; resolving `f₁²` needs
  `ε ≪ f₁²`, i.e. the required resolution scales with the **unknown atom** ⇒ **no uniform lower bound**.
- **(C) Capacity / equilibrium are blind** `[T]`: `cap([s₁,S])=(S−s₁)/4`, equilibrium = arcsine,
  **independent of `f₁²`**; a thin atom at the floor leaves the support unchanged ⇒ no information.
- **(D) Energy / equilibrium-defect**: `I(ν)=+∞` for **any** atom (`f₁²>0`) and `0` for none — a
  *qualitative* detector but gives **no quantitative bound** (it is `+∞` regardless of `f₁²`), and knowing
  `I(ν)=∞` already assumes the atom ⇒ circular.

## §3 · Witness (`w7/tool_potential.py`)
Measure `ν = f₂δ_{s₁} + (1−f₂)·Unif[s₁,S]`, `s₁=1, S=4`.
- **Reading** works to `1e-6` at `ε≈1e-6` (e.g. `f₂=1e-4 → 0.000101`) — but complex-`z` = A3.
- **Coarse `ε≈1e-2`**: `err/f₂ = 0.20, 2.04, 20.44` for `f₂=1e-2,1e-3,1e-4` — a thin atom is **swamped**;
  the precision must scale with `f₂` (unknown).
- **Capacity** `cap([1,4])=0.75`, equilibrium arcsine — **independent of `f₂`** (blind).
- **Acid**: `f₂=0 →` extraction `1.3e-6≈0` (does not misfire on the atomless continuum), but a thin atom is
  also read `≈0` ⇒ cannot force `f₁²>0`.

## §4 · Verdict
- **[NG]** for "potential theory forces a lower bound on a thin atom."
- **Real yield (the clean *why*)**: potential theory is the **Stieltjes/A3 object in disguise** — its one
  atom-sensitive datum (the log-coefficient) is the residue, its global invariants (capacity/equilibrium)
  are **support-only and blind to the mass**, and its complex-point evaluations are integrals in the
  OPE-RHO class. It **sharpens** the established "atom is usc / finite data lsc" picture: *even the
  analytic-continuation potential cannot force a thin atom; it can only read it (A3).*
- **No QED; Clay 0%.** A clean negative that closes the "potential-theory / equilibrium / Fekete" tool
  class: not forcing, and its reading route is the same A3 wall.

*Files:* `research/unconditional/w7/tool_potential.py` (runs clean, `py -3.12 -I`). No commit (parent stages).
