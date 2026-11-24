# TOOL-CARATHEODORY — Carathéodory/Toeplitz/Christoffel positivity of the Stieltjes transform gives an UPPER bound on the atom only; the Taylor-at-`z₀` data are integrals (OPE-RHO class) and the LP lower end is exactly 0 — **no forcing lower bound, no escape**

> **VERDICT (no QED).** The stated hope was a *forcing, non-integral* lower bound on `f₁²` from the Taylor
> coefficients of `Π(z)=∫dν(s)/(s+z)` at a point `z₀`, using the positive-real (Carathéodory/Toeplitz)
> structure. It **fails on both counts**. **(i) The data are integrals, not non-integral.** The Taylor
> coefficients are `Π^{(k)}(z₀)/k! = (−1)^k∫dν(s)/(s+z₀)^{k+1}` — integrals with continuous bounded
> kernels, so they sit **directly in the OPE-RHO weak-continuous class**; positivity over them is blind to
> a thin atom. **(ii) Positivity is an upper-bound machine.** The Christoffel/Carathéodory bound gives the
> **maximum** atom at `s₁`; the **minimum is 0** (witness LP range `[0, 0.5954]` for a true atom `0.5`; a
> continuum-only measure matches the *exact* Taylor data to `2.8×10⁻¹⁷`). A lower bound on the atom is not
> a positivity consequence — it is the opposite direction (atom is usc, positivity is lsc).

---

## §0 · Labels

**[T]** 真定理 · **[C]** 条件定理 · **[O]** 精确缺失步 · **[NG]** no-go for the named class · **[V-num]** new numeric. Acid: `≤0`/`=0` on free massless *and discriminatingly* on `:φ²:`.

---

## §1 · The construction

**Stieltjes transform and its Taylor data.** For the 0⁺⁺ KL measure `ν` on `[m₁²,∞)`,
`Π(z)=∫ dν(s)/(s+z)`, analytic and positive for `z>0`. Its normalised Taylor coefficients at `z₀>0` are

```
a_k := Π^{(k)}(z₀)/k! = (−1)^k ∫ dν(s) / (s+z₀)^{k+1},     k = 0,1,2,…
```

**Key structural fact `[T]`.** Each `a_k` is an **integral** of `dν` against the kernel `φ_k(s)=(−1)^k/(s+z₀)^{k+1}`,
which is **continuous and bounded on `[m₁²,∞)`**. Hence the Taylor-at-`z₀` data are **weak-continuous
integrals `∫ρ φ_k`** — **exactly the OPE-RHO class**. They are *not* "non-integral": derivatives of the
Stieltjes transform are still integrals. (The genuinely non-integral datum — the pole *location* and
*residue* — requires analytic continuation; that is A3, read-not-force, §3.)

**Carathéodory / positive-real / Christoffel.** `Π` is positive-real (`Re Π>0` off the negative axis, `Π(z)>0`
for `z>0`). Carathéodory–Fejér / Toeplitz positivity constrains the coefficients: for a positive-real
function the Toeplitz matrices are PSD and `|c_k|≤2 c_0`. The classical **Christoffel function** gives the
**maximum** atom at a point consistent with the data:
`ν({s₁}) ≤ λ_n(s₁) = 1/K_n(s₁)`, where `K_n` is the Christoffel–Darboux kernel of the truncated problem.
**This is an upper bound.** The minimum consistent atom is `0` (a density, or a narrow bump, reproduces the
same integrals).

---

## §2 · The test (witness `w7/tool_caratheodory.py`, `py -3.12 -I`, runs clean)

True measure: atom `0.5` at `s₁=1`, continuum weight `0.04` at `s∈{1.05,…,8}`; `z₀=5`, `K=8`.

| quantity | value |
|---|---|
| Taylor coeffs `a_0..a_7` | `1.464e-1, −2.258e-2, 3.558e-3, −5.686e-4, …` |
| atom contribution to `a_0,a_1,a_2,a_7` | `8.33e-2, 1.39e-2, 2.31e-3, 2.98e-7` (decays fast; dominated by low `k`) |
| **(1)** continuum-only `nnls` fit (no atom at `s₁`) | matches the **exact** Taylor data to `max|Mw−a| = 2.8×10⁻¹⁷` |
| **(2)** LP range of atom at `s₁` consistent with the data | `[0.0, 0.5954]` (true `0.5`) |
| **(3)** `Re Π(z)>0` for atom vs no-atom measures | `True, True` — Carathéodory bound satisfied by **both** |

**Reading.**
- **(1)** A measure with **no atom at `s₁`** reproduces the *same* Taylor coefficients to **machine
  precision** — a narrow continuum bump near `s₁` is indistinguishable from the `δ` under smooth-kernel
  integrals. This is the OPE-RHO thin-atom blindness, now for the `Π`-Taylor data.
- **(2)** The LP lower end is **exactly `0`**; the upper end `0.5954` is the (Christoffel-type) maximum.
  **Positivity gives an upper bound only; there is no positive lower bound.**
- **(3)** Both the atom and the no-atom measures are positive-real, so the Carathéodory coefficient bound is
  a constraint they *both* satisfy — it cannot discriminate.

---

## §3 · Does it force? Is it non-blind? Does it escape OPE-RHO?

- **Force a positive lower bound `f₁²>0`? NO `[NG]`.** Positivity (Carathéodory/Toeplitz/Christoffel) is an
  **upper-bound** machine: it bounds how *large* the atom can be, never how large it *must* be. The lower
  end is `0` (LP + nnls, exactly). This is the structural statement "atom is usc, positivity is lsc"
  (`PROFILE`): a positivity/lower-bound tool cannot print a lower bound on an usc object.
- **Non-blind (atom+continuum)? NO `[NG]`.** The atom's Taylor contribution decays as `(1+z₀)^{−k−1}` and
  is mimicked by a continuum bump near `s₁` (residual `2.8×10⁻¹⁷`). It is blind to the atom even at finite
  weight, because the data are integrals against smooth kernels.
- **Escape OPE-RHO? NO `[T]`.** The Taylor-at-`z₀` coefficients are **weak-continuous integrals** with
  continuous bounded kernels — *verbatim* the OPE-RHO class. (Contrast `WI`, which was a *relational*
  invariance and escaped the integral class; the Carathéodory coefficients are *absolute integrals*, so
  they do not.) The only non-integral content of `Π` — the **pole location/residue** — is recovered by
  analytic continuation, i.e. **A3, read-not-force = `H1`**.

**Acid test.** Free massless field: `ν` has no atom, `f₁²=0`; a gapped-atomless `:φ²:` also `f₁²=0`. The
Carathéodory bound is satisfied by both (it does not fire) — so it is **acid-vacuous**: it does not
discriminate the atom from `:φ²:`, exactly the failure criterion.

---

## §4 · Bottom line

- **The tool fails.** Carathéodory/Toeplitz/Christoffel positivity of the Stieltjes transform is a **positivity
  (upper-bound) tool on integral data**; it gives `ν({s₁}) ≤ λ_n(s₁)` and **nothing on the lower side**
  (LP min `= 0`, continuum-only matching `2.8×10⁻¹⁷`). *`[NG]` for "forcing lower bound via Carathéodory."*
- **Not an escape.** The Taylor data are integrals (OPE-RHO class); the non-integral datum is the pole
  residue = A3 = `H1`. So the Carathéodory route re-lands on `A3/OPE-RHO/H1` — the **fifth or sixth**
  independent confirmation that "positivity cannot force the atom."
- **Honest yield:** a clean statement that the Carathéodory coefficient bound is **equivalent to an upper
  bound on the atom** and is **acid-vacuous** for the KL atom — a negative result that closes the
  "positivity-from-Taylor" idea, not a new lower-bound tool.

**Clay 0%. No QED.**

---

*Files:* `research/unconditional/TOOL-CARATHEODORY.md` (this), `research/unconditional/w7/tool_caratheodory.py`
(runs clean, `py -3.12 -I`, exit 0). **Not committed** (parent stages).
*Refs:* `PNT-OPE-RHO.md` (OPE-RHO theorem), `PROFILE.md` (atom usc / positivity lsc), `WI-ESCAPE.md`
(relational vs integral), A3 (`PROOF-F1SQ` / Stieltjes residue).
