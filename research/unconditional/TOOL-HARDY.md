# TOOL-HARDY — the Hardy/H² norm of the resolvent boundary value is a genuine NON-INTEGRAL atom detector (escapes the OPE-RHO integral no-go, robust to any continuum), but its quantitative weight is A3 (needs the pole location + the L→∞ limit) = H1

> **VERDICT (no QED, Clay 0%).** Filling the "non-integral atom detector" part of the tool blind spot: the
> L² (Hardy H²) norm of the Stieltjes resolvent boundary value `G(x+i0)` **diverges iff an atom is present
> and is finite iff atomless** (`[T]`); its divergence rate is `2f₁⁴`, so it **recovers f₁²** (witness:
> 0.49992 vs true 0.5) and is **robust to any continuum weight** (coeff 0.49984 for continuum mass 0, 0.25,
> 1.0). Because it is **quadratic with a singular kernel** it **distinguishes an atom from an equal-mass
> smooth bump** (bump: finite `N→14.04`; atom: `N·ε→0.5`) — i.e. it **escapes the OPE-RHO integral
> no-go**, which is blind to exactly this. **BUT** using it in YM needs (a) the **pole location `m₁²`**
> (locating the pole = the gap) and (b) the **ε→0 limit** (A3, read-not-force); and the atom-vs-continuum
> distinction is really the **L→∞** limit (finite-L is all poles ⇒ divergence trivial). So it is a genuine
> new non-integral tool, but the L→∞ content is **isolation = H1**.

---

## §0 · Labels

**[T]** 真定理 · **[C]** 条件定理 · **[O]** 精确缺失步 · **[NG]** no-go · **[V-num]** new numeric.
Acid: `≤0`/`=0` on free massless *and discriminatingly* on `:φ²:`.

---

## §1 · The object (and why it is non-integral)

KL measure `ν` on `[0,∞)`, atom `f₁²` at `m₁²`. Resolvent / Stieltjes transform
`G(z)=∫ dν(s)/(s−z)`, boundary value `G(x+i0)=PV∫dν(s)/(s−x) + iπρ(x)`.
The atom is a **pole** of `G` at `m₁²`: near it `G(x) ≈ f₁²/(m₁²−x) + G_cont(x)`.

**Functional:** the truncated Hardy-type L² norm
`N(ε) = ∫_{ε<|x−m₁²|<δ} |G(x+i0)|² dx`.
This is **quadratic** in the measure and has a **singular kernel** (`|1/(·−x)|²`), so it is **not** of the
OPE-RHO form `∫ρφ` with bounded `φ`. It is the natural Hardy H² (`∫|G|²`) truncated around the pole.

---

## §2 · The theorem `[T]`

With `G = f₁²/(m₁²−x) + G_cont` (`G_cont` smooth near `m₁²`), on a **symmetric** annulus:

`N(ε) = 2 f₁⁴ (1/ε − 1/δ) + O(1)`  ⇒  `ε·N(ε)/2 = f₁⁴ (1 − ε/δ) + O(ε) → f₁⁴`.

So:
1. **Atom present (`f₁²>0`): `N(ε) ~ 2f₁⁴/ε → ∞`.** The divergence coefficient `lim_{ε→0} εN(ε)/2 = f₁⁴`, hence `f₁²`.
2. **Atomless (`f₁²=0`): `N(ε) → ∫|G_cont|² < ∞`** — no divergence. **Acid-clean** (free massless ⇒ finite).
3. **Smooth bump of mass `f₁²`, width `γ`:** `G=f₁²/(m₁²−x−iγ)+…`, `|G|²=f₁⁴/((m₁²−x)²+γ²)`,
   `N(ε)→πf₁⁴/γ<∞` — **finite, no divergence** (for any `γ>0`).

⇒ The norm **detects the atom** (divergence) and **distinguishes it from an equal-mass smooth bump.**

---

## §3 · Witness `[V-num]` (`w7/tool_hardy_witness.py`, exact quadrature)

Model: atom `f₁²=0.5` at `m₁=1` + flat continuum `ρ₀` on `[4,8]`, `δ=0.3`.

| ε | ATOM `N·ε` | BUMP(γ=0.05) `N·ε` | ATOMLESS `N·ε` |
|---|---|---|---|
| 1e-2 | 0.48332 | 0.12087 | 0.00026 |
| 1e-3 | 0.49833 | 0.01386 | 0.00003 |
| 1e-4 | **0.49984** | **0.00140** | 0.00000 |

- Atom `N·ε → 0.5 = 2f₁⁴`; **recovered `f₁² = 0.49992`** (true 0.5). Bump `N→14.04` finite; atomless `N→0.027` finite.
- **Continuum robustness:** `N·ε = 0.49984 / 0.49984 / 0.49987` for continuum mass `0 / 0.25 / 1.0` — the atom
  divergence coefficient is **independent of the continuum**. (Contrast the WI complex-zero `𝒜`, destroyed by
  0.01% continuum.)

---

## §4 · It escapes the OPE-RHO integral no-go `[T]`

OPE-RHO: every finite measurable **integral** `∫ρφ` (bounded `φ`) is weak-continuous (lsc), so a thin atom and
an equal-mass smooth bump give the **same** value — the atom is unreadable. `N(ε)` is **quadratic with a
singular kernel**: atom ⇒ `∞`, equal-mass bump ⇒ finite. So `N` **separates them** — it is **outside** the
OPE-RHO class. This is the genuine escape: **the no-go is specific to bounded-kernel integrals; the L²/singular
route is not blind.**

---

## §5 · The price (why it is still = H1) `[O]`

1. **Pole location.** `N(ε)` must be centered on `m₁²`. Locating the pole `m₁²` **is** locating the gap = `H1`.
2. **The limit.** Exact `f₁²` needs `ε→0` (the `O(ε)` error is a strict *underestimate*, `f₁⁴ ≥ εN(ε)/2`, tight
   in the limit) — an A3-style read-not-force.
3. **Finite-L triviality.** On a finite lattice the spectrum is **discrete** ⇒ `G` is a finite sum of poles ⇒
   `N` **always diverges** (every finite-L state is an "atom"). The atom-vs-continuum distinction is really the
   **L→∞ limit** (does an isolated pole survive? = **isolation = H1**).

---

## §6 · Nehari / Hankel view `[T]/[NG]`

Cayley transform `θ=(G−i)/(G+i): ℂ⁺→𝔻`. The atom ⇒ `θ→1` at `m₁²` (a singular-inner/Blaschke boundary point).
**Nehari:** `dist_{H∞}(θ)=‖H_θ‖` (Hankel norm) — measures the atom's **presence** (the singular-factor degree =
a **non-integral count**), but `‖H_θ‖≤1` in magnitude does **not** give the residue. The **weight** needs the
`L²` divergence rate `2f₁⁴` = the **Stieltjes residue** = **A3 (read-not-force)**. So Nehari certifies *that* an
atom is there (non-integral), not *how big* — and locating it is the gap.

---

## §7 · Honest verdict

- **Genuine new tool `[T]`:** a **non-integral atom detector** exists — the Hardy/H² norm of the resolvent
  boundary value — that **escapes the OPE-RHO integral no-go** (atom vs equal-mass bump separable), is
  **acid-clean** (finite iff atomless), **robust to any continuum**, and recovers `f₁²` from its divergence rate.
  This fills the **non-integral atom detector** part of the blind spot.
- **Still `H1` `[O]`:** the detector needs the **pole location** (locating the pole = the gap), the **ε→0 limit**
  (A3), and the real content is the **L→∞ survival of an isolated pole = isolation = H1**.
- **Sharpening:** the OPE-RHO no-go is **specific to bounded-kernel integrals**; the singular/`L²` route is not
  blind. That is a real structural clarification — but it does **not** yield a forcing lower bound free of the
  isolation problem. **Clay = 0%; no QED.**

*Refs:* `PNT-OPE-RHO.md`, `WI-ESCAPE.md`, `WI-MODIFIED-FUNCTIONAL.md`, `A3`/`PROOF-F1SQ`; witness
`w7/tool_hardy_witness.py`.
