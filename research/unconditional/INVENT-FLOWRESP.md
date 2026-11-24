# INVENT-FLOWRESP — the Lüscher gradient-flow *response* is a renormalon-free **gap detector**, but it is **atom-blind** (fires on the `:φ²:` continuum) and its response is a **moment** (homogeneous) — a new handle, the old wall

> **VERDICT (no QED).** The Yang–Mills gradient flow gives the single best *renormalon-free* handle
> the campaign has found: the flowed energy `E(t)` is finite by construction (Lüscher), so any
> flow-time response is UV-finite with no `u=2` ambiguity. **This genuinely addresses the strut's
> renormalon obstruction.** But the flow-time response is **not a gap certificate**: (i) `E(t)` is the
> **Laplace transform** of the spectral measure, `E(t)=∫ s\,ν(ds)\,e^{−2st}`, so it detects the
> **continuum edge** (a scale) but is **blind to the atom** `f₁²`; the atomless `:φ²:` continuum gives a
> **finite freeze time**, so the response fires on it — **acid FAIL** (it is a gap detector, not an atom
> certifier); (ii) the "response" `∂E/∂κ` is a **moment** (linear functional of `ν`) — homogeneous in the
> source, i.e. the RHS of the fixpoint, **not** the `κ`-independent source the base round requires.
> **New handle (renormalon-free), but it re-encodes the gap scale, it does not bound it.** Clay 0%.

---

## §0 · Labels, iron laws, acid test

**[T]** 真定理 · **[C]** 条件定理 · **[O]** 精确缺失步 · **[NG]** no-go for the named class · **[V-num]** numeric
witness. **Acid**: the tool must be `≤0`/`=0` on the free massless field, **and discriminatingly** on the
gapped-atomless `:φ²:` composite (`f₁²=0`). **Iron laws**: (IL-i) non-positivity; (IL-ii) nonabelian;
(IL-iii) no small parameter.

---

## §1 · The species, defined

Lüscher flow: `∂_t B_μ = −D_ν F_{νμ}(B)`, `B_0 = A`, `t` of mass dimension `−2`. The flowed energy
`E(t) = ⟨tr G_{μν}(t)G_{μν}(t)⟩` is **finite to all orders** and nonperturbatively well-defined (no UV
divergence — the flow provides a gauge-invariant, renormalized probe). The flow damps an intermediate state
of invariant mass-squared `s=p²` by `e^{−2p²t}=e^{−2st}`, so structurally

```
E(t) = ∫ s · ν(ds) · e^{−2 s t}      (Laplace transform of the 0⁺⁺ spectral measure ν).
```

**The invented response.** Insert an auxiliary coarse mass `κ` (a shift of the spectral threshold / the
operator's effective mass) and define the **flow-time response** `∂E(t)/∂κ`, and the **freeze time**
`t_* := 1 / (−d log E/dt)` at large `t` (the inverse of the late-time decay rate).

---

## §2 · Acid test — the witness (`w7/invent_flowresp_witness.py`)

Three spectral measures, exact quadrature:

| measure | large-`t` rate `−dlogE/dt` | freeze `t_*` | acid |
|---|---|---|---|
| **free massless** `ν~s\,ds` | `3/t → 0` (power law) | `∞` (no freeze) | **correct** |
| **atomless `:φ²:`** continuum from `θ=(2m)²=4` | `→ 2θ = 8` | `0.123` | **FIRES ⇒ FAIL** |
| **gapped atom** `f₁²δ(s−m₁²)` + continuum | `→ 2m₁² = 2` | `0.500` | — |

**The acid failure is decisive:** the atomless `:φ²:` theory (`f₁²=0`) is gapped (continuum edge at
`2m`), and the flow response **fires on it** (finite `t_*=0.123`, positive `1/√t_*`). So the response
**cannot distinguish the glueball atom (`f₁²>0`) from the `:φ²:` two-particle continuum (`f₁²=0`)** — it is
a **gap detector, not an atom certifier**. This is exactly the base round's `gap-detector ≠ atom-certifier`
fact (`:φ²:`), re-derived in the flow coordinate.

**Why, structurally:** `E(t)` is the **Laplace transform** of `ν`. The Laplace transform is injective in
principle (all `t` determines `ν`), but (a) the inverse is **discontinuous** — a finite `t`-window /
finite precision cannot stably resolve an atom (the F4 "computable data" wall), and (b) the response's
large-`t` rate `→ 2×(continuum edge)` in **both** the atomless and atom cases; the atom's weight `f₁²`
never enters the rate. **The flow sees the scale, never the atom.**

---

## §3 · The response is a moment (homogeneous), not the source

`∂E/∂κ` (with `κ` shifting the threshold `θ`) is, from §1,

```
∂E/∂θ = −∫_θ^∞ s e^{−2 s t} ds      (witness: ratio ∂E/∂κ : E(t) = O(t), both ∝ moments).
```

This is a **linear functional of `ν`** — a connected correlator / susceptibility. It is **homogeneous in
the deformation**: it scales with `E` and vanishes with the measure. It is therefore the **RHS** of a
fixpoint `κ = F(κ)`, **not** the `κ`-independent inhomogeneous source `𝒮` that the base round proved is
required (the source must be nonzero at `κ=0`, a constant added to a homogeneous relation). The flow
response has the wrong homogeneity.

---

## §4 · Does it bound the base strut (`χ_lattice ≲ a²` / `σ_w(2.5)>c₀` / `R_vac`)? No.

The flow response **measures** the scale: `t_* ~ 1/m₁²` presupposes the spectrum `ν`. To turn `t_*` into a
**lower** bound on `m₁` you must already know the spectral shape (which is the theory). It gives no
inequality `χ_lattice ≲ a²` etc. **The renormalon-free advantage is real, but it buys a cleaner
*measurement*, not a *certificate*.**

---

## §5 · Verdict — new handle, old wall

- **(ii) renormalon-free: YES — genuine.** `E(t)` finite by construction ⇒ no `u=2` ambiguity. This is the
  best renormalon cure found in the campaign (better than the raw `⟨tr F²⟩` of `ROOF-CONDENSATE-FLOOR`).
- **(i) acid-clean: NO.** Atom-blind; fires on the atomless `:φ²:` continuum; a gap detector, not an atom
  certifier.
- **(iii) bounds the strut: NO.** A measurement of the scale, not a lower bound.
- **New or re-encoding?** A **genuinely new handle** (renormalon-free flow-time observable) that
  **re-encodes the gap scale** (`t_* ~ 1/m₁²`) rather than bounding it — the old `H1` wall in the flow
  coordinate.
- **Exact failing step `[O]`:** `E(t)` is the Laplace transform of `ν` ⇒ atom-blind (the rate is set by the
  continuum edge, not `f₁²`); and `∂E/∂κ` is a moment (homogeneous) ⇒ not the `κ`-independent source.

**What survives as a tool:** a clean **renormalon-free *measurement*** of the `0⁺⁺` scale via `t_*` — useful
for numerics, **not** a gap certificate. **Clay 0%; no QED.**

*Refs:* `ROOF-CONDENSATE-FLOOR.md` (renormalon obstruction this address), `BASE-FIELDMEASURE.md`
(Theorem BE, homogeneity), `BASE-NONLOCAL-BOUND.md` (1/√χ), `PROFILE.md` (P1 data-mode: spectral
data are `{blind, circular}`); witness `w7/invent_flowresp_witness.py`.
