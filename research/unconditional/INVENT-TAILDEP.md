# INVENT-TAILDEP — the extreme-value / tail-dependence functional is a *new type* but **not a new handle**: it is the correlation length (`κ`) relabelled, atom-blind, and fires on the free-massless field; as a transition detector it is killed by B2

> **VERDICT (no QED).** Bare-handed species: an extreme-value-theory (EVT) functional of the gauge
> measure — the **tail index `ξ`** and the **extremal index `θ_EV`** (clustering of `tr F²` exceedances
> across coarse blocks), equivalently the **tail-dependence / extremal coefficient**. It is formally a
> new *type* (not an order parameter, not a Markov-semigroup functional, not a topological index). But it
> **is not a new handle**: `ξ≈0` (light tail) for every adversary (no discrimination); `θ_EV` is a
> **monotone function of the correlation strength = `κ`** (homogeneous in `κ`, not a source); it cannot
> tell gapped-YM from gapped-atomless `:φ²:` (**atom-blind**); it fires on the **free-massless** field
> (long-range ⇒ low `θ`) on the **wrong side** of the acid test. And as an EVT/large-deviation
> **transition detector** it is killed by **B2** (zero-T 4D YM is single-phase). **Exact failing step: the
> functional reduces to `H1` in the `κ` direction, and is atom-blind (it is not a KL-atom reader). Clay
> 0%.**

---

## §0 · Labels, iron laws, acid test

**[T]** 真定理 · **[C]** 条件定理 · **[O]** 精确缺失步 · **[NG]** no-go for the named class · **[V-num]** new numeric.
**Acid test** (carried): the functional must be `>0` on the interacting YM measure, `≤0` on the
**free-massless** field (`f₁²=0`, `G″=0`), and **discriminatingly `≤0` on the gapped-atomless `:φ²:`
composite** (gapped, `f₁²=0`). — The last clause forces **atom-sensitivity** (a mere gap detector fires
on `:φ²:` and is rejected).
**B2** (`ROOF-TRANSFER`, `[T]`): zero-T 4D pure YM is a **single phase** (crossover, no transition).

---

## §1 · The species (bare-handed definition)

For the 0⁺⁺ local operator `O(x)=tr F²(x)` (or `1−N⁻¹Re tr U_□` on the lattice), viewed as a stationary
random field under `μ`:

- **Tail index `ξ`**: `P(O>t) ~ t^{−1/ξ}` (Fréchet `ξ>0`; Gumbel `ξ=0`, light; Weibull `ξ<0`, bounded).
- **Extremal index `θ_EV ∈ [0,1]`**: `P(max_{block} O > t) ≈ θ_EV · P(O>t)` for large `t` — clustering of
  exceedances over a coarse block. `θ_EV=1` ⇔ exceedances essentially independent (fast correlation decay);
  `θ_EV<1` ⇔ clustering (slow/long-range correlation).
- **Tail-dependence / extremal coefficient `λ_T`**: `lim_{t} P(O_i>t | O_j>t)` for two separated blocks.

The hope (bare-handed): `ξ` or `θ_EV` is a **phase-sensitive, non-order-parameter, non-Markov,
inhomogeneous-in-`κ`** number. It is a genuinely different *type*. We test whether it is a different
*handle*.

---

## §2 · Acid test — FAILS (witness `w7/invent_taildep_witness.py`, clean)

AR(1) Gaussian proxy: `ρ→0` = short-range/gapped, `ρ→1` = long-range/massless; `:φ²:` proxy = the squared
series (Wick/χ²).

| `ρ` (correlation) | `ξ` | `θ_EV` | role |
|---|---|---|---|
| 0.00 | −0.12 | **0.974** | gapped proxy |
| 0.30 | −0.13 | 0.921 | |
| 0.60 | −0.13 | 0.768 | |
| 0.90 | −0.10 | 0.364 | |
| 0.99 | −0.15 | **0.102** | massless (long-range) |

`:φ²:` composite (squared series): `ξ≈0.01–0.04`, `θ_EV = 0.973, 0.935, 0.810, 0.402` at `ρ=0,0.3,0.6,0.9`
— **the same `θ_EV` as the Gaussian adversary** (within run-to-run noise).

Three exact consequences:

1. **`ξ≈0` everywhere (light/Gumbel); on the lattice `O` is bounded (Weibull `ξ<0`).** The tail index does
   **not** discriminate — and in the continuum the tail is **UV/local** (short-distance fluctuations, which
   are universal by asymptotic freedom), hence **floor-blind** (the same F4 obstruction that kills the
   moments).
2. **`θ_EV` decreases monotonically with correlation strength** (0.974→0.102 as `ρ`:0→0.99): it is a
   **monotone function of the correlation length = `κ`** — `κ`-relabelled, **homogeneous in `κ`** (fails
   S3). It is not a `κ`-independent source.
3. **Atom-blindness:** gapped-YM-proxy and gapped-atomless `:φ²:` give the **same `θ_EV` (~1)**. EVT cannot
   distinguish a gapped theory *with* a KL atom from a gapped theory *without* one — it is not an atom
   reader.

**Sign impossibility:** with only `θ_EV(gapped YM) ≈ θ_EV(:φ²:) ≈ 1 > θ_EV(massless)`,
- `F=θ_EV−1`: `≤0` on massless ✔, `=0` on `:φ²:` ✔, but `=0` on YM (**not `>0`**) ✘;
- `F=1−θ_EV`: `=0` on YM (**not `>0`**) ✘, `=0` on `:φ²:` ✔, but `>0` on massless ✘.

**No sign convention passes the acid test**, because `θ_EV` cannot be positive *only* on the atom.

---

## §3 · Inhomogeneity / `κ`-homogeneity — FAILS S3

`ROOF-UNIFORM-CHI` demands a `κ`-**independent** source (`κ=F(κ)+𝒮`). But `θ_EV` is determined by the
correlation function (Gaussian: by the covariance), which **is** the mass `κ`; `θ_EV(κ)→1` as `κ→0⁺`
(gapped), `<1` as `κ→0` (massless). Hence `θ_EV` is the **κ-dependent RHS**, exactly like the
Poincaré/LSI/Wasserstein functionals killed by **Theorem BE** — not an additive source. It supplies **no**
new bound on `χ_lattice ≲ a²`, `σ_w(2.5)>c₀`, or `R_vac`.

---

## §4 · EVT/LDP is a transition detector — killed by B2

The natural EVT companion, the **large-deviation rate function** `I(x)` for the block-averaged operator, is
a **phase-transition detector**: a kink/non-analyticity in `I` marks a transition. But **B2** (`[T]`):
zero-T 4D pure YM is a **single phase** (crossover, no transition) ⇒ `I` is analytic ⇒ **nothing to
detect**. Likewise `θ_EV`'s clustering is the EVT shadow of long-range dependence, i.e. the *gapless* side —
the opposite of what the acid test wants. The EVT/LDP class is therefore a **third** casualty (after the
order-parameter and Markov-functional paradigms) of the same single-phase fact.

---

## §5 · Verdict and exact failing step

- **New type, not a new handle.** EVT/tail-dependence is formally none of {order parameter, Markov
  functional, topological index} — a genuine bare-handed invention. But it is **equivalent to the
  correlation length `κ`** (homogeneous, `κ`-relabelled), **atom-blind**, and **fires on the wrong side**
  (free-massless).
- **Exact failing step `[O]`:** `θ_EV` reduces to `H1` along the `κ` direction, and the tail index is a
  UV/local (floor-blind) datum; the class cannot read the KL atom. It also dies as a transition detector
  via **B2**.
- **Genuine (negative) yield:** it **closes the EVT/large-deviation class** — a new named category ruled out
  for the 0⁺⁺ gap, by the same single-phase theorem (B2) that killed order parameters. It does **not**
  provide the strut.
- **Clay = 0%; no QED.**

*Witness:* `w7/invent_taildep_witness.py` (numpy+scipy, `py -3.12 -I`, clean).
*Refs:* `ROOF-SYNTHESIS.md`, `ROOF-TRANSFER.md` (B2), `BASE-FIELDMEASURE.md` (Theorem BE),
`ROOF-UNIFORM-CHI.md` (S3 / `κ`-homogeneity), `PROFILE.md`.
