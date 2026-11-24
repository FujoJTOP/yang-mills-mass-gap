# KAPPA-LARGEN — large-N does NOT compute κ: it fixes only the 1/N² scaling (κ is N‑independent at LO, κ_∞ ≈ band bottom); the value is the large‑N gap = H1

> **VERDICT (no QED).** κ (the rigid‑string bending rigidity, `κ∈[0.138,0.193]`, the single missing
> input of the flux‑tube/string family = `H1` in string coordinates) **is not computed by the large‑N
> attack**. The `1/N` expansion controls the **N‑scaling** of κ (κ is **N‑independent at leading 't Hooft
> order**, with `1/N²` corrections) but **not its value**: `κ_∞ = (m_∞/√σ)²/(8π²)` is the large‑N
> **glueball mass** — the gap itself. Two further obstructions: the physical `κ(N=3)` is a *finite‑N*
> object distinct from `κ_∞` (κ *decreases* with N; even a successful large‑N value would need the `1/N²`
> correction), and the flux tube is **not thin** at large N (width `~1/√σ`, N‑independent), so no
> thin‑tube derivative expansion supplies κ. **κ = H1 survives; the large‑N lever on κ is closed.**

---

## §0 · Labels, acid test

**[T]** 真定理 · **[C]** 条件定理 · **[O]** 精确缺失步 · **[NG]** no‑go for the named class · **[V-num]** new numeric.
**Acid test:** the whole rigid‑string construction needs a **confining string** (`σ>0`); the free
massless field and the gapped‑atomless `:φ²:` have **no flux tube** (`σ=0`), so the construction is
**absent** there — vacuously clean (it does not even fire on the gapless / atomless adversaries).

---

## §1 · What κ is, and the rigid‑string identity

The neutral `0⁺⁺` as a closed flux tube with **bending rigidity** (energy cost per unit length of
curvature `E_bend=(κ/2)∮K²ds`, κ dimensionless) has the rigid‑string energy
```
E(R) = 2πσR + πκ/R ,      K=1/R on a circle  ⇒  (κ/2)·(2π/R) = πκ/R .
```
Interior minimum `R* = √(κ/2σ)`, ground state `m_{0⁺⁺}=2π√(2σκ)` ⇒
```
   κ = (m_{0⁺⁺}/√σ)² / (8π²) .                                             (★)
```
`(★)` is exact **within the model**. It is the corpus's established relation
(`CONSTRUCT-MASTERFIELD`), and it is also the reason every route to κ is circular: **(★) says κ is the
gap rescaled** — `κ∈[0.138,0.193] ⟺ m_{0⁺⁺}/√σ∈[3.30,3.90]` (witness §1).

---

## §2 · The large‑N scaling law `[C]`

Under 't Hooft scaling (`λ=g²N` fixed) the dynamically generated scale `Λ` is **N‑independent**, so at
leading order `m_{0⁺⁺}=c_m Λ` and `√σ=c_σ Λ` are both **N‑independent**; the pure‑gauge `1/N`
expansion is in **`1/N²`** (orientable surfaces; odd powers correspond to non‑orientable surfaces, absent
for unitary `SU(N)`). Hence, via `(★)`,
```
   κ(N) = κ_∞ + c/N² + O(1/N⁴) ,     κ_∞ = (m_∞/√σ)²/(8π²) .               (2.1)
```
So the `1/N` expansion **fully controls the N‑dependence** of κ (a clean `1/N²` series) — a genuine
`[C]` structural result. What it does **not** control is the leading constant `κ_∞`.

---

## §3 · Why the long‑N attack does not compute `κ_∞` `[NG]` (three independent reasons)

**(a) κ_∞ is the gap in disguise `[T]`.** By `(★)`, `κ_∞=(m_∞/√σ)²/(8π²)`; `dκ/d(m/√σ)=(m/√σ)/(4π²)≠0`.
Knowing `κ(N)` for all N is *equivalent* to knowing `m_{0⁺⁺}(N)/√σ` for all N. So "compute κ at large N"
= "compute the large‑N `0⁺⁺` glueball mass" — the **same wall `H1`**, now in the N‑scaling coordinate. The
large‑N master field gives **σ** (the 1‑form / loop‑equation side, `PNT-LARGEN-LOOP`, closed in `W`) but
**not** `m_{0⁺⁺}`: the `0⁺⁺` mass is the linearization of the loop‑equation fluctuation around the **unknown**
master field — the 1‑form/0‑form channel wall.

**(b) The tube is not thin at large N `[T]`.** The flux‑tube width is `w²(R)=(D−2)/(2πσ)ln(R/r₀)+w₀²` with
intrinsic core `w₀~1/√σ`; under 't Hooft scaling `σ~Λ²` is N‑independent ⇒ **`w₀` is N‑independent**. The
tube does **not** shrink as `N→∞` (the large‑N simplification is that the string does not split/emit
mesons, *not* that it is thin). So there is **no** `1/(σw₀²)→0` small parameter, hence **no** thin‑tube
derivative expansion that would compute κ.

**(c) Pure Nambu–Goto is unstable `[T]`.** The free (NG) string determinant gives `κ=0` (the closed string
shrinks, `E₀′(ℓ)=σ+π/(3ℓ²)>0`): the large‑N *free* string does **not even possess a stable `κ>0` closed
state**. A physical `κ>0` is a non‑NG, non‑perturbative input (`CONSTRUCT-WORLDSHEET-CFT`: the `(∂²X)²`
coupling is an **irrelevant dim‑4 LEC**, invisible to Virasoro/modular/RP).

**Holographic loophole closed.** The only context where `κ_∞` *could* be a controlled tree‑level number is
a large‑N **string dual** with a **derived** IR wall (`CONSTRUCT-HOLO-FLUXTUBE`: Witten cigar gives a genuine
`0⁺⁺` atom, `m/√σ≈3.5–4`). But for pure YM the wall is **assumed** (put in by hand) and equals the
confinement scale — so it is `H1` in the bulk‑geometry coordinate, not a large‑N computation.

---

## §4 · Quantitative: κ_∞ lands at the **bottom** of the band, and κ *decreases* with N `[V-num]`

Using `(2.1)` fit to representative large‑N data (`m/√σ`: 3.79 at `SU(2)`, 3.33 at `SU(8)`;
`[V-num]`, inherited from `CONSTRUCT-LARGEN-MASTER`, illustrative):

| N | m/√σ | κ |
|---|---|---|
| 2 | 3.790 | 0.182 |
| 3 | 3.517 | 0.157 |
| 8 | 3.330 | 0.140 |
| ∞ | 3.299 | **0.138** (band bottom) |

So **κ is N‑independent at LO and equals the band bottom `κ_∞≈0.138`; `1/N²` corrections raise it at
finite N**. Two honest caveats:
- The physical `κ(N=3)` is a **finite‑N** object **distinct** from `κ_∞`; even a successful large‑N value
  would require the `1/N²` correction. (And the direct `SU(3)` lattice value `m/√σ≈3.9 ⇒ κ≈0.193` sits at
  the **top**, **above** the naive `1/N²` extrapolation `≃0.157` — a real `~0.036` finite‑N/lattice
  systematic, flagged not resolved.)
- The number `κ_∞≈0.138` is **not an independent confirmation** of the band: the band was *defined* as
  `m/√σ∈[3.30,3.90]`, so `κ_∞` at `m/√σ=3.30` is the same datum. It is a **consistency** statement (κ has a
  clean `1/N²` scaling and the large‑N limit sits at the band edge), not a derivation.

---

## §5 · Verdict and exact missing step

- **κ is not computed at large N.** `1/N` fixes the **N‑dependence** (`κ(N)=κ_∞+c/N²`, a genuine `[C]`),
  but the value `κ_∞=(m_∞/√σ)²/(8π²)` is the large‑N **gap**. `[NG]` for "large‑N computes κ".
- **Three independent reasons κ_∞ is not a softer object:** κ_∞ is the gap in disguise; the tube is not
  thin (`w₀~1/√σ`); pure NG is unstable (`κ=0`). The holographic loophole needs an assumed IR wall (= H1).
- **Exact missing step `[O]`:** the large‑N `0⁺⁺` closed‑string mass `m_∞/√σ` — equivalently the
  linearization of the Makeenko–Migdal fluctuation around the unknown master field — `= H1` in the
  N‑scaling coordinate.
- **What is genuinely gained:** κ has a clean, universal **`1/N²` scaling** (N‑independent at LO), and the
  large‑N attack is **closed** — it is not a lever on κ. The physical value remains the finite‑N gap.

**Clay 0%; no QED.**

*Files:* `KAPPA-LARGEN.md`; witness `w7/kappa_largen.py` (numpy/math only, runs clean).
*Refs:* `CONSTRUCT-MASTERFIELD.md`, `CONSTRUCT-KAPPA-DERIVE.md`, `CONSTRUCT-WIDTH-KAPPA-LATTICE.md`,
`CONSTRUCT-WORLDSHEET-CFT.md`, `CONSTRUCT-HOLO-FLUXTUBE.md`, `PNT-LARGEN-LOOP.md`,
`CONSTRUCT-LARGEN-MASTER.md`, `CONSTRUCT-LARGEN-BOOTSTRAP.md`.
