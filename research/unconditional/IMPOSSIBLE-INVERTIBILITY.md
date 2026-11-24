# IMPOSSIBLE-INVERTIBILITY — the gap as operator INVERTIBILITY / nonzero DEGREE, not spectral positivity of a measure: a clean reformulation, but the degree is only `0`/`+1` (Krein–Rutman) and no index is conserved across β (rigidity) ⇒ **no new lever, reduces to the gap**

> **VERDICT (no QED, Clay 0%).** Jumping out of the failed bases (measure positivity, weak-* topology, order
> parameters, `A3`). **Impossible base tried:** the mass gap = **invertibility of `(I − T_phys)`** on the excitation
> space, proved by a **degree / topological-index / Leray–Schauder** argument on a **nonlinear operator**, not by
> spectral positivity of `ν`. Result: this is a genuine **reformulation** (frames the gap as an operator property,
> not a measure property) and it is **acid-clean in form** (`deg = 0` ⟺ gapless; `deg = +1` ⟺ gapped), **but it
> yields no new lever**: for a positive compact transfer operator the degree is only `0` or `+1` (Krein–Rutman),
> i.e. it **is** the gap's binary indicator; the **homotopy-invariance** hope (carry the gap from `β<1/48` to
> `β=2.5`) fails because **no integer index is conserved** across the crossover (the gap is a non-analytic
> `β`-jump — rigidity `B3`); and the nonlinear `Leray–Schauder` route needs the self-energy map to be **compact**,
> which it is not (infinite tower), and any truncation loses the gap (admits the gapless fixed point). So the
> "impossible base" **lands on the gap itself** — consistent with `ONLY-EXIT`.

---

## §0 · Labels

**[T]** 真定理 · **[C]** 条件定理 · **[O]** 精确缺失步 · **[NG]** no-go · **[V-num]** new numeric. Acid: distinguishes free-massless from gapped.

---

## §1 · The base: gap = invertibility on the excitation space

On the physical Hilbert space `H_phys` (gauge-invariant, BRST cohomology), the transfer operator
`T = e^{−aH}` is a **positive** operator with Perron eigenvalue `λ₀ = 1` (the vacuum `Ω`). The gap is
`λ₁ < 1`, equivalently:

> **`0 ∉ spec(H|_{Ω^⊥})` ⟺ `(I − T)` is invertible on `Ω^⊥` ⟺ `(I − T/λ₀)` has a nonzero degree.**

This is a reformulation: "the gap" = "an operator is invertible on the excitation space". It does **not**
mention the KL measure `ν`, positivity cones, or weak-* topology — so it is genuinely **off every failed base**.

---

## §2 · The linear / finite-dimensional construction: the degree is only `0` or `+1`

Normalize `λ₀ = 1` (divide by the Perron eigenvalue). On the excitation space,

`I − T/λ₀` is diagonal in the energy basis with entries `1 − e^{−a(E_i − E₀)}`, so

**`det(I − T/λ₀) = ∏_{i≥1} (1 − e^{−a(E_i − E₀)}) ≥ 0`, and it is `0` exactly when `E₁ = E₀` (gapless).**

Hence the degree `deg(I − T/λ₀) = sign det ∈ {0, +1}`: **`0` ⟺ gapless, `+1` ⟺ gapped.** This is `[T]` and
it **does** distinguish the free massless field (eigenvalue `E₁ = E₀ = 0` ⇒ `det = 0`) from a gapped theory
(`det > 0`) — **acid-clean in form**.

**But this is not a richer invariant — it is the gap's binary indicator.** By Krein–Rutman, a positive compact
operator has a simple Perron eigenvalue and all others `|λ_i| < λ₀`; so `det(I − T/λ₀)` is a **single sign**
whose value is exactly "is `λ₁ < λ₀`" — i.e. exactly the gap. There is no winding number, no rich index that
could carry extra information (unlike a Hopf-degree problem); the "degree argument" computes `∏(1 − e^{−aΔ_i})`,
which presupposes the `Δ_i` = the spectrum = the gap. **Witness (A):** the sign-det is `+1` iff gapped, `0` at
`Δ→0`; it is a `0/1` detector, not an independent handle.

---

## §3 · The homotopy hope — and why it fails (no conserved index)

The only way a degree argument could **force** the gap (rather than restate it) is **homotopy invariance**:
deform `β` continuously from the **rigorously gapped** strong-coupling end (`β < 1/48`, OS'78/SZZ) to
`β = 2.5`, and if an integer degree were **conserved** and nonzero at `β<1/48`, the gap would persist.

**This fails, structurally.** The transfer-operator degree is `0`/`+1` and it **drops to `0` at any level
crossing** `E₁(β) = E₀(β)` — the witness (B) drives `E₁(s) = s − ½`, so the indicator is `+1` for `s≠½` and
`0` at `s = ½`: **there is no nonzero integer invariant that is conserved through the crossing.** This is the
degree-theoretic face of the **rigidity no-go (`B3`)**: the gap is a **non-analytic `β`-jump**, not a
protected topological constant; and of `BASE-INDEX`: the whole index class (instantons, spectral flow, eta,
Callias) is **floor-blind** — homotopy invariants are constant on connected `β`-components and cannot detect a
jump. **No symmetry protects the gap** (pure YM has no protecting index above the vacuum), so there is nothing
for a degree argument to conserve. `[NG]`.

---

## §4 · The nonlinear / Leray–Schauder angle — compactness fails

The directive names a **nonlinear** operator (Brouwer / Leray–Schauder / Conley). The natural nonlinear map
whose zero/fixed set encodes the gap is the **self-energy / DSE / RG gap equation** `κ = F(κ)` (`κ` the coarse
mass). `κ = 0` is always a fixed point (the gapless solution, `SLOT-DSE`/`SLOT-RICCATI`); the gap = a
**nonzero, isolated** fixed point. A degree argument `deg(I − F, 𝒜, 0) ≠ 0` on an annulus `𝒜` avoiding `0`
would **force** such a nonzero fixed point.

**It fails on compactness — the exact step:**

1. **Leray–Schauder needs `F` compact.** The physical `F` is the **infinite DSE/ERG tower** (1-contour sources
   2-contour, …) — not compact on any controlled space; the LS degree is undefined. `[NG]`
2. **Truncation loses the gap.** Any finite truncation `F_N` is compact and has a degree, but the truncated
   equation admits the **gapless** fixed point `κ=0` and its nonzero fixed point depends on the truncation —
   degree `0` or a truncation artifact, not the gap (the campaign's `SLOT-DSE`/`SLOT-BOOTSTRAP`: truncation
   ⇒ `Δ_min=0`). `[NG]`
3. **Toy is circular.** `F(κ) = λκ/(1+κ²)`: `F(0)=0`; `deg(I−F)` forces a nonzero fixed point **iff `λ > 1`**
   — and `λ>1` **is** the `k→0` tachyon = the gap. So the degree "forces" the gap only after assuming the
   tachyon, i.e. after assuming the gap. `[O]`

---

## §5 · Acid test

**Passes in form, fails in content.** A degree/invertibility argument distinguishes the free massless field
(`E₁=E₀ ⇒ det(I−T)=0`, singular) from a gapped theory (`det>0`, invertible) — so it is **acid-clean in
form**. But "`det ≠ 0`" **is** `∏(1 − e^{−aΔ_i}) ≠ 0`, which presupposes the `Δ_i` = the gap. The invariant
carries no information beyond the sign of the gap. It is a clean **reformulation** (operator-invertibility), not
a **handle**.

---

## §6 · Verdict

- **Clean new framing `[T]`:** gap = `(I − T_phys)` invertible on `Ω^⊥` = nonzero degree — genuinely off the
  measure/positivity bases.
- **No new lever `[NG]`:**
  - the degree is only `0`/`+1` (Krein–Rutman) — it **is** the gap's indicator;
  - **no index is conserved across `β`** (rigidity `B3`; `BASE-INDEX` floor-blind) — the homotopy hope fails;
  - the **nonlinear LS route needs compactness**, which the infinite DSE/ERG tower does not have; truncation
    loses the gap; the toy is circular (`λ>1` = the tachyon = the gap).
- **Exact missing step `[O]`:** a **compact, gap-encoding nonlinear operator with a conserved, computable,
  nonzero degree invariant across the crossover** — none exists (no protecting symmetry; the tower is not
  compact). The base lands on the gap.

**Net:** the "impossible base" is a real, acid-clean **reformulation** of the gap as an operator-invertibility
condition, but it is **not a path around `H1`** — the degree is trivial (`0/+1`) and unprotected (rigidity).
This is the same wall, in the language of operator degree instead of measure positivity. **Clay 0%; no QED.**

*Refs:* `BASE-INDEX.md` (rigidity, floor-blind index), `SLOT-DSE-FIXPOINT.md` / `SLOT-RICCATI-FLOW.md`
(gapless fixed point always exists), `ONLY-EXIT-SYNTHESIS.md`, `lean/YMStatement/Reduction.lean`.
Witness: `w7/impossible_invertibility.py`.
