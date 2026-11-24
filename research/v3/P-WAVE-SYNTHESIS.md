# P-WAVE-SYNTHESIS — the **periodic** constant mode is the physical case, and it **saturates**: exact global-holonomy anchor `0.205115` [T], periodic `k=0` floor `c=0.222` [V], transfer gap `Δ=0.6444`; the center sector does not change it; the 4D floor remains `H1`, and a curvature floor is never an atom

> **VERDICT (no QED, Clay 0%).** Six routes attacked the periodic (untwisted, no fixed boundary)
> constant mode — the physical decompactification case, which is neither the twisted strip nor the
> fixed-`b=0` strip. **In every tractable model it saturates at a positive `L`-independent floor.**
> The exact global constant-holonomy slice is flat (`S_W≡0`) in every center sector, so its marginal
> is pure Haar with `G″=1/(2·Cin(2π))=0.205115` exactly [T]. The physical `k=0` Fourier mode has an
> exact flat direction in `S_W`; its curvature is `0.2077→0.2165` for `L=2..8` at `β=2.5`, with an
> exponential fit `c=0.22196±0.00011` (`q=0.849`) and a gapped periodic transfer matrix
> (`Δ=0.6444`, `ξ=1.55`). The center sector does not change the floor. **The 4D statement is still
> `H1`**, and the floor is a susceptibility ceiling, never an atom.

---

## §0 · The six verdicts

| Route | Verdict | Result |
|---|---|---|
| **P1** exact global holonomy | `[T]`/`[V]` | On the periodic torus the constant-holonomy slice has every plaquette `=I`, so `S_W≡0`; the marginal is pure Haar and `G″_anch(0)=−ρ_Haar″/ρ_Haar=1/(2·Cin(2π))=0.20511530` (SU(2)), `L`-independent. The acid test (free `U(1)`: `−ρ″/ρ=0` vs the spurious `1/Var=3/π²`) proves the exact object is `−ρ″/ρ`, **not** `1/Var`; `1/Var=0.26883` is the Gaussian/CLT value. |
| **P2** periodic `k=0` mode | `[T]`/`[V]` | The uniform Cartan direction is an exact flat direction of `S_W`; exact `Tr[M_k^L]` gives `G″(0)=0.21146,0.21304,0.21439,0.21553,0.21649` (`L=4..8`), fit `c+Aq^L` with `c=0.22196±0.00011`, `q=0.849`, `rms=3.3×10⁻⁶`; `Var(a)` plateaus `4.39→4.45`, `G″·Var→0.963`; `β=0` reproduces D1's exact convolution. |
| **P3** Polyakov/center | `[T]`/`[V]` | The constant temporal mode is the Polyakov holonomy. The flat slice is flat in **every** center sector (`I`, `−I`, `(−I,−I)`); the periodic measure is exactly `θ`-independent (`∫ρ da=26.0980322` for all `θ`, spread `2.7×10⁻⁸`); center symmetry is unbroken at finite `L` (`⟨tr P⟩=0`, gapped kernel `Δ=0.644/3.176/1.096`). The floor is the Haar-carried P2 value `≈0.222`, not a center-sector floor. |
| **P4** periodic transfer matrix | `[T]`/`[V]` | The periodic bulk kernel has `Δ=−log(λ₁/λ₀)=0.6444`, `ξ=1.5519`, `L`-independent (identical to the D3 `θ=0` bulk kernel). Curvature `0.2077→0.2165` (`L=2..8`), fits `c=0.2126±0.0018` (`L=2..8`) / `0.2142±0.0016` (`L=4..8`). The periodic Wilson part is a bounded negative `−0.0319…−0.0451`, fit `c=−0.0461±0.0003` (`p=2.05`). |
| **P5** adversarial universality | `[T]`/`[V]`/erratum | Route B (saturation) survives; Route A (decay) survives only for the free `U(1)` field (`Var≍L`, `G″≍1/L²`). P5's apparent `SU(2)` `Var(a)` growth is a **witness bug**: at `β=0` the exact identity is `Var(a)=3.7198` for every `L`, which that moment formula violates from `L=4`. Its independent `G″` quadrature agrees with P2/P4. |
| **P6** floor → H1 logic | `[T]`/`[C]`/`[NG]` | A positive floor is a curvature floor (and conditionally a `k=0` susceptibility ceiling `χ_{k=0}≤1/c₀`), **never an atom**: at `β=0`, `G″=0.2688` while `f₁²=0`. It gives no bound on the non-constant tower, and even a proven `H1` does not force the bottom of the Källén–Lehmann support (thin-atom; `H1_does_not_imply_atom`). |

---

## §1 · The boundary-class table (the question was ill-posed until now)

| class | uniform-mode action | `Var(a)` as `L↑` | `G″` | verdict |
|---|---|---|---|---|
| fixed `θ=π` | end plaquettes pin, `O(1/L)` | plateaus | `→0.41` | saturates, **actionful 't Hooft flux** (`S_W=2β`) |
| fixed `θ=0` / `b=0` | end plaquettes pin, sign depends on `a₀` | grows / ambiguous | `1/L` fit gives `c<0` | not protected |
| **periodic** | **`S_W≡0` exactly, every `β,L`, every center sector** | **plateaus (`4.45` at `β=2.5`; exact `3.7198` at `β=0`)** | **`→0.222`** | **saturates — the physical case** |
| free `U(1)` | `≡0`, but `J=1` (no Jacobian) | **grows `≍L`** | `≍1/L²→0` | **decays (acid test)** |

The periodic class is a **third** class. On the open strip the end plaquettes `P_1=X_1X_2`,
`P_L=BX_L^{-1}` are not functions of the differences, so a uniform shift costs `O(1/L)` and the twist
`B` decides the sign; `θ=π` (`B=iσ₃`, `B²=−I`) is an **actionful** sector with `S_W=2β` per twist
plaquette, i.e. an explicitly center-broken boundary. On the periodic chain every plaquette is a pure
difference, the uniform shift is an exact symmetry of `S_W`, and the only pinning is the Haar measure.
Hence the physical periodic floor is the **Haar-carried** one, and the fixed-boundary `≈0.48` is a
boundary artifact, not the decompactification value.

---

## §2 · The exact anchors

**Global holonomy [T].** Every plaquette of a constant commuting gauge field is the identity:
`U·I·U†·I=I` (one direction constant) or `U·U·U†·U†=I` (all links equal), verified to `10⁻¹⁶`. The
marginal is pure Haar on one modulus, and

```
   G″_anch(0) = −ρ₁″(0)/ρ₁(0) = 1/(2·Cin(2π)) = 0.20511530  (SU(2)),   L-independent.
```

**Why not `1/Var` [T].** `Var_Haar(ξ₃)=(4π²−6)/9=3.719824` gives `1/Var=0.268830`; the two differ by
`0.064`, because the Haar marginal is bounded and non-Gaussian. The acid test decides: for the free
`U(1)` field `−ρ″/ρ=0` exactly, while `1/Var=3/π²=0.304` would be a positive floor for a free field.
So `1/Var` is the Gaussian/CLT limit, and every statement of the form `G″=1/Var` in the corpus is a
Gaussian-limit statement `[C]`, not an identity.

**Object A vs object B.** The physical soft mode is the `k=0` Fourier mode (object B), whose `β=0`
marginal is the `L`-fold convolution and whose curvature rises `0.2053→0.2688` as `L→∞` [T]:
`0.20532, 0.23966, 0.24940, 0.25431, 0.25918, 0.26160, 0.26402, 0.26522, 0.26703` for
`L=1,2,3,4,6,8,12,16,32`. The global holonomy (object A) is an exact anchor, not the soft-mode
marginal.

---

## §3 · The periodic `k=0` result (the physical number)

At `β=2.5`, two independent quadratures of `Tr[M_k^L]` agree:

| `L` | 2 | 3 | 4 | 5 | 6 | 7 | 8 |
|---|---|---|---|---|---|---|---|
| `G″(0)` | 0.2077 | 0.2097 | 0.2115 | 0.2130 | 0.2144 | 0.2155 | 0.2165 |
| `G″(β=0)` | 0.2395 | 0.2493 | 0.2542 | 0.2572 | 0.2591 | 0.2605 | 0.2616 |
| Wilson part | −0.0319 | −0.0397 | −0.0428 | −0.0441 | −0.0447 | −0.0450 | −0.0451 |

Fits: `c+Aq^L` gives `c=0.22196±0.00011`, `q=0.849`, `rms=3.3×10⁻⁶`; `c+A/L` gives `0.2126±0.0018`
(`L=2..8`) and `0.2142±0.0016` (`L=4..8`); the Wilson part fits a bounded negative
`c=−0.0461±0.0003` (`p=2.05`). The amplitude variance plateaus (`4.391→4.449`) and
`R=G″·Var=0.929→0.963`, so the invariant statement is `G″→1/χ>0`, `χ≈4.45`. The transfer matrix is
gapped: top eigenvalues `228.174, 119.788, 63.818, 33.401`, `Δ=0.6444`, `ξ=1.55`, `L`-independent.

**Why the periodic Wilson part is negative and harmless.** The flat direction has no end pinning, so
the Wilson action can act on the constant-mode marginal only through the measure-mediated correlation
of non-constant fluctuations. That correlation is negative (`−0.045`) but `O(1)`; it subtracts a
fixed amount from the Haar floor and never drives the total to zero. The total still rises toward
`≈0.22`.

---

## §4 · The center sector does not change the floor

The constant temporal mode is the Polyakov holonomy. `Z(N)` acts on temporal links; the Wilson action
is center-invariant (`z z̄=1` per temporal plaquette); the Polyakov loop has `N`-ality 1, while every
gauge-invariant local operator (`tr F²`, glueballs, the `0⁺⁺` channel) has `N`-ality 0 and is
center-blind. On the constant-holonomy slice every center sector is flat (`S_W≡0`), so every sector's
marginal is pure Haar and `G″=0.2051`. The full-3D `L=1` periodic measure is **exactly**
`θ`-independent (`∫ρ da=26.0980322` for `θ∈{0,π/2,π,3π/2,2π}`, spread `2.7×10⁻⁸`, right-Haar); the
flat-chart spread `0.730` vs `0.013` is a chart artifact (`J(2π)=0` at the `−I` edge). Center symmetry
is unbroken at finite `L` (`⟨tr P⟩=0`; gapped kernels `Δ=0.644, 3.176, 1.096`). So the physical
periodic floor is the center-neutral Haar-carried value `≈0.222`, bracketed by P1's exact `0.2051` and
the CLT `0.2688`; it is **not** a center-sector floor and **not** the `≈0.48` explicitly-broken twist.

---

## §5 · Corrections and errata

1. **`1/Var` is not the exact curvature [T].** The exact object is `−ρ″/ρ`; `1/Var` is its
   Gaussian/CLT limit. The free-field acid test (`0` vs `3/π²`) settles it. Numerically unchanged;
   logical status sharpened to `[C]` where `1/Var` is used.
2. **P5's `Var(a)` moment formula is invalid for `L≥4`.** At `β=0`, `Var(a)=3.7198` exactly for every
   `L` (independent links), which the formula violates from `L=4`; its apparent `SU(2)` growth is an
   artifact. Its independent `G″` quadrature and the P2/P4 solvers agree.
3. **`L=1` twist values in the flat chart are chart-dependent.** On the periodic minimal torus the
   measure is exactly `θ`-independent (right-Haar); the `0.73` vs `0.013` spread is the `J(2π)=0`
   edge. D1's `L≥2` twist fits are unaffected; the `L=1` flat-chart `θ`-spread should not be read as
   a physical twist effect.
4. **Boundary-class language.** D5's "saturation is twist-specific; at `b=0` it decays" is correct
   for the **fixed-boundary** classes, but the physical periodic class is a third class and it
   saturates. The correct statement is: *fixed-boundary strips can lose the floor; the periodic torus
   does not.*

---

## §6 · What this changes, and what it does not

**Changed:**
- The decompactification objection **in the periodic tractable model is removed**: the constant mode
  saturates at a positive `L`-independent floor, with an exact Haar anchor and a gapped transfer
  matrix. The earlier `≈1/L` was a fixed-boundary transient, not the periodic asymptote.
- The physical number is `≈0.222` (P2), bracketed by the exact `0.2051` (P1) and the CLT `0.2688`.
- The exact curvature object is `−ρ″/ρ`, not `1/Var`.

**Not changed:**
- **4D remains `H1` `[O]`.** The `k→0` susceptibility's uniform boundedness in 4D is the mass-gap
  input; the periodic strip proves the 1D analogue only.
- **The floor is not the gap.** A curvature floor is a susceptibility ceiling; at `β=0` it is
  positive while `f₁²=0`. The atom `f₁²≥r₀` with isolation remains the missing half (P6, `[T]`/`[NG]`).
- **No QED. Clay 0%.**

---

## §7 · Cleanest statements

> **[T · periodic anchor]** On the periodic torus the global constant-holonomy slice has `S_W≡0` in
> every center sector; its marginal is pure Haar and `G″=1/(2·Cin(2π))=0.20511530` (SU(2)),
> `L`-independent.

> **[T · flat direction]** On the periodic chain every plaquette is a pure difference, so the uniform
> Cartan shift is an exact symmetry of `S_W` for every `β,L`; the only pinning is the Haar measure.

> **[T]+[V · periodic saturation]** The periodic `k=0` curvature is positive and converges to
> `1/χ>0`, with `χ=Σ_d C(d)<∞` by the gapped transfer matrix (`Δ=0.6444`); the measured
> `c=0.22196±0.00011` (`q=0.849`) and `R=G″·Var→0.963` confirm the Gaussian/CLT identification.

> **[C · 4D]** If the 4D `k→0` susceptibility is bounded uniformly in `(L,a)` — which is `H1` —
> then the periodic 4D constant-mode curvature saturates at a positive floor. The converse is not
> claimed; and saturation does not give the atom.

**One line.** *The periodic constant mode — the physical decompactification case — saturates: the
global holonomy is exactly Haar with `G″=0.205115` [T], the `k=0` mode has `c=0.222` with a gapped
transfer matrix `Δ=0.6444`, the center sector does not change it, and the exact object is `−ρ″/ρ`,
not `1/Var`. 4D is still `H1`, and a curvature floor is never an atom. No QED. Clay 0%.*

*Refs:* `P1-PERIODIC-HOLONOMY.md`, `P2-PERIODIC-K0.md`, `P3-POLYAKOV-CENTER.md`,
`P4-PERIODIC-TRANSFER.md`, `P5-PERIODIC-DECISION.md`, `P6-PERIODIC-TO-H1.md`;
`D-WAVE-SYNTHESIS.md`, `SM-WAVE-SYNTHESIS.md`, `STAR2-CLASSIFICATION.md`; witnesses in `w7/`.
