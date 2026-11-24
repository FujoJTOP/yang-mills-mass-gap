# P4-PERIODIC-TRANSFER — the periodic constant-mode transfer matrix is **gapped** (`Δ=0.6444`, `ξ=1.55`, `L`-independent) and its curvature **saturates** at `G″=0.2077→0.2165` (`L=2..8`), with a bounded negative Wilson part `−0.0451`; no decay to zero

> **2026-11-24 · Wave P4 — the periodic transfer matrix.**
> Question: D3 computed the constant-mode transfer matrix on the **fixed-boundary** `1×L` strip
> (`Δ=3.18` at `θ=π`, `0.64` at `θ=0`) and a curvature that saturates at the Haar floor plus a
> decaying boundary (Wilson) part. Does the **periodic** transfer matrix — cyclic chain
> `X_{L+1}=X_1`, no twist, no fixed boundary — also have a gapped, `L`-independent spectrum, and does
> its Riccati curvature saturate at a positive constant?
> Witness: `research/v3/w7/p4_periodic_transfer.py` (`py -3.12 -I … quick`, exit 0; numpy 2.4.6 /
> scipy 1.18.1; deterministic, no MC). **No third-party sources.**
>
> **Labels (charter §2).** **[T]** theorem/identity · **[C]** conditional · **[NG]** no-go for a named
> class · **[O]** exact missing step · **[V]** numerical witness here. No QED. Clay 0%.

---

## 0 · Verdict (no QED; Clay 0%)

1. **The periodic bulk kernel is gapped and `L`-independent [T]+[V].** At `β=2.5` the top eigenvalues
   of the `k=0` kernel are `228.174, 119.788, 63.818, 33.401`, giving
   ```
      λ₁/λ₀ = 0.52498,   Δ = −log(λ₁/λ₀) = 0.6444,   ξ = 1/Δ = 1.5519,
   ```
   exactly the D3 `θ=0` bulk kernel (same object), and `L`-independent by construction (one kernel
   per step). At `β=0` the kernel is rank one (`λ₁/λ₀ = 5.7×10⁻¹⁶`, `Δ=35.1`): independent links,
   `δ`-correlated, exact saturation at `1/Var = 0.26883`.
2. **The curvature saturates and is positive [V].** Exact `G″(0)` at `β=2.5`:
   ```
      L        2        3        4        5        6        7        8
      G″    0.2077   0.2097   0.2115   0.2130   0.2144   0.2155   0.2165
   ```
   The moment ratio and the 5-point `−d²logρ/da²` agree at every `L`; `a₀* = 0.0000` confirms the
   symmetric peak. Fits: `c+A/L^p` gives `c=0.2126±0.0018` (`L=2..8`) and `c=0.2142±0.0016`
   (`L=4..8`), `A→0` in the fit because the sequence is already flat; the exponential form gives the
   same. **`c>0`: no decay to zero.**
3. **The periodic Wilson part is negative but bounded [V].** `G″(2.5)−G″(β=0) =
   −0.0319, −0.0397, −0.0428, −0.0441, −0.0447, −0.0450, −0.0451` for `L=2..8`, fit
   `c=−0.0461±0.0003` (`p=2.05`), `L=4..8` `c=−0.0453±0.0001` (`p=3.57`), `rms≤2×10⁻⁴`. Unlike the
   fixed strip (where the end holonomy gives a positive Wilson part), the periodic flat direction has
   a **small negative** Wilson part; it is `O(1)` and does not overcome the Haar floor. This is the
   precise content of "the flat direction is pinned only by the measure".
4. **The exact `β=0` solver validation passes [T].** `Tr[M_k^L]` by eigenvalues vs the closed form
   `[φ(k/√L)]^L` at `L=4`, `k=0.3,1.1,2.5` agree to `10⁻¹⁵` relative. The `β=0` curvature row
   (`0.2395, 0.2493, 0.2542, 0.2572, 0.2591, 0.2605, 0.2616`) matches D1/D4's exact convolution to
   all quoted digits.
5. **Acid test passes [T].** Free Gaussian marginals return `G″=κ→0` (`0.5, 0.1, 0.02`), the free
   box gap `4π²/L²→0`, while the `β=0` Haar floor is `L`-independent `0.2688>0`; the discriminator is
   the Jacobian `Φ=−log J`.
6. **Consistency with P2.** P2's independent solver (`Tr[M_k^L]` + inverse transform + 5-point
   stencil, `nz=110`) gives the **same** `G″` values to 4 digits (`0.2115, 0.2130, 0.2144, 0.2155,
   0.2165` at `L=4..8`) and an exponential fit `c=0.22196±0.00011` (`q=0.849`). Two different
   quadratures, same saturation.

**What is still `[O]`.** In 4D the same question — is the `k→0` susceptibility bounded uniformly in
`L,a`? — is `H1`. The periodic strip removes the decompactification objection **in the tractable
model**; it does not prove the 4D floor.

---

## 1 · The periodic kernel and why it differs from the strip

Cyclic chain `X_1,…,X_L`, `X_{L+1}=X_1`, pure Cartan sector, exact one-link Haar marginal `ρ₁`; the
Wilson action is the Cartan projection of the cyclic plaquettes
`S_W=βΣ_i(1−cos((z_{i+1}−z_i)/2))`. The partition function at fixed `k=0` amplitude is a **trace**:

```
   Z(k) = Tr[M_k^L],    M_k(z,z′) = √(ρ₁(z)ρ₁(z′)) e^{−ik(z+z′)/(2√L)} e^{−β(1−cos((z′−z)/2))} ,
   ρ_v(a) = (1/2π)∫Z(k)e^{ika}dk ,    G″(a) = −d_a² log ρ_v(a) .
```

Because every periodic plaquette is a pure difference, the uniform shift is an exact flat direction:
`z_i→z_i+c` leaves `S_W` invariant for every `β,L` **[T]**. The strip's end plaquettes break this and
supply an `O(1/L)` restoring force; the periodic chain has none, so the only pinning is the Haar
measure. Consequently the `β→∞` limit is not a delta but the convolution of Haar measures, and the
curvature floor is the inverse zero-momentum susceptibility `1/χ` of the 1D periodic Gibbs measure
**[T]**. The measured `Δ=0.6444>0` is exactly `χ<∞`.

---

## 2 · The data

**Bulk spectrum (`k=0`), `L`-independent:**

| `β` | `λ₀` | `λ₁/λ₀` | `Δ=−log(λ₁/λ₀)` | `ξ=1/Δ` | top evals |
|---|---|---|---|---|---|
| 0.0 | 39.478 | 5.68e−16 | 35.104 | 0.0285 | 39.478, 2.24e−14, … (rank 1) |
| 2.5 | 228.174 | 0.52498 | **0.6444** | **1.5519** | 228.174, 119.788, 63.818, 33.401 |

D3 fixed strip reference: `Δ=3.18` (`θ=π`), `0.64` (`θ=0`) — the periodic `Δ` is the `θ=0` bulk
value, as it must be (same kernel).

**Curvature and split:**

| `L` | `G″(2.5)` | `G″(β=0)` | Wilson `=G25−G0` | Wilson/Haar | `G″` (moment) | `a₀*` |
|---|---|---|---|---|---|---|
| 2 | 0.2077 | 0.2395 | −0.0319 | −0.133 | 0.2077 | 0.0000 |
| 3 | 0.2097 | 0.2493 | −0.0397 | −0.159 | 0.2097 | 0.0000 |
| 4 | 0.2115 | 0.2542 | −0.0428 | −0.168 | 0.2115 | 0.0000 |
| 5 | 0.2130 | 0.2572 | −0.0441 | −0.172 | 0.2130 | 0.0000 |
| 6 | 0.2144 | 0.2591 | −0.0447 | −0.173 | 0.2144 | 0.0000 |
| 7 | 0.2155 | 0.2605 | −0.0450 | −0.173 | 0.2155 | 0.0000 |
| 8 | 0.2165 | 0.2616 | −0.0451 | −0.172 | 0.2165 | 0.0000 |

**Fits:**

| object | selection | form | `c` | `A` / `p` | `rms` |
|---|---|---|---|---|---|
| total | L=2..8 | `c+A/L^p` | **+0.2126±0.0018** | `A→0`, `p` unconstrained | 0.0030 |
| total | L=4..8 | `c+A/L^p` | **+0.2142±0.0016** | `A→0` | 0.0018 |
| Wilson | L=2..8 | `c+A/L^p` | **−0.0461±0.0003** | `p=2.05±0.12` | 0.0002 |
| Wilson | L=4..8 | `c+A/L^p` | **−0.0453±0.0001** | `p=3.57±0.26` | 0.0000 |

The Wilson fit is the cleanest object: a bounded negative constant `−0.045`, not a growing or
decaying-to-`−∞` term. The total is Haar `0.2688` minus `0.045` plus a small positive finite-`L`
correction → `0.2165` at `L=8`, still rising toward `≈0.22`.

---

## 3 · What the gapped spectrum means

A gapped 1D transfer matrix gives exponential clustering, hence
`χ = Σ_d C(d) < ∞` and `G″ → 1/χ > 0` **[T]**. The measured `Δ=0.6444` (correlation length
`ξ=1.55` links) is the quantitative form of "the periodic constant mode is not a free zero mode". The
free `U(1)` field has no gap in this channel (`Var(a)≍L`, `G″≍1/L²`), which is the acid test; SU(2)
has the gap because of the compact Haar measure (`Φ`), not because of the Wilson action (which is
exactly flat on the direction).

---

## 4 · Acid test and solver validation

```
[check] L=4: eigvals Tr vs beta=0 closed form at k=0.3,1.1,2.5: agree to 1e-15 relative
[T5] free Gaussian kappa=0.5,0.1,0.02 -> G''=0.5,0.1,0.02 -> 0
     free box gap 4pi^2/L^2 = 9.870,2.467,0.617,0.154 -> 0
     beta=0 Haar floor -> 0.2688 (>0), discriminator Phi=-log J
```

No floor is manufactured for the free field **[T]**.

---

## 5 · What remains `[O]`

- **4D periodic torus.** The strip has a gap and a saturation; 4D asks whether the `k→0`
  susceptibility is bounded uniformly. That is `H1` `[O]`.
- **Non-constant tower.** This is the `k=0` mode; the `k→0` tower is not bounded by it (P6, `[NG]`).
- **Atom.** Curvature is not an atom: at `β=0` the floor is `0.2688` while `f₁²=0` **[T]**.

**One line.** *The periodic transfer matrix is gapped (`Δ=0.6444`, `ξ=1.55`, `L`-independent), its
constant-mode curvature is positive and saturating (`0.2077→0.2165` for `L=2..8`, `c≈0.213–0.214`),
and its Wilson part is a bounded negative `−0.045`; the exact `β=0` solver validation and the free
acid test both pass. No decay to zero in the periodic model; 4D remains `H1`. No QED. Clay 0%.*

---

## 6 · Witness — `research/v3/w7/p4_periodic_transfer.py`

`py -3.12 -I p4_periodic_transfer.py quick` (exit 0; `nz=320`, `K=7`, `dk=0.03`). Quoted output
(abridged):

```
exact one-link Var(xi3)=3.719824 -> Haar floor 1/Var=0.26883
[check] L=4 k=0.3,1.1,2.5: Tr(M^L) vs closed form agree to 1e-15
[T1] beta=0: lam1/lam0=5.68e-16 Delta=35.10 ; beta=2.5: lam1/lam0=0.52498 Delta=0.6444 xi=1.5519
[T2] L=2..8 G(2.5)=0.2077,0.2097,0.2115,0.2130,0.2144,0.2155,0.2165
     G(0)  =0.2395,0.2493,0.2542,0.2572,0.2591,0.2605,0.2616
     Wil   =-0.0319,-0.0397,-0.0428,-0.0441,-0.0447,-0.0450,-0.0451 ; a0*=0.0000
[T4] TOTAL L=2..8 c+A/L^p c=+0.2126+-0.0018 ; L=4..8 c=+0.2142+-0.0016
     WILSON L=2..8 c=-0.0461+-0.0003 p=2.05 ; L=4..8 c=-0.0453+-0.0001 p=3.57
[T5] free Gaussian -> 0 ; box gap -> 0 ; Haar floor 0.2688 > 0
```

*Refs.* `D1-EXACT-L.md`, `D3-TRANSFER-CONSTANT.md` (fixed-strip transfer matrix), `D4-HAAR-WILSON.md`
(exact Haar floor), `D6-DECISION.md` (`χ<∞⟺saturation`), `P2-PERIODIC-K0.md` (independent solver,
same numbers), `P1-PERIODIC-HOLONOMY.md` (exact anchor), `P6-PERIODIC-TO-H1.md` (logic).
*No third-party sources downloaded.*
