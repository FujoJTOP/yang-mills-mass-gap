# P2-PERIODIC-K0 — the periodic `k=0` fluctuation mode **saturates** at a positive `L`-independent floor: exact `Tr[M_k^L]` gives `G″(0)=0.2115→0.2165` for `L=4..8` at `β=2.5`, fit `c=0.22196±0.00011` with `q=0.849`, `Var(a)` plateau `4.45`, `G″·Var→0.963`; the fixed-boundary decay is **not** the periodic behaviour

> **2026-11-24 · Wave P2 — the periodic constant mode, physical case.**
> Question: on a periodic torus (no fixed external holonomy), does the curvature of the `k=0` Fourier
> mode of the gauge-field fluctuations saturate (`c>0`) or decay (`c=0`) as `L` grows? D1 computed the
> **fixed-boundary** strip and got `c=0.248` (`θ=0`), `0.41` (`θ=π`); D5 reported that saturation is
> twist-specific and that at fixed `b=0` the variance grows. The periodic chain is a third class.
> Witness: `research/v3/w7/p2_periodic_k0.py` (`py -3.12 -I … quick`, 28.7 s, exact/deterministic, no
> MC). **No third-party sources.**
>
> **Labels (charter §2).** **[T]** theorem/identity · **[C]** conditional · **[NG]** no-go for a named
> class · **[O]** exact missing step · **[V]** numerical witness here. No QED. Clay 0%.

---

## 0 · Verdict (no QED; Clay 0%)

**The periodic `k=0` mode saturates at a positive floor [T]+[V].**

1. **The uniform Cartan direction is an exact zero mode of `S_W` on the periodic torus [T].** Every
   plaquette is a pure difference, `P_i = X_i^{-1}X_{i+1}` with `X_{L+1}=X_1`; the uniform shift
   `z_i → z_i + c` leaves every `Δz_i` invariant, so `S_W` is **flat** along the `k=0` direction for
   every `β`, `L`. The only pinning is the Haar measure. This is the structural difference from the
   fixed-boundary strip, where the end plaquettes `P_1=X_1X_2`, `P_L=BX_L^{-1}` carry an `O(1/L)`
   restoring force.
2. **The floor is positive and `L`-independent [V].** At `β=2.5`, exact `Tr[M_k^L]`:
   ```
      L        4        5        6        7        8
      G″(0)  0.21146  0.21304  0.21439  0.21553  0.21649
   ```
   Every fit form gives `c>0`: `c+A/L` → `c=0.22126`, `c+A/L²` → `0.21765`, `c+Aq^L` →
   **`c=0.22196±0.00011`**, `q=0.849`, `rms=3.3×10⁻⁶`. The approach is exponential, not a power.
3. **The marginal Gaussianises [V].** `Var(a) = 4.391, 4.422, 4.437, 4.445, 4.449` for `L=4..8`
   (plateau), and `R := G″·Var = 0.929, 0.942, 0.951, 0.958, 0.963 → 1`. So the invariant statement is
   `G″→1/χ > 0` with `χ ≈ 4.45`, exactly the D5/D6 dichotomy, now in the periodic class.
4. **The `β=0` validation is exact.** With `M_0` rank one, `Z_0(k)=[φ(k/√L)]^L`, and the witness
   reproduces D1's exact Haar values `0.23955, 0.25425, 0.25914, 0.26157` (`L=2,4,6,8`) **to the
   digit** — the transfer-matrix solver is validated against the known exact convolution.
5. **The fixed-boundary classes are different [V].** The same extraction on the D1 strip gives
   `θ=0`: `c=−0.033` under `c+A/L` (negative intercept, i.e. no positive floor in that fit) and
   `c=+0.235` under `c+A/L²` with `rms 2.8×10⁻³`; `θ=π`: `c=+0.148` / `0.308` with `rms ≈9×10⁻²`.
   The periodic class is the only one whose exponential fit is clean to `10⁻⁶` and whose `Var(a)`
   plateaus. So **the periodic physical case saturates**, and the fixed-boundary `b=0` decay found in
   D5 is a boundary artefact, not the periodic behaviour.
6. **Acid test passes [T].** For a free massive field, `k=0` is the Laplacian zero mode,
   `Var(a)=1/m²` independently of `β,L`, and `G″=m²→0` as `m→0` (witness: `1.00, 0.090, 0.010,
   0.0001` for `m=1,0.3,0.1,0.01`). The tool manufactures no floor for the free field.

**What is still `[O]`.** This is a **1D-strip** theorem/number. In 4D the same question is whether the
`k→0` susceptibility is bounded uniformly in `L` — which is `H1`. The periodic strip result removes the
"decompactification kills the floor" worry in the tractable model; it does not prove the 4D floor.

---

## 1 · The model and why the periodic chain is flat

Pure-Cartan sector of SU(2), flat chart `U=e^{izσ₃/2}`, per-link **exact** one-component Haar weight
`ρ₁(z)=π(Cin(2π)−Cin|z|)` (transverse components integrated exactly), Wilson action the Cartan
projection of the cyclic plaquettes:

```
   S_W = β Σ_{i=1}^{L} (1 − cos((z_{i+1}−z_i)/2)),   z_{L+1} = z_1 .
```

For the open strip (D1) the first and last plaquettes are `P_1=X_1X_2`, `P_L=BX_L^{-1}`; these are
**not** functions of `z_{i+1}−z_i`, so the uniform shift is pinned at the ends with an `O(1/L)` force.
For the periodic chain every plaquette is a difference, so

```
   z_i → z_i + c  ⇒  S_W invariant,  for every β, L.      [T]
```

The `k=0` mode is `a = L^{-1/2}Σ_i z_i`; the partition function at fixed `a` is
`Z(k)=Tr[M_k^L]` with

```
   M_k(z,z′) = √(ρ₁(z)ρ₁(z′)) e^{−ik(z+z′)/(2√L)} e^{−β(1−cos((z′−z)/2))} ,
```

and `ρ_v(a)=(1/2π)∫Z(k)e^{ika}dk`, `G″(a)=−d_a² log ρ_v(a)`. By `z→−z` symmetry `a₀*=0` and
`G″(0)` is the exact moment ratio. The computation is deterministic (Gauss–Legendre nodes `nz=110`
in quick mode; a 5-point stencil on `log ρ`), no Monte Carlo.

---

## 2 · The data

**Periodic, `β=2.5`:**

| `L` | `G″(0)` | `ρ(0)` | `G″_{β=0}` (Haar) | `Var(a)` | `1/Var` | `R=G″·Var` |
|---|---|---|---|---|---|---|
| 4 | 0.21146 | 2.490e4 | 0.25425 | 4.3913 | 0.22772 | 0.9286 |
| 5 | 0.21304 | 4.488e5 | 0.25718 | 4.4217 | 0.22616 | 0.9420 |
| 6 | 0.21439 | 8.248e6 | 0.25914 | 4.4368 | 0.22539 | 0.9512 |
| 7 | 0.21553 | 1.531e8 | 0.26053 | 4.4445 | 0.22500 | 0.9579 |
| 8 | 0.21649 | 2.854e9 | 0.26157 | 4.4485 | 0.22479 | 0.9631 |

`ρ` is even and peaked at `a=0` (`argmax a=0.000`; `ρ(±0.6)/ρ(0)=0.974` at `L=4`, `0.973` at `L=8`).
The robustness scan over `(kw,δ,K,nz)` gives **zero spread** at fixed `L` (`rel 0.00%`), so the
numbers are quadrature-stable.

**Fits (periodic, `β=2.5`):**

| form | `c` | `A` | `q` | `rms` |
|---|---|---|---|---|
| `c+A/L` | **+0.22126** | −0.0400 | — | 1.9e−4 |
| `c+A/L²` | **+0.21765** | −0.1041 | — | 3.6e−4 |
| `c+A/L³` | **+0.21647** | −0.3450 | — | 5.2e−4 |
| `c+Aq^L` | **+0.22196±0.00011** | −0.0202 | 0.849 | **3.3e−6** |

The exponential form wins by two orders of magnitude in `rms`; the approach to the floor is
`q^L` with `q≈0.85` (correlation length `ξ≈1/|log q|≈6.2` links), not a power law. This is the
signature of a **gapped** 1D transfer matrix.

**Fixed-boundary replica (same extraction):**

| `L` | `θ=0` | `θ=π` |
|---|---|---|
| 4 | 0.74401 | 0.73513 |
| 5 | 0.55534 | 0.37086 |
| 6 | 0.45774 | 0.45947 |
| 7 | 0.40116 | 0.45672 |
| 8 | 0.36529 | 0.44098 |

Fits: `θ=0` `c+A/L` → `c=−0.033` (`rms 1.6×10⁻²`), `c+A/L²` → `c=+0.235` (`rms 2.8×10⁻³`);
`θ=π` `c=+0.148` (`rms 9.1×10⁻²`) / `c=+0.308` (`rms 8.4×10⁻²`). These are the D1 numbers; the
periodic class is the clean one.

---

## 3 · Interpretation

1. **The flat direction is pinned only by the measure.** Because `S_W` is exactly flat on the `k=0`
   direction, the `β→∞` limit of the marginal is not a delta but the convolution of Haar measures;
   the curvature floor is `1/χ`, where `χ` is the zero-momentum susceptibility of the 1D periodic
   Gibbs measure. The 1D periodic measure has a spectral gap (`q=0.849<1`), so `χ<∞` and the floor is
   positive **[T]** (Perron–Frobenius/1D mixing).
2. **The Wilson part is negative but bounded.** Comparing to the Haar column, `G″(2.5)−G″(β=0) =
   −0.043, −0.044, −0.045, −0.045, −0.045`: the periodic Wilson part is **negative** (unlike the
   strip, where the end holonomy gives a positive Wilson part), but it is `O(1)` and does **not**
   overcome the Haar floor. The total still saturates at `≈0.22 > 0`.
3. **This resolves the D5 caveat for the periodic class.** D5's statement "saturation is
   twist-specific; at `b=0` the variance grows" refers to the fixed-boundary `b=0` strip, whose
   variance does grow (the end pinning is removed but the measure is not periodic). The **periodic**
   chain is the third class, and it saturates. The physical decompactification model is periodic, so
   the correct caveat is: *the fixed-boundary strip decays; the periodic torus saturates.*
4. **Not an atom.** A positive curvature floor says `χ` is finite; it says nothing about `f₁²`. At
   `β=0` the periodic floor is `0.2688` while `f₁²=0`. P6 makes this precise.

---

## 4 · Acid test

Free massive field, constant-mean constraint, `k=0` is the Laplacian zero mode:
`Var(a)=1/m²` independently of `β` and `L`, `G″=m²`. Witness:

```
   m = 1.00, 0.30, 0.10, 0.01  ->  G″(free) = 1.0000, 0.0900, 0.0100, 0.0001  -> 0.
```

The tool returns `0` as `m→0`: no floor is manufactured for the free field. **[T]**

---

## 5 · What remains `[O]`

- **4D.** The periodic strip has a gap and a positive floor; the 4D periodic torus asks whether the
  `k→0` susceptibility is bounded uniformly in `L,a`. That is `H1` `[O]`.
- **The non-constant tower.** This computation is the `k=0` mode only; the `k→0` soft tower is not
  bounded by it (P6, `[NG]`).
- **The atom.** Curvature ≠ atom, at every `β` (the `β=0` counterexample).

**One line.** *On the periodic chain the `k=0` mode is an exact flat direction of `S_W`, pinned only by
the Haar measure; exact `Tr[M_k^L]` gives `G″(0)=0.2115→0.2165` for `L=4..8`, exponential saturation
`c=0.22196±0.00011` (`q=0.849`), `Var(a)` plateau `4.45`, `G″·Var→0.963`; the fixed-boundary decay is
a boundary artefact. The 4D version of the same statement is `H1`. No QED. Clay 0%.*

---

## 6 · Witness — `research/v3/w7/p2_periodic_k0.py`

`py -3.12 -I p2_periodic_k0.py quick` (exit 0, 28.7 s; `nz=110`). Quoted output (abridged):

```
[0] Var(xi3)=3.719824 -> Haar floor 1/Var=0.268830
[V0] beta=0 periodic validation vs D1: L=2,4,6,8 -> 0.23955,0.25425,0.25914,0.26157 (exact)
[P] beta=2.5 periodic: L=4..8 -> G''=0.21146,0.21304,0.21439,0.21553,0.21649
[U] rho even, argmax a=0.000; rho(+-0.6)/rho(0)=0.974 (L=4), 0.973 (L=8)
[F] fixed boundary: theta=0 -> 0.744,0.555,0.458,0.401,0.365 ; theta=pi -> 0.735,0.371,0.459,0.457,0.441
[X] Var(a)=4.391..4.449 ; 1/Var=0.2277..0.2248 ; R=G''Var=0.929..0.963
[R] robustness: zero spread at fixed L (rel 0.00%)
[FIT] c+A/L: c=+0.22126 rms1.9e-4 ; c+A/L^2: c=+0.21765 ; c+A/L^3: c=+0.21647 ;
      c+Aq^L: c=+0.22196+-0.00011 A=-0.0202 q=0.849 rms=3.3e-6
[ACID] free field: m=1,.3,.1,.01 -> G''=1,.09,.01,.0001 -> 0
```

*Refs.* `D1-EXACT-L.md` (fixed-boundary strip, exact Haar validation), `D5-UNIVERSALITY.md`
(boundary-class dependence), `D6-DECISION.md` (dichotomy `χ<∞⟺saturation`), `P1-PERIODIC-HOLONOMY.md`
(exact Haar anchor, exact object `−ρ″/ρ`), `P6-PERIODIC-TO-H1.md` (logical strength).
*No third-party sources downloaded.*
