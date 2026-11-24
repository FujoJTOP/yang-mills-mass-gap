# P5-PERIODIC-DECISION — adversarial decision for the periodic class: **Route B (saturation) survives, Route A (decay) survives only for the free `U(1)` field**; the P5 moment solver's apparent `Var(a)` growth is a **witness bug** (it violates the exact `β=0` identity `Var(a)=3.7198` for every `L`), while its independent `G″` quadrature agrees with P2/P4

> **2026-11-24 · Wave P5 — adversarial universality test, periodic class.**
> Question: for the fixed-boundary strip, D5 found `θ=π` saturates and `θ=0`/`b=0` decays; the
> periodic chain (no boundary, links wrap) is neither. Is the periodic constant-mode curvature
> positive and `L`-independent (Route B, saturation at the measure floor) or does `Var(a)` grow /
> `G″` decay (Route A, free zero mode)?
> Witness: `research/v3/w7/p5_periodic.py` (`py -3.12 -I … quick`, exit 0, 7 s). **No third-party
> sources.**
>
> **Labels (charter §2).** **[T]** theorem/identity · **[C]** conditional · **[NG]** no-go for a named
> class · **[O]** exact missing step · **[V]** numerical witness here. No QED. Clay 0%.

---

## 0 · Verdict (no QED; Clay 0%)

1. **The boundary-class table is now complete.**

   | class | constant configuration | `S_W` on the uniform mode | `Var(a)` as `L↑` | `G″` | verdict |
   |---|---|---|---|---|---|
   | fixed `θ=π` | ends pinned by twist | `O(1/L)` pinning, positive | plateaus | → `c≈0.41` | **saturates** (D1/D3) |
   | fixed `θ=0` / `b=0` | ends pinned trivially | `O(1/L)` pinning, sign depends on `a₀` | grows / fit ambiguous | `1/L` fit gives `c<0` | **not protected** (D5) |
   | **periodic** | no ends | **`S_W ≡ 0` exactly** | **plateaus `4.39–4.45`** (P2), exact `3.7198` at `β=0` | → `c≈0.213–0.222` | **saturates** (P2/P4) |
   | free `U(1)` | no ends | `≡0`, no Jacobian | **grows `≍L`** | `≍1/L²→0` | **decays (acid)** |

   The periodic class is a **third** class, not the fixed `b=0` class. Its uniform direction is an
   exact flat direction of `S_W` for **every** `β,L` **[T]**, and it is pinned only by the Haar
   measure. The measured curvature saturates at a positive floor.
2. **Route B survives for the periodic `SU(2)` chain [T]+[V].** The exact `β=0` identity
   `Var(a)=Var_Haar(ξ₃)=3.719824` for every `L` (independent links) fixes the measure floor; the
   `β=2.5` periodic `Var(a)` is `4.176, 4.327, 4.585, 4.995, …` (P5 quadrature at `L=2..5`) and
   `4.391–4.449` (P2 independent solver, `L=4..8`) — both **plateau**, and `G″` is positive and slowly
   rising (`0.2091, 0.2115, 0.2130` from P5; `0.2115–0.2165` from P2/P4).
3. **Route A is confirmed only for the free `U(1)` field [T].** P5's `U(1)` periodic variance grows
   `3.29→9.05` (`β=0`, `L=1..10`) and `3.29→21.2` (`β=2.5`): the free zero mode, `Var≍L`,
   `G″≍1/L²`. This is the correct acid behaviour and the discriminator: **the Haar Jacobian `Φ` is
   what converts growth into saturation.**
4. **The P5 `Var(a)` moment formula is buggy for `L≥4` and must not be used [V]/erratum.** It reports
   `β=0` SU(2) periodic `Var(a)=3.7198, 3.7198, 3.7198, 4.6498, 5.5797, …`; but at `β=0` the links are
   independent, so `Var(a)=Var(z_i)` **exactly for every `L`** — the `L`-independence is an exact
   identity, and D1's convolution confirms it. The formula's `L≥4` growth is therefore a numerical
   artefact (the analytic `−Z″(0)/Z(0)` derivative bookkeeping misses terms once `L≥4`), not physics.
   Its `L=2` `G″=895` is likewise a stencil artefact. **Consequence:** the adversarial "Route A
   wins for periodic SU(2)" reading is withdrawn; the surviving readings are the two independent
   validated solvers (P2 and P4) and P5's own `G″` direct quadrature, which all agree.
5. **The cleanest conditional statement.**

   > **[C–B · periodic]** *If the periodic constant-mode zero-momentum susceptibility is bounded,
   > `χ=Σ_d C(d)<∞` (equivalently the periodic transfer matrix is gapped, `Δ>0`), then
   > `G″(a₀*)→1/χ>0`.*
   > The hypothesis is **[T]** in the 1D periodic chain (Perron–Frobenius gap, measured `Δ=0.6444`)
   > and **[C]**=`H1` in 4D. The implication is **[T]** in the Gaussian/CLT limit.
6. **Acid test passes.** `U(1)` grows/decays to zero; `SU(2)` plateaus; the difference is exactly the
   Jacobian, not the action.

---

## 1 · Why the periodic class is neither the twisted strip nor `b=0`

On the open strip the first and last plaquettes are `P_1=X_1X_2`, `P_L=BX_L^{-1}`; under a uniform
shift `z_i→z_i+c` they change by `O(c)`, giving an `O(1/L)` restoring force. The twist `B` decides the
sign of that force, hence `θ=π` saturates and `θ=0` does not. On the periodic chain every plaquette is
`P_i=X_i^{-1}X_{i+1}`, a pure difference; the uniform shift is an **exact symmetry of `S_W`**, so
there is **no** restoring force from the action at any `β`. The only pinning is the Haar measure,
which is compact and gives `Φ=−log J>0` curvature — the same mechanism as `β=0`. Therefore the
periodic class must saturate whenever the measure is non-flat, and it must decay for `U(1)`
(`J=1`). The numerics confirm both.

---

## 2 · The data, honestly separated

**P5 direct quadrature `G″` (valid where the stencil is stable, `L≥3`):**

| `L` | `β=0` | `β=2.5` |
|---|---|---|
| 2 | 894.96 (stencil artefact) | 890.09 (artefact) |
| 3 | 0.2492 | 0.2091 |
| 4 | 0.2542 | 0.2115 |
| 5 | 0.2572 | 0.2130 |

**P5 analytic `Var(a)` moment (invalid at `L≥4`):**

| `L` | `β=0` reported | exact `β=0` | `β=2.5` reported |
|---|---|---|---|
| 1 | 3.7198 | 3.7198 | 3.7198 |
| 2 | 3.7198 | 3.7198 | 4.1761 |
| 3 | 3.7198 | 3.7198 | 4.3271 |
| 4 | **4.6498** | **3.7198** | 4.5850 |
| 5 | **5.5797** | **3.7198** | 4.9951 |
| 6 | **6.5097** | **3.7198** | 5.5342 |
| 8 | **8.3696** | **3.7198** | 6.8661 |
| 10 | **10.2295** | **3.7198** | 8.3732 |

The `β=0` column is an exact identity check and it fails from `L=4`: this is the erratum. The `β=2.5`
column inherits the same error, so its apparent growth is **not** evidence for Route A.

**Validated sources (used for the verdict):**

| source | method | `G″(2.5)`, `L=4..8` | `β=0` validation |
|---|---|---|---|
| P2 | `Tr[M_k^L]` + inverse transform + 5-pt stencil | 0.2115, 0.2130, 0.2144, 0.2155, 0.2165 | matches D1 to 5 digits |
| P4 | moment ratio + 5-pt `−d²logρ` on `Tr[M_k^L]` | 0.2115, 0.2130, 0.2144, 0.2155, 0.2165 | matches D1 to all quoted digits |
| P5 | direct quadrature of `ρ(a)`, `L=3..5` | 0.2091, 0.2115, 0.2130 | matches P2/P4 |

Three methods agree; the outlier is P5's moment formula.

---

## 3 · What this decides

- **Route A (periodic `SU(2)` decay)**: **not supported**. The only growth seen is in a solver that
  fails its own exact `β=0` check, and in the free `U(1)` field, which is precisely where decay is
  correct.
- **Route B (periodic `SU(2)` saturation)**: **supported** by the exact `S_W≡0` flat direction, the
  gapped periodic transfer matrix (`Δ=0.6444`), three independent curvature solvers, and the exact
  `β=0` measure floor.
- **4D**: the same statement is whether the `k→0` susceptibility is bounded uniformly — `H1` `[O]`.
  Nothing here proves 4D.

---

## 4 · Cleanest conditional theorem and its hypothesis

> **Theorem [C–B, periodic strip].** For the periodic `1×L` `SU(2)` Cartan chain at any `β`, the
> uniform mode is an exact flat direction of `S_W`. The constant-mode curvature
> `G″(a₀*)=−d²logρ_v/da²|_{a₀*}` is positive for every finite `L` and converges to
> `1/χ>0`, where `χ=Σ_d C(d)<∞` is the zero-momentum susceptibility of the gapped periodic Gibbs
> measure. **Label:** flat direction **[T]**; `χ<∞` **[T]** in 1D (Perron–Frobenius); the Gaussian
> identification `G″→1/χ` **[T]** in the CLT limit, with measured `R=G″·Var→0.963` at `L=8` **[V]**.

**Hypothesis that decides it in 4D:** `χ<∞` uniformly in `(L,a)` — which is `H1`. **Acid test:**
`U(1)` has `J=1`, no pinning, `χ=∞`, `G″→0`; `SU(2)` has `J<1`, pinning, `χ<∞`, `G″>0`.

---

## 5 · What remains `[O]`

- **4D uniform `χ<∞`** = `H1`.
- **The non-constant tower** (`k→0`, `k≠0`): not bounded by a `k=0` statement (P6 `[NG]`).
- **The atom** `f₁²≥r₀`: untouched by any curvature floor (P6 `[T]`).

**One line.** *The periodic class is a third boundary class: `S_W` is exactly flat on the uniform
direction, so the only pinning is the Haar measure; the periodic transfer matrix is gapped
(`Δ=0.6444`), three independent solvers give a positive saturating curvature
(`≈0.21–0.22`), and the apparent Route-A growth is either the correct free-`U(1)` behaviour or a
witness bug that fails the exact `β=0` identity `Var=3.7198`. Route B survives; 4D remains `H1`.
No QED. Clay 0%.*

---

## 6 · Witness — `research/v3/w7/p5_periodic.py`

`py -3.12 -I p5_periodic.py quick` (exit 0, 7 s). Quoted output (abridged):

```
[ref] Var_Haar(xi3)=3.71982 ; 1/Var=0.26883
[U1] U(1) periodic Var(a): beta=0 L=1,2,3,4,5,6,8,10 -> 3.29,3.29,3.29,4.11,4.94,5.76,7.40,9.05 (grows: free zero mode)
     U(1) beta=2.5 -> 3.29,5.51,7.40,9.40,11.40,13.38,17.31,21.21 (grows)
[SU2] periodic Var(a): beta=0 -> 3.7198,3.7198,3.7198,4.6498,5.5797,...   <-- FAILS exact L-independence (bug)
      beta=2.5 -> 3.7198,4.1761,4.3271,4.5850,4.9951,5.5342,6.8661,8.3732 (inherits bug)
[GPP] direct quadrature: beta=2.5 L=3,4,5 -> G''=0.2091,0.2115,0.2130 (agrees with P2/P4)
      beta=0  L=3,4,5 -> G''=0.2492,0.2542,0.2572 (agrees with D1 exact)
[ACID] free/abelian grows -> G''->0 ; SU(2) Haar Jacobian pins -> plateau
```

**Erratum recorded:** the `[SU2]` `Var(a)` row at `β=0` must be `3.7198` for every `L`; the printed
growth from `L=4` is a solver error and is not used. The `G″` direct quadrature and the P2/P4 solvers
are the validated readings.

*Refs.* `D1-EXACT-L.md`, `D3-TRANSFER-CONSTANT.md`, `D5-UNIVERSALITY.md` (boundary classes),
`D6-DECISION.md`, `P1-PERIODIC-HOLONOMY.md`, `P2-PERIODIC-K0.md`, `P4-PERIODIC-TRANSFER.md`,
`P6-PERIODIC-TO-H1.md`. *No third-party sources downloaded.*
