# Q2-LOOP-BOOTSTRAP — the exact finite-`N` loop equation, `χ = χ_G + χ_conn`, and why **reflection positivity + UV asymptotics do not close a bootstrap on `χ`**: the connected 4-point enters with the *wrong sign* (`G″ = S − D`, `D ≥ 0`), 2D closes exactly (factorization) and the free field diverges (acid pass); the missing term is an **upper** bound on `D` (the Savvidy rank-one part) `= H1`

> 2026-11-24 · Wave **Q2** (loop equations + moment bootstrap). Attack on the last open input of the
> campaign, `H1`: `χ = Var_μ(ξ_{k→0}) ≤ 1/λ` uniformly in `(L,a)` at `β=2.5`, equivalently
> `G″_{k=0} ≥ λ > 0` (`paper-v2/ch02-squeeze.tex` eq. (H1); `paper-v3/main.tex` §6).
> Reads: `P-WAVE-SYNTHESIS.md`, `D6-DECISION.md`, `H1-RG-DESCENT.md`, `H2-DEFECT-BOUND.md`,
> `C2-BALABAN-UNIFORM.md`, `docs/paper-v2/ch02-squeeze.tex`, `docs/paper-v3/main.tex` §6; and the
> corpus's prior loop/bootstrap faces `PNT-LARGEN-LOOP`, `E3-LOOPSPACE`, `CONSTRUCT-LARGEN-BOOTSTRAP`,
> `A4-SPECTRAL-BOOTSTRAP`, `SLOT-BOOTSTRAP`, `SLOT-DSE-FIXPOINT`, `GLUEBALL-BOOTSTRAP`,
> `FORGE5-BOOTSTRAP-EXCLUDE`, `BROWNIAN-LOOP-SOUP` (all in `research/unconditional/`).
> Witness: `research/v3/w7/q2_loop_bootstrap.py` (`py -3.12 -I`, numpy 2.4.6 + scipy 1.18.1; output
> quoted §7). **No third-party sources downloaded.**
>
> **Labels (charter §2).** **[T]** true theorem/identity · **[C]** conditional · **[NG]** no-go for a
> named class · **[O]** exact missing step · **[V]** numerical witness here. No QED. Clay 0%.

---

## 0 · Verdict (no QED; Clay 0%)

**The loop equation is exact and, at finite `N`, is *closed in `W(C)`* up to a connected 4-point
term of order `1/N²`; it decomposes `χ` into a Gaussian/quadratic term plus a connected 4-point term
(equivalently `G″(k=0) = S − D`). Reflection positivity plus the known UV asymptotics do NOT close a
bootstrap inequality on `χ` at `β=2.5` — and the failure is structural, not a missing estimate.**

Precisely:

1. **The exact finite-`N` loop equation `[T]`.** The lattice Schwinger–Dyson identity
   `0 = ∫DU D^a_l[f e^{−S}]` gives the loop equation; in the continuum it is the
   Makeenko–Migdal equation with the `1/N²` connected correction (`§2`):
   ```
   ∂_μ (δ/δσ_{μν}(x)) W(C) = λ ∮_C dy_ν δ⁴(x−y) [ W(C_{xy})W(C_{yx}) + (1/N²) Ĝ(C_{xy},C_{yx}) ],
   λ = g²N ,   Ĝ = O(1) .
   ```
   At `N=∞` it is closed (RHS a product of the same functional); at finite `N` the `O(1/N²)` `Ĝ` is
   the connected two-loop correlator. **Linearizing** gives a loop-space operator whose spectrum
   equals the gauge-invariant field spectrum (Giles `[T]`) — so the loop object carries **no new
   spectrum** and sits on the same `0`-form/`1`-form wall (`E3`, `PNT-LARGEN-LOOP`).

2. **The decomposition `χ = χ_G + χ_conn` `[T]`.** The Dyson–Schwinger / Bethe–Salpeter structure of
   `⟨O O⟩_c`, `O = tr F²`, is
   ```
   (𝒦 ⟨OO⟩_c)(k) = C(k) + (𝒦 Π_B ⟨OO⟩_c)(k) ,      χ(0) = χ_G(0) + χ_conn(0) ,
   ```
   with `C` the contact (Wick) term = the **Gaussian/quadratic** contribution and `Π_B` the
   amputated connected 4-point. In the Legendre-dual curvature coordinate this is the exact
   variance identity
   ```
   G″(k=0) = S − D ,   S = E[Hess V] (Gaussian term) ,   D = Var(force) (connected 4-point) ,
   ```
   verified to six digits on a smooth two-variable toy (`§7 (1b)`). `D ≥ 0` by reflection positivity.

3. **RP + UV does NOT close it `[NG]`.** The loop equation plus positivity yields the inequality
   ```
   G″(0) ≤ S   (⇔  χ ≥ 1/S)        [T, from D≥0] ,
   ```
   i.e. an **upper** bound on the curvature and a **lower** bound on `χ`. The wanted statement
   `G″ ≥ λ` (⇔ `χ ≤ 1/λ`) needs `D ≤ S − λ`, an **upper** bound on the connected 4-point — which is
   the opposite of what a variance/positivity constraint gives. UV asymptotics fix only the
   large-`s` moments `M_k=∫s^k dν` and `M_{−2}` (anomaly); the gap `inf supp ν` is not a moment
   (`A4 §2`, `SLOT`). The free massless field is RP- and UV-consistent with `χ=∞`. So the class
   {exact loop equation + RP + UV asymptotics + large-`N` factorization} is **gap-blind**.

4. **Exact-test calibration `[T]/[V]`.**
   - **2D Yang–Mills** closes *exactly*: plaquettes factorize (Migdal), `χ_O = Var_{1pl}(O)`,
     `L`-independent, cross-plaquette connected `=0` → only the Gaussian/one-plaquette term survives.
   - **One-plaquette** (`SU(2)`, plaquette trace `x=cos(θ/2)`, `ρ∝√(1−x²)e^{βx}`) is exact:
     `χ=Var(x)=0.134119` at `β=2.5`, `κ₄=+1.66×10⁻² ≠ 0`, and the **operator** marginal is
     log-concave (`G″>0`). (The Savvidy non-log-concavity is the **gauge-field fibre**, a different
     variable — `STAR2`.)
   - **Free-field acid test** passes *discriminatingly*: free massless gives `χ → ∞` (exact zero
     mode; regularized `1/k̂_min² ~ L²` in 4D, `G″=1/χ→0`), massive gives `χ→1/m²`; the connected
     4-point of a Gaussian is `0`. The tool returns divergence, **no spurious floor**.

5. **The exact missing input `[O]`.** An **upper** bound on the connected 4-point `D(k=0)`, uniform
   in `(L,a)`: `D(k=0) ≤ S(k=0) − λ`. In the exact rank-one split `D = D₊ + D_Sav` (`D₊≤0` by
   Brascamp–Lieb, `D_Sav≥0` the Savvidy negative mode, `H2 §1`), this is the **softened
   sub-criticality** `D_Sav < S_bare − λ + |D₊|` (`H2 §2.4`) — i.e. `H1` in the loop-equation
   coordinate. It is *not* a sign reflection positivity can supply; it is the sign of the connected
   correlator, which crossing makes indefinite.

**One line.** *The finite-`N` loop equation is exact and closed-in-`W` up to an `O(1/N²)` connected
4-point; it writes `χ = χ_G + χ_conn` (equivalently `G″ = S − D`, `D≥0`). Reflection positivity +
UV asymptotics give only `χ ≥ 1/S` (the wrong side of `χ ≤ 1/λ`); 2D closes exactly and the free
field correctly diverges, so the machinery is honest but gap-blind in 4D. What is missing is an
**upper** bound on `D` — the Savvidy rank-one connected 4-point — i.e. `H1`. **No QED. Clay 0%.***

---

## 1 · Set-up, conventions, and the two `χ`'s

**Lattice gauge theory.** Hypercubic lattice, gauge group `SU(N)`, links `U_l∈SU(N)`, Haar measure.
Wilson action in the corpus normalization
```
   S_W[U] = β Σ_p ( 1 − (1/N) Re tr U_p ) ,      β = 2N/g² ,
```
so that `β=2.5` corresponds to `g²=2N/β = 1.6` at `N=2`. Wilson loop
`W(C) = ⟨(1/N) tr U(C)⟩`, `U(C)` the ordered product around `C`.

**Two distinct susceptibilities.** This note is careful about an identification the corpus usually
leaves implicit.
- **Field susceptibility (H1).** `χ_field = Var_μ(ξ_{k→0})`, `ξ = ⟨v, log U⟩` the `k=0` coarse mode
  of the (log-)link field (`ch02-squeeze.tex` eq. (H1)). This is the object `H1` bounds.
- **Operator susceptibility (loop-equation native).** `χ_O = ∫d⁴x ⟨O(x)O(0)⟩_c`, `O = tr F²` (or
  the plaquette trace `(1/N)Re tr U_p` as its lattice proxy). This is what the loop equation
  *directly* computes, because `O` is the natural small-loop/field-strength operator.

The bridge is two-fold and **neither leg is an identity in the non-Gaussian case**:
`G″ · Var(ξ) = 1` **only for a Gaussian marginal** (`ch02 §G″ vs 1/χ`; measured `2–5%` gap at
`β=2.5`, `[C]`), and the trace anomaly `T^μ_μ = (β(g)/2g) tr F²` ties `χ_O` to the anomaly
correlator (`[T]`, but it constrains the `M_{−2}` moment, an upper-bound machine — §4). **So the
loop equation controls `χ_O`, not `χ_field`; transporting `χ_O` to `H1` is itself `[O]`.** This
channel distinction is the first structural point of the Q2 wave.

---

## 2 · The exact finite-`N` loop equation

### 2.1 The lattice Schwinger–Dyson identity `[T]`

For any function `f` of the links and the left-invariant derivative `D^a_l` acting on link `l`
(`D^a_l f(U) = d/dt f(e^{tT^a}U_l)|_{t=0}`), Haar invariance of the compact measure gives
```
   0 = ∫ DU  D^a_l [ f(U) e^{−S_W(U)} ]  =  ⟨D^a_l f⟩ − ⟨f D^a_l S_W⟩ .              (2.1)
```
(The boundary term vanishes on the compact group.) Take `f = (1/N) tr U(C)`. Writing the loop as the
word `U(C) = A U_l B`, cyclicity gives `tr U(C) = tr(U_l · BA)`, so `D^a_l f = (1/N) tr(T^a U_l BA)`:
the derivative "inserts a generator" — the loop becomes an open path ending in `T^a`. The RHS act of
`D^a_l` on `S_W` hits each of the `2(d−1)` plaquettes through `l`, giving the loop **with one extra
plaquette spliced in** (an "ear"). Hence the exact lattice loop equation:
```
   ⟨ loop with generator ⟩  =  β Σ_{p∋l} ⟨ (plaquette-dressed loop) ⟩ .               (2.2)
```
This is an exact identity, but it **sources loops with one more plaquette** — it does not close on a
finite set of loop shapes. `[T]`

### 2.2 The continuum Makeenko–Migdal equation at finite `N` `[T]`

The continuum limit of (2.1)–(2.2) is the Migdal–Makeenko equation. With
`σ_{μν}` the loop-surface (area) derivative and `λ = g²N`, the **exact** equation (any `N`, before
factorization) is
```
   ∂_μ (δ/δσ_{μν}(x)) W(C)
     = λ ∮_C dy_ν δ⁴(x−y) [ W(C_{xy}) W(C_{yx}) + (1/N²) Ĝ(C_{xy},C_{yx}) ] ,       (2.3)
```
where at the self-intersection `x=y` the loop `C` splits into `C_{xy}, C_{yx}`, and
`Ĝ = N² Cov(W,W) = O(1)` is the **connected** two-loop correlator. The counting is standard:
`W(C)=⟨(1/N)tr U(C)⟩ = O(1)`, and `⟨(1/N)tr U(C_1)(1/N)tr U(C_2)⟩ = W(C_1)W(C_2) + O(1/N²)`. Thus:
- **`N=∞`:** (2.3) is **closed** in the single functional `W(C)` (RHS `= λ W·W`), *no* `n`-point
  tower. This is the genuine structural advantage over the DSE tower (`PNT-LARGEN-LOOP §1`). `[T]`
- **finite `N`:** the `O(1/N²)` connected term is present and is exactly the connected 4-point/data
  that no large-`N` factorization removes. `[T]`

### 2.3 Linearization and the channel wall `[T]/[NG]`

Linearizing (2.3) around the (unknown) vacuum `W₀` gives a **loop-space fluctuation operator** `ℳ`;
its inverse is the connected two-loop correlator, and its lightest mode is the `0⁺⁺` glueball.
But `{W(C)}` separates gauge orbits (**Giles**, PRD 24 (1981) 2160), so the generated `σ`-algebra is
the full gauge-invariant one and
```
   spectrum(ℳ) = spectrum( L_field |_{gauge-invariant} ) .                              (2.4)
```
Hence the loop equation carries **no spectrum the field theory does not already carry**: the "loop
gap" is either the physical gap itself (if `ℳ` is the physical transfer) or a fictitious-time rate
(if `ℳ` is a Langevin generator) — `E3-LOOPSPACE` maps this exhaustively. Additionally, the loop
equation's *large-loop* content is the string tension `σ` (`1`-form / confinement), whereas `χ` and
`m_{0⁺⁺}` are `0`-form objects. **`[T]` the identity (2.4); `[NG]` for "the loop equation yields an
independent `0⁺⁺` gap".**

---

## 3 · `χ` as a Gaussian term plus a connected 4-point term

### 3.1 The Dyson–Schwinger / Bethe–Salpeter decomposition `[T]`

For the composite `O=tr F²`, the exact SD equation is a Bethe–Salpeter equation. In momentum space,
with `𝒦` the free (composite) inverse propagator, `C(k)` the contact/Wick term, and `Π_B` the
amputated connected 4-point irreducible kernel,
```
   (𝒦 ⟨OO⟩_c)(k) = C(k) + (𝒦 Π_B ⟨OO⟩_c)(k) .                                          (3.1)
```
Solving, `⟨OO⟩_c(k) = (𝒦 − 𝒦Π_B)^{-1} C(k)`, and at `k=0`
```
   χ ≡ ⟨OO⟩_c(0) = χ_G(0) + χ_conn(0) ,   χ_G = C(0)/𝒦(0) (Gaussian/quadratic term),
                                           χ_conn = (connected 4-point contribution).     (3.2)
```
`[T]` — this is just the SD equation rewritten. `χ_G` is the self-contraction (tadpole/one-loop)
piece, controlled by the free propagators and hence by **UV asymptotics**; `χ_conn` is the
interaction piece, controlled by the connected 4-point.

### 3.2 The Legendre-dual curvature form and the exact variance identity `[T]`

The susceptibility is the inverse of the second variation of the effective action. In the curvature
coordinate the object is `G″(v;a₀) = E[∂_v²V|a₀] − Var(∂_vV|a₀)` (`ch02`; `INVENT-VARIANCE-INEQ`;
`H2 (1.1)`). For a single coarse coordinate this is an **exact** identity, proved in three lines:

> **Lemma Q2.1 (variance decomposition).** Let `ρ(a)=∫db e^{−V(a,b)}` be a marginal of the
> Boltzmann measure `e^{−V(a,b)}`. Then
> ```
>    G″(a) := −(log ρ)″(a) = E[∂_a²V | a] − Var(∂_aV | a) ,
> ```
> with `E, Var` over the conditional measure `e^{−V(a,·)}/ρ(a)` restricted to `b`.
> *Proof.* `ρ'=∫(−∂_aV)e^{−V}` and `ρ''=∫[(∂_aV)²−∂_a²V]e^{−V}`, so
> `−ρ''/ρ= E[∂_a²V]−E[(∂_aV)²]`; adding `(ρ'/ρ)²=(E∂_aV)²` gives `G″=E[∂_a²V]−Var(∂_aV)`. ∎

The corpus's **multi-direction** version adds the cross-coupling term `⟨K⁻¹c,c⟩` (`H2 (1.2)`,
`H1-RG-DESCENT (1.2)`); the scalar-coordinate identity above is the clean backbone. Witness `§7 (1b)`
confirms `G″ = E[∂²V] − Var(∂V)` to six digits. Write
```
   S := E[∂_a²V] (= E[Hess V], the Gaussian/quadratic term) ,   D := Var(∂_aV) (connected 4-pt) ,
   G″(a₀) = S − D .                                                                       (3.3)
```
`D` is the connected 4-point function of the force `f=∂_aV` — the variance that appears as the
subtraction. **This is the same "connected 4-point term" as `χ_conn` in (3.2), read in the curvature
coordinate.** `[T]`

### 3.3 Relating the two coordinates, and the sign

`χ` and `G″` are reciprocal only in the Gaussian-marginal limit: `G″·Var(ξ)=1` (`ch02 §"G″ versus
1/χ"`, `[C]` at `β=2.5`). Hence in the Gaussian/CLT approximation
```
   χ ≈ 1/G″ = 1/(S − D)  (Gaussian) ,      exact:  χ = Var(ξ) ,   G″ = S − D .              (3.4)
```
The decisive sign is in (3.3): **`D = Var(force) ≥ 0`, so the connected 4-point *subtracts* from the
Gaussian stiffness `S`.** Reflection positivity / Källén–Lehmann guarantee `D ≥ 0` (a variance) and
`χ ≥ 0`. Both are *lower* bounds in the direction that makes `G″` **smaller**, `χ` **larger**.

---

## 4 · Does reflection positivity + UV asymptotics close a bootstrap on `χ`? **No `[NG]`**

### 4.1 The exact inequality one gets, and the one that is wanted

From (3.3) and `D≥0`:
```
   G″(0) = S − D ≤ S      [T]     (the connected 4-point can only LOWER G″) .
```
(In the Gaussian regime with `S>0` this reads `χ ≥ 1/S` — a lower bound on `χ`, the opposite of what
`H1` needs; the *direction* conclusion below is independent of the sign of `S`.)
The wanted `H1` statement is `G″(0) ≥ λ > 0` (⇔ `χ ≤ 1/λ`). It requires
```
   D(0) ≤ S(0) − λ      (an UPPER bound on the connected 4-point)   [O] .
```
Reflection positivity gives `D ≥ 0` — a **lower** bound on the very quantity that must be bounded
**above**. This is the recurrent "upper-bound machine, lower-bound target" mismatch of the campaign
(`ch02 §convex-order`; `A4`; `SLOT`), now read in the loop-equation coordinate: **positivity controls
the wrong end of `D`.** `[T]` for the inequality `G″≤S`; **`[NG]` for closing `H1` with it.**

### 4.2 Dispersion (Kramers–Kronig): the UV fixes the tail, the IR fixes the value `[T]`

Write the connected self-energy `Π(k)` of the `O`-channel; analyticity + RP give a spectral
representation `Im Π(s) ≥ 0`, and the value at `k=0` is the dispersion integral
```
   Re Π(0) = (2/π) ∫_0^∞ Im Π(s) / s  ds .                                                  (4.1)
```
Reflection positivity supplies `Im Π(s) ≥ 0` (a **lower** bound on the integrand); asymptotic freedom
supplies `Im Π(s)` at **large `s`** (the perturbative/OPE tail, `[T]`). But the integral (4.1) is
**IR-dominated** (weight `1/s`, so small `s` dominates): `Re Π(0)` is fixed by `Im Π` at small `s`,
where the theory is nonperturbative and `Im Π` is unknown. So RP + UV leave `Re Π(0)` — hence `χ` —
unbounded from above. `[T]`

### 4.3 UV asymptotics are moments, and moments are gap-blind `[T]`

"Asymptotic freedom / OPE" makes precise the following data: the KL moments `M_k=∫s^k dν(s)` for
`k≥0` (from the UV expansion of `⟨OO⟩_c`), and `M_{−1}=χ`, `M_{−2}∝Λ⁴` (the anomaly, `SLOT §1`). All
`M_k` are **Taylor-germ data at `s=∞`**; the gap is `inf supp ν`, a **pole at finite `s`**. Every
moment/Markov/Christoffel constraint bounds `inf supp` from **above** (`M_0/M_{−1}`, `M_{−1}/M_{−2}`,
…) and admits the free-field gapless measure (`ν` absolutely continuous, `f_1²=0`) at residual
`<1.5×10^{-11}` for every truncation order (`A4 §2–3`, `SLOT §2`). **`[T]` Finitely many exact
moments — which is exactly what UV asymptotics furnishes — never determine `inf supp`.** So UV
asymptotics is *gap-blind by construction*, not for lack of data. **The moment bootstrap for `χ`
itself** (set `χ=M_{−1}`, supply `M_k` from UV, ask feasibility) therefore has an
**unbounded-above** feasible set — the free field (`χ=∞`) is feasible at every truncation
(`SLOT §2`: `Δ_min(K)=0` for all `K`) — so it yields **no** upper bound `χ ≤ 1/λ`; this is the
`A4`/`SLOT` verdict read at `χ=M_{−1}`.

### 4.4 The large-`N` factorization error is the wrong smallness `[T]/[O]`

At `N=∞` the connected term in (2.3)/(3.2) drops and the loop equation closes. The error is exactly
`O(1/N²)`:
```
   N=2 : 1/N² = 0.25 ,      N=3 : 1/N² = 0.11      (not small) .                            (4.2)
```
At the physical `N=2,3` the "factorization error bound" is `11–25%`, the same order as the effect.
More seriously, a **uniform** bound on the connected term is *equivalent to* `⟨OO⟩_c(0) < ∞`, i.e.
`χ < ∞` = `H1` itself. So large-`N` factorization does not furnish an independent bound; using it
uniformly is **circular** (`[O]`). `[T]` that the error is `O(1/N²)`; `[O]` that it cannot be
controlled uniformly without `H1`.

### 4.5 Verdict on the class

> **`[NG]` (for the named class).** No argument combining **(i)** the exact finite-`N` loop equation,
> **(ii)** reflection positivity / Källén–Lehmann, **(iii)** the known UV asymptotics (moments,
> anomaly), and **(iv)** large-`N` factorization, can prove `χ ≤ 1/λ` uniformly at `β=2.5`. The
> obstruction is the **sign/direction of the connected 4-point**: RP gives `D≥0` (subtracts from
> `G″`), UV gives the large-`s` moments (gap-blind), and factorization is `O(1/N²)` and circular.

The *only* way this class could carry a gap would be if one of (i)–(iv) secretly supplied an **upper**
bound on `D`. None does: (i) is an identity (free field satisfies it), (ii) is one-sided the wrong
way, (iii) is gap-blind, (iv) is circular. This is the loop-equation restatement of
`A4`/`SLOT`/`FORGE5`'s "positivity is Schur-concave / cannot force the edge."

---

## 5 · Exact tests: 2D Yang–Mills, the one-plaquette model, the free field

### 5.1 2D Yang–Mills: the loop bootstrap closes **exactly** `[T]`

In `D=2` pure YM the plaquette variables are the independent degrees of freedom: on the torus the
map `links → plaquettes` has unit Jacobian, and the Wilson measure factorizes into a product of
single-plaquette measures (Migdal; the exact solvability underlying `W(C)=e^{−λ C₂ \mathrm{Area}/2}`,
`PNT-LARGEN-LOOP §2`). Hence any per-plaquette operator `O_p = (1/N)Re tr U_p = x_p` is **i.i.d.**:
```
   ⟨O_p O_{p'}⟩_c = δ_{p p'} Var_{1pl}(O) ,   χ_{2D}(O) = Σ_{p'} ⟨O_p O_{p'}⟩_c = Var_{1pl}(O) ,   (5.1)
```
**`L`-independent**, and the **connected 4-point across distinct plaquettes vanishes** (it is a pure
contact at `p=p'`). Therefore in `D=2` the loop/moment bootstrap **closes exactly**: only the
Gaussian/one-plaquette term survives, `χ_conn` (propagating) `=0`. Witness `§7 (2)`: `χ_{2D}(β=2.5)
= 0.134119` for any `L`, cross-plaquette connected `=0`. **`[T]`**

*Honest caveat.* 2D closes **trivially** — it does *not* test the hard part (the propagating
connected 4-point in `D=4`) because that part is absent. Its value is calibration: it shows the
decomposition (3.2)–(3.3) is exact when the measure is known, and isolates that the 4D difficulty is
*entirely* the non-contact connected 4-point.

### 5.2 The one-plaquette model: exact `χ`, non-Gaussianity, and the channel distinction `[T]/[V]`

`SU(2)` one plaquette, `U=e^{iθ n·σ/2}`. In the exponential chart the class-angle marginal is
`∝ sin²(θ/2)dθ`; in the **plaquette-trace** variable `x=(1/N)Re tr U = cos(θ/2)∈[−1,1]` this becomes
the exact marginal
```
   ρ(x) ∝ √(1−x²) e^{βx} ,   Z = π I_1(β)/β ,   ⟨x⟩ = I_0(β)/I_1(β) − 2/β .                  (5.2)
```
The exact susceptibility and connected 4-point (`§7 (1)`):

| `β` | `⟨x⟩` | `χ=Var(x)` | `1/Var` | mode `x₀` | `G″(x₀)=(1+x²)/(1−x²)²` | `κ₄` (conn 4-pt) |
|---|---|---|---|---|---|---|
| 0 | 0 | **0.250000** | 4.000 | 0 | 1.000 | −6.25×10⁻² |
| 2.5 | +0.507 | **0.134119** | 7.456 | +0.820 | 15.55 | **+1.66×10⁻²** |
| 4.30661 | +0.679 | 0.065359 | 15.30 | +0.891 | 41.93 | +1.18×10⁻² |
| 6 | +0.763 | 0.036886 | 27.11 | +0.920 | 78.52 | +4.82×10⁻³ |

Two exact facts. **(a)** In the plaquette-trace (`tr F²`) variable, `G″(x₀)=(1+x²)/(1−x²)²>0`
**always** ⇒ this marginal is **log-concave**, `χ` finite. **(b)** `κ₄≠0` ⇒ the marginal is
non-Gaussian, so the connected 4-point term is genuinely nonzero — the bootstrap does **not** reduce
to a Gaussian even in the solvable model; it is exactly computable only because there is a single
variable. **`[T]`**

*Channel distinction (important).* The Savvidy non-log-concavity of the campaign is a fact about the
**gauge-field fibre** (`STAR2 §1`: `G″=0.20511530−0.01105929β²<0` for `β>4.307`), a **different**
variable from the plaquette-trace `x`. The loop equation's operator `χ_O` is the tame (log-concave)
one above; `H1`'s field `χ_field` is the Savvidy one. **They are not the same object** (§1) — the
one-plaquette model makes the difference explicit and shows the loop-equation machinery is
*optimistically* well-behaved on `χ_O` while `H1` lives on the harder `χ_field`.

### 5.3 Free-field acid test: divergent `χ` `[T]/[V]`

A free massless field must give `χ=∞` (equivalently `G″=1/χ=0`); a tool that returns a positive floor
on the free field is vacuous (charter §1). The zero-momentum susceptibility is the `k→0` propagator.
The exact zero mode gives `χ=∞` identically (no restoring force); regularized on a periodic `L⁴`
lattice, `χ(L)=1/k̂_min²` with `k̂_min²=16 sin²(π/L) → 0`, so `χ ≍ L²`:
```
   L:      4        8       12       16       24       32
   χ:    0.1250   0.4268  0.9330   1.6421   3.6685   6.5054     (∝ L² → ∞)
   G″:   8.000    2.343   1.072    0.609    0.273    0.154      (1/χ → 0)
```
Free massive (`m²=0.5`): `χ → 1/m² = 2.000` (§7: `1.53→1.86` at `L=32,64`, converging). The
connected 4-point of a Gaussian field is `0` (Wick), so the divergence is **purely the Gaussian
(zero) mode** — there is no connected term to blame. **The tool returns divergence, not a spurious
floor. `[T]+[V]`** This is the acid test, passed *discriminatingly*: the same machinery that gives
`χ=∞` on the free field gives `χ<∞` on the gapped one-plaquette/2D models.

---

## 6 · The exact missing input `[O]`

Assembling §3–§5:

> **`[O]` (exact missing term).** An **upper** bound on the connected 4-point function at `k=0`,
> uniform in `(L,a)`:
> ```
>    D(k=0) ≤ S(k=0) − λ        (⇔  G″(k=0) ≥ λ  ⇔  χ ≤ 1/λ) ,
> ```
> where `D = Var(force)` (equivalently `χ_conn`, the amputated connected 4-point). In the exact
> rank-one split `D=D₊+D_Sav` (`D₊≤0` Brascamp–Lieb; `D_Sav = E|⟨φ₋|c⟩|²/(g|B|) ≥ 0` the Savvidy
> negative mode, `H2 §1`), this is the **softened sub-criticality**
> ```
>    D_Sav  <  S_bare − λ + |D₊|      at  β≈2.5 ,     [H2 §2.4] ,
> ```
> which is `H1` in the loop-equation coordinate.

The missing term is, concretely, **a sign-definite upper bound on the Savvidy rank-one contribution to
the connected 4-point** — the `D_Sav` that `H2` shows is `≥0`, saturates `~S` (`D/S≈0.97`), and is
**not** bounded above by Brascamp–Lieb (which gives only `D₊≤0`). Reflection positivity cannot supply
it: RP is the lower edge `D≥0`; the needed statement is the upper edge. Crossing symmetry does not
help either — it makes the *t*-channel connected 4-point **sign-indefinite** by construction (the
`s`-channel positivity of `GLUEBALL-BOOTSTRAP`/`FORGE5` is the wrong channel; `D` lives in the
*t*-channel force-force contraction). Three equivalent disguises of the same input:
1. the **upper bound on `D`** (`= D_Sav < S_bare − λ`), loop-equation coordinate;
2. the **large-`N` factorization error bound with the right sign at `N=3`** (an *upper* bound on the
   `O(1/N²)` connected term, uniform — §4.4; circular unless independent);
3. the **IR spectral-weight bound** `ν([0,ε]) ≤ Cε^p` (a non-moment, edge-sensitive datum,
   `A4 §6`/`SLOT §6`).

All three are the single open slot `H1 = I2 = H*` (`ch02 §correct sufficient condition`). It is **not**
a statement UV asymptotics or reflection positivity can close — those are moments-and-positivity (upper
bounds / free-field-blind / the wrong sign of `D`).

---

## 7 · Witness — `research/v3/w7/q2_loop_bootstrap.py`

`py -3.12 -I q2_loop_bootstrap.py` (numpy 2.4.6, scipy 1.18.1; exit 0). Quoted output (abridged):

```
(1) ONE-PLAQUETTE SU(2): x=cos(theta/2), rho~sqrt(1-x^2)e^{beta x}
    beta    <x>      chi=Var(x)  1/Var    mode x0   G''(x0)   kappa4(conn4)
    0.000 +0.00000  0.250000    4.000000  +0.00000   1.0000   -6.25e-02
    2.500 +0.50720  0.134119    7.456065  +0.81980  15.5495   +1.66e-02
    4.307 +0.67930  0.065359   15.300137  +0.89062  41.9293   +1.18e-02
    6.000 +0.76273  0.036886   27.110634  +0.92013  78.5208   +4.82e-03
(1b) EXACT variance identity  G''=E[d^2V]-Var(dV):  a=0.3,0.8,1.4 -> 1.627479,1.544411,1.394681
     match to 6 digits;  D=Var(dV)=0.0098,0.0641,0.1548 (>=0).  Fisher: E[G'']=I=1.5625.
(2) 2D YM (Migdal): chi_2D = one-plaquette Var = 0.134119 (any L); cross-plaquette connected = 0.
(3) FREE massless: chi=1/khat_min^2 ~ L^2 -> inf;  G''=1/chi -> 0.
      L=4,8,16,32: chi=0.1250,0.4268,1.6421,6.5054 ;  G''=8.00,2.34,0.61,0.15.
    free massive m^2=0.5: chi -> 1/m^2=2.000 (L=64: 1.857).  Gaussian connected 4-pt = 0.
(4) DIRECTION: G''=S-D, D=Var(force)>=0 subtracts.  RP gives D>=0 (lower); need D<=S-lam (upper).
```
Readings: **(1)** exact one-plaquette `χ`, non-Gaussian (`κ₄≠0`), log-concave in the operator
variable; **(1b)** the variance identity `G″=S−D`, `D≥0`, exact to six digits, with the Fisher
identity on a clean marginal; **(2)** 2D factorizes ⇒ bootstrap closes exactly; **(3)** free-field
divergence `∝L²` (acid pass), massive `→1/m²`; **(4)** the sign/direction obstruction. **Not a proof
— a witness for the exact decomposition, the two solvable closures, and the direction of the missing
term.**

---

## 8 · Status ledger

| item | statement | label |
|---|---|---|
| lattice SD identity `0=∫D^a_l[fe^{−S}]`; loop equation sources an extra plaquette | (2.1)/(2.2) | **[T]** |
| finite-`N` MM `∂_μ(δ/δσ)W = λ∮δ[W·W+(1/N²)Ĝ]`; closed at `N=∞` | (2.3) | **[T]** |
| `spectrum(loop) = spectrum(field\|gauge-inv)` (Giles); loop `0⁺⁺` = field `0⁺⁺` or fictitious time | (2.4), E3 | **[T]** / **[NG]** |
| `χ = χ_G + χ_conn` (SD/BSE); `G″=S−D`, `S=E[Hess]` (Gaussian), `D=Var(force)` (conn 4-pt) | (3.1)–(3.3), Lem Q2.1 | **[T]** |
| `χ≈1/G″` Gaussian only; exact `χ=Var`, `G″·Var=1` only Gaussian (2–5% at β=2.5) | (3.4), ch02 | **[T]** / **[C]** |
| RP/positivity ⇒ `D≥0`, `χ≥0`; hence `G″≤S`, `χ≥1/S` — the **wrong** side of `χ≤1/λ` | §4.1 | **[T]** |
| needed inequality `D≤S−λ` (upper bound on the connected 4-pt) | §4.1 | **[O]** |
| Kramers–Kronig: `ReΠ(0)=(2/π)∫ImΠ(s)/s ds` is IR-dominated; RP gives `ImΠ≥0`, UV the large-`s` tail | (4.1) | **[T]** |
| UV asymptotics = moments `M_k`, `M_{−2}`; moments are `s=∞` Taylor germ, gap-blind (`δ₀+δ₂` vs `2δ₁`) | §4.3, A4/SLOT | **[T]** |
| large-`N` factorization error `O(1/N²)` = 0.25 (N=2), 0.11 (N=3); uniform control is circular | §4.4 | **[T]** / **[O]** |
| **loop eq + RP + UV + factorization cannot prove `χ≤1/λ` at β=2.5** | §4.5 | **[NG]** |
| 2D YM: plaquettes factorize ⇒ `χ=Var_{1pl}`, L-indep, conn cross-plaquette `=0`, bootstrap closes exactly | (5.1) | **[T]** |
| one-plaquette `SU(2)`: exact `ρ∝√(1−x²)e^{βx}`, `χ=Var(x)=0.134` (β=2.5), `κ₄≠0`, log-concave in operator var | (5.2), §7 | **[T]**/**[V]** |
| operator `χ_O` (log-concave) ≠ field `χ_field` (Savvidy, `G″=0.2051−0.0111β²`) | §1, §5.2, STAR2 | **[T]** (distinct objects) |
| free-field acid: `χ→∞` (zero mode; `∝L²`), `G″→0`; massive `→1/m²`; Gaussian conn 4-pt `=0` | (5.3), §7 | **[T]**/**[V]** |
| missing term `= D_Sav` upper bound `= 2g²/λ₄ < S_bare−λ` (softened H1) `= H1 = I2 = H*` | §6, H2 §2.4 | **[O]** |
| **4D `χ≤1/λ` at β=2.5** | — | **[O] = H1** |
| **Clay** | unbuilt | **0%** |

---

## 9 · Sharpest line

**The finite-`N` loop equation is exact and, at `N=∞`, closed in the one functional `W(C)` — but the
`1/N²` connected 4-point it leaves is precisely the term that decides `χ`. In the curvature
coordinate it reads `G″=S−D`, `S=E[Hess]` (Gaussian) and `D=Var(force)` (connected 4-pt), and
reflection positivity fixes the sign of the *subtraction*: `D≥0`, so positivity, Källén–Lehmann, the
UV moments (`M_k`, the anomaly `M_{−2}`) and large-`N` factorization all give the **lower** bound
`χ≥1/S` — never the upper `χ≤1/λ`. 2D closes exactly (plaquettes factorize, `χ=Var_{1pl}`,
`L`-independent), the one-plaquette model is exact (`χ=0.134` at `β=2.5`, non-Gaussian but
log-concave in the operator variable), and the free massless field correctly diverges (`∝L²`), so the
machinery is honest — and, being honest, gap-blind in 4D. The exact missing input is an **upper** bound
on the connected 4-point `D`, equivalently on the Savvidy rank-one part `D_Sav`, uniformly in `(L,a)`:
the softened sub-criticality `D_Sav < S_bare − λ` — `H1` in the loop-equation coordinate, and the sign
of a connected correlator that no positivity can supply. **No QED. Clay 0%.**

---

*Refs (internal, this session).* `P-WAVE-SYNTHESIS.md`, `D6-DECISION.md`, `H1-RG-DESCENT.md`,
`H2-DEFECT-BOUND.md` (§1 rank-one split, §2.4 softened sub-criticality), `C2-BALABAN-UNIFORM.md`;
`docs/paper-v2/ch02-squeeze.tex` (H1, thin atom, squeeze, `G″` vs `1/χ`), `docs/paper-v3/main.tex`
§6 (defect, Fisher, descent, soft mode); `research/unconditional/`: `PNT-LARGEN-LOOP.md` (MM closed in
`W`, 2D exact, `0`-form/`1`-form wall), `E3-LOOPSPACE.md` (Giles, loop spectrum = field spectrum),
`CONSTRUCT-LARGEN-BOOTSTRAP.md` (D=4 truncation deficit, two functionals), `A4-SPECTRAL-BOOTSTRAP.md`
(moments = `∞` Taylor germ, gap-blind), `SLOT-BOOTSTRAP.md` (feasible floor `=0` at every `K`),
`SLOT-DSE-FIXPOINT.md` (DSE hierarchy, seagull certifies `m²=0`, vertex-pole residue `=H1`),
`GLUEBALL-BOOTSTRAP.md`/`FORGE5-BOOTSTRAP-EXCLUDE.md` (crossing/positivity Schur-concave, wrong
direction), `BROWNIAN-LOOP-SOUP.md` (no positive 4D loop measure); `STAR2-CLASSIFICATION.md` (Savvidy
one-plaquette `G″=0.2051−0.0111β²`).
*External, standard, from the corpus (not re-fetched):* Makeenko–Migdal, Phys. Lett. B88 (1979) 135;
Migdal, Sov. Phys. JETP 42 (1976) 413; Giles, PRD 24 (1981) 2160; Brascamp–Lieb (1976), Prékopa
(1971); Källén (1952)/Lehmann (1954); Kramers–Kronig (standard dispersion). Witness
`research/v3/w7/q2_loop_bootstrap.py` (`py -3.12 -I`, numpy 2.4.6 + scipy 1.18.1). **No third-party
sources downloaded. Clay 0%. No QED.**
