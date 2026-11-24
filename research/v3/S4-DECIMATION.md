# S4-DECIMATION — the decimation step: what `(★)` and `(★dec)` each assert, and their exact logical relation

> 2026-11-24 · Attack on the UV RG step (L3 node), angle: **the decimation step**.
> Reads: `v3/C2-BALABAN-UNIFORM.md` §2.4, `docs/paper-v2/ch03-existence.tex` §G1.2 (L3),
> `unconditional/PROOF-COARSE-POINCARE.md`, `G1-UNIFORM-ZA.md`, `G1-L3-TOWER-POINCARE.md`,
> `L3-band-decay.md`, `PROOF-L3.md`, `PROVE-Gpp-positive.md`, `PROVE-H3.md`, `PROVE-H2-uvgap.md`,
> `PROOF-H2-GLOBAL.md`, `BLOCKING-DESIGN.md`.
> Witness: `research/v3/w7/s4_decimation.py` (deterministic; output quoted in §8).
>
> **Labels.** **[T]** true theorem/identity · **[C]** conditional · **[NG]** no-go for a named class ·
> **[O]** exact missing step / open · **[V]** verified here (numerics). No QED.

---

## 0 · Verdict (no QED)

**`(★)` and `(★dec)` are logically *independent*. Neither implies the other `[T]`.**

They are the two halves of **one matrix** — the curvature `G″ = E[H] − Cov(s,s) + W` — read on the diagonal
and off the diagonal respectively:

| | statement | what it controls |
|---|---|---|
| `(★)` | `Var_μ(⟨s_B,v⟩) ≤ E_μ⟨H_B v,v⟩ + W_B(v) − δ_W` | **diagonal** of `G″`: a curvature *comparison* (variance **vs** mean Hessian) |
| `(★dec)` | `\|Cov_μ(s_B,s_{B′})\| ≤ C e^{−κ dist(B,B′)/a_j}` | **off-diagonal** of `G″`: a *locality/decay* statement |

Consequences, sharp:
1. `(★)` ⇏ `(★dec)` **[T]**: a diagonal bound gives only `|Cov| ≤ sqrt(Var_B Var_B′) ≤ const` — **no decay**.
   To pass from a variance bound to exponential covariance decay one needs a **mixing/Poincaré hypothesis**;
   the only one that works is the *global* fiber LSI (H2-global, Zegarlinski `C(β)<1`), which **fails at the
   head** `β≈2.5` (`C≈15`). So the missing hypothesis is named, and it is the same Savvidy wall (§3.5).
2. `(★dec)` ⇏ `(★)` **[T]**: at `B=B′` it gives only `Var(s_B) ≤ C` (boundedness). `(★)` additionally needs
   the **Hessian comparison** `E[H_B]+W_B ≥ C+δ_W`, and `E[H_B]` can be *negative* (Savvidy floor `−1.1β`).
   Boundedness is not curvature. The two meet only under that explicit quantitative condition.
3. `(★) ∧ (★dec) ⇒ G″ ⪰ ½δ_W > 0` (operator form) **[T]** by **Gershgorin / diagonal dominance**, with the
   explicit interface `κ > (2c_d C/δ_W)^{1/d}`, `c_4=12π²≈118` (§3.4, §4).
4. **Decimation *preserves and strengthens* the decay rate `κ`** (`κ′ ≥ Lκ` per coarse block) but **not the
   prefactor `C`** (`C′ ≤ (block multiplicity)·C`) **[T]** (§4.3). So the "multiplicative structure" helps
   `κ`, hurts `C`; hence `(★dec)` must be stated with `C` **per unit volume** (the same Λ-tempering that
   already freezes Lemma B).
5. **Free-field acid test** (§5): the *inequality* `(★)` **fails** on the free massless field (`G″(0)=0` at
   the `k→0` mode, so `G″ ≥ δ_W` is false), while `(★dec)` **holds trivially** (the Gaussian score covariance
   is the coarse precision, a *finite-range* kernel). They are different objects, not one inequality renamed.

**Net.** The UV step stays `[C]` on `(★) ∧ (★dec)` — but the two are now *separated*: `(★)` is a curvature
comparison at one block, `(★dec)` is a mixing statement across blocks, and **only their conjunction** (via
Gershgorin) feeds the decimated Hessian. The bridge between them is H2-global, i.e. exactly the head wall.

---

## 1 · The decimation map, exactly

### 1.1 Fine measure and block map `P`

Fine lattice: spacing `a` (shell index `j`, so `a≡a_j`), periodic torus `Λ`, link set `E`, gauge group
`G=SU(N)`, configuration `U∈G^E`. Wilson measure
```
dμ(U) = Z^{-1} e^{−S_W(U)} Π_e dHaar(U_e),   S_W(U)=β Σ_p (1 − (1/N) Re tr U_p).
```
**Block map (campaign's `P`, `PROVE-H3` §"P, precisely").** Partition sites into `2^d`-blocks (factor
`L`); for each coarse link `e′` crossing a block face,
```
U′_{e′} = P(U)_{e′} = proj_G exp( Σ_f w_f Ad(T_f) log U_f ),
```
`w` a Gaussian kernel at the block scale, `T_f` parallel transports to the block edge. Gauge covariance:
`U→gU ⇒ U′→g′U′`, `g′=g(center)`. Linearised at `U=1`, `dP` is convolution with an **even low-pass symbol**
`p(k)`: `p(0)≠0` (DC passed), `p(π/a)→0` (UV killed), surjective onto the coarse DOF. On tangents at least,
`dP` is a genuine (co)isometry: `dP dP† = I` on coarse, `Π := I − dP†(dPdP†)^{-1}dP` the orthogonal
projector onto the fiber `ker dP` (`PROVE-H3(a)`, holds).

### 1.2 Coarse variables, fiber coordinates, and the induced measure

Write `U ↔ (U′, ξ)`, `U′ = P(U)`, `ξ` the fiber remainder; this is a local diffeomorphism, with Jacobian
`J(U′,ξ)` (`log J` = Haar-volume/Weyl factor). The **indigenous coarse (decimated) measure** is the
pushforward
```
ν := P_* μ,    S_eff(U′) := −log(dν/dU′).
```
The **fiber conditional** at fixed `U′` is
```
dμ(·|U′) ∝ e^{−S_f(U′,ξ)} dξ,   S_f(U′,ξ) := S_W(U) − log J(U′,ξ).
```
**This conditional — not `ν`, not `μ` — is the measure in `(★)` and `(★dec)`.**

### 1.3 The score `s_B`: which observable

Define the **coarse score in block `B`**
```
s_B(U′,ξ) := ∂_{U′_B}( −S_W + log J ) = −∂_{U′_B} S_f(U′,ξ).
```
This is the campaign's `s_B` (`L3-band-decay`: `s_B=∂_{U′_B}(−S_W+log J)`). Three facts, all **[T]**:
- `s_B` is a **local** fiber observable: `∂_{U′_B}` at fixed `ξ` touches only the `O(1)` fine plaquettes of the
  coarse link's block, so `supp s_B ⊆ W(B)`, a width-`O(a_j)` neighbourhood (`L3-band-decay` fact 1,
  `PROVE-Gpp §1`).
- `∂_{U′_B} S_eff(U′) = E_{μ(·|U′)}[s_B]` (differentiate `Z(U′)=∫e^{−S_f}`).
- `Σ_B s_B ≡ 0` (uniform coarse shift = block gauge orbit), so `Σ_{B′}Cov(s_B,s_{B′}) ≡ 0`; only the
  **absolute** row-sum carries content (`L3-band-decay` fact 4).

`(★)`'s LHS `Var_μ(⟨Pᵀ∇S_j,v⟩)` is `Var_μ(⟨s_B,v⟩)` summed over the coarse direction `v` supported in `B`;
`E_μ⟨PᵀH_jP v,v⟩` is the coarse-restricted mean Hessian. So `(★)` is a statement about `s_B` under
`μ(·|U′)` — the same `s_B` that appears in `(★dec)`. That is the whole reason the two are confused.

### 1.4 The exact curvature identity (law of total variance for `−log Z`)

Let `H := ∂²_{U′U′} S_W` (the campaign's force Hessian, `H=Hess S_W`: positive at small field, Savvidy-indefinite
on strong backgrounds) and `W := −∂²_{U′U′}(log J)` (so `W≥0` on the geodesic ball, `PROVE-Gpp-positive §2`).
A two-line computation (`∂_Z`, then `∂²_Z`) gives the **exact** identity **[T]** (it is the Hessian form of
`Var=ΣVar_conditional + Var_mean`):
```
G″_{BB′} := ∂²_{U′_B U′_{B′}} S_eff
          =  E_{μ(·|U′)}[ H_{BB′} ]  −  Cov_{μ(·|U′)}( s_B , s_{B′} )  +  W_{BB′}.
```
Equivalently, in matrix form, with `D_B := E[H_{BB}]+W_{BB}−Var(s_B)` and `[Cov]_{BB′}:=Cov(s_B,s_{B′})`:
```
G″  =  D  −  (Cov − diag Cov).      (diagonal already has −Var; off-diagonal is −Cov)
```
So:
- **`(★)` is exactly `D_B ≥ δ_W` for every block `B`** — a statement about the **diagonal** of `G″`.
- **`(★dec)` is exactly a bound on the off-diagonal `Cov(s_B,s_{B′})`** — a statement about the **off-diagonal**
  of `G″`.

This is the structural core: the two inequalities are the diagonal and the off-diagonal of one and the same
matrix. (Sign convention: the campaign writes `G″=E[H]−Var+W`; the derivation above fixes the sign of the
off-diagonal term as `−Cov`, matching `G1-L3-TOWER-POINCARE` §5's Gershgorin row-sum over
`|Cov(s_B,s_{B′})|`.)

### 1.5 What each statement is about (summary)

| | `(★)` | `(★dec)` |
|---|---|---|
| observables | `s_B = ∂_{U′_B}(−S_W+\log J)` | `s_B , s_{B′}` |
| measure | `μ(·\|U′)` (fiber conditional, one shell `j`) | `μ(·\|U′)` (same) |
| object | `Var(s_B)` **vs** `E[H_{BB}]+W_{BB}` (a comparison) | `Cov(s_B,s_{B′})` (a correlation) |
| distance | none (single block) | `dist(B,B′)` in shell-`j` lattice units, i.e. `dist/a_j` |
| matrix position | **diagonal** of `G″` | **off-diagonal** of `G″` |
| uniform in | `j, U′, a` | `j, U′, a`; `C,κ` **per unit volume** |

**"Cross-shell" caveat.** As written by the campaign, `(★dec)` is *same-shell, cross-block*: one decimation
step `P_j`, one conditional `μ(·|U′)`, two blocks at that scale. The genuinely **cross-shell** requirement is
*uniformity of `(C,κ)` in `j`* and the tower composition (§4.4) — **not** a different inequality. We adopt the
campaign's wording but flag the distinction.

---

## 2 · Sharp logical content

From §1.4 both statements are consequences/properties of the single matrix `G″`:

- `(★)` ⟺ `∀B,v: Var_μ(⟨s_B,v⟩) ≤ E_μ(H_B v,v) + W_B(v) − δ_W`.   (diagonal curvature comparison)
- `(★dec)` ⟹ `∀B: Var_μ(s_B) = |Cov(s_B,s_B)| ≤ C`.                (diagonal *boundedness*, from `dist=0`)

They share only the weak consequence "`Var(s_B)` is bounded"; each carries an **extra, different** clause:

```
(★)      =  [Var bounded]  ∧  [Var ≤ E[H]+W−δ_W]      (a comparison with the mean Hessian)
(★dec)   =  [Var bounded]  ∧  [Cov(s_B,s_B′) decays]  (a decay across blocks)
```
The comparison clause and the decay clause are logically disjoint, hence:

> **Independence [T].** `(★)` and `(★dec)` are logically independent: their intersection is
> "`Var(s_B)` bounded"; `(★)` adds the Hessian comparison, `(★dec)` adds the covariance decay, and neither
> extra clause implies the other.

---

## 3 · The implications, proved or disproved

### 3.1 `(★) ⇒ (★dec)`? — **No `[T]`.**

*Proof (general).* `(★)` gives `Var(s_B) ≤ E[H_B]+W_B−δ_W ≤: V_max` (finite per block). By Cauchy–Schwarz,
`|Cov(s_B,s_{B′})| ≤ (Var_B Var_{B′})^{1/2} ≤ V_max`. This is a **uniform** (`r`-independent) bound; it
cannot produce `≤ C e^{−κ r/a_j}` for any `κ>0` (take `r→∞`). No amount of diagonal information bounds an
off-diagonal entry below the diagonal. ∎

*Counterexample (Gaussian force field, block lattice; witness §8.2).* Take a probability measure on the block
variables under which the force field `(s_B)` has mean zero and covariance `Σ(r)`; the diagonal
`Var(s_B)=Σ(0)=A` is fixed, so `(★)` can be met by positing the (free) curvature `E[H]+W = A+δ_W`. Set
`Σ(r) = A/(1+r)^{α}`, `α = d−2+ε` (the massless Green's function; positive-definite). Then
`Cov(s_B,s_{B′})` decays as a **power** `r^{−α}`, so `(★dec)` fails for **every** `κ>0`, while `(★)` holds.
This is realisable as a Gaussian measure on the block variables; it is the log-density of a **critical
(gap-1)** field. Witness prints the log–log slope `≈−α` (no exponential) and the row-sum growing with the
torus size. **Scoping note `[O]`:** for the *actual Wilson fiber* `S_W` is finite-range, so the *free* part of
`Cov(s,s)` is the coarse precision (local, `§5`); the power-law piece is the **interacting (entropic)**
correction, which is exactly the object `L3-band-decay` says is open. Hence: the non-implication is `[T]` at
the level of *measures satisfying (★)*; for the Wilson fiber it says `(★) alone cannot certify (★dec)` — the
certificate is the open band-decay lemma.

### 3.2 `(★dec) ⇒ (★)`? — **No `[T]`.**

*Proof.* `(★dec)` at `B=B′` gives `Var(s_B) ≤ C`. `(★)` is `Var(s_B) ≤ E[H_B]+W_B−δ_W`. These agree only if
`E[H_B]+W_B−δ_W ≥ C`, i.e. `E[H_B]+W_B ≥ C+δ_W`. `(★dec)` says nothing about `E[H_B]`. Since `E[H_B]` can be
**negative** — Savvidy floor `E[H]≥−β C_N #□`, numerically `≈−1.1β` at `β=2.5` (`PROVE-Gpp-positive §3`) —
the comparison can fail while `(★dec)` holds. ∎

*Counterexample (witness §8.3).* An exponential kernel (`(★dec)` holds, `κ>0`) with the diagonal variance
`C` and the *indefinite* mean Hessian `E[H]+W = −1.1β + W < C+δ_W` gives `D_B < δ_W`: `(★)` fails, `(★dec)`
holds. Boundedness of a variance is not curvature of the action.

### 3.3 The one true *conditional* implication `[T]`

The exact content of 3.1–3.2, collected:
```
(★dec) with diagonal C   +   [ E[H_B]+W_B ≥ C + δ_W  ∀B ]   ⇒   (★).          [T, trivial]
```
So `(★)` and `(★dec)` are not "independent" in the sense of unrelated — they are two clauses on one matrix —
but **neither implies the other**; the bridge is the explicit quantitative comparison, which at the head is
`E[H]+W ≥ C+δ_W` and is exactly where the Savvidy negativity bites.

### 3.4 `(★) ∧ (★dec) ⇒ G″ ⪰ ½δ_W I` — the two-line UV step `[T]`

**Proposition.** If `(★)` holds with margin `δ_W` and `(★dec)` holds with `(C,κ)` satisfying
`R := Σ_{B′≠B}|Cov(s_B,s_{B′})| ≤ ½δ_W` (uniform in `B`), then `G″ ⪰ ½δ_W I` on the coarse space.

*Proof (two lines).* By §1.4, `G″ = D − (Cov−diag)`. By `(★)`, `D ⪰ δ_W I`. By Gershgorin / diagonal
dominance, `λ_min(G″) ≥ min_B (D_B − R) ≥ δ_W − ½δ_W = ½δ_W > 0`. ∎

This is the precise sense in which the UV step needs **both**: `(★)` supplies the diagonal `δ_W`, `(★dec)`
keeps the off-diagonal row-sum `R` below `δ_W`. It is *not* "`(★)` is the step, `(★dec)` is a technicality" —
without `(★dec)`, `G″` can have a negative eigenvalue purely from off-diagonal covariance (witness §8.1: a
long-range kernel with bounded diagonal drives `λ_min` negative).

### 3.5 The *only* bridge `(★) ⇒ (★dec)`: a uniform **global** fiber LSI `[C]`

`(★)` gives a variance bound; exponential covariance decay needs mixing. The standard mechanism
(Zegarlinski/Guionnet–Zegarlinski, `PROOF-H2-GLOBAL`):
```
[uniform global fiber LSI with constant 1/λ_f]  +  [finite-range interaction]
      ⇒  |Cov(s_B,s_{B′})| ≤ C₀ e^{−κ r/a_j},   κ ≍ √λ_f · a_j / (block scale).          [C]
```
But the **proven** fiber constant is only the block-local Holley–Stroock `λ_f = Λ_UV = e^{−2β_j N_p}/C_G`
(`PROVE-H2-uvgap`) — it collapses into the UV (`G1-L3-TOWER-POINCARE` §2) — and upgrading block-local to
**global** is the Zegarlinski contraction `C(β)=sup_k Σ_{j∼k}δ_{kj}`, which at `β=2.5` is `C≈15≫1`
(`PROOF-H2-GLOBAL`). So:
- the *rate* from the proven block gap is `κ ≲ √(e^{−2β_j N_p}) = e^{−β_j N_p} → 0`, and the row-sum
  `≳ c_d e^{β_j N_p d}` diverges — **no** `(★dec)` from proven inputs;
- a *uniform global* LSI (H2-global) **would** give `(★dec)` with `κ≍√λ_f`, but that hypothesis **is** the
  head wall (`C(β)<1` fails at `β≈2.5`), the same wall as W6/H2-globalization/BREAK-IN-1.

Hence `(★) ⇒ (★dec)` holds **only** under H2-global `[C]`, not from `(★)` itself, and H2-global is open at the
head. This is *not* a new no-go: it is the identification of the exact missing hypothesis.

---

## 4 · Constants: `κ` vs `δ_W`, and the multiplicative (decimation) structure

### 4.1 The row-sum

`(★dec)` with `n = dist/a_j` the block separation. The exact lattice row-sum is
`S_d(κ) := Σ_{r∈Z^d\setminus0} e^{−κ|r|}`, whose small-`κ` asymptotics is the polar integral
```
S_d(κ)  ~  c_d / κ^d,     c_d = |S^{d−1}|·Γ(d) = 2 π^{d/2} Γ(d) / Γ(d/2)
        [d=2: c_2 = 2π ≈ 6.28;   d=4: c_4 = 12π² ≈ 118.4].
```
So the row-sum `R ≤ C·S_d(κ) ≍ C c_d/κ^d` is **finite for every `κ>0` and volume-free** — the same Λ-tempering
role that Lemma B plays for large-field entropy (`C2-BALABAN-UNIFORM §2.2`): a *tempered per-unit-volume*
quantity, not a total. (For a **power-law** kernel `S` diverges — §3.1, witness §8.2 — so the *exponential*
form is what makes the row-sum volume-free.)

### 4.2 The interface condition `[T]`

`R ≤ ½δ_W` holds in the small-`κ` regime if
```
κ  >  κ*(δ_W,C) := ( 2 c_d C / δ_W )^{1/d}.        [d=4:  κ* = ( 24π² C / δ_W )^{1/4} ≈ (237 C/δ_W)^{1/4} ]
```
(`c_d` = polar constant above; the factor `2` makes `R≤½δ_W`.) This is the quantitative coupling of the decay
rate `κ` and the curvature margin `δ_W`. Witness §8.1 sweeps `κ`: on a positive covariance the row-sum bound is
**exact** (`λ_min(G″)=δ_W−R` by Perron–Frobenius on the positive circulant), giving a sign flip at a **finite**
`κ*≈2.52` in `d=2` (`C=1,δ_W=½`); the small-`κ` polar estimate `3.55` is a slight **overestimate**, as expected
(the lattice is sparser than the continuum at `r=1`). Two regimes:
- **deep UV** (`g_j→0`): `δ_W=O(1)` (Lemma B, scale-independent) and the fiber is massive in lattice units,
  so `κ=O(1)` in block units — `κ>κ*` is **easy**;
- **head** (`β≈2.5`): `κ` is set by the *box floor* `c_L>0` only if `(★dec)` is granted; the proven block
  constant gives `κ≈e^{−β_jN_p}`, and the Savvidy negativity pushes `δ_W` *down* — `κ>κ*` is **hard**. Same
  UV-easy / head-hard split as the campaign's `(★)`.

### 4.3 Decimation preserves `κ`, not `C` `[T]`

Let `A` be an operator on the shell-`j` block lattice with `|A(B,B′)|≤C e^{−κ|B−B′|}` (block distance), and
`P′` the next block averaging (factor `L`, kernel supported in `O(1)` blocks). Then
```
|(P′AP′ᵀ)(B,B′)|  ≤  C′ e^{−κ′ |B−B′|},   κ′ ≥ Lκ − o(1) ≥ κ,   C′ ≤ (block multiplicity)·C.
```
*Reason.* Each coarse entry is a sum of `O(1)` shell-`j` entries at fine block-distance `≈ L·`(coarse
distance); the exponential survives with the **same physical decay length**, i.e. **stronger** `κ` per coarse
block, while the prefactor picks up the `O(1)`-per-step multiplicity. So the "multiplicative structure of the
RG map" **helps `κ` and hurts `C`**:
```
κ is preserved (indeed ≥ Lκ per step);   C may grow as (a/a_j)^d = L^{Jd} over J steps.
```
**Consequence:** `(★dec)` must be stated with `C` **per unit volume** (tempered), or the accumulated
prefactor `C·L^{Jd}` will beat `δ_W` at the bottom of the tower. This is the same "per-unit-volume vs total"
issue Lemma B settles for the large-field channel — here it recurs in the covariance prefactor.

### 4.4 Cross-shell composition is a **min**, not a product `[T]`

The decimated Hessian at shell `j+1` is `G″^{(j+1)} = D^{(j+1)} − Cov^{(j+1)}` where both entries are
`P′`-transforms of shell-`j` objects; `G″⪰½δ_W` at every shell requires `(★)` and `(★dec)` **at every shell**
with `κ>κ*` and `R<½δ_W` **uniform in `j`**. The tower constant is the `min` over shells
(`G1-L3-TOWER-POINCARE` §5), so it is enough that `(C_j,κ_j,δ_W)` be `j`-independent; by §4.3
`κ` is fine and only the **tempered `C_j`** must be checked. This is the precise "cross-shell" content.

---

## 5 · Free-field acid test (it must bite the *inequality*, not the implication)

Free massless scalar on a ring/lattice, `K=Q^{−1}`, `Q=−Δ` (massless: `Q(k)=4Σ_μ sin²(k_μ/2)`, `Q(0)=0`).
With a sharp passband `P` commuting with `Q`:
- **`G″ = (P K Pᵀ)^{−1}`** is the **coarse free Laplacian**: `G″(k)=Q(k)` on the passband, hence
  `G″(0)=Q(0)=0`. So `(★)` — `G″ ≥ δ_W` — **fails** at the `k→0` mode. **[V]**
- **`Cov(s,s) = PQPᵀ − G″`** (the identity `G″=E[H]−Cov`, with `E[H]=PQPᵀ`), and each term is a **coarse
  local operator** — `PQPᵀ` the coarse Laplacian, `G″` the coarse precision — so the kernel is **short-range**
  (witness §8.4: geometric decay, `~0.17` per block). So `(★dec)` **holds**, even though `(★)` fails. **[V]**

(The linear-Gaussian statements above are exact for a **sharp** passband `P` commuting with `Q`; the witness
uses the *physical* block-averaging `P`, where `G″(0)=m²` (→0 as `m→0`) and `Cov(s,s)` decays geometrically
`~0.17`/block — the same qualitative facts, §8.4.) This is the cleanest statement that the two are **not one
inequality**: the free field **fails `(★)`** and **passes `(★dec)` (trivially)**. The acid test is satisfied where it must be
— on `(★)` — and `(★dec)` is seen to be a different (locality) object. (Do not confuse this with the L3
*operator-norm* statement `sup_a Z_a^{−1}‖H_a‖_D<∞`, which *passes* on the free field; the *inequality*
`(★)`, a positivity statement, fails. `C2-BALABAN-UNIFORM §5`.)

---

## 6 · Implication diagram

```
                 (★)                                  (★dec)
        [diagonal of G″: D_B ≥ δ_W]          [off-diagonal: |Cov(s_B,s_B′)| ≤ C e^{−κr/a_j}]
            "curvature comparison"                       "mixing / locality"
                  │  ╲                                   ╱  │
   Var ≤ V_max ──►│   ╲   ⇒ only Var ≤ V_max          ╱   │  ◄── needs H2-global C(β)<1
     (boundedness) │    ╲  (Cauchy–Schwarz)          ╱    │      [C], fails at head (C≈15)
                   │     ╲                          ╱     │
                   │      ╲    ✗ no implication    ╱      │
                   │       ╲                      ╱       │
                   │        ▼                    ▼        │
                   │        (★) ∧ (★dec)   +   [E[H]+W ≥ C+δ_W]   ── (§3.3, [T])
                   │                    │
                   │                    │  Gershgorin / diagonal dominance (§3.4, [T])
                   │                    ▼
                   │             G″ ⪰ ½δ_W I  (L3 operator positivity)
                   │                    │
                   │                    │  + tower uniformity: κ, C/vol, δ_W  j-independent (§4)
                   │                    ▼
                   │             block-Poincaré leg of L3      [C]
                   │
                   └──── (★dec) at B=B′ gives Var ≤ C;  but (★) needs the Hessian half
                         (E[H] may be Savidy-negative)  ⇒  (★dec) ⇏ (★)   (§3.2, [T])
```
Green arrows (implications that hold): `(★dec) ⇒ Var bounded`; `(★)∧(★dec) ⇒ G″⪰½δ_W`;
`(★)∧[E[H]+W≥C+δ_W]∧(★dec) ⇒ (★)`. Red (implications that fail): `(★)⇒(★dec)`, `(★dec)⇒(★)`.
Dashed: `(★)+[H2-global] ⇒ (★dec)` is `[C]` and open at the head.

---

## 7 · Status ledger

| item | statement | label |
|---|---|---|
| decimation map `P`, fiber conditional `μ(·\|U′)`, score `s_B` | §1.1–1.3, campaign's `PROVE-H3`/`L3-band-decay` | **[T]** |
| exact identity `G″=E[H]−Cov(s,s)+W` | §1.4, Hessian form of law of total variance | **[T]** |
| `(★)` = diagonal `D_B≥δ_W`; `(★dec)` = off-diagonal decay | §1.4–2 | **[T]** |
| `(★) ⇏ (★dec)` | §3.1; Gaussian power-law score; Wilson version `[O]` | **[T]** (abstract), **[O]** (Wilson) |
| `(★dec) ⇏ (★)` | §3.2; Savvidy-negative `E[H]` vs variance bound `C` | **[T]** |
| `(★dec)+[E[H]+W≥C+δ_W] ⇒ (★)` | §3.3 | **[T]** |
| `(★)∧(★dec) ⇒ G″⪰½δ_W` | §3.4, Gershgorin | **[T]** |
| `(★) ⇒ (★dec)` needs uniform global fiber LSI | §3.5; = H2-global `C(β)<1`, `C≈15` at head | **[C]**, open at head |
| interface `κ > (2c_d C/δ_W)^{1/d}`, `c_4=12π²≈118` | §4.2 | **[T]** |
| decimation preserves `κ` (`κ′≥Lκ`), not `C` | §4.3 | **[T]** |
| cross-shell = min over shells; need `C_j` tempered | §4.4 | **[T]** |
| free field: `(★)` fails (`G″(0)=0`), `(★dec)` holds (`κ=∞`) | §5 | **[V]** |
| **the UV step** | still `[C]` on `(★)∧(★dec)`; now *separated*: diagonal comparison vs off-diagonal mixing | **[C]** / **[O]** |

**Sharpest line.** `(★)` and `(★dec)` are the **diagonal** and **off-diagonal** of one matrix `G″`; they are
*independent* (`[T]`), their intersection is only "variance bounded," and the UV step `G″⪰½δ_W` is their
**conjunction** via Gershgorin with the explicit interface `κ>(24π²C/δ_W)^{1/4}`. Decimation preserves the
decay rate but inflates the prefactor, so `(★dec)` must be tempered per unit volume — the same fix Lemma B
applied to the large-field channel. The free field fails `(★)` and passes `(★dec)`: they are not one
inequality renamed. Clay 0%.

---

## 8 · Witness (`research/v3/w7/s4_decimation.py`)

Deterministic, `py -3.12`, numpy 2.4.6; all four sections run (exit 0). Quoted output (abridged):

```
[S4.1] Gershgorin interface  (d=2 torus, C=1, delta_W=0.5)
  volume-freeness of the row-sum R(M) = sum_{r != 0} C exp(-kappa |r|):
    kappa    M=16        M=32        M=64        M=96      (R -> R_inf?)
    0.30     51.0310    66.6220    68.8562    68.8814
    0.50     22.5676    24.2042    24.2465    24.2465
    0.80      8.9220     8.9979     8.9981     8.9981
    1.20      3.6280     3.6298     3.6298     3.6298
  spectrum on M=24:  G'' = delta_W I - (Cov - diag);  lambda_min = delta_W - lambda_max(off)
    kappa    rowsum R     lambda_min     Gershgorin lb   (star dec)->pos?
    0.20    116.0815   -115.5815     -115.5815       NO
    0.80      8.9943     -8.4943       -8.4943       NO
    1.60      1.8017     -1.3017       -1.3017       NO
  empirical kappa* (lambda_min crosses 0, M=240) = 2.519
  small-kappa continuum estimate sqrt(c_d C/delta_W) = 3.545  (c_2=2*pi; an UPPER heuristic)
  => (star) alone (diagonal = delta_W) does NOT give positivity; sign flip at kappa*=2.52.
[S4.2] (star) NOT=> (star dec): bounded diagonal, power-law covariance
  same diagonal A=Var(s_B)=1, kernel A/(1+r)^alpha, alpha=0.50;  log-log slope = -0.447
    r=2:0.5774  r=4:0.4472  r=8:0.3333  r=16:0.2425  r=32:0.1741  (power law: NO exp decay)
  row-sum R(M):  16:100.2  32:299.7  64:874.1  96:1623.5  128:2513.7  (diverges with volume)
  => diagonal bounded (Var=A=1) yet (star dec) fails for every kappa>0.
[S4.3] (star dec) NOT=> (star): exponential kernel C=1,kappa=0.8 (decay HOLDS)
  mean Hessian E[H]+W = -1.1*beta = -2.75 ; Var=1.00 ; delta_W=0.50
  diagonal D_B = -3.75 < delta_W  =>  (star) FAILS  (boundedness != curvature).
[S4.4] free-field acid test (1D ring N=128, factor-2 block averaging)
  mass m     min eig G'' (=G''(k=0)=m^2)   max|Cov_s|   max|Cov_s(r>4 blocks)|
    0.50      2.500000e-01                 1.577e-01    3.013e-04
    0.10      1.000000e-02                 1.710e-01    4.147e-04
    0.02      4.000000e-04                 1.715e-01    4.203e-04
  Cov_s kernel decay profile (m=0.10): s=1:1.5e-02 4:2.4e-03 6:7.1e-05 8:2.1e-06 11:1.0e-08
  => G''(k=0)=m^2 -> 0 (star FAILS on free massless field); Cov_s short-range (star dec HOLDS).
```

`S4.1` confirms §4.2 (row-sum volume-free, exact sign flip `κ*≈2.52` for a positive covariance); `S4.2`
confirms §3.1 (power-law ⇒ no decay, row-sum diverges); `S4.3` confirms §3.2 (decay ⇏ curvature); `S4.4`
confirms §5 (`G″(k=0)=m²→0` ⇒ `(★)` fails; short-range score covariance ⇒ `(★dec)` holds).

*Refs.* `C2-BALABAN-UNIFORM §2.2,§2.4,§5`; `ch03-existence.tex §G1.2`; `G1-UNIFORM-ZA §1,§3`; `PROOF-L3`;
`L3-scale-uniform`; `L3-band-decay`; `G1-L3-TOWER-POINCARE §2,§4,§5`; `PROVE-Gpp-positive §1,§3,§4`;
`PROVE-H3` (P, score, (a),(c)); `PROVE-H2-uvgap` (`Λ_UV`); `PROOF-H2-GLOBAL` (`C(β)`, Zegarlinski);
`PROOF-COARSE-POINCARE (c)`; `BLOCKING-DESIGN` (score = mixed Hessian); `D1-BLOCK-MARTINGALE` (block
decomposition, martingale identity); Shen–Zhu–Zhu CMP 400 (2023) arXiv:2204.12737 (β<1/48, off the tower);
Guionnet–Zegarlinski (uniform LSI ⇒ clustering); Holley–Stroock Ann. Probab. 15 (1987).
