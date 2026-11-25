# W3-SMALL-VOLUME — exact small-volume alignment `w₁=f₁²/M₀` at SU(2), `β=2.5`: **every exactly solvable small volume gives `w₁=0` or a degenerate atom — small volume does *not* already deliver positivity; the only exact positive small-volume anchor is strong coupling `w₁=s/2` (L-uniform), and the full 4D `L=2,3` at `β=2.5` is `[O]`**

> **2026-11-24 · Wave W3 — the exact small-volume end of the weight problem.**
> The W-charter (§1 W3) asks: compute `w₁` **exactly** on `L=2,3` at `β=2.5` by transfer matrix /
> character expansion (no Monte Carlo); give the exact table, the finite-size law, and decide whether
> `w₁` is **already positive at small volume**. The master identity (`U1`) is
> `w₁=|<ψ₁|o−1>|²/Var(o)`, `o=O/⟨O⟩`, `O=¼tr F²` (the local source); `ψ₁` the first excited transfer
> eigenvector. **Reads.** `W-WAVE-CHARTER.md` §1–§2; `U1-HARNACK-BOUND.md` (master identity, NG3),
> `U3-LARGE-N-BINDING.md` (`w₁=s/2`), `U4-CONTINUUM-ATOM.md` (γ=0), `R5-EXACT-ANCHORS.md` (2D/one-plaquette/
> β=0), `P4-PERIODIC-TRANSFER.md` (transfer kernel), `D1-EXACT-L.md`; `PROOF-P3-UNCONDITIONAL.md`.
> **Witness.** `research/v3/w7/w3_small_volume.py` (`py -3.12 -I`, numpy/scipy, deterministic, no MC,
> ~1 s, **10/10 PASS**; output quoted §7). No third-party downloads.
> **Labels (charter §2).** **[T]** true theorem/identity · **[C]** conditional/model · **[NG]** no-go for a
> named class · **[O]** exact missing step · **[V]** numerical witness here. **No QED. Clay 0%.**

---

## 0 · Verdict (up front; no QED, Clay 0%)

1. **At `β=2.5`, every small volume where `w₁` is *exactly defined* returns `w₁=0` or a degenerate atom — none is positive `[T]+[V]`.** The exact anchors are: (i) **2D** `tr F²=H` (conserved energy) ⇒ the connected source annihilates the vacuum ⇒ `f₁²=M₀=0`, `w₁=0`; (ii) the **periodic 1D Cartan transfer matrix** at `β=2.5`: the first excited eigenvector `ψ₁` is **odd** while the local source is **even** ⇒ `⟨ψ₁|o−1⟩=0` exactly ⇒ `w₁=0` at every `L`; (iii) **`β=0`** Haar: rank-one kernel ⇒ no time correlation ⇒ `w₁=0`; (iv) **free field**: the composite is number-even, the one-particle state number-odd ⇒ orthogonal ⇒ `w₁=0` (acid); (v) **one plaquette / 4D `L=1` single hypercube**: a *single time slice*, no first excited state — the atom is **degenerate/undefined**, not positive.

2. **The one exact positive anchor is strong coupling, and it is `L`-uniform, not a small-volume accident `[T]`.** At `β<1/48` the P3 tube residue gives `a₀=(2N)^{-2}s`, `M₀=Var((1/N)Re\,tr U)=1/(2N²)`, hence `w₁=s/2=O(N⁰)`, **uniform in `L`** (Witness [F]). So **small volume *can* be positive — but only when the dynamics *binds*** a four-wall `0⁺⁺` tube that the source reaches; that is a strong-coupling binding fact, absent in every solvable `β=2.5` small volume.

3. **The finite-size law points the *wrong* way for a floor `[V]/[T]`.** The solvable symmetric cases have `w₁(L)=0` **exactly, `L`-independent**. The only small-volume mechanism that gives a local source weight on a light mode is **delocalisation**, and it gives the level-1 **fraction** `=n_lev/(L−1)=2/(L−1)` for `L≥4` (`=1` at `L=2,3`), which **decreases** with `L` (Witness [G]); it is not a uniform floor and it is **not** the `tr F²` composite atom (which is `0` by parity/conservation). So small volume is where a generic local source has its *largest* level-1 weight, yet the *atom fraction of the composite* is still zero — the two facts coexist and neither certifies a floor.

4. **The full 4D `L=2,3` lattice at `β=2.5` is `[O]`.** The exact character expansion is the only MC-free route, and at `β=2.5` it is **outside its convergence radius** (`uK≫1`; activity `≈1.01`, Ursell `≈4.6×10³`), so a truncated character sum is a controlled *model* (a character-truncated lattice), not the exact theory. Its 4D `L=2,3` value is exactly known only in the two solvable limits — **strong coupling `s/2`** and **`β=0` zero** — and at `β=2.5` remains the open number. The obstruction is stated precisely in §5.

**One line.** *At `β=2.5` the exact small-volume atom fraction is zero or degenerate in every case a deterministic MC-free method reaches (2D conserved/center, periodic-strip parity, `β=0` contact, free parity, no-time one-plaquette/4D-`L=1`); the only exact positive small-volume anchor is the strong-coupling `w₁=s/2`, which is `L`-uniform and is *binding*, not small volume; the finite-size trend from delocalisation is `2/(L−1)→0`, never a floor; and the full 4D `L=2,3` at `β=2.5` sits beyond the character-series radius, hence `[O]`. **Small volume does not already give `w₁>0` at `β=2.5`; positivity is binding.** No QED. Clay 0%.*

---

## 1 · The object and what "exact at `L=2,3`" can mean

`w₁=f₁²/M₀` with `f₁²=|<ψ₁|O|Ω>|²`, `M₀=Var(O)`. By the U1 master identity it is the squared
direction cosine of `o−1` onto `ψ₁`. "Exact small-volume" therefore needs **three** things the finite
lattice supplies only in special models: (a) the transfer spectrum `{λ_n,ψ_n}`, (b) the source matrix
elements `⟨ψ_n|O|Ω⟩`, (c) `Var(O)`. The four tools below are exactly solvable; each computes `w₁` with
**no approximation and no Monte Carlo**.

| tool | what it fixes exactly | `w₁` at `β=2.5` | label |
|---|---|---|---|
| 2D continuum YM | `tr F²=tr E²=H` (conserved density); spectrum `C₂(j)` | **0** (source annihilates vacuum) | **[T]** |
| periodic 1D Cartan transfer (P4 kernel) | spectrum `λ₀>λ₁>…`, eigenvectors `ψ_n(z)` | **0** (`ψ₁` odd, source even) | **[T]+[V]** |
| `β=0` Haar | rank-one transfer kernel | **0** (t=0 contact) | **[T]+[V]** |
| free field | `:φ²:` / `tr F²` number-even vs 1-particle odd | **0** (parity) | **[T]** |
| one plaquette / 4D `L=1` | one time slice, no `ψ₁` | **degenerate** | **[T]** |
| strong coupling `β<1/48` | P3 four-wall tube residue | **s/2** (>0, `L`-uniform) | **[T]** |

The 4D `L=2,3` *primary* lattice is **not** in this table: no deterministic exact method reaches it at
`β=2.5` (§5). The table is the honest exact answer to *which* small volumes can and cannot certify a
positive `w₁`.

---

## 2 · 2D SU(2): exact `w₁=0`, by conservation **and** by representation theory `[T]`

Two independent exact statements, both verified in Witness [A].

**(a) Continuum operator — conservation.** In `1+1`D the field strength has one component, `F₀₁=E`, so
`tr F²=tr E²=H` *is* the energy (R5 §1.1, standard 2D solvability). The connected source is
`H−⟨H⟩`, and `H|Ω⟩=0`, so `(H−⟨H⟩)|Ω⟩=0`: the source **annihilates the vacuum**, the connected two-point
function is a pure `t=0` contact, `f₁²=0` and `M₀=0`, and `w₁=0` (the atom is absent). This is the 2D acid
test and it holds at **every** `L`.

**(b) Character selection rule.** `∫dU\,χ_jχ_k=δ_{jk}` checks to `1.7×10⁻¹⁵`; the SU(2) CG rule
`χ_{1/2}χ_j=χ_{|j−1/2|}+χ_{j+1/2}` (j>0), `χ_{1/2}χ_0=χ_{1/2}`, checks to quadrature precision. The
lightest `0⁺⁺` (center-even, integer-`j`) level above the vacuum is `j=1` (`C₂=2`), and
```
   ⟨χ_1|χ_{1/2}|χ_0⟩ = ∫χ_1χ_{1/2}χ_0 = 0     (tensor rule sends j=0 → j=1/2 only).
```
So the fundamental character has **zero** vacuum→lightest-`0⁺⁺` matrix element (Witness [A]).
**[C] caveat:** the identification of the *lattice plaquette* with `χ_{1/2}` of the holonomy is a
model statement (the exact 2D result is the continuum operator (a)); we record (b) as the group-theory
selection rule behind the vanishing, and label the *operator* identification `[C]`, consistent with
`R5`'s "2D `tr F²` atom `=0`". Either way, `w₁(2D)=0` is `[T]` for the continuum `tr F²`.

---

## 3 · Periodic 1D Cartan transfer matrix at `β=2.5`: `ψ₁` odd ⇒ `w₁=0`, exactly, all `L` `[T]+[V]`

The P4/D3 periodic bulk kernel is
```
   K(z,z') = √(ρ₁(z)) e^{β cos((z−z')/2)} √(ρ₁(z')) ,     ρ₁(z)=π(Cin(2π)−Cin|z|) ,
```
`β=2.5`, `z∈(−2π,2π)` on an exact Gauss–Legendre grid. Diagonalising (Witness [B]):
```
   λ₀=228.1738, λ₁=119.7875, λ₂=63.8178 ;  λ₁/λ₀=0.524984,  Δ=0.6444,  ξ=1.5519
```
(matches P4's `0.52498 / 0.6444 / 1.5519` to `10⁻⁵`). The kernel is parity-symmetric
(`K(−z,−z')=K(z,z')`), so eigenvectors carry definite parity; the **ground state is even**
(oddness `4×10⁻³¹`) and the **first excited state is odd** (oddness `4.00`). The local source
`o(z)=1−cos(z/2)` is **even**. Hence
```
   ⟨ψ₁|o−1⟩ = 0   ⇒   w₁ = 0     (machine zero: 8.4×10⁻³⁴) .
```
The full even-source weight sits in the **next even level**: `w=(0.119, 0, 0.576, 0, 0.237, 0, 0.057, …)`
(Witness [B]). Because the bulk kernel is translation invariant, `w₁` is **`L`-independent** — the
periodic strip gives `w₁(L)=0` for `L=2,3,4,…` The mechanism is the cleanest possible instance of the
charter's small-volume obstruction: **a symmetry that `T` respects makes `ψ₁` orthogonal to the source,
no matter how small the volume.**

---

## 4 · `β=0`, free field, and the degenerate one-time-slice cases

- **`β=0` Haar `[T]+[V]`.** The kernel degenerates to rank one (`λ₁/λ₀=4.9×10⁻¹⁶`, Witness [C]):
  links are i.i.d., the `tr F²` two-point function is a `t=0` contact, the KL measure has no atom,
  `w₁=0`. (This is the `R5 §1.3` anchor: compactness gives the *ceiling* `χ=3.7198`, **not** the atom.)
- **Free field `[T]`.** The composite `:φ²:` / `tr F²` is number-even; the one-particle state is
  number-odd, so `⟨1|O|0⟩=0` (Witness [D]: `0.00e+00`). The composite reaches only the two-particle cut
  ⇒ `f₁²=0`, `w₁=0`. This is the acid test: any small-volume "theorem" producing `w₁>0` for the free
  composite is refuted here.
- **One plaquette / 4D `L=1` single hypercube `[T]`.** A **single time slice** has a `1×1` transfer
  matrix (vacuum only): there is **no first excited state**, so `w₁` is **degenerate/undefined**, not
  "positive". The exact one-plaquette curvature at `β=2.5` is `G″(0)=0.729924` (Witness [E]) — a real,
  positive **ceiling** `χ=1/G″`, but a curvature is not an atom (at `β=0` the same model has `G″=0.205`
  while `f₁²=0`). The 4D `L=1` periodic single hypercube (4 links, 6 plaquettes, the minimal lattice
  carrying the Savvidy commutator) is the same `T=1` case: one Euclidean time slice, hence degenerate.

This is exactly the charter's acid test #3: *"no time direction; the atom is degenerate; the route must
say so, not invent one."* We say so.

---

## 5 · The full 4D `L=2,3` lattice at `β=2.5`: why it is `[O]`, and what *is* exact

The exact MC-free route on the 4D torus is the **character expansion**,
`Z[O]=∫∏_ℓ dU_ℓ ∏_p[Σ_j c_j(β)χ_j(U_p)]·O`, with the link integrals done by SU(2) recoupling. Two exact
limits, one open number:

- **`β=0`:** all `c_{j>0}=0`; i.i.d. links; `w₁=0` (§3). Exact at every `L`.
- **`β<1/48` (strong coupling):** the P3 four-wall tube is the unique leading connected polymer; its
  residue is a positive modulus-square (`f⊗f̄→1`), and the β=0 normaliser is exact ⇒ `w₁=s/2>0`,
  `L`-uniform (Witness [F]; `R5 §1.4`, `U3 §4`). **This is an exact 4D `L=2,3` value — at strong
  coupling.**
- **`β=2.5`:** the same series has left its radius. The character coefficients do *not* have a small
  parameter at `β=2.5` (activity `≈1.01`, Ursell `≈4.6×10³`, `≫1`; `T2 §5.1`, `R5 §1.4`), so a
  truncated sum is a *character-truncated model*, not the exact SU(2) theory. The 4D `L=2` lattice has
  64 links / 96 plaquettes and 6-valent link vertices; an exact full recoupling contraction at `β=2.5`
  is outside both the convergence radius and practical deterministic reach. **The exact 4D `L=2,3`
  `w₁` at `β=2.5` is `[O]`.**

**What the exact limits already settle:** the answer to "is `w₁` positive at small volume?" is
**coupling-dependent, and at `β=2.5` the exact solvable cases say no.** Strong coupling *does* make
small-volume `w₁` positive (`s/2`), because binding fuses the flux into a single reachable `0⁺⁺` tube; at
`β=2.5` no solvable small volume binds one. The 4D `L=2,3` number at `β=2.5` would be a *lattice mass
gap / binding* datum, not a small-volume triviality.

---

## 6 · The finite-size law, stated exactly

| case | `w₁(L)` | behaviour | label |
|---|---|---|---|
| 2D continuum `tr F²=H` | `0` | `L`-independent zero (conserved) | **[T]** |
| periodic 1D Cartan transfer | `0` | `L`-independent zero (parity) | **[T]+[V]** |
| `β=0` Haar | `0` | `L`-independent zero (rank-1 / contact) | **[T]+[V]** |
| free composite | `0` | `L`-independent zero (parity) | **[T]** |
| one plaquette / 4D `L=1` | degenerate | no `ψ₁` (one time slice) | **[T]** |
| strong coupling `β<1/48` | `s/2` | `L`-**uniform** positive floor (binding) | **[T]** |
| local point source on a delocalised mode | `n_lev/(L−1)` (`=2/(L−1)`, `L≥4`; `=1` at `2,3`) | **decreases** with `L`; not a floor; not the composite atom | **[V]/[T]** |

**Reading.** The solvable small-volume cases give a *flat zero* (`w₁(L)=0`), not a floor. The only
finite-size *growth* mechanism (delocalisation of a point source) gives a weight that shrinks with `L`
and is not the composite atom. The one genuine `L`-uniform positive floor (`s/2`) is *coupling*-driven
(strong-coupling binding), not volume-driven. **So there is no small-volume window, at `β=2.5`, that
certifies `w₁≥w₀>0`; and the finite-size trend available to a local source is decreasing, the opposite
of a uniform floor.** This is the W3 end of the U1-NG2/NG3 picture, now computed.

---

## 7 · Witness — `research/v3/w7/w3_small_volume.py` (`py -3.12 -I`, 10/10 PASS)

Deterministic, numpy/scipy, no MC, ~1 s. Quoted output (abridged):

```
[A] SU(2) character orthonormality max|G-I| = 1.7e-15 ; CG chi_{1/2}*chi_j residual 1.3e-08
    <chi_1|chi_{1/2}|chi_0> = 0.00e+00  (vacuum->lightest 0++ zero)   [C operator / T continuum]
    2D continuum trF^2=H: source annihilates vacuum => f1^2=0, M0=0, w1=0   (ACID 2D)  [T]

[B] periodic Cartan kernel beta=2.5: lam1/lam0=0.524984 Delta=0.6444 xi=1.5519  (P4 match)  [T]+[V]
    psi0 oddness 4e-31 (even) ; psi1 oddness 4.00 (ODD)
    even source o(z)=1-cos(z/2): w_n = [0.119, 0.000, 0.576, 0.000, 0.237, 0.000, 0.057, ...]
    w1 (first excited) = 8.4e-34  => w1 = 0 exactly (parity), every L                    [T]+[V]

[C] beta=0: lam1/lam0 = 4.9e-16  => rank-one kernel, t=0 contact => f1^2=0, w1=0     [T]+[V]

[D] free field: <1|:phi^2:|0> = 0.00e+00  => composite orthogonal to one-particle => w1=0 (ACID)  [T]

[E] one plaquette: G''(2.5;theta=0)=0.729924 (ceiling chi=1/G''); one time slice => atom DEGENERATE  [T]

[F] strong coupling beta<1/48: a0=(2N)^-2 s, M0=1/(2N^2) => w1 = s/2 = 0.500, L-uniform  [T]

[G] point source on delocalised free cycle: level-1 fraction
      L=2:1.000  L=3:1.000  L=4:0.6667  L=8:0.2857  L=16:0.1333  L=32:0.0645  (=2/(L-1), L>=4)  [V]

[H] VERDICT TABLE (see §6).  10/10 checks PASS.  No QED. Clay 0%.
```

---

## 8 · Ledger

| item | statement | label |
|---|---|---|
| master identity | `w₁=<ψ₁\|o−1>²/Var(o)`, `Σ w_n=1` (U1) | **[T]** |
| 2D continuum | `tr F²=H`, source annihilates vacuum ⇒ `f₁²=M₀=0`, `w₁=0`, all `L` | **[T]** (acid 2D) |
| 2D selection rule | `⟨χ_1\|χ_{1/2}\|χ_0⟩=0` (verified `1e-15`); operator identification | **[T] group theory / [C] operator** |
| periodic Cartan transfer | `Δ=0.6444`, `ξ=1.552` (matches P4); `ψ₁` odd, source even ⇒ `w₁=0` exactly, all `L` | **[T]+[V]** |
| `β=0` Haar | rank-one kernel, `t=0` contact ⇒ `w₁=0` | **[T]+[V]** |
| free field | composite number-even vs 1-particle odd ⇒ `w₁=0` | **[T]** (acid free) |
| one plaquette / 4D `L=1` | one time slice, no `ψ₁` ⇒ atom degenerate (not positive) | **[T]** (acid one-plaquette) |
| one-plaquette ceiling | `G″(2.5;θ=0)=0.729924` (curvature, `χ=1/G″`, **not** atom) | **[T]+[V]** |
| strong coupling | `w₁=s/2=O(N⁰)`, **`L`-uniform**, `β<1/48` (the one exact positive anchor) | **[T]** (acid strong-coupling) |
| finite-size law (symmetric) | `w₁(L)=0` exactly, `L`-independent | **[T]+[V]** |
| finite-size law (delocalised) | point-source level-1 fraction `=2/(L−1)` (`L≥4`), **decreasing**; not a floor | **[V]/[T]** |
| full 4D `L=2,3` at `β=2.5` | exact value unknown; character series outside radius (`uK≫1`) | **[O]** |
| **W3 answer @ β=2.5** | **small volume does NOT already give `w₁>0`: solvable small volumes give 0 or degenerate** | **[T]+[V]** |
| **W3 answer @ strong coupling** | **`w₁=s/2>0` at `L=2,3` (L-uniform)** — positivity needs **binding**, not volume | **[T]** |
| **Clay mass gap** | unbuilt (this is the `W` weight input only; `W⊥I`, U6) | **[O] · 0%** |

---

## 9 · Sharpest line

**Exact small-volume `w₁` at `β=2.5` is a table of zeros and one degeneracy — and that is the finding.**
Every small volume where the MC-free transfer/character machinery reaches `w₁` gives `0` by an exact
reason: 2D by `tr F²=H` conservation (the source annihilates the vacuum), the periodic Cartan transfer by
**parity** (`ψ₁` odd, source even — checked to `8×10⁻³⁴`), `β=0` by the rank-one kernel (a `t=0`
contact), the free field by number parity, and the one-plaquette / 4D-`L=1` hypercube because a **single
time slice has no first excited state** (degenerate, and we say so rather than invent one). The **only**
exact positive small-volume anchor is the strong-coupling `w₁=s/2`, which is `L`-**uniform** and comes
from **binding** a four-wall `0⁺⁺` tube that the source reaches — a coupling fact, not a volume effect;
and the sole finite-size *growth* mechanism, delocalisation of a local source, gives `2/(L−1)`, which
**decreases** with `L` and is not the composite atom. So at `β=2.5` **small volume does not already
deliver `w₁>0`**, and the finite-size trend available to a source is the opposite of a uniform floor.
The full 4D `L=2,3` number at `β=2.5` lies beyond the character-series radius (`uK≫1`), so it is `[O]` —
but the exact limits (strong coupling `s/2`, `β=0` zero) already bracket the physics: **positivity is
binding, not size.** No QED. Clay 0%.

---

*Refs.* `W-WAVE-CHARTER.md` §1–§2; `U1-HARNACK-BOUND.md` (master identity §1, NG2/NG3 §5, free parity §7D);
`U3-LARGE-N-BINDING.md` §4 (`w₁=s/2`); `U4-CONTINUUM-ATOM.md` §0 (`γ_{w₁}=0`); `R5-EXACT-ANCHORS.md`
§1.1 (2D `tr F²` atom `=0`), §1.2 (one-plaquette `G″=0.730`), §1.3 (`β=0` `f₁²=0`), §1.4 (P3, `β<1/48`);
`P4-PERIODIC-TRANSFER.md` (`λ₁/λ₀=0.525`, `Δ=0.644`, `ξ=1.55`); `D1-EXACT-L.md` (transfer kernel);
`unconditional/PROOF-P3-UNCONDITIONAL.md`; `SM1-SOFT-EXACT.md` (4D `L=1` single hypercube, Savvidy
commutator). **Standard inputs (textbook, not re-fetched):** SU(2) character orthogonality / CG
(Weyl, Peter–Weyl); Migdal–Witten 2D solvability (Migdal, Sov. Phys. JETP 42 (1976) 413; Witten, CMP 141
(1991) 153); Källén–Lehmann. **Witness:** `research/v3/w7/w3_small_volume.py` (`py -3.12 -I`, numpy/scipy,
deterministic, no MC, 10/10 PASS, ~1 s). **No third-party downloads; no `[C]` upgraded to `[T]`; no QED;
Clay 0%.**
