# FORGE-SQUEEZE-BOGOLIUBOV — squeezing / hyperbolic canonical transformations against the indefinite Savvidy form: the squeeze is a congruence, so it preserves the tachyon; the surviving invariant *is* the mass (tautological)

> 2026-11-23 · Forge the squeezing (Bogoliubov–Valatin hyperbolic canonical) transformation as a *lower*-bound engine for the physical mass, across the indefinite Savvidy Hessian. **判詞 (verdict): NO 5/5 — a clean no-go + one structural identity, no independent lower bound.** *Does the squeeze's cost reintroduce the problem?* **Yes, exactly.** Every Bogoliubov transform is a **congruence** `Ω²↦SᵀΩ²S`; by **Sylvester's law of inertia** the signature `(N₋,N₀,N₊)` is invariant, so no real squeeze makes the physical Hessian positive — any such `S` is a singular limit `‖S‖→∞`, and the hyperbolic amplitude `e^{2r}` *equals* the negative-mode growth `e^{μt}`: it is the tachyon, gauged, not regularized. *Is there a squeezing-invariant lower bound on the mass?* **Only the tautological one:** the (Gaussian) gap is the **Williamson symplectic floor**, a *complete* symplectic invariant; any squeezing-invariant lower bound *equals* it (circular), and indefiniteness makes it imaginary (`iμ`). The exact stuck step is again **(H)=H1** — the lift must be a **non-congruence** (nonperturbative measure back-reaction), not a squeeze.

## 1 Setup
*(Distinct from `PROOF-BOGOLIUBOV.md`, which is the Kubo–Mori correlation **inequality** — an upper bound, a different object. Here `Bogoliubov` = the hyperbolic canonical **transformation**.)* On the gauge-reduced phase space `z=(A,E)` (positions = transverse 1-forms; physical sector `𝒜/𝒢`, orbit-orthogonal) the background sits in the **potential block** `Ω²≡−D²+2adF` ("frequency-squared") — the exact Savvidy operator (`PROOF-SAVVIDY-EXACT-SIGN`): `λ_min=−g|B|`, `N_neg=c_N(gB)²V/8π²`. A squeezing is a linear field redefinition `A=SÃ` (resp. real symplectic `S` on `z`). The whole question: does some `S` send indefinite `Ω²` to a positive form?

## 2 Lemma SQ-I (inertia no-go): the squeeze cannot remove the tachyon
`A=SÃ` is the **congruence** `Ω²↦SᵀΩ²S`; **Sylvester** gives `sign(SᵀΩ²S)=sign(Ω²)`. The physical form is genuinely indefinite: the transverse-magnetic Landau rung `λ_min=−g|B|<0` is orbit-transverse (not pure-gauge), and after the exact `Q⊥` quotient **7.1% stay negative, min `−0.82β`** at SU(2) β=2.5 (`PROOF-PHYSICAL-SUBSPACE-CONVEXITY`; no BRST choice deletes them, `D-BBS-INTEGRATE-NEG` §a). So `N₋>0` is congruence-invariant: **no real squeezing yields `Ω²⪰0`.** *(Bankable, exact.)* Hyperbolic check: the inverted mode `½(P²−μ²X²)` is canonically reachable only as the **saddle** `{+μ,−μ}` (real `S` with `cosh/sinh`, amplification `e^{2r}`); the stable form `½(P²+ν²X²)` is a different inertia class, unreachable. The cost `e^{2r}` **is** the tachyon growth and diverges as `μ→0`. So *the squeeze's cost reintroduces the problem* quantitatively: amplification = instability. One-loop, same fact: complex Savvidy determinant, `Im Γ=(π/2)N_neg` (`D-BBS` §b). The "rotation by `ln i`" that flips sign is purely hyperbolic, i.e. physical decay — never a stable vacuum.

## 3 Lemma SQ-W (invariant floor): the surviving object is the mass itself
For positive-definite physical `M`, **Williamson** gives `M=Sᵀdiag(σ₁,σ₁,…,σₙ,σₙ)S`, `S` symplectic, `{σ_i}` **unique** — a complete symplectic invariant; the Gaussian gap is `Δ=min_i σ_i`. Since `SᵀMS` shares `D`, **every Bogoliubov transform preserves the floor**: squeezing is passive relabeling of normal coordinates. Two horns:
- **Benign:** any *squeezing-invariant* lower bound **equals `min σ_i` = the gap** — circular (the BLIND-LOWER §1 Temple/Weinstein disease, now linear-algebraic).
- **Indefinite:** the symplectic spectrum is `σ=iμ` and, by SQ-I, stays imaginary under every real `S` — **no real floor.**

## 4 The one-loop complex form and the exact stuck step
The only non-congruence change is `δ²(½Tr log K)`. `D-BBS`: exact Fubini = vacuous (Prékopa needs *joint* log-concavity, which fails); one-loop = **complex**, `Re⪰Σ_fast>0` (Lemma NS, bankable), `Im=(π/2)N_neg≠0`. A real congruence preserves the **Hermitian inertia** (Sylvester) *and* the imaginary (antisymmetric) part; **no real squeeze removes `Im`**. The physical quadratic form is Hermitian-but-not-real, with complex symplectic bottom; the only operation lifting it to a *real positive* form uniformly in `(L,a)` is nonperturbative dynamics, not a linear canonical map:

```
squeeze-invariant: inertia(Ω²) preserved  ⇒ tachyon immovable       (SQ-I)
gap = min_i σ_i(M) = symplectic floor       ⇒ invariant ⇒ circular   (SQ-W)
lift must be non-congruence, real, > g|B| on 𝒜/𝒢, uniform (L,a):
   E_μ[PᵀHP|U′] − Var_μ(⟨Pᵀ∇S,v⟩|U′) ≥ λ‖v‖²  ⟺  G″ ≥ λ  ⟺  (H)=H1.
```

## 5 Acid test and `(L,a)`-uniformity — why it is not a qualifying tool
Free massless field: `M` positive, `σ_min=0` ⇒ the floor reports `Δ=0`, failing to manufacture positivity — but *symmetrically* (an identity, not an inequality), so it does **not** pass the acid test (vacuous on the free field, nothing beyond the gap on YM). `(L,a)`-uniformity is not saved: the finite-box floor is the geometric `(2π/L)²→0` (`TRANSPORT-MACHINE` §4), which squeezing preserves. The residue is no better: `⟨0|O|1⟩` is linear in `sinh r`, but the unsqueezed vacuum has `r=0`, so no nonzero *signed amplitude* is manufactured from nothing — the lightest-glueball profile is the gap eigenvector, again H1.

| leg | status |
|---|---|
| SQ-I: congruence preserves inertia ⇒ no squeeze removes tachyon; cost `e^{2r}`=growth | **真定理**（Sylvester; free field vacuous） |
| SQ-W: Gaussian gap = Williamson floor = complete symplectic invariant | **真定理** |
| Squeezing-invariant lower bound ⇒ equals `min σ_i` | **真定理 = 循環**（tautology） |
| Real squeeze removes `Im Γ` | **NO-GO**（congruence ⟹ Re&Im invariant） |
| `Δ≥m` crossing Savvidy, `(L,a)`-uniform, fails free field | **= H1**（未證） |

**Sharpest line.** The squeezing transform is exactly the wrong species for a *lower* bound: a **congruence** is blind to what we must change (the inertia/tachyon, by Sylvester) and *identical* to what we must bound (the symplectic floor) — immovable where it would help, circular where defined. The hyperbolic amplitude `e^{2r}` is not a technicality: it *is* the negative-mode growth `e^{μt}`, divergent at criticality `μ→0`; **the squeeze's cost reintroduces the instability exactly.** The one-loop physical form is Hermitian-but-complex (`Re⪰Σ_fast>0`, `Im=(π/2)N_neg`); no real `S` removes `Im`. Remaining: **(H)=H1**, a *non-congruence*, real, positive lift `>g|B|` on `𝒜/𝒢`, uniform in `(L,a)` — the measure back-reaction the squeezing algebra cannot supply. **Not 5/5; Clay 0%.** Salvage is structural: the (Gaussian) mass *is* a squeezing invariant, so the wall is re-expressed in canonical normal-mode coordinates — sharpened, not weakened.

*Refs.* Internal: PROOF-SAVVIDY-EXACT-SIGN (`−D²+2adF`, `λ_min=−g|B|`, `N_neg`); PROOF-PHYSICAL-SUBSPACE-CONVEXITY (orbit-transverse, 7.1% neg after `Q⊥`, min −0.82β); D-BBS-INTEGRATE-NEG (§a negatives, §b one-loop complex `i^{N_neg}`, Lemma NS `Σ_fast⪰+g|B|`); BLIND-LOWER §1/§6 (Temple circularity; (H)); PROOF-COARSE-POINCARE (exact `G″`, H1); TRANSPORT-MACHINE §4 (`(2π/L)²→0`); MS-RG-ODE §4 (`m₁` RG-invariant). External: Williamson, *Amer. J. Math.* **58** (1936) 141; N. N. Bogoliubov, *Nuovo Cim.* **7** (1958) 794; Valatin, *Nuovo Cim.* **7** (1958) 843; Simon–Høegh-Krohn–Wehlen (Gaussian/symplectic invariance); Savvidy, *PLB* **71** (1977) 133; Nielsen–Olesen, *NPB* **144** (1978) 376.
