# ANALYTIC MICROLOCAL (WF_A) — the analytic wavefront *set* is mass-blind (identical light cone for massless and massive); the mass lives only in the FBI decay *rate* `Λ = m₁`, which by Paley–Wiener **is** the KL edge; and propagation-of-singularities runs the wrong way (it bounds the singular set from *above*, the gap needs a singularity-free *strip*). **Reduction to the KL edge; no lever.**

> 2026-11-23 · NEW SPECIES, distinct from the C∞/subelliptic species in `HORMANDER.md`. That file analyzed the *principal symbol* / Fefferman–Phong subelliptic estimate (`1/r` derivative gain at `k→∞`); this file analyzes the **analytic** category — Sato–Kawai–Kashiwara microfunctions, Sjöstrand's `WF_A`, the Wightman tube, Jost–Lehmann–Dyson, edge-of-the-wedge. The two share the verdict "principal-symbol/sign-blind regularity data," but the object and the exact failure differ, and the present species yields two *sharper* no-gos (set-blindness; singular-set monotonicity) plus one clean exact dictionary (`w = m₁`). **判词 (verdict). The reduction is exact and nonperturbative, and it is a relabeling: "distance to the nearest analytic singularity" `=` "tube width" `=` "FBI decay rate" `=` the KL mass `m₁` `=` H1. There is no microlocal lever, because (1) the `WF_A` *set* is literally identical for the gapless free field and the gapped one — it carries no mass information at all; (2) propagation of singularities is a lower-bound-on-*regularity* theorem, while a gap is an upper-bound-on-*singular-set* statement, and nonlinear propagation can only *add* singularities. The atom content (pole vs cut) is already banked as the Fredholm determinant; the mass content is already banked as H1; nothing new closes.**

---

## §1 Setup — the Wightman tube, `WF_A`, and the KL / decay / width dictionary

Take the gauge-invariant Hermitian local field `O = ¼ tr F_{μν}F^{μν}` (the scalar density of the corpus) or `s = (β/2g)\,tr F²`; the choice is immaterial — the operator must be local, gauge-invariant, and have a nontrivial `0⁺⁺` projection. Let

```
W_O(ζ) = ⟨0| O(x) O(y) |0⟩ ,   ζ = x − y ∈ ℝ⁴ ,
```

the Wightman two-point function. Three standard facts fix the analytic landscape.

**(F1) Tube theorem / boundary value.** `W_O(ζ)` is the boundary value, as `Im ζ → 0` from inside `−V_+` (`V_+` = open forward cone), of a function `W_O(z)` holomorphic in the forward tube `T = ℝ⁴ − iV_+`. (Streater–Wightman; Reed–Simon II §IX.) Equivalently, in momentum space `W̃_O(p)` is the boundary value of a function analytic for `Im p ∈ V_+`.

**(F2) Källén–Lehmann.** Reflection positivity gives `W̃_O(p) = ∫₀^∞ dν(s)/(p² − s + i0)` with `ν ≥ 0` a positive measure, `supp ν ⊂ [m₁², ∞)`, `m₁ = inf supp ν`. The **atom** is `f₁² = ν({m₁²}) = |⟨0|O|G⟩|²`; `f₁² > 0` means an isolated pole, `f₁² = 0` means a pure cut from possibly `m₁ = 0`. The **mass** is `m₁`; the corpus's `H1` is `Var_μ(a_{k→0}) = χ = ∫ s⁻¹ dν ≤ 1/λ`, and the squeeze `f₁²/χ ≤ m₁² ≤ M₀/χ` ties `m₁` to `(f₁², χ)`.

**(F3) Paley–Wiener / decay duality.** An analytic continuation of `W_O` into a tube of half-width `w` in the difference variable is equivalent to exponential decay `|W_O(ζ)| ≤ C e^{−w|ζ|}` at large `|ζ|`. Sending `W̃_O(p)`'s analyticity: `W̃_O` is analytic in `p` off the real cut `[m₁², ∞)`, so the distance from the physical point `p² = 0` to the nearest singularity (the branch point `p² = m₁²`) is `m₁²`, and in the `|p|` metric the "tube width" is `m₁`.

Gathering (F1)–(F3):

```
w  :=  sup{ half-width of the analyticity tube of W_O in ζ }        (x-space)
    =  dist(|p|, nearest singularity of W̃_O above p²=0)             (p-space)
    =  − lim_{|ζ|→∞} (1/|ζ|) log|W_O(ζ)|                            (decay rate)
    =  m₁  =  inf supp ν .                                          (KL)
```

All four readings are the **same number**, and they are `m₁` *by the KL representation itself*. This identity is exact, nonperturbative, and the entire positive content of the species. Call it **(W)**.

**(F4) Analytic wavefront set.** `(x₀, ξ₀) ∉ WF_A(W_O)` iff `W_O` is real-analytic in a neighborhood of `x₀` (analytic `WF_A` is empty in the fiber over any real-analytic point; Sjöstrand's FBI definition: `|T^λ W_O(x,ξ)| ≤ C e^{−λN}` for all `N` in a neighborhood). The *complement of the projection* of `WF_A` to `x`-space is the **analytic singular support** `sing_A W_O`.

**(F5) The microlocal spectrum condition (precise form of the free-field statement).** For **any** Hadamard two-point function — free massive, free massless, and (BFK) the interacting case as a *condition* — the `C∞` wavefront set is exactly

```
WF(W_O) = { (x,ξ; y,−η) : (x,ξ) ~ (y,η) on a null bicharacteristic, ξ future-directed null } ,
```

i.e. **the light-cone bicharacteristic set, with no mass parameter in it** (Radzikowski; Brunetti–Fredenhagen–Köhler). The analytic version has `sing_A W_O ⊂ {ζ : ζ² = 0}` on the light cone. **This is the first decisive fact: the singular set is mass-independent.** It is identical for `m = 0` and `m > 0`.

**What a propagation-of-singularities theorem would need to give a LOWER bound on `m`.** A theorem of the form `m₁ ≥ m₀ > 0` via propagation would need three ingredients, and each is unavailable:

- **(i) A seed.** An open region where `WF_A(W_O) ≠ ∅` at a momentum shell `p² = μ² ≥ m₀²` of *known* location — a known massive state. The only structural seed is the canonical light cone `p² = 0`, which is mass-blind (F5) and sits at `w = 0`.
- **(ii) A flow that preserves separation.** The Hamiltonian flow `H_{□}` of the wave operator moves singularities along *null* bicharacteristics at *fixed* `p²`; it preserves the mass shell and cannot *create* a gap between a seed and the physical point.
- **(iii) A regularity conclusion.** That `WF_A` is *empty* in the strip `{0 < p² < m₀²}`. But propagation gives `WF_A(W_O) ⊂ WF_A(P W_O) ∪ Char(P)`: it bounds the singular set from **above** (it is a lower-bound-on-*regularity*/upper-bound-on-*singular-set* statement only relative to `Char P`), and `Char(□)` = the massless light cone is *always* present. Emptiness in a strip is a **regularity** input, never a propagation output.

So a propagation theorem yielding `m₁ ≥ m₀` must already *be given* an isolated massive singularity at `≥ m₀` — i.e. `f₁² > 0`. It is circular at (i), inert at (ii), and direction-reversed at (iii).

---

## §2 The strongest concrete candidate — stated precisely, then adjudicated

I build the strongest candidate the analytic category admits, state it with explicit hypotheses, and then exhibit exactly why it is a **reduction, not a lever**.

**Candidate M1 (FBI decay-rate lever / singularity-free strip).**
*Hypotheses:* Wightman axioms (spectral condition, positivity `ν ≥ 0`, locality, Hermiticity), `O` local gauge-invariant. *Object:* the Sjöstrand FBI exponential rate
```
Λ(ξ) := − lim_{λ→∞} (1/λ) log |T_λ W_O(0, ξ)| ,     T_λ the FBI transform .
```
*Claim to be built:* a lower bound `Λ(ξ) ≥ m₀ > 0` uniformly for `|ξ| = 1` (equivalently: `W_O` is analytic in the tube `{|Im ζ| < m₀}`, a **singularity-free strip** of width `m₀` around the physical domain), obtained *from the microlocal axioms alone*, uniform in `(L,a)`.

**Adjudication — it is the KL edge (a reduction), by the following exact chain.**

1. By Paley–Wiener (F3), `Λ(ξ) ≥ m₀` ⟺ tube width `w ≥ m₀` ⟺ decay rate `≥ m₀`.
2. By KL (F2/(W)), `w = m₁ = inf supp ν`. So `Λ ≥ m₀` ⟺ `m₁ ≥ m₀` ⟺ **the mass gap** (equivalently `H1` with the residue floor, since `m₁² ≥ f₁²/χ` and `Λ = m₁`).
3. **What the axioms supply:** the spectral condition gives `supp ν ⊂ [0,∞)` hence `Λ ≥ 0`; positivity `ν ≥ 0` gives nothing more on the *bottom* of `supp ν` (the free massless `ν = δ₀` is positive and has `m₁ = 0`). The axioms are **consistent with `Λ = 0`** (a gapless interacting CFT is a legitimate Wightman theory). So `Λ ≥ m₀ > 0` is *not implied* by the microlocal axioms; it *is* the excluded alternative `¬(iii)` of the reduction.

Hence **M1 reduces exactly to the KL edge `= H1`**: the "lower bound on the singularity distance" **is** the gap. This is not a defect of a particular candidate; §3 shows it is structural, and §4 (Condition 4) shows the set-level object cannot distinguish the two regimes at all.

**The one "positive" microlocal reading, and why it is also banked.** The *strength* of an isolated singularity is the residue `f₁²`; "`f₁² > 0`" in microlocal language is "`WF_A(W_O)` has an **isolated** singularity — a pole, not only a cut — at `p² = m₁²`." The acid test is clean here: the free field's `ν` is a pure cut, no pole, `f₁² = 0`. **But this is exactly the already-banked Fredholm determinant identity** `f₁² = −λ₁ Res G(z)` (`INVENT-NEW-IDENTITY`; `docs/new-tools-registry.md` A3; refined in `CONTINUUM-FREDHOLM.md`, which shows the determinant zero *coalesces into the KL pole* at `a→0`). The pole-vs-cut / atom-detection half of the microlocal story is therefore **not new** — it is the Fredholm zero-set, restated. What would be new is the *distance* (mass) half, and that is the KL edge.

So the species splits cleanly:
- **atom half** (pole vs cut, `f₁² > 0`) = Fredholm determinant — **already banked**;
- **mass half** (singularity distance / decay rate / `w = m₁`) = KL spectral edge = `H1` — **already banked**;
- **the gap between "has an isolated singularity somewhere" and "the nearest singularity is at distance `m₀ > 0`"** = exactly `H1`, and microlocal analysis supplies no independent handle on that distance.

---

## §3 Three-regime test (honest)

| Regime | `ν` | `WF_A` structure (x-space set) | tube width `w` | what the tool returns | verdict |
|---|---|---|---|---|---|
| **FREE MASSLESS** (`O = :φ²:` or `F²` on Maxwell) | `δ₀` (cut from 0) | canonical light cone, reaching `ζ = 0`; `sing_A = {ζ²=0}` | `0` | `w = 0`, `Λ = 0` — **correct**, gapless, no spurious gap | **acid-NEUTRAL as an identity** (`w=m₁` holds trivially), **acid-FAIL as a lever** (free-field-blind: the set is *identical* to the gapped case, see below) |
| **FREE MASSIVE** (`ν = δ_{m²}`) | `δ_{m²}` | **the same light cone** `sing_A = {ζ²=0}` — a pole appears only in *momentum* space at `p²=m²`, i.e. in the *rate* | `m` | `w = m` — **correct**, but read off the pole `δ_{m²}` that was handed in | **CIRCULAR**: `w = m` is `dist(p², pole)`; recovering `m` from `w` just reads the input spectral function. (Decisive: the x-space `WF_A` *set* is literally the same as massless — `K₁(m√{−ζ²})` is analytic off `ζ²=0`; the mass lives only in `e^{−m|ζ|}`, the rate, not the set.) |
| **INTERACTING YM** (`β≈2.5`) | `ν`, `supp ⊂ [m₁²,∞)`, atom `f₁² > 0` | light cone **+** pole at `p²=m₁²` (the atom) **+** cut above (continuum from `4m₁²`/excited) | `m₁` | needs `ν`: it gives `Λ ≥ 0` only; locating `m₁` requires knowing the bottom of `supp ν` | **BREAKS here**: to get `m₁ > 0` one must exclude the strip `(0, m₁²)`, a *regularity* statement the tool assumes, never proves. The microlocal axioms are *consistent with `m₁ = 0`* (gapless CFT). The sign of `w` is the gap. |

**The decisive row is FREE MASSIVE.** Its x-space `WF_A` *set* is the light cone — **identical to the free massless set** — yet `m = m* ≠ 0`. Therefore the `WF_A` *set* is **mass-blind**: two theories with `m = 0` and `m = m*` have the *same* analytic singular set. The mass is not a point of `WF_A`; it is the *exponential rate* `Λ`, and the rate is `m₁` by Paley–Wiener = KL. This is why the species cannot be a sign-controlling lever: the very object it computes (the set) does not carry the number that matters, and the object that carries it (the rate) is, by an identity, the answer.

---

## §4 Five-condition scorecard

| # | Condition | Score | Reason |
|---|---|---|---|
| 1 | **Nonperturbative** | **PASS** | (F1)–(F5)/(W) are exact, no loop expansion; even the failure is nonperturbative (it is an identity + a direction obstruction, not a weak-coupling artifact). |
| 2 | **Sign-controlling LOWER bound** | **FAIL** | `WF_A` is a *binary set*; propagation yields `WF ⊂ WF(source) ∪ Char` — an **upper bound on the singular set / lower bound on regularity**. The gap is the *opposite* (singularity-free strip = upper bound on singular set near the physical point). Nonlinear propagation can only *add* singularities (Bony: `WF_A(uv) ⊂ WF_A(u) ∪ WF_A(v) ∪ {ξ=ξ₁+ξ₂}`). No lower bound on `m` is produced; (W) is an identity, not an inequality. |
| 3 | **Non-circular** | **FAIL** | `w = m₁ = inf supp ν` **by definition** (KL). Reading `w` reads the spectrum; a "lower bound on the singularity distance" **is** the gap. Jost–Lehmann–Dyson and the edge-of-the-wedge/BEG analytic completion all *take the absorptive part = `ρ`* as input; the FBI threshold *is* the rate. M1's premise `Λ ≥ m₀` is the conclusion `m₁ ≥ m₀`. |
| 4 | **Non-SUSY + fails free massless** | **FAIL (max)** | Non-SUSY yes. But the free *massless* and free *massive* and interacting Hadamard 2-pt functions share the **same** `WF` set (F5) — the canonical light cone, mass-independent. The tool is not merely free-field-blind; the set is *literally invariant* under changing the mass, so it is structurally incapable of carrying the gap. (A tool that two theories with and without a gap cannot separate has "failed" the acid test in the strongest sense: it returns the same value on both.) |
| 5 | **Uniform `(L,a)`** | **FAIL** | `w = m₁` *inherits* the uniformity question with no independent source. On `T⁴_L` the free value `w = 2π/L → 0`; on `ℝ⁴` the gapless free field has `w = 0`. The uniform positivity of `w` is exactly the open `H1`; microlocal analysis supplies no uniformity lever (and `Char(□)` is scale-free, giving `w=0` at every `L,a`). |

**Score: 1 / 5.** The single point is the exact nonperturbative dictionary (W) and the clean, acid-consistent *reduction*. Conditions 2–5 — every condition that could turn a reduction into a proof of the gap — fail. This matches and sharpens the C∞-species verdict in `HORMANDER.md` (which scored 1/5 on the analogous subelliptic statement) while contributing two facts that file did not state: **set-blindness across the mass** (§3, FREE MASSIVE row) and **singular-set monotonicity / direction-reversal of propagation** (§1 (iii), §4 #2).

---

## Exact stuck estimate / where it breaks

```
Wightman tube (F1):  W_O(z) holomorphic in T = ℝ⁴ − iV_+ ;  W̃_O(p) = limbdy ∫ dν(s)/(p²−s+i0).
Paley–Wiener (F3):   analytic tube half-width w  ⟺  |W_O(ζ)| ≤ C e^{−w|ζ|}.
KL (F2):             w = m₁ = inf supp ν.                                   (W) — the identity

The "lower bound on the singularity distance":
   Λ(ξ) := −lim_λ (1/λ) log|T_λ W_O(0,ξ)|  ≥ m₀
   ⟺ (PW)  tube width w ≥ m₀  ⟺ (KL)  m₁ ≥ m₀   ←—— this IS the mass gap / H1.

Axioms supply only:   spectral condition ⟹ Λ ≥ 0 ;  positivity ⟹ nothing on inf supp ν.
Gapless interacting CFT:  ν ⊂ [0,∞) continuous at 0  ⟹  Λ = 0  — a legitimate Wightman theory.
                       The axioms are CONSISTENT with Λ = 0 ;  they cannot imply Λ ≥ m₀ > 0.

MASS-BLINDNESS (decisive, §3):
   free massless:  sing_A W = {ζ²=0} ;  w = 0.
   free massive:   sing_A W = {ζ²=0} ;  w = m.        ←—— SAME set, different w.
   ⟹ the WF_A *set* carries no mass information; the mass is the *rate*, = KL = the answer.

DIRECTION (the wall):
   propagation:  WF_A(W) ⊂ WF_A(PW) ∪ Char(P)   (upper bound on the singular set,
                                                   = lower bound on regularity).
   gap needs:    WF_A(W) ∩ {0 < p² < m₀²} = ∅    (upper bound on the singular set
                                                   near the physical point).
   nonlinearity: Bony — WF_A(uv) can only ADD singularities to WF_A(u) ∪ WF_A(v).
   the only structural singularity is Char(□) = the massless light cone, at w = 0, mass-blind.
   ⟹ a propagation theorem giving m₁ ≥ m₀ needs an isolated massive seed at ≥ m₀ as INPUT
      = f₁² > 0  (circular); inert in the flow; direction-reversed at the conclusion.
```

**Sharpest line.** Analytic microlocal analysis hands you the one singularity that is **mass-blind** — the canonical light-cone `Char(□)`, present at `w = 0`, *identical* for the gapless free field and the gapped one, so the `WF_A` *set* literally does not know whether there is a mass — and it is silent on the one object that **carries** the mass, the exponential *rate* `Λ`, which by Paley–Wiener *is* the KL edge `m₁ = inf supp ν`. The "lower bound on the singularity distance" a propagation theorem could offer is therefore either the light cone (`= 0`, gapless, and free-field-fitting) or the KL edge (`= m₁`, circular) — there is no third option. The atom content (pole vs cut) is already the banked Fredholm determinant; the distance content (mass) is already `H1`; **analytic microlocal analysis re-derives the KL spectral edge, and the lower bound on the singularity distance is the gap itself.**

**Exact missing step (if one insists on a lever).** A theorem `m₁ ≥ m₀ > 0` in this species must supply a **singularity-free strip** `{0 < p² < m₀²}` from structural axioms — an *upper bound on the singular set* — whereas the entire analytic-microlocal toolkit produces *lower bounds on the singular set* (propagation) and *identities* (KL/Paley–Wiener). Bridging that direction gap is `H1`; it is not a microlocal statement, and no propagation, FBI, JLD, or edge-of-the-wedge theorem closes it, because each presupposes the spectral bottom as input.

---

*Refs.* L. Hörmander, *The Analysis of Linear Partial Differential Operators* I (Springer, 1983/2003; §7 Paley–Wiener, §8 wave-front sets/FBI) and IV (1985; Ch. XXV–XXVI analytic wave-front sets, propagation of singularities Thm 26.1.1); M. Sato, T. Kawai, M. Kashiwara, *Microfunctions and pseudo-differential equations*, Lecture Notes in Math. **287** (1973) 265; J. Sjöstrand, *Singularités analytiques microlocales*, Astérisque **95** (1982); J.-M. Bony, *Calcul symbolique et propagation des singularités … non linéaires*, Ann. Sci. ÉNS **14** (1981) 209; M. Radzikowski, *Micro-local approach to the Hadamard condition …*, Comm. Math. Phys. **179** (1996) 529; R. Brunetti, K. Fredenhagen, M. Köhler, *The microlocal spectrum condition and Wick polynomials …*, Comm. Math. Phys. **180** (1996) 633; R. Streater, A. Wightman, *PCT, Spin and Statistics, and All That* (1964) (tube theorem, `W̃` boundary values, KL); N. Bogoliubov, B. Medvedev, M. Polivanov, *Questions in the Theory of Dispersion Relations* (1958) (edge of the wedge); R. Jost, H. Lehmann, Nuovo Cim. **5** (1957) 1598, and F. Dyson, Phys. Rev. **110** (1958) 1460 (Jost–Lehmann–Dyson representation, takes `ρ` as input); M. Reed, B. Simon, *Methods of Modern Mathematical Physics* II (1975) §IX (Wightman axioms, tube domains). Internal: `HORMANDER.md` (the C∞/subelliptic sibling species — principal symbol, `1/r`; same 1/5, complementary failure); `docs/YANG-MILLS-CONDITIONAL-REDUCTION.md` §1.2, §3 (KL measure `ν`, `f₁²`, `H1`, spectral squeeze); `INVENT-NEW-IDENTITY.md` + `CONTINUUM-FREDHOLM.md` + `docs/new-tools-registry.md` A3 (Fredholm zero-set = the pole-vs-cut / atom half, and its coalescence into the KL pole at `a→0` — i.e. the atom content of this species is already banked); `PROOF-IR-H1.md` (the `∇_⊥` vs `∇_∥` two-coordinate structure of `H1`); `REDUCTION-AUDIT.md` (`H1` necessary-not-sufficient; Schur-convex master obstruction → sign-blind regularity data).
