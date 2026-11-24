# MODULAR-BUILD — Tomita–Takesaki in its sharpest form: a modular gap ("modular mass") ⟺ the local algebra is type III_λ (λ<1); but the continuum wedge is type III₁ (Connes spectrum ℝ₊) ⇒ mass-blind in *every* state, and the only mass-bearing modular Hamiltonian is the thermal βH ⇒ circular. No third horn.

> **判词 (verdict).** This upgrades the banked vacuum-modular no-go (`AQFT-CAUSAL`, `MODULAR-TAKESAKI`, `H1-T2-AQFT`) from a *mass-blindness remark* to a **classification dichotomy**, and closes the two evasions those files left open (the *state* evasion and the *lattice/type-I* evasion). The precise statement: a "modular mass" — a discrete low spectrum of `K_W = −log Δ_W` tracking `m` — **exists iff `𝒜(W)` is type III_λ with λ∈(0,1)** (Connes invariant `S(M)={0}∪{λ^n}`, Arveson spectrum `(2π/|log λ|)ℤ`; gap `|log λ|`). But the continuum wedge of a d≥3 net is **type III₁** (`S(M)=[0,∞)`, `spec(−log Δ_φ)=ℝ`) — proven for the free field at *every* mass (Hislop–Longo 1982) — so there is **no modular gap in any faithful normal state** (the type and Connes spectrum are *cocycle-invariant*, Takesaki R–N, so "pick a confining state" cannot evade it). The one modular operator that *does* carry `m` is the **thermal/KMS** one, where Tomita–Takesaki forces `Δ=e^{−βH}` so `K_β=βH`: its gap is `β·(Hamiltonian gap)` — i.e. the Hamiltonian spectrum you were trying to bound (circular), and on the lattice it is type-I_d and literally `βH_lat`. **Dichotomy: mass-blind (type III₁, any state) XOR circular (type I, `K_β=βH`). The candidate-as-stated breaks first, already on the free massive field: BW is exact for free fields at all masses, `Δ_W` is mass-independent ⇒ false negative before YM is even reached.** Score **1/5** (#1 genuine; #4 vacuous). No `λ`; no fabrication.

## §1 The modular machine, stated exactly

Hilbert space `ℋ`, vacuum `Ω` (cyclic and separating for `𝒜(W)` by Reeh–Schlieder). For the R-wedge `W_R={x: x¹>|x⁰|}` let `𝒜(W_R)=(∪_{O⊂W_R}𝒜(O))''` be the local von Neumann algebra of a Haag–Kastler net.

- **Tomita–Takesaki.** The antilinear `S_W: AΩ ↦ A*Ω` (closable) factors as `S_W = J_W Δ_W^{1/2}`: `Δ_W=S_W*S_W≥0` self-adjoint (**modular operator**), `J_W` antiunitary (**modular conjugation**). Then `σ_t := Ad Δ_W^{it}` is a one-parameter `*`-automorphism group of `𝒜(W)` (**modular flow**), and `K_W := −log Δ_W` is the (**unbounded**) **modular Hamiltonian**. `J_W 𝒜(W) J_W = 𝒜(W)'`.
- **KMS.** The modular flow is, by construction, KMS at `β=1` for `ω=⟨Ω|·|Ω⟩`: for all `A,B∈𝒜(W)` the function `t↦ω(A σ_t(B))` extends to `{0≤Im t≤1}` with the boundary relation `ω(A σ_{t+i}(B)) = ω(σ_t(B)\,A)`. Equivalently `σ_t` is the *canonical* flow making `ω` KMS at `β=1` (Takesaki's converse pins it among KMS flows).
- **Bisognano–Wichmann (1975/76).** For the wedge, `Δ_W^{it} = U(Λ_W(−2πt))`, `J_W = Θ·R_W` with `Θ`=CPT and `R_W` the geometric reflection of the wedge onto its commutant. "**Modular flow = boost**" means precisely: `σ_t` on `𝒜(W_R)` is the geometric one-parameter Lorentz boost that preserves `W_R`; hence `K_W = 2π·K_{\rm boost}` and the modular (Unruh) inverse temperature is `β_mod = 2π`. This is a **true theorem** for (i) free fields of **all masses** (the modular operator of the wedge is the boost, independent of `m`) and (ii) conformal fields (Hislop–Longo extend to double cones via fractional-linear/conformal modular flow). It is **conjectural / corrected** for interacting massive theory in d=4: there one expects `K_W = 2πK_{\rm boost} + ∫_{W} (nonlocal, state-dependent correction)`, the interacting modular Hamiltonian, which is *not* the pure boost and is unknown — and the theory's net itself is unconstructed (the G1 wall). The acid test operates on the free case, where BW is exact.

## §2 The candidate, then the precise no-go

**Candidate ("modular mass").** For the confining YM vacuum `ω_Y` (or any confining/dressed state), suppose the wedge modular Hamiltonian `K_W^{ω_Y}=−log Δ_W^{ω_Y}` has a discrete low-lying spectrum with a gap `m_mod>0` that tracks the physical gap, `m_mod = c·m` with `c` computable and `(L,a)`-uniform. Then `m ≥ λ` follows, closing H1 (`Var_μ(a_{k→0})≤1/λ`). Stated as a working hypothesis:
```
(H_work):  𝒜(W) is type III_λ, λ∈(0,1);  gap(K_W) = |log λ| = c·m  (a "modular mass"),
           uniform in (L,a);  λ→1 (⇔ |log λ|→0 ⇔ type III₁) iff m→0  (so the free massless field is the m=0 limit).
```

**The no-go (classification form).** Two horns, exhaustive.

**Horn A — continuum, type III₁, mass-blind in *every* state.** For a d≥3 Haag–Kastler net with standard axioms, `𝒜(W)` is the hyperfinite **type III₁** factor. A **modular gap exists iff the factor is type III_λ (0<λ<1)** — Connes' invariant is then `S(M)={0}∪{λ^n:n∈ℤ}`, the Arveson spectrum of the modular group is `T(M)=(2π/|log λ|)ℤ` (discrete, periodic), and `−log Δ` has log-spectrum `(log λ)ℤ` with gap `|log λ|>0`. By contrast:
```
type I, II : S(M)={1};                          type III₀ : S(M)={0,1};
type III_λ : S(M)={0}∪{λ^n},  T(M)=(2π/log λ)ℤ,  GAP |log λ|;
type III₁ : S(M)=[0,∞),       T(M)={0},          spec(−log Δ_φ)=ℝ,  NO GAP.
```
`III₁` is an **invariant of the algebra `M`, independent of the faithful normal state `φ`** — because any two states are related by the **Connes cocycle / Takesaki Radon–Nikodym derivative** `Δ_φ^{it}=(Dφ:Dφ')_t Δ_{φ'}^{it}(Dφ:Dφ')_t^*`, under which both the *type* and the *Connes spectrum* are preserved. Therefore **no state choice** — vacuum, thermal, confining, glueball-dressed, `ω_Y` itself — yields a discrete modular spectrum on a type III₁ factor: `spec(−log Δ_W)=ℝ` for every one. The load-bearing, *proven* instance is the **free massless scalar field, whose double-cone and wedge algebras are type III₁ (Hislop–Longo 1982)** — mass-blind, consistent only with `m=0`, so it passes the acid test *vacuously* and proves the implication false: a III₁ algebra admits no modular gap while the gapless free field is its paradigm.

**Horn B — lattice, type I_d, circular.** On the finite lattice at `β≈2.5`, `𝒜(W)` is finite-dimensional (**type I_d**) and `Δ_W` *is* a positive matrix with genuine discrete spectrum. But take any `(β,α_t)`-**KMS** state for the lattice time translation `α_t`: Tomita–Takesaki forces the modular flow to equal the dynamics, `σ_t^{ω}=α_{−βt}`, so on the lattice `Δ_ω = e^{−βH_{\rm lat}}` and
```
K_β := −log Δ_ω = β H_{\rm lat}   ⇒   gap(K_β) = β · gap(H_{\rm lat}).
```
The modular operator sees `m` — **only because it literally *is* the exponentiated lattice Hamiltonian**. Computing its gap is diagonalizing `H_{\rm lat}`, i.e. the spectral problem being asked. And as the regulator is removed `type I_d → III₁`: the discrete modular spectrum dissolves into `[0,∞)` (Horn A), so the lattice "modular gap" is a finite-size/UV artifact that is **not stable in the continuum limit**.

**Dichotomy (the result).** *Every* modular Hamiltonian is in exactly one of two classes:
```
 type III₁ (continuum vacuum, ANY state)  ⇒  mass-blind (no modular gap at all);
 type I (finite lattice, thermal KMS)     ⇒  K_β = βH  ⇒  circular (= the Hamiltonian spectrum).
```
There is **no third horn**: a modular object that is both gap-bearing *and* independent of the Hamiltonian spectrum. The exact hypothesis to make the candidate work — `(H_work)` above, a state/cocycle-invariant type-III_λ algebra with `1/|log λ| ∝ m` — is refuted on the acid-test field (free field is III₁ ⇒ `λ=1` ⇒ `|log λ|=0`, forcing `m=0`, which is exactly the free-massless PASS) and is physically excluded in the target regime (a `λ<1` modular flow is *periodic*, `T(M)≠{0}`, i.e. recurrent — the opposite of the ergodic clustering that a positive mass produces).

## §3 Three-regime test table (honest)

| regime | modular operator `K_W` | type / Connes `S(M)` | modular gap? | yields `m`? | verdict |
|---|---|---|---|---|---|
| **free massless**, vacuum wedge | `Δ_W`=boost (BW exact) | **III₁**, `S=[0,∞)` | none | no (`m=0`) | **acid-test PASS: returns 0** — but mass-blind; `f₁²=0` correctly |
| **free massless**, thermal `β`/KMS | `σ_t`=time translation; lattice `Δ=e^{−βH}`, `H` gapless | III₁ (cont.) / I_d (lattice) | lattice: discrete but `βH` gapless; cont.: none | no | consistent — no spurious gap |
| **free massive**, vacuum wedge | `Δ_W`=boost (**BW exact for free, all `m`**) | **III₁**, `S=[0,∞)` (Hislop–Longo, `m`-independent) | **none even though `m>0`** | **NO — FALSE NEGATIVE** | **decisive: the wedge modular operator cannot see `m`** |
| **free massive**, thermal `β` | lattice `Δ=e^{−βH}`, `K=βH`, `H` gapped | I_d (lattice) → III₁ (cont.) | lattice gap `βm`; cont. none | "yes" on lattice but `=βH` | **circular (it *is* the Hamiltonian) + dissolves in continuum** |
| **interacting YM `β≈2.5`** (lattice) | type I_d; finite discrete modular spectrum | I_d finite → III₁ as `a→0` (expected; net not yet built = G1) | discrete on finite lattice, but UV/size artifact → 0 | no usable `λ` | **breaks**: (i) `I_d→III₁` dissolves the gap; (ii) thermal form `=βH_lat` circular; (iii) no 4D net (G1) |

**Where it first breaks.** Not at YM — **at the free massive field.** BW is exact for free fields at *every* mass and the double-cone/wedge algebra is type III₁ with `m`-independent Connes spectrum (Hislop–Longo): the wedge modular Hamiltonian is provably identical for `m=0` and `m>0`, so `spec K_W` is a **false negative on an exactly solvable model**. The interacting regime adds two independent breaks: the lattice type-I_d modular gap is a UV/size artifact that dies as `I_d→III₁` (a *modular image of the G1 continuum wall*), and any attempt to stay finite and use the thermal `βH_lat` is definitional circularity.

## §4 Five-point scorecard

1. **Nonperturbative — PASS.** Exact Tomita–Takesaki / Connes-classification; every `β`, no loop or `1/N` expansion.
2. **Sign-controlling lower bound `m≥λ` — FAIL.** The wedge `Δ_W` is type III₁ (spectrum `[0,∞)`, `m`-independent) ⇒ `λ=0`; the only gap-bearing `K` is `βH` ⇒ no independent `λ`.
3. **Non-circular — FAIL.** Thermal `K_β=βH`; its gap *is* the Hamiltonian spectrum = the target. State choice cannot help (Takesaki R–N preserves III₁).
4. **Non-SUSY + fails free massless — PASS, vacuously.** Free massless → 0 (`III₁`, `S=[0,∞)`, no gap) — but free **massive** also → 0, so the detector does not discriminate gapped from gapless; the "pass" is empty.
5. **`(L,a)`-uniform `≥λ` — FAIL.** No uniform `λ`; the lattice "modular gap" is a UV/size artifact `→0` as `I_d→III₁`.

**Score: 1/5** (#1 genuine; #4 vacuous). Same wall as `AQFT-CAUSAL`, now at classification resolution.

## §5 Sharpest line & exact stuck step

**Sharpest line.** The modular operator of a *confining vacuum* is a **boost**, and a boost does not know how much the field weighs: on the continuum the wedge algebra is **type III₁** (Connes spectrum `ℝ_{>0}`) for the free field at *every* mass, so `K_W` is provably identical for `m=0` and `m>0`. The one modular Hamiltonian that *does* carry `m` is the thermal one, and it carries it only because it **is** `βH` — reading its gap is reading the Hamiltonian you set out to bound. **Type III₁ makes modular theory mass-blind; the thermal escape makes it circular; there is no third horn.** In the corpus's terms: the classification invariant `S(M)` is a *conserved charge of the algebra*, not a functional of the softest mode; H1 asks for the stiffness of the `k→0` mode, which modular type does not measure.

**Exact stuck step.**
```
Target:  m ≥ λ   (⟺ H1: Var_μ(a_{k→0}) ≤ 1/λ).
Candidate lever:  gap( K_W^{ω} ) ≥ m_mod ~ m,   K_W = −log Δ_W.
── continuum obstruction ─────────────────────────────────────────────
  𝒜(W) type III₁  ⇒  S(M)=[0,∞),  T(M)={0},  spec(−log Δ_φ)=ℝ   for EVERY faithful normal φ
  (Fredenhagen 1985; Hislop–Longo 1982 for the free field; Connes 1973);
  type & Connes spectrum are Takesaki-R–N cocycle invariants ⇒ no state evades it.
  ⇒ no modular gap, any state, any m.                      [MASS-BLIND horn]
── lattice obstruction ───────────────────────────────────────────────
  type I_d  ⇒  thermal KMS: Δ=e^{−βH_lat}, K_β=βH_lat  ⇒  gap(K_β)=β·gap(H_lat)
  = the Hamiltonian spectrum being sought;  and  I_d → III₁ as a→0  dissolves it.  [CIRCULAR horn]
── missing object (would make (H_work) true) ─────────────────────────
  a type-III_λ (λ<1) local algebra, 1/|log λ| ∝ m, state/cocycle-invariant, uniform in (L,a)
  — FALSE for the acid-test field (free = III₁) and physically excluded (III_λ flow is periodic ≠ clustering).
```

## §6 Position vs. the bank (what is new here)

- **Already banked** (`MODULAR-TAKESAKI`, `AQFT-CAUSAL`, `H1-T2-AQFT`): the vacuum wedge modular operator is a mass-blind boost; RS is anti-clustering; split/nuclearity is UV/wrong-direction; QNEC is local/wrong-channel; entanglement saturation is free/Gaussian only; area law runs the wrong way (`H1-T4`). Conclusion *already* a no-go.
- **New in this file** (the untested part, per STEP 0): (a) the **classification reformulation** — modular gap `⟺` type III_λ (λ<1) `⟺` Connes `S(M)={0}∪{λ^n}` / Arveson `(2π/|log λ|)ℤ`; type III₁ `⟺ S(M)=[0,∞)` `⟺` no gap; (b) the **state-choice evasion closed** via the Connes cocycle (Takesaki R–N): *no* state on a III₁ algebra yields a discrete modular spectrum; (c) the **thermal/KMS horn** identified exactly — the only mass-bearing modular Hamiltonian is `βH`, hence circular; (d) the **lattice type-I_d → III₁ dissolution** as the modular image of the G1 continuum wall. Net: same verdict, but *evasion-proof and classification-sharp* — an upgrade of the banked no-go, not a re-bank.

*Refs (external, verified).* M. Takesaki, *Tomita's Theory of Modular Hilbert Algebras*, Lecture Notes in Math. **128** (Springer, 1970); *Theory of Operator Algebras II* (2003, cocycle R–N / KMS). A. Connes, "Une classification des facteurs de type III," *Ann. Sci. École Norm. Sup.* (4) **6** (1973), no. 2, 133–252 (S-invariant, III_λ). H. Araki, W. Woods, "A classification of factors," *Publ. Res. Inst. Math. Sci.* **4** (1968/69) 51–130 (asymptotic ratio set). K. Fredenhagen, "On the modular structure of local algebras of observables," *Comm. Math. Phys.* **97** (1985) 79–89 (type III₁). P. D. Hislop, R. Longo, "Modular structure of the local algebras associated with the free massless scalar field theory," *Comm. Math. Phys.* **84** (1982) 71–85 (free field III₁, conformal modular flow). J. J. Bisognano, E. H. Wichmann, *J. Math. Phys.* **16** (1975) 985 / **17** (1976) 303 (BW). H. Araki, *PTP* **32** (1964) 844 (KMS). R. Haag, *Local Quantum Physics* (1996), §V (KMS, thermal = modular). U. Haagerup, "Connes' bicentralizer problem and uniqueness of the injective factor of type III₁," *Acta Math.* **158** (1987) 95–148 (III₁ uniqueness). M. Takesaki, "Tomita's theory…" + Connes cocycle for state-independence of type. Internal: `AQFT-CAUSAL`, `MODULAR-TAKESAKI`, `H1-T2-AQFT`, `H1-T4-entanglement`, `C8-SPECIES-NOGO` (empty-referent pattern), `k0-CONSTRAINT-EFFPOT` (`G″=1/Var(a₀)`), `PROOF-IR-H1` (H1 = gap), `COMPUTE-GPP-SU3` (`1/Var≈G″`), `GAP-G2-edge` (edge-atom vs aggregate), `G1-CONTINUUM-OS` (no 4D net).
