# INVERSE-SCATTERING-IST — the Gel'fand–Levitan inverse problem is a bijection that takes the bound-state atom as *input*; its only rigorous lower bounds (Temple) need `E₂`, its count (Levinson) is an integer frozen by the rigidity no-go, its one real theorem is 1+1D integrable, and its new content is A3's Fredholm determinant

> **2026-11-23 · NEW SPECIES — inverse scattering / Gel'fand–Levitan–Marchenko (IST/GLM).** 判词 (verdict). Forge: cast the `0⁺⁺` coarse operator as a half-line scattering problem `L=−d²/dr²+V_eff`; read the KL atom `f₁²=ν({m₁²})` as a bound state `(κ₁,c₁²)`; use the inverse problem to certify the atom / the gap `G″(k→0)≥λ ⟺ E₁(L_{0⁺⁺})≥λ`. **Verdict: exact dictionary, identity-not-inequality, circular, one counting theorem frozen by the rigidity no-go.** GLM is a *bijection* `{potentials}↔{spectral measures}`: the atom enters `F(ξ)=Σc_j²e^{−κ_jξ}+∫R e^{ikξ}` as **input**, the potential `V=−2∂_xK(x,x)` is **output** — certifying the atom runs backward (WIENER–HOPF §2, second alphabet). The only rigorous *lower* bounds on a Schrödinger ground state (Temple/Weinstein/Kato) consume a bound on `E₂`, the gap one rung up (GAP-G2-edge regress); min–max is upper-only; Levinson's count `N=(1/π)ΔargT∈ℤ` ⇒ rigidity no-go (B3) pins it to the free value `0`. The one real theorem — 1+1D integrable TBA — is dimension-locked, already banked (WH §5). New content = A3: `det(1+K_GLM)=1/T(k)`, zeros = atoms = `det(1−zT_O)`. **#1 PASS (character); #2,#3,#5 FAIL; #4 PASS (vacuous). No `λ`.**

## §1 Dictionary — atom `=` bound state, gap `=` floor on `E₁`
`L_{0⁺⁺}=−d²/dr²+ℓ(ℓ+1)/r²+V(r)` on the half-line (ℓ=0 radial). Scattering data `R(k),T(k)`, bound states `k=iκ_j` (poles of `T`), normalizations `c_j²`. The KL measure `ν=Σf_i²δ(s−m_i²)+ρ(s)` *is* the spectral measure of `L`: atom `f₁²` ↔ point mass `=` `c₁²`; `m₁²=−κ₁²`; threshold `s_th` ↔ a.c. branch point. Target `G″(k→0)≥λ ⟺ Var(a)≤1/λ ⟺ m_eff²≥λ` is exactly a **uniform floor on `E₁(L)`** — bound state separated from threshold by `≥√λ`. The dictionary is exact — and *symmetric*, which is the tell.

## §2 The atom is GLM *input*
`K(x,y)+F(x+y)+∫_x^∞K(x,s)F(s+y)ds=0`, `F⊃Σc_j²e^{−κ_jξ}+(1/2π)∫R e^{ikξ}dk`, `V=−2∂_xK(x,x)`. The bound state `(κ_j,c_j)` sits *inside* `F`; GLM reconstructs `V` **from** it. To certify the atom you must already supply `(κ₁,c₁²)` — the asymptotic glueball and its overlap, which *is* the atom. Worse, the scattering data presupposes an asymptotic Fock space: Haag–Ruelle/LSZ take the isolated shell (**the gap**) as hypothesis (WH §2). Two inversions, one circle.

## §3 No lower bound escapes `E₂` — the wall is intrinsic to spectral theory
- **Min–max.** `E₁=min_ψR(ψ)` ⇒ any fixed trial gives `R(ψ)≥E₁` = **upper** (direction wall, GAP-G2-edge).
- **Temple/Weinstein/Kato.** Each lower bound **consumes a certified bound on `E₂`** — the excited `0⁺⁺`/multi-glueball threshold, *the gap one rung up* (GAP-G2-edge "edge ≡ two-scale separation"). Provably circular.
- **Comparison.** `V≥V_ref ⇒ E₁(V)≥E₁(V_ref)` — needs a known reference well and pointwise domination; Savvidy's indefinite sign falsifies `V≥V_ref`, and the reference well's gap is the unknown.

The sole counting law, **Levinson** `N=(1/π)ΔargT`, is an **integer**: locally constant along the connected β-family ⇒ rigidity no-go (B3) freezes it to `0`. The gap is a β-non-analytic jump, not a topological index.

## §4 Acid test — letter pass, continuously vacuous
Free massless `V=0`: `R=0`, no bound state, `dρ` continuous from `0`, gap `0` — IST returns `0` (**letter PASS**). Gapped free `V=m²`: threshold `m²`, still **no bound state**, `f₁²=0` (reproduces PROOF-F1SQ) — IST sees threshold ≠ atom sharply, but reads it off the `dρ` it cannot obtain; and `κ₁(q)` is **continuous** with `κ₁(0)=0`, so no discrete jump from free to interacting: it cannot *force* a floor (F2/QFB vacuous pass).

## §5 What is genuinely true — A3, and only 1+1D integrable
`det(1+K_GLM)=1/T(k)`; zeros in `ℂ⁺` are atoms — **exactly A3's** `D_O(z)=det(1−zT_O)`, `f₁²=−λ₁Res`. IST is A3's 1D-scattering realization; dictionary enrichment `c₁²↔f₁²`, `κ₁↔m₁`. Nonperturbative only in 1+1D integrable (sine-Gordon, O(3) NLσM): integrability supplies the exact `S`-matrix, TBA/GLM derives the gap (Hasenfratz–Niedermayer) — real, dimension-locked, the datum 4D SU(N) withholds (WH §5).

## Five-point scorecard
1. **Nonperturbative — PASS (character).** Exact GLM/TBA, no loop/`1/N`; concretely 1+1D only.
2. **Sign-controlling `G″≥λ` — FAIL.** Bijection is an identity; Temple needs `E₂`; min–max upper; Levinson integer.
3. **Non-circular — FAIL.** Atom is GLM input; Haag–Ruelle assumes the shell (= the gap).
4. **Non-SUSY + fails free field — PASS (letter, vacuous).** Free ⇒ no atom ⇒ `0`; `κ₁(q)` continuous ⇒ no forced floor.
5. **Uniform `≥λ` — FAIL.** No `(V,a,β)` floor; bound energy continuous in the unknown potential.

## Exact failing estimate
```
Target: G″(k→0)≥λ ⟺ Var(a)≤1/λ ⟺ E₁(L_{0⁺⁺})≥λ       [uniform floor on the lowest bound state]
dictionary: f₁² ↔ c₁² (normalization),  m₁²=−κ₁² (bound energy),  s_th ↔ a.c. branch point  [exact]
GLM: K+F+∫KF=0, F⊃c₁²e^{−κ₁ξ}, V=−2∂_xK(x,x)  ⇒ atom INPUT, potential OUTPUT          [circular, WH §2]
min–max: E₁=min_ψR(ψ) ⇒ trial gives R(ψ)≥E₁ (UPPER)                                    [direction wall]
Temple/Kato: E₁ lower bound consumes E₂ — gap one rung up (GAP-G2-edge)                 [circular]
Levinson: N=(1/π)ΔargT∈ℤ ⇒ rigidity no-go (B3) ⇒ frozen to free value 0                [index=0]
det(1+K_GLM)=1/T(k), zeros=atoms ≡ A3 D_O(z)=det(1−zT_O)                                [= A3]
free massless: V=0, R=0, no bound state, dρ from 0 ⇒ 0   [acid PASS; κ₁(q) continuous ⇒ vacuous]
```

**Sharpest line.** IST/GLM supplies the cleanest dictionary yet — the KL atom is literally a bound-state normalization `c₁²` at `−κ₁²`, and pole-vs-threshold is exactly what scattering data separates natively. But the dictionary is symmetric, and the symmetry is the disease: GLM reconstructs the potential *from* the atom, so certifying the atom certifies the input. The one subfield that *engineers* lower bounds on the lowest Schrödinger eigenvalue — Temple/Weinstein/Kato — is provably circular here (it consumes `E₂`, the gap one rung up); min–max gives only upper bounds (GAP-G2-edge's direction wall); and the one exact counting law, Levinson, is an integer the rigidity no-go freezes to zero across the connected β-family. So IST is A3's Fredholm determinant in a 1D-scattering alphabet, nonperturbative only where an exact `S`-matrix exists (1+1D), silent exactly at 4D SU(N). The missing object is a uniform floor `E₁≥λ` that no reconstruction identity can manufacture. **#1,4 pass; #2,3,5 fail. No fabrication.**

*Refs.* (external) Gel'fand–Levitan, *Izv. Akad. Nauk SSSR* **15** (1951) 309; Marchenko, *Dokl. Akad. Nauk SSSR* **104** (1955) 695; Krein, *Dokl. Akad. Nauk SSSR* **76** (1951) 21; Levinson, *Kgl. Danske Vid. Selsk. Mat.-Fys. Medd.* **25** (1949) 9; Temple, *Proc. Roy. Soc. A* **119** (1928) 276; Kato, *J. Phys. Soc. Jpn.* **4** (1949) 334; Deift–Trubowitz, *CPAM* **32** (1979) 121; Bargmann, *PNAS* **38** (1952) 961; Hasenfratz–Niedermayer, *Nucl. Phys. B* **337** (1990) 233; Zamolodchikov–Zamolodchikov, *Ann. Phys.* **120** (1979) 253. (internal) WIENER-HOPF (§2 Haag–Ruelle circularity, §5 1+1D TBA); PROOF-DISSOLUTION-SPECTRAL-DENSITY (`ν`, `s_th`, atom); GAP-G2-edge (moment-upper, residue–level-gap regress); C8-SPECIES-NOGO (empty-referent / false-implication pattern); INVENT-NEW-IDENTITY (A3 Fredholm determinant, zeros=atoms); INVENT-NEW-INDEX (rigidity no-go, index=0); PROOF-F1SQ (gapped free `f₁²=0`); PROOF-COARSE-POINCARE / PROOF-IR-H1 (H1 = gap); F2-KREIN-RUTMAN-QUANT (vacuous acid-test pass); QUANTUM-FISHER-BURES (continuous/circular pattern); PROOF-P3-UNCONDITIONAL (`f₁²>0` at `β<1/48`). No fabricated bound.
