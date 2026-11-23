# GIRSANOV — change of measure to a gapped reference: `Var_μ(a)≤C_ν e^{osc W}` is **right-sign** (first in the E-family) but carries the extensive relative entropy `H(μ|ν)∼βV`, so "uniform" ⟺ "gapped reference with `O(1)` RN" ⟺ the gap

> **判词 (verdict).** `dμ∝e^{−S_YM}`, `dν∝e^{−S_ref}`, `W=S_YM−S_ref`. **Girsanov/Holley–Stroock is exact and sign-correct:** `ν` PI`(C_ν)` (gapped, `C_ν=1/m_ref²`) ⇒ `μ` PI`(C_ν e^{osc W})`, i.e. `Var_μ(a)≤C_ν e^{osc W}E_μ|∇a|²` — an **upper** bound on `Var_μ`, hence a **lower** bound on `S=1/Var` and the mass: the **first right-sign transfer** in the E-family (E2 was sign-wrong). **But it is a relabelling.** *Global:* `osc W=O(βV)` ⇒ `C_μ≲e^{−cV}`, vacuous. *Transport:* `Var_μ(a)≤2/m_ref²+O(H(μ|ν))`, `H(μ|ν)∼βV` ⇒ no uniform floor; `O(1)` entropy ⟺ bounded perturbation of a gapped `ν` ⟺ gap. *Curvature:* Bakry–Émery needs `Hess S_YM⪰λ`, killed by the physical Savvidy `−2g|F_B|`. Only non-circular cell = Dobrushin ⇒ `β<1/48` = **P3 re-derived**. #1 PASS; #2 direction-PASS/floor-FAIL; #3 vehicle-PASS/close-FAIL; #4 PASS (letter); #5 FAIL.

## 1. The transfer (exact)
`S=1/Var(a)=Γ^{(2)}_T(k)`, gapped ⟺ `S→m_eff²>0` (CONTINUUM-VARIANCE-DEFECT §3). **Holley–Stroock:** `dμ∝e^{−W}dν`, `ν` PI`(C_ν)` ⇒ `μ` PI`(C_ν e^{osc W})`. An **upper** bound on `Var(a)` **is** a **lower** bound on `m_eff²` — the direction every information species lacked (E2 bounded `Var` from *below*); the sign is inside Girsanov, free.

## 2. Why it doesn't cash: `osc W` is the entropy
`W=βΣ_p(1/N)Re tr U_p−½⟨A,(−Δ+m_ref²)A⟩` sums over `≈4V` links ⇒ `osc W=O(βV)`, `C_μ≲e^{−cV}`: the global transfer says nothing as `V→∞`. The local repair (Dobrushin one-site influence `<1`, Zegarlinski uniform LSI) is genuinely non-circular — bounded one-site oscillations, **not** `Var(a)` — but holds only below the uniqueness threshold, the strong-coupling cell `β<1/48` already owned by Lemma SC / P3. At `β=2.5` (Dobrushin condition badly violated — strong-coupling/scaling region, not high temperature) no gapped reference has bounded influence; the reference one *needs* is itself confining, so its "known gap" is the gap.

## 3. Transport face = the same wall
Otto–Villani `W₂²(μ,ν)≤K C_ν H(μ|ν)` with `Var_μ(a)≤2Var_ν(a)+2Lip²W₂²` gives `Var_μ(a)≤2/m_ref²+O(H(μ|ν))`. `Var_ν(a)=1/m_ref²` is finite **by construction**; all rests on `H(μ|ν)∼βV` (C7 §c). Uniform `Var_μ` ⟺ `H(μ|ν)=O(1)` ⟺ gapped `ν` with bounded RN ⟺ gap: an *identity* shape, not an undone estimate.

## 4. Curvature face dies at Savvidy
Bakry–Émery: gap if `Hess S_YM⪰λ`. `Hess S_YM|_mag=−D²+2adF_B` has transverse-magnetic `−2g|F_B|` (W4: 7.1% negative, physical, survives the quotient — PROOF-SAVVIDY-WGEO-DOMINATION). A reference `Hess S_ref=m_ref²` only moves the bar: need `Hess W⪰λ−m_ref²`, defeated on the exact soft mode. Girsanov cannot launder an indefinite Hessian.

## 5. Free-field acid test (letter PASS, spirit tautology)
Free massless = the **gapless** reference: `λ_ref=0`, `1/λ_ref=∞`, transfer returns `0` — no spurious gap. And `H(massless|massive)=½Σ_k[m²/k²+log(k²/(k²+m²))]→+∞` at `k→0`: a gapless measure **cannot** be dominated by a gapped one with finite entropy. Passes *in letter* — but the divergence `∼Var(a)∼m²L²` **is** `Var(a)`: E2's spirit-failure, sign repaired.

## Five-point scorecard
1. **Nonperturbative — PASS.** HS/Girsanov/Otto–Villani exact, all β, finite `V`; no loops.
2. **`G″(k→0)≥λ` — direction PASS, floor FAIL.** First right-sign transfer, but `e^{−cV}`/`βV` entropy ⇒ no uniform `λ`.
3. **Non-circular — vehicle PASS, closing FAIL.** Dobrushin uses `osc W`, not `Var(a)`; non-vacuous ⟺ gapped reference with `O(1)` RN ⟺ gap.
4. **Non-SUSY + fails free field — PASS (letter).** Massless `λ_ref=0`⇒0, `H=∞`; spirit-tautology caveat.
5. **Uniform `≥λ` — FAIL.** `osc W=O(βV)`, `H∼βV`; no `V`-floor.

## Exact failing estimate
```
Girsanov/HS : Var_μ(a) ≤ C_ν e^{osc W} E_μ|∇a|² , C_ν=1/m_ref² ; S=1/Var ⇒ m_eff² ≥ (C_ν e^{osc W})^{-1}
  global : osc W=O(βV) ⇒ e^{−cV} → 0                         [vacuous]
  T2/OV  : Var_μ(a) ≤ 2/m_ref² + O(H(μ|ν)) , H(μ|ν)~βV        [extensive ⇒ no uniform floor]
  ⟺ Var uniform ⟺ H(μ|ν)=O(1) ⟺ bounded-entropy perturbation of gapped ν ⟺ H1/gap
  BE cur.: need Hess S_YM ⪰ λ ; Savvidy −D²+2adF has −2g|F_B| ⇒ FAIL at soft mode
  Dobrushin : influence<1 ⇒ uniform LSI ⇒ β<1/48 (P3), NOT β=2.5
  acid : massless λ_ref=0 ⇒ 0 ; H(massless|massive)=½Σ_k[m²/k²+log(k²/(k²+m²))]→+∞ (=Var)
  sign↔gap : S≥λ is H1 ; G″=S−D≥λ also needs D<S−λ (D/S≈2%, still H1)
```

**Sharpest line.** Girsanov is the **first mechanism in the information/transport family to get the sign right**: a gapped reference upper-bounds `Var_μ(a)` (lower-bounds the mass) with **no `Var_μ` on the RHS** — Holley–Stroock is exact, a real advance over E2's sign-wrong I-MMSE. But the species *is* the wall: the transferred constant `e^{osc W}`/`H(μ|ν)` is exactly the extensive `∼βV` entropy C7 pinned; demanding it `O(1)` is demanding `μ` be a bounded-entropy perturbation of a gapped measure — H1 as an entropy budget. The only non-circular window (Dobrushin) is the strong-coupling cell `β<1/48` = P3; the Bakry–Émery face dies at the physical Savvidy negative mode. The residual object, unchanged, is now an entropy: **a gapped reference with `O(1)` RN oscillation at `β=2.5`**. No fabrication.

*Refs.* Holley–Stroock, *J. Stat. Phys.* **46** (1987) 1159 (PI/LSI `C_ν e^{osc W}`); Girsanov, *Theory Probab. Appl.* **5** (1960) 285; Bakry–Émery, *LNM* **1123** (1985); Bakry–Gentil–Ledoux, *Analysis and Geometry of Markov Diffusion Operators* (Springer, 2014); Otto–Villani, *J. Funct. Anal.* **173** (2000) 361 (T₂/entropy); Villani, *Optimal Transport: Old and New* (2009); Dobrushin, *Theory Probab. Appl.* **13** (1968) 197, **15** (1970) 458; Zegarlinski (uniform log-Sobolev, lattice spin systems, early 1990s). House: E2-INFORMATION (`G″Var(a)=1`, `S=1/Var`, sign-wrong contrast); C7-SPECIES-CONSTRUCT §(c) (`H(μ_β|μ_0)∼βV`, HWI); CONTINUUM-VARIANCE-DEFECT §3 (`S≥λ`=H1, `D/S≈2%`); BANACH-FIXED-POINT (Savvidy breaks contraction); C8-SPECIES-NOGO; ABC-NEW-INEQ (Ward `E[af]=1`); INVENT-VARIANCE-INEQ (`G″=S−D`); PROOF-P3 / PROOF-COARSE-POINCARE (Lemma SC, `β<1/48`); PROOF-SAVVIDY-WGEO-DOMINATION (`−2g|F_B|` physical); W4-7.1%; YANG-MILLS-MASTER §2.1/§2.5.
