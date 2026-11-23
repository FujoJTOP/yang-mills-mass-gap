# QUANTUM-FISHER-BURES — species QFB: SLD / Bures / quantum Cramér–Rao on the coarse q-mode. A genuinely noncommutative metric and a new *coherence* diagnostic — but the QCRB bounds variance from below, the free vacuum is QFI-extremal, and "diverges-free / finite-YM" is again `Var(a)=1/G″`. Same wall, noncommutative alphabet.

> **判词 (verdict).** Forge: treat the `k→0` transverse mode `a=⟨v,X⟩` as a **quantum state family** `{ρ_θ}`. The SLD quantum Fisher information `F_Q=Tr(ρ L_θ²)`, `½{L_θ,ρ}=∂_θρ`, gives the Bures metric `ds²=¼F_Q dθ²` and the quantum Cramér–Rao bound `Var(θ̂)≥1/(νF_Q)`. Target `G″(k→0)≥λ ⟺ Var(a)≤1/λ`. **Verdict: real new *diagnostic*, cleanest free-field statement in the family — but sign-wrong and circular like E2.** (Banked) `F_Q` is the *optimal-over-all-POVM* Fisher (Helstrom); its gap above a fixed homodyne Fisher, `F_Q−F_{cl}`, is a **coherence measure** of the mode — positive for the mixed non-Gaussian YM state (`β=2.5`, `D/S≈2%`), zero for the *pure Gaussian* free vacuum. (Why it dies) The QCRB bounds estimator variance from **below**; the target needs the *prior* `Var(a)` from **above**; and for any Gaussian/pure mode `F_Q∝Var(a)`, so "finite QFI" ⟺ `Var(a)<∞` ⟺ `S=1/Var≥λ` ⟺ **H1**. The Gaussian is the QFI-extremizer as it was the entropy extremizer in E2; the free vacuum has `F_Q∝1/k²→∞` — the tool is **loudest on the free field**, not failing on it. **#1,#4 pass (letter); #2,#3,#5 fail. Not 5/5.**

## 1 Forge and dictionary
`ρ` = reduced state of the mode; family `ρ_θ=e^{−iθĜ}ρ e^{iθĜ}` with generator `Ĝ` conjugate to `a`. Pure-state SLD QFI `F_Q=4Var(Ĝ)`; with `Ĝ=a`, `F_Q=4Var(a)` — the quantum face of E2's score `J=Var_p(E[f|a])`, **optimized over all measurements** instead of fixed on the `a`-marginal. `G″≈1/Var(a)` (COMPUTE-GPP-SU3; AQFT-CAUSAL §1), so `G″·(F_Q/4)≈1` iff Gaussian — the noncommutative generalization of `G″Var(a)=1`.

## 2 Sign — QCRB points the wrong way
`Var(θ̂)≥1/(νF_Q)` is a **lower** bound; the quantum de Bruijn identity makes `F_Q` the entropy-derivative under noise, same monotonicity as classical. The target needs `Var(a)≤1/λ` from **above**. No QCRB, Bayesian/van Trees, or Bures speed limit bounds a *prior* variance from above — every quantum information quantity bounds it from below. E2's diagnosis verbatim, with the SLD metric.

## 3 Novelty — a diagnostic, not a lever
Helstrom's `F_Q=max_{POVM}F_{cl}` makes the **gap** `F_Q−F_{cl}` meaningful: the mode's quantum advantage, positive iff `[ρ,∂_θρ]≠0` (mixed/non-Gaussian), zero for the pure Gaussian free vacuum. A **new alphabet** distinct from the classical defect `D` — a *coherence* reading of the Wall. But it is `O(D)` and `∝Var(a)`; bounding it gives no upper bound on `Var`. The SLD relaxation gap of the block generator *is* the dissipative gap `G″` (cf. F2's Bushell rate) — circular.

## 4 Acid test — letter met, spirit vacuous
Free massless `k→0`: `ρ` is the pure Gaussian vacuum, `Var(a)∝c/k²→∞`, `F_Q∝1/k²→∞`; YM finite `∝1/m²`. "Diverges free / finite YM" — **letter** met, and here more cleanly than E2 (free mode *pure Gaussian*, QFI explicit). But the divergence **is** `Var(a)`: a bijection with the free covariance, and the free field is QFI-**maximal**. The charter wants the tool to *fail* at `G″=0` free; it instead peaks there — so the discriminator is vacuous, not discriminating.

## 5 Five-point scorecard
1. **Nonperturbative — PASS.** SLD/Bures/QCRB exact, all β, no loops / `1/N` / `g`.
2. **Sign-controlling `m_eff²(k→0)≥λ` — FAIL.** QCRB bounds variance *below*; target needs the prior from *above*.
3. **Non-circular — FAIL.** `F_Q∝Var(a)`; `G″(F_Q/4)=1` iff Gaussian; "finite QFI" ≡ `S≥λ` ≡ H1; SLD relaxation gap = `G″` (F2).
4. **Non-SUSY + fails free field — PASS (letter), FAIL (spirit).** Free pure-Gaussian QFI diverges `∝1/k²` (cleanest free statement) — but free is QFI-*maximal* ⇒ vacuous.
5. **Uniform `≥λ` — FAIL.** No positive floor without already containing `1/Var`.

## Exact failing estimate
```
SLD QFI: F_Q=Tr(ρL²), ½{L,ρ}=∂_θρ ;  Bures ds²=¼F_Q dθ² ;  QCRB Var(θ̂)≥1/(νF_Q)
coarse q-mode: k→0 transverse a=⟨v,X⟩ ;  target G″(k→0)≥λ ⟺ Var(a)≤1/λ
Gaussian/pure: F_Q ∝ Var(a)  ⇒  G″·(F_Q/4) ≈ 1           [generalizes G″Var(a)=1]
free massless: pure Gaussian vacuum ⇒ F_Q ∝ 1/k² → ∞ ; YM finite ∝ 1/m²   [letter met]
QCRB: LOWER bound on estimator variance ⇒ cannot give Var(a) ≤ 1/λ (upper)   [sign FAIL]
F_Q = max_{POVM} F_cl (Helstrom) ⇒ F_Q − F_cl = coherence excess, O(D), ∝ Var [new diagnostic]
finite QFI ⟺ Var(a)<∞ ⟺ S=1/Var ≥ λ ⟺ H1 ;  SLD relaxation gap = dissipative gap = G″  [circular]
```

**Sharpest line.** QFB is the most *quantum* language yet for the coarse mode — the SLD/Bures metric and the coherence excess `F_Q−F_cl` see off-diagonal information classically invisible — and it sharpens the free-field test (the free mode is a *pure Gaussian vacuum*, QFI-extremal, `∝1/k²`). But it is the same wall in noncommutative dress: the QCRB is a *lower* bound (wrong sign for `Var(a)≤1/λ`), `F_Q` is a covariance function `∝Var(a)` (finite QFI ≡ H1), and the SLD gap of the relaxation *is* the dissipative gap `G″` (F2). The bankable residue is a **diagnostic, not a lever**: the coherence excess `F_Q−F_cl=O(D/S)` as an independent noncommutativity reading of the `k→0` mode. The missing step is unchanged — an upper bound on `Var(a)=1/G″` at `k→0` not derived from `Var(a)` — H1, now in a third information alphabet.

*Refs.* Helstrom, *Quantum Detection and Estimation Theory*, Academic Press (1976); Braunstein–Caves, *Phys. Rev. Lett.* **72** (1994) 3439; Bures, *Trans. Amer. Math. Soc.* **135** (1969) 199; Uhlmann, *Rep. Math. Phys.* **9** (1976) 273; Holevo, *Probabilistic and Statistical Aspects of Quantum Theory*, North-Holland (1982) (SLD≤RLD, monotonicity); Petz, *Linear Algebra Appl.* **244** (1996) 81 (monotone metrics); Paris, *Int. J. Quantum Inf.* **7** (2009) 125; Giovannetti–Lloyd–Maccone, *Nat. Photonics* **5** (2011) 222; Pinel–Jian–Treps–Fabre–Braun, *Phys. Rev. A* **88** (2013) 040102 (single-mode Gaussian QFI); Mandelstam–Tamm, *J. Phys. USSR* **9** (1945) 249. Internal: E2-INFORMATION (classical Fisher score = Ward force `J=Var(E[f|a])`, Gaussian extremality, sign-wrong); ABC-NEW-INEQ (`G″Var(a)=1` iff Gaussian); CONTINUUM-VARIANCE-DEFECT / INVENT-VARIANCE-INEQ (`G″=S−D`, `D/S≈2%`); COMPUTE-GPP-SU3 (`1/Var=0.40/0.38`); F2-KREIN-RUTMAN-QUANT (rate = `G″` = H1, sign-correct-but-circular); AQFT-CAUSAL (mass-blind, `G″≈1/Var`); C8-SPECIES-NOGO (empty-referent pattern).
