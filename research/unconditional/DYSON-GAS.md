# DYSON-GAS — the transfer spectrum is **not** a Vandermonde log-gas; its low edge is a soft (Airy) continuum edge, and the only exact log-gas corners (2D YM, Eguchi–Kawai) make the **effective action**, not the mass

> 2026-11-23 · NEW SPECIES: Dyson-gas / random-matrix transfer operator. Treat the low eigenvalues of `T=e^{−aH}` as a Coulomb gas (log-gas) sourced by the Haar/Vandermonde measure, and bound the low edge `G″(k→0)≥λ` by equilibrium-measure / edge-repulsion methods (`G″(k→0)=1/Var(a₀)=m_eff²(k→0)`, INVENT-VARIANCE-INEQ; PROOF-COARSE-POINCARE). **判词: in 4D the transfer spectrum is not a Vandermonde log-gas — the Vandermonde lives on one link's holonomy, not on the energies. Where a log-gas is exact (2D YM = GWW/Douglas–Kazakov; large-N Eguchi–Kawai) its edge is the free energy = effective action = H1. The energy edge, where it exists, is a *soft* `ρ~√(E−E₀)` continuum edge — vanishing density = no atom = anti-mass.** Five-point: #1 PASS-in-principle(unbuilt); #2,#3,#5 FAIL; #4 PASS narrow, vacuous.

**§1 Where the Vandermonde actually is.** Weyl's integration formula makes a *single* Haar link a circular Dyson gas at β=2: `dg` pushes to torus eigenvalues with density `∏_{i<j}sin²((θ_i−θ_j)/2)` (CUE). The repulsion is on **holonomy eigenvalue spacings** — a property of the *gauge field*, not of `T`. `T=e^{−aH}` is a *fixed* positive operator for each `(g,a,L)`; its eigenvalues `e^{−aE_n}` carry no Haar pushforward and no `|Δ(E)|²`. There is no Vandermonde on the energies.

**§2 The exact log-gas corners do not contain the 4D mass.** Two honest identifications exist, neither 4D. *(2D YM)* is a unitary matrix model `Z=Σ_R(dim R)²e^{−g²C₂(R)A/2N}`, a Gross–Witten–Wadia log-gas with a Douglas–Kazakov edge. *(Large-N reduction)* Eguchi–Kawai makes 4D YM a single-site multi-matrix model = a (multi-)matrix Coulomb gas. Both are `N→∞` statements whose "edge" is the **free energy** `−log Z`, i.e. the effective action, H1; and EK is center-unstable off the quenched/twisted patch. The 4D master field that would define the gas is unbuilt (H1-T8; free probability proven only 0D/2D).

**§3 The would-be bound, and its exact failure.** If the energies were a β-ensemble with potential `W=ΣW(E_n)`, the equilibrium measure `μ_eq=argmin(I[μ]−∫W dμ)` would give `supp μ_eq=[E₀,E₁]` with `N^{−2/3}` rigidity about `E₀` (Tracy–Widom). For `G″≥λ` we would need `E₀≥λ` with `W` known *without* the gap — but `W=aH`, so `E₀=inf spec H` **is** the gap. Repulsion controls **spacing**, never **location**. Worse, the energy-edge density is `ρ(E)~√(E−E₀)` (Airy soft edge; Bleher–Its, Deift): *absolutely continuous, vanishing density at the edge* — the categorical opposite of the mass atom `f₁²=ν({m₁²})>0` (INEQ-R1; PROOF-F1SQ; GAP-G2-edge). A soft edge is a **no-atom theorem**.

**§4 Conditioning is again H1.** A gap in a Coulomb gas is a large-deviation event, rate `N²I(δ)` (Ben Arous–Guionnet): a *probability* of a low eigenvalue over an ensemble, not a deterministic floor, with tilt rate = effective action — the same LD-tilt = H1 structure as E4.

**Five-point scorecard.**
1. **Nonperturbative — PASS in principle, unbuilt.** Edge/equilibrium-measure methods are exact `N→∞`; the 4D gas is the unbuilt master field.
2. **Sign LOWER bound — FAIL.** Repulsion is spacing, not location; `ρ~√(E−E₀)` is a vanishing-density continuum edge (no `f₁²>0`).
3. **Non-circular — FAIL.** Edge location = potential `W=aH` = `inf spec H` = gap; exact corners give `−log Z` = H1.
4. **Non-SUSY + fails free field — PASS (vacuous).** Non-SUSY ✓; Gaussian transfer has no Vandermonde, `[0,∞)` edge at 0, no gap ✓ — but the free field *is not a log-gas*, so the machinery never bites.
5. **Uniform `λ` — FAIL.** `N^{−2/3}` rigidity about a potential-set location; no `λ` uniform over `β∈[1/48,2.5]`, `V→∞`; zero soft-edge density.

**Exact failing estimate / missing object.**
```
Vandermonde: on holonomy torus of ONE link (CUE, β=2), NOT on energies E_n of T=e^{−aH}.
T fixed by (g,a,L): no ensemble, no |Δ(E)|² ⇒ no log-gas over the gauge measure.
Exact corners: 2D YM = GWW/Douglas–Kazakov log-gas ⇒ edge = −log Z = eff. action (H1).
               Eguchi–Kawai large-N = (multi-)matrix Coulomb gas ⇒ same; center-unstable off patch.
Energy edge:  SOFT Airy, ρ(E)~√(E−E₀) ⇒ vanishing density ⇒ NO atom (f₁²=0) = anti-mass.
Location:     E₀ = inf spec H = confining potential W=aH = the gap ⇒ circular.
Gap prob.:    P(min E>λ) ≲ e^{−N²I} (Ben Arous–Guionnet) ⇒ probabilistic, rate=eff.action=H1.
Free field:   Gaussian, no Vandermonde, ρ on [0,∞), edge 0, gapless ✓ — vacuous pass.
```

**Sharpest line.** The Vandermonde is real, but it repels **holonomy eigenvalues of a link**, not the **energies of the 4D transfer operator**; the transfer's own low edge is a soft `√(E−E₀)` continuum edge — zero density, hence **no mass atom**, the exact opposite of what the gap needs. The only exact log-gas edges (2D YM, Eguchi–Kawai) compute `−log Z` = the effective action = H1, at `N→∞` on an unbuilt/unstable patch. Dyson repulsion controls **spacing**, never **location**; the location is the confining potential, which is the Hamiltonian, which is the gap. The missing object is not a Coulomb-gas theorem — it is the identification of a 4D energy edge with a **hard, atom-bearing** log-gas edge whose potential is known without the gap: exactly the absent 4D master field, renamed.

**Regime pin / Record.** `G″(k→0)≥λ` not reached. CANNOT exist as an energy-spectrum log-gas (no Vandermonde on `E_n`) or as a repulsion floor (spacing ≠ location; soft edge has zero density). COULD exist — conditional only — inside a 2D YM / EK large-N log-gas, where the "edge" is the effective action (H1) and the model is center-unstable off the quenched/twisted patch. Genuinely absent in the open 4D cell.

*Refs.* Dyson, *J. Math. Phys.* **3** (1962) 140; Mehta, *Random Matrices* (3rd ed., 2004); Weyl, *The Classical Groups* (1939); Deift–Kriecherbauer–McLaughlin, *J. Approx. Theory* **95** (1998) 388; Deift, *Orthogonal Polynomials and Random Matrices* (AMS 1999); Bleher–Its, *Ann. Math.* **150** (1999) 185; Tracy–Widom, *CMP* **159** (1994) 151; Ben Arous–Guionnet, *CMP* **218** (2001) 177; Gross–Witten, *PRD* **21** (1980) 446; Wadia, *PLB* **93** (1980) 403; Douglas–Kazakov, *PLB* **319** (1993) 219; Eguchi–Kawai, *PRL* **48** (1982) 1063; 't Hooft, *NPB* **72** (1974) 461; Osterwalder–Seiler, *Ann. Phys.* **110** (1978) 440. Internal: E4-RMT-EDGE, k0-FREEPROB-AVERAGE, PROOF-GAUSS-INVARIANT-PF, C8-SPECIES-NOGO, k0-SPECTRAL-REDUCTION, GAP-G2-edge, PROOF-F1SQ, INEQ-R1, H1-T8-free-prob, PROOF-COARSE-POINCARE, INVENT-VARIANCE-INEQ.
