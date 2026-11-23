# G1-CONTINUUM-SCALING — the λ_a covariance is free; "uniform G″_a≥λ_a + AF ⇒ m_phys>0" is the gap, not G1

> **判词 (verdict).** Split three ways. **(S) Scaling is free.** Along the critical line the coarse curvature renormalises **multiplicatively** through the trace anomaly, `λ(μ)=Z_O²λ(μ₀)`, `Z_O=(β/2g)(μ₀)/(β/2g)(μ)>0` — the *same* transport that carries `f₁²`, exact, a square (positive for free), `O(1)` in `V`, acid-passing (`0→0`). The required engineering law is fixed: the dimensionless lattice curvature `G″_a=a²m_phys²` must run as `a²Δ₀²` (the lattice gap `m_lat=√G″_a~aΔ₀`); a bound **uniform in `a`** is the free-box scale, not a particle. **(T) The implication is not a theorem given G1.** G1's uniform `Z_a` is the `∇_⊥` UV-fiber bound = existence, explicitly `⊥` gap (G1-EXISTENCE, G1-UNIFORM-ZA §4). The claim needs the `∇_∥` `k→0` positivity (H1), `L`-uniformity, and uniqueness/`(CP)` — none in G1. **(W) It is a distinct wall.** Read literally (`λ_a` uniform-positive in `a`) it gives `m_phys=+∞` (lattice artifact); read correctly (`λ_a≳a²Δ₀²`) the premise *is* H1 = the gap. **Net: covariance = free (banked, conditional on `(CP)`); base positivity = H1 = G2; the "bridge" is the gap restated.**

## 1 · The exact scaling of λ_a
`G″_a` is the Hessian of the *dimensionless* Wilson action, so it is `a²` times a physical mass²: `G″_a=a²m_phys²`, lattice gap `m_lat=√G″_a=a·m_phys`, continuum `Δ=lim m_lat/a=m_phys`. Hence `Δ≥Δ₀>0` requires
```
G″_a ≥ a²Δ₀²   ⟺   m_lat=√G″_a ≥ aΔ₀.       (the "λ_a~a·Δ_phys" law is √G″, not G″)
```
A bound **uniform in `a`** (`G″_a≥λ₀>0`) yields `m_phys=lim√λ₀/a=+∞`: it is the free-box/UV gap `c_L=4π²/L²` (P1-L0), a lattice artifact, *not* a mass. So "uniform fixed-`a` lower bound" is dimensionally wrong for the gap; the right object is a lower bound on the **ratio** `G″_a/a²` — the asymptotic-scaling rate (string-tension template `a²σ_lat↔σ_phys`, confinement L7). A physical gap means the lattice curvature *vanishes like `a²`*, not stays positive.

## 2 · The covariance is free — and carries zero
`S=1/Var(a)=m_eff²(k→0)` and `f₁²` are both marginal `Tr F²`-insertions, so both transport multiplicatively through `θ^μ_μ=(β/2g)Tr F²`:
```
λ(μ)=Z_O(μ,μ₀)²λ(μ₀),   γ_O=−d ln(β/2g)/d ln μ=2b₀g²+O(g⁴).
```
`Z_O²` is a **square** and `O(1)` in `V` (never `e^{−cV}`); this is literally INVENT-TRANSPORT-INEQ applied to the curvature. It passes the acid test **by construction**: a homomorphism sends `0→0`, so the free field (`G″=0`, `f₁²=0`) rides `0↦0`. The covariance is exact and cheap *given a seed* — but it transports only the seed's **value**, never its positivity or size; it cannot manufacture `Δ₀>0`.

## 3 · Why G1 does not supply the implication
G1 gives fixed-`L` existence of a **subsequential** `μ_L` (OS0/1/3), conditional on L3. Its uniform `Z_a=sup_a Z_a⁻¹‖H_a‖_D<∞` is the **fiber** (`∇_⊥`) moment bound; the ledger calls L3 "`⊥` gap." The claim's `G″_a≥λ_a` at `k→0` lives on the orthogonal `∇_∥` coarse side (PROOF-IR-H1). Two readings, neither free:
- **(A) `λ_a`=L3's per-shell `δ_W`.** Then `G″_a≥λ_a` *is* L3, whose entire payoff is **existence, never gap** — the sign is BPHZ-invisible (G1-UNIFORM-ZA §4). The implication fails exactly there.
- **(B) `λ_a`=`m_eff²(k→0)`.** Then the premise *is* H1; the implication is circular.

AF fixes only the **rate** at which a would-be gap approaches the Gaussian point (`g→0`); the fixed point is free and gapless (PROOF-STEP-SCALING). AF forces no positivity.

## 4 · What the continuum needs beyond uniform Z_a
1. **Subspace** — the bound must sit on `∇_∥` (`k→0`), not the `∇_⊥` fiber G1 controls.
2. **Volume** — `L`-uniform (G1's `4π²/L²→0`; the `L→∞` limit is G2 / stage-2 UVS-L).
3. **Rate** — the `a²`-engineering law `G″_a/a²≥Δ₀²`: covariance free via `Z_O²`, but the constant `Δ₀>0` is H1.
4. **Seed preservation** — the subsequence must keep the edge atom (Kato isolation/uniqueness; G1 is subsequential).
5. **Phase continuity** — `(CP)`: strong- and weak-coupling lattices one `T=0` phase (INVENT-TRANSPORT-INEQ).

## 5 · Honest boundary
The scaling is genuinely free and banked: *given* H1 at one scale and `(CP)`, the physical gap follows by exact multiplicative transport (`Z_O²`, explicit, acid-safe). But the proposed bridge is **not** cheap from G1: G1's uniform `Z_a` is an existence bound on the wrong (`⊥`) half with the wrong (`a`-independent) scaling. The remainder factorises cleanly — **base positivity = H1 = G2; volume/uniqueness = G2; phase = `(CP)`** — and none is delivered by AF, `Z_a`, or covariance. **The wall is not the transport; it is the nonzero seed the transport multiplies.**

*Refs.* G1-EXISTENCE-THEOREM; G1-UNIFORM-ZA §4; G1-L3; CONTINUUM-VARIANCE-DEFECT; INVENT-TRANSPORT-INEQ; PROOF-GAP-DIFFERENTIAL-INEQ; PROOF-STEP-SCALING; PROOF-G1-BPHZ; PROOF-IR-H1; P1-L0-proof; GAP-G2-mass; CONDITIONAL-LATTICE-GAP; A1-twosided-UV; confinement L7.
