# PROPAGATION-OF-CHAOS — the McKean–Vlasov particle dual is genuine and acid-passing, but propagation of chaos transfers the **mean-field** gap (`=G″_MF`, dynamically `≈0`, only the Haar floor `N/12` survives); the finite-`z` correction is **order-reducing** (mean field over-orders above `β_c`) and **log-marginal at `d=4` (`d_u=4`)** — i.e. `O(1)` at `z=6,N=2,3`, worst exactly at `β=2.5`

> **判词 (verdict).** Forge the self-consistent (McKean–Vlasov) stochastic quantization of the links: field replicas whose empirical measure `μ_n` converges (propagation of chaos) to the nonlinear Fokker–Planck limit `μ_∞`, fluctuation field `η=√z(μ_n−μ_∞)` given by the Bolthausen–Sznitman CLT. This yields an exact dictionary `Var(η_{k=0})=1/λ_edge=1/G″` (Ward `G″·Var(a)=1`): the lattice gap *is* the inverse slowest-mode variance of the McKean–Vlasov operator. **Genuine vehicle:** non-vacuous, cleanly reproduces `λ=0` for free massless (centered Gaussian ⇒ free `k=0` OU; acid PASS). **But it transfers no dynamical gap.** The limit gap is `λ_∞=G″_MF`, which by equipartition is the *non-dynamical* Haar floor (`E[Φ″]≥N/12`, `Δ_dyn≈0`; CAVITY-MEAN-FIELD §2). The `1/z` correction bridging to the lattice is **order-reducing** (Weiss over-orders just above `β_c≈2.47`, and `β=2.5` sits there) and **marginal at `d=4`** (`d_u=4`, log-enhanced), so at `z=2(d−1)=6`, `N=2,3` it is `O(1)` and *is* H1; every PoC→gap theorem needs a convexity at the fixed point that Savvidy falsifies. **#1 PARTIAL, #2,3,5 FAIL, #4 PASS — one clean pass.**

## §1 The forge — self-consistent stochastic quantization is McKean–Vlasov
YM is the invariant law of link Langevin `dU_l=−∇_l S_W dt+√2\,dW_l`, a *short-range* particle system with coordination `z=6`. Averaging the environment (Weiss) — each link sees `m_∞=E_μ[U]`, not its `6` neighbours — gives the nonlinear generator `L_∞` with drift `b(x;μ)=−∇(S_0+Φ+βΣ\,Re\,tr(x\,m_∞))`; equilibrium = cavity/Weiss measure. Propagation of chaos (Sznitman; Dawson; Dobrushin) gives `μ_n→μ_∞`, replica tensorization, and a Gaussian fluctuation (Bolthausen) with `Var(η_{k=0})=1/λ_edge=1/G″` — the dictionary we would need.

## §2 What the transfer delivers — the Haar floor, not a mass
`λ_edge|_{k=0}=G″_MF:=E[H+Φ″|0]` (CAVITY §1). In `G″=Δ_dyn+E[Φ″]+C` the dynamical part self-cancels by equipartition (`Δ_dyn=E[H]−Var(G)≈0`; B4 2.01→2.05), leaving only the **non-dynamical** Haar curvature `E[Φ″]≥N/12` (Lemma B) — present at `β=0`, a Jacobian, not a mass. The dual reproduces the *vacuous* floor, as the cavity did.

## §3 Wrong sign — and worst at `β=2.5`
A lattice lower bound would read `G″_lattice ≥ λ_∞ − (correction)`. **(i) Sign.** Mean-field `T_c` exceeds the true one, so just above `β_c` Weiss *over*-orders: `G″_MF > G″_lattice`, i.e. the correction is **order-reducing** — the wrong direction for a lower bound. `β=2.5` lies within `0.03` of `β_c≈2.47` (frontier), where Weiss is least reliable. **(ii) Size.** The correction is `O(1/z)` only for `d>4`; `d=4` is the upper critical dimension (`d_u=4`, Ginzburg), where it is **log-marginal** — concretely `1/z≈0.17`, `1/N_c=O(1)` at `N=2,3`. There is no large parameter (replicas `n` are auxiliary; physical coupling is local): the transfer does not approach H1, it *equals* it.

## §4 Acid test — PASS, structurally clean
Free massless: `m_∞=0`, `k=0` mode is a free OU (rate `0`) ⇒ `λ_∞=0`, and `Var(η_{k=0})=1/λ_∞=∞` ↔ `Var(a)=1/G″→∞` (Ward). The dual manufactures **no** gap on the free field. **Acid PASS.**

## §5 Five-point scorecard
1. **Nonperturbative — PARTIAL.** Exact stochastic-quantization + FP machinery (no `β`-expansion in the transfer), but the YM↔Weiss bridge is a `1/z` mean-field expansion failing at `z=6` (cf. CAVITY §1).
2. **Sign lower bound `G″≥λ` — FAIL.** Limit gap `=G″_MF` dynamically `≈0`; the correction is order-reducing and `O(1)` at physical `z,N_c`.
3. **Non-circular — FAIL.** `λ_∞=G″_MF=G″` (CAVITY §1); `Var(η_{k=0})=1/λ_edge`, Ward `Var(a)=1/G″`; proving the particle gap `≥λ` at `β=2.5` *is* the indefinite-Hessian/Savvidy problem = H1, zero slack.
4. **Non-SUSY + fails free field — PASS (clean).** Massless `λ_∞=0`; YM non-SUSY.
5. **Uniform — FAIL.** `d_u=4` ⇒ log-marginal; `z=6,N=2,3` ⇒ `O(1)`; Savvidy is bulk (equipartition, W6-BE 21–42% negative).

## Exact failing estimate
```
G″_lattice ≥ λ_∞ − C(β,d)/z^{α},  λ_∞=G″_MF=Δ_dyn+E[Φ″](+C)
  limit:  Δ_dyn≈0 (equipartition) ⇒ λ_∞=E[Φ″]≥N/12     [non-dynamical Haar floor]
  sign:   β=2.5>β_c≈2.47 ⇒ G″_MF≥G″_lattice             [order-reducing, wrong sign]
  size:   d=4=d_u ⇒ log-marginal, not O(1/z); z=6, 1/N_c=O(1) ⇒ correction O(1)
  ⇒ the "correction" is H1: G″≥λ ⟺ Var(a)≤1/λ
PoC→gap (Guillin–Liu–Wu; Malrieu; Cattiaux–Guillin–Malrieu) need convexity of the
  McKean–Vlasov potential — Savvidy H⪰̸0 (PᵀHP min≈−0.3)
```

**Sharpest line.** The dual is real: self-consistent stochastic quantization, a genuine McKean–Vlasov system whose Bolthausen–Sznitman fluctuations encode `G″` exactly as `1/Var(η_{k=0})`, and it correctly gives `λ=0` on the massless free field. But propagation of chaos converges to the *mean-field* measure, whose dynamical gap is zero by equipartition and whose only curvature is the non-dynamical Haar floor `N/12`; the `1/z` correction that would make it dynamical is order-reducing (Weiss over-orders near `β_c≈2.47`, where `β=2.5` sits) and log-marginal at `d=4=d_u`, hence `O(1)` at `z=6,N=2,3`. What the pipeline transfers at `k→0` is not a mass but `Var(a)≤1/λ` — H1 in fluctuation-field clothing. **#1 PARTIAL; #2,3,5 fail on the identical object; #4 cleanest.** *Not* 5/5.

*Refs.* McKean, PNAS **56** (1966) 1907; Sznitman, *Topics in Propagation of Chaos*, LNM **1464** (1991) 165; Dawson, J. Stat. Phys. **31** (1983) 29; Dobrushin, Theory Probab. Appl. **15** (1970) 458; Bolthausen, PTRF **72** (1986) 305; Ginzburg, Sov. Phys. Solid State **2** (1960) 1824 (`d_u=4`); Malrieu, Ann. Appl. Probab. **13** (2003) 540; Cattiaux–Guillin–Malrieu, PTRF **140** (2008) 19; Guillin–Liu–Wu, Ann. Probab. **50** (2022) 1590; 't Hooft, Nucl. Phys. B **72** (1974) 461; Savvidy, Phys. Lett. B **71** (1977) 133. Internal: CAVITY-MEAN-FIELD §1–3; D2-COUPLING-DOEBLIN §4; G1-UNIFORM-ZA §2; C8-SPECIES-NOGO; B4-VARIANCE-MAG; PROOF-JACOBIAN-SAVES; ABC-NEW-INEQ §1; PROOF-COARSE-POINCARE; W6-BE-at-2.5; docs/yang-mills-frontier (`β_c≈2.47`).
