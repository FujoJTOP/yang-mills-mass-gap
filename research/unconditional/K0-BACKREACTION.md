# K0-BACKREACTION — does the fiber back-reaction dominate at k=0? The conditioning shift is real but SMALL (+0.05); the sign is the Haar term `W_geo` (+0.3), and the literal target `E[H|0]+W > Var+|A_neg|` is just `G″>|A_neg|≈0.3` — marginal, not proven

> **2026-11-23 · Verdict (判词).** The angle is right about the *architecture* and wrong about the *carrier*. The exact identity is `G″(0)=E[H|0]−Var(∇S|0)+W` (flat chart; `H`=action Hessian, `W=W_geo+X`). At the softest accessible k (L=8, `k=π/4`) the **action bracket `E[H|0]−Var(∇S|0)` is statistically ZERO** (`+0.019±0.073` SU3, `+0.039±0.036` SU2) and **`W_geo` carries 91–111 % of `G″`** (`+0.319±0.005`→`G″+0.31`; `+0.252±0.004`→`+0.19`). The proposed "conditioning raise" `E[H|0]≫A_neg` is a **category error**: on its own footing **`E[H|0]=E[H]` to three digits** (1.989=1.989; 1.404=1.404) — conditioning does *not* raise the mean plane-wave Hessian. Its only action effect is a hair lower force variance (`Var(G)` 2.02→1.97 SU3): a **+0.05 conditioning shift**, not a +1.7 raise. The `−0.3→+1.4` story compares the *bare projected soft-mode eigenvalue* with the *conditioned mean of a different mode*. So back-reaction is a `+0.05` effect; the sign is the measure-geometric `W_geo`, present with or without conditioning. Two exact theorems pin this down (reflection ⇒ zero conditioned score; full equipartition ⇒ the entire `G″` is the *conditioning defect* of an otherwise-balanced IBP), but their sign is the Savvidy-indefinite fiber = H1. **Not 5/5; exact stuck step = (◆_k0) below.**

## (1) The exact split and the measured back-reaction
`G″(0)=bracket+W`, `bracket=E[H|0]−Var(G|0)`, `W=W_geo+X`. At L=8 `k=π/4`, delete-block jackknife (COMPUTE-GPP-K0-CONTINUUM):

| group | `G″` | `bracket` | `W_geo` | `X` | `E[H]` | `E[H|0]` | `Var(G)` |
|---|---|---|---|---|---|---|---|
| SU3 β=5.7 | +0.307±0.055 | +0.037±0.099 | **+0.319±0.005** | +0.035 | 1.989 | 1.989 | 1.951 |
| SU2 β=2.4 | +0.192±0.058 | +0.050±0.067 | **+0.254±0.004** | +0.025 | 1.404 | 1.404 | 1.354 |

`E[H]=E[H|0]` and `Var(G)≈E[G²]=E[H]` (equipartition, unconditioned `E[H]−E[G²] = −0.03/+0.02`). The conditioning changes the action bracket by only `+0.05`(SU3)/`+0.02`(SU2). **The sign lives in `W_geo`, not in the back-reaction.**

## (2) Steelman — three structural legs specific to k=0
**(a) Reflection theorem (exact, any k).** `a→−a` (`U→U†`) is an exact symmetry of the Wilson+Haar measure; `a=⟨v,log U⟩` flips sign, `S_W,Φ` are even. Hence `p(a)=p(−a)`, the score `d/da log p|_0=0`, so **`E[G|0]=E[Φ′|0]=0`** exactly (verified: `E[G]=−0.01/−0.02`, skew`(a)≈0.05`). Therefore `Var(G|0)=E[G²|0]`, `Var(Φ′|0)=E[Φ′²|0]`: the bracket is a *pure* variance-vs-Hessian defect, no leaking mean force.
**(b) Full-equipartition theorem ⇒ gap = conditioning defect.** The full-measure IBP gives `E[H+Φ″]=E[(G+Φ′)²]` exactly (compact group, no boundary) — the IBP is *balanced*. Since `p` is even, `G″(0)=E[H+Φ″|0]−E[(G+Φ′)²|0]` is the **difference between the slice defect and the full defect**: *the whole curvature is, verbatim, the a=0 conditioning defect.* "Back-reaction = 100 % of the gap" is a theorem (SUB1, ABC2-SYNTHESIS §27).
**(c) k=0 specificity.** The constant/uniform mode is where the Haar volume element is probed most coherently (`corr(a,Φ′)≈+0.85`), and `W_geo` is **β- and k-independent to 0.5 %** (`0.25` SU2 across β=0→2.4; `0.318/0.319/0.319/0.320` SU3 at (L4,kπ/2),(L6,kπ/3),(L8,kπ/4),(L8,kπ/2)). Physically, the Savvidy negative direction is the *delocalized lowest-Landau* mode at `k≈0`, so conditioning on the `k=0` amplitude directly targets the most negative direction — the soundest steelman for *why k=0*.

## (3) Attack — where it breaks
**(i) Wrong carrier.** `E[H|0]=E[H]` (§1): conditioning is a no-op on the mean plane-wave Hessian; the `+1.7` "raise" subtracts two different objects (bare soft eigenvalue `PᵀHP≈−0.3` vs conditioned mean of a *generic* mode `≈+1.4`). The special soft mode does raise (`E[H|0]=0.22>λ_E=0.04`, Part B), but by `+0.18` in a normalization where `W_geo` then has no such margin — a chart artifact, not dominance.
**(ii) Equipartition consumes the raise.** Even granting a raise, equipartition `E[f²]=E[H]` forces `E[H|0]≈Var(G|0)`; the bracket cannot stay positive by "raise dominates Var" — they track. Measured bracket ≈ 0, not `+`(raise).
**(iii) `W_geo` magnitude unproven.** Prékopa gives `W_geo≥0`; the floor `N/12` is a β=0 statement. At β=2.5 the measured `+0.32` has *no* proven lower bound beyond 0; keeping only `W_geo≥0` and `bracket+X` statistically, 95 % one-sided `−0.098` (SU3), `−0.035` (SU2) — **0 not excluded**.
**(iv) The literal target is `G″>|A_neg|`.** `E[H|0]+W>Var+|A|` ⟺ `G″>|A|`; with `G″≈+0.3`, `|A|≈0.3` that is `+0.30 vs +0.30` — *equality within errors*, a numerical coincidence, not a margin.

## (4) Forged tool — RASER (Reflection-Antisymmetric Score + Equipartition Reduction) + acid test
Two exact lemmas — **R** (`E[G|0]=E[Φ′|0]=0`, reflection) and **E** (`G″(0)=` slice−full IBP defect) — collapse the competition to one sign: `sign G″(0)=sign{E[H+Φ″|0]−E[(G+Φ′)²|0]}`. **Acid test (free massless):** `Φ=0`, `V=½∫(∂A)²` Gaussian ⇒ IBP balanced (defect 0) *and* the `k=0` zero-mode is flat (`Var(a)=∞` ⇒ `G″=0`); RASER returns `0`, not `>0` — non-vacuous. But its sign **is** the indefinite fiber: `E[H+Φ″|0]−E[(G+Φ′)²|0]` is `S−D` in SCBL (`D>0 ⟺ K` Savvidy-indefinite). RASER is a clean rename, not a closure.

## (5) Exact stuck step + five-point verdict
```
(◆_k0)  E[H+Φ''|a=0] − E[(G+Φ')²|a=0] ≥ λ>0  at k→0, uniform (L,U')
        ⟺ Var(a_{k=0}) ≤ 1/λ  ⟺  m_eff²(k→0) ≥ λ   = H1.
```
1. **Nonperturbative PASS** (exact identity + Lemmas R,E + Prékopa). 2. **Sign-controlling FAIL** (only `W_geo≥0` banked; action defect sign open). 3. **Non-circular FAIL** (closing = sign of the conditioning defect = H1). 4. **Free-field PASS** (RASER returns 0). 5. **Uniform `m_eff²≥λ` FAIL.**

**Sharpest line.** Back-reaction is *literally* 100 % of `G″` (the gap is the a=0 conditioning defect of a balanced IBP — a theorem), but its measurable action part is only `+0.05`: `E[H|0]=E[H]`, the raise is consumed by `Var(G)`. The sign is the Haar term `W_geo` (+0.3, β-robust, provably `≥0`, magnitude unproven). The target `E[H|0]+W>Var+|A|` is `G″>|A_neg|≈0.3`, satisfied marginally (`+0.30` vs `+0.30`) — a 1–2σ numerical near-miss, not a structural dominance. **Exact stuck step: `(◆_k0)` = the sign of the conditioned full-V IBP defect = the indefinite-fiber `D>0` = H1.**

*Refs.* k0-CONSTRAINT-EFFPOT (`V_eff″=G″`); PROOF-JACOBIAN-SAVES (`G″=Δ_dyn+W_geo+X`, `W_geo≥N/12`); INVENT-VARIANCE-INEQ (`G″=S−D`, `D>0⇔K` indefinite); B4-VARIANCE-MAG (`VarG≈2` equipartition-pinned); FORGE-BLIND (FB/TCI); COMPUTE-GPP-K0-CONTINUUM (L=8 table, −0.098/−0.035); OTHER-GPP-KZERO ((◆)=H1); ABC2-SYNTHESIS §27 (conditioning = whole gap); PROOF-COARSE-POINCARE (H1, `PᵀHP` min −0.3); C7-SPECIES-CONSTRUCT (`D>0`=indefinite `K`); MARGINAL (`PᵀHP` −0.3); W6-BE-at-2.5 (90–97 % fiber negatives); SUB1-COV-BOUND; data `k0_su3_b5.7_L8_n2_results.json`, `k0_su2_b2.4_L8_n2_results.json`; Prékopa 1971; Brascamp–Lieb 1976; Savvidy 1977.
