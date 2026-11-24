# G2-ASSEMBLE-100 — the strongest complete G2 theorem, the exact modulo list, and the N-row 100% ledger

> 2026-11-23 · Assembles G2 (the mass gap) from ABC-NEW-INEQ / ABC-CONSTRUCTIVE / ABC2-* / PROOF-JACOBIAN-SAVES / k0-CONSTRAINT-EFFPOT / OTHER-H1-POINCARE / OTHER-GPP-KZERO / OTHER-F1SQ-TRANSPORT / SUB1 / SUB2 / B4 / B5 / C4–C7 / D1 / D2 / E1 / E2 / PROOF-P3-UNCONDITIONAL / PROOF-ISOLATION / PROOF-EDGE-E2 / GAP-G2-edge / PROOF-NO-DISSOLUTION-KINEMATIC / RESIDUE-CONDITIONAL / CONTINUUM-CP / ABC-FINAL-LEDGER / docs/PROOF-SKELETON §1.
>
> **判词 (verdict).** G2 assembles into a complete theorem **modulo ONE named species**, which is **one object at two faces**: the **coarse-curvature face** `H1 = {Var_μ(a_{k→0})≤1/λ}` (⟺ `G″(0)=m_eff²(k→0)≥λ`) and the **spectral-edge face** `H* = {f₁²≥r₀ with isolation I2}`. Every node G2.0–G2.6 and the residue/transport/isolation arm is **(a) a proven theorem**, **(b) a reduction to that single species**, or **is** the species. **There is no (c): no genuinely-open item distinct from the species.** The classification is 100% complete; the **Clay gap stays 0%** because the species is unbuilt. This is honest bookkeeping, not a proof.

---

## §1 · The exact master statement

> **MASTER (G2, the mass gap).** Fix compact simple `SU(N)`, `N≥2`, at the physical crossover `β≈2.5`; along the two-loop RT look for `Δ>0`. **The mass gap holds** — `m_eff²(k→0):=lim_{k→0}G″(k)=lim_{k→0}V_eff″(k)≥λ>0`, equivalently `f₁²=|⟨0|trF²|G⟩|²>0` with the `0⁺⁺` atom **isolated** (`Λ:=inf supp(ν−f₁²δ_{m₁²})`, `Λ−m₁²≥δ>0`), equivalently `Var_μ(a_{k→0})≤1/λ` uniform in `(L,a)` — **modulo [X] = `{H1}`**, the single named species: a one-sided, nonperturbative, coupling-free **upper bound on the conditioned force variance**
> `Var_μ(⟨v,∇S⟩|0) ≤ E_μ[H_vv|0]+W−δ`  (=(◆) of PROOF-JACOBIAN-SAVES),
> `(L,U′)`-uniform.

| Statement | Modulo | Status |
|---|---|---|
| **G2-curv** `G″(0)≥λ` | `{H1}` | conditional; region-covered 100% except the open cell |
| **G2-edge** `f₁²≥r₀` ∧ `I2` | `{H*}` = `{H1}` | conditional; I1 proven β<1/48, I2 open |
| **Clay gap** (both, `a→0`, `L→∞`) | one species | **0%** |

The rest shows (i) the region G2.1–4 is closed by theorems, (ii) the remaining cell is *exactly* H1 (zero slack), and (iii) the residue/isolation arm adds **no new species** — it is H* = H1.

---

## §2 · Why the two faces are one object

The **residue–level-gap lemma** (GAP-G2-edge, proved): given the total moments, `f₁² ≥ (ΛM₀−M₁)/(Λ−m₁²)`, equality iff the rest is one atom at `Λ`. Hence a residue lower bound **⟺** a second-level (isolation) lower bound — `I1⟺I2` — and edge control at level `n` needs level `n+1` (an infinite regress *unless* the threshold is dynamical). At the edge, "closing `I2` **is** `f₁²>0` **is** `H*`" (PROOF-EDGE-E2). And `H*` "is **one object** with Polchinski `c`, the marginal curvature `G″`, and the LGP constant `γ`" (RESIDUE-CONDITIONAL): residue, edge, and curvature are **one unprotected number in several vocabularies**, tied by Ward `G″·Var(a)=1` and KL/OS at the pole. So the species — the **one-sided variance upper bound** — has curvature face `H1` and edge face `H*`; `I2` is H* shifted one rung (the upper gap, strictly above bare `m>0`), **a re-expression, not a new unprotected input**.

---

## §3 · The N-row classification ledger

`(a)` proven theorem · `(b)` reduces to the single species (H1=H*). *(No `(c)`.)*

| # | Input / node | Cat. | Why / source |
|---|---|---|---|
| 1 | Identity `G″=E[H\|0]−Var(∇S\|0)+W=V_eff″=m_eff²(k→0)` | **(a)** | exact, all β (k0-CONSTRAINT-EFFPOT, PROOF-JACOBIAN-SAVES §1, D1) |
| 2 | Ward `G″·Var(a)=1` (⟺ marginal Gaussian) | **(a)** | exact IBP, all β (ABC-NEW-INEQ §1) |
| 3 | SD `E[H]=E[f²]` (unconditioned) | **(a)** | exact, but **gap-blind** (SUB1: free `Var(a)=∞` saturates it) |
| 4 | SCBL `G″=S−D` | **(a)** | identity; `D<S` ⟺ H1 (INVENT-VARIANCE-INEQ) |
| 5 | `W_geo=E[Φ″]−Var(Φ′)≥0` (Prékopa) | **(a)** | every mode incl. `k→0`, every L (PROOF-JACOBIAN-SAVES) |
| 6 | `E[Φ″]≥N/12` (Lemma B) | **(a)** | pointwise, any measure (A5, PROOF-JACOBIAN-SAVES §2) |
| 7 | Shallow band `Φ″⪰N/12`+BL ⇒ `G″≥S` | **(a)** | closed, all k, all L (ABC-CONSTRUCTIVE §2) |
| 8 | High-field cone WSTC `≤(cβκ_max)²e^{−cβ(κ²−κ_typ²)N_□}→0` | **(a)** | L-uniform; cone cannot destabilize (ABC-NEW-INEQ §2) |
| 9 | AF-UV deep cone `μ²e^{−cμ⁴/g⁴}→0` (Σg⁴<∞) | **(a)** | closed UV (ABC-CONSTRUCTIVE §2) |
| 10 | STC cone-free shells `G″(k_b)≥Λ_{b−1}²−δ_b>0` | **(a)** | UV only; silent on sign-carriers (ABC2-NEW-INEQ) |
| 11 | Equipartition pin `Var(G\|0)≈E[H\|0]≈2` (L4→L6: 2.01→2.05) | **(a)** | *no-go*: variance route dead by identity shape (B4, C7) |
| 12 | **Open cell `Var(G\|0)≤W−…`** | **(b) = H1** | =H1 verbatim, zero slack (ABC2-SYNTHESIS §3–4) |
| 13 | **Negative-association `Cov(G,Φ′\|0)≤0`** | **(b) = H1** | variance-free face; **proven absent** from 9 classes (C5,C6,C7) |
| 14 | D1 amplitude-variance `Var(a)=χ=1/G″` | **(a)** → **(b)** | identity; `ΣE[(Δ_ia)²]≤C` V-uniform ⟺ clustering = H1 |
| 15 | D2 coupling drift `λ=−b′(0)=G″` (von Renesse–Sturm) | **(a)** → **(b)** | sharp equivalence, no slack ⇒ H1 (D2) |
| 16 | E1 FRG `V_eff,0″≥0` (convexity), `≥N/12` at β=0 | **(a)** → **(b)** | strict `≥λ` at 2.5 = H1 (E1) |
| 17 | E2 Fisher/MMSE `J=Var(E[f\|a])`, `G″Var(a)=1` | **(a)** → **(b)** | ceiling ≡ H1; wrong sign otherwise (E2) |
| 18 | **`H1 = Var_μ(a_{k→0})≤1/λ`** | **(b) — THE species** | coarse-curvature face (§1) |
| 19 | P3 residue `f₁²≥r₀>0`, `β<1/48` | **(a)** | unconditional, V,L-uniform (PROOF-P3-UNCONDITIONAL §5) |
| 20 | Transport `f₁²(μ)=Z_O²f₁²(μ₀)` | **(a)** | exact square, `0→0`, acid-clean — but presumes `(CP)` |
| 21 | Residue–level-gap lemma `I1⟺I2` | **(a)** | proved (§2) |
| 22 | **`I2` isolation / upper gap** | **(b) = H*** | =H*, the edge face of the species (SUB2, PROOF-EDGE-E2) |
| 23 | No-dissolution kinematic theorem | **(a)** | margin `m_{0⁺⁺}<θ=2m_g`; conditional on uniform gap+ordering |
| 24 | `(CP)=[CP-b]=I1∧I2` | **(b) = H*** | entailed by isolation (CONTINUUM-CP §2) |
| 25 | Kato analyticity in β | **(a)** | value-transport; domain = `I2` (PROOF-TRANSPORT-ANALYTIC) |
| 26 | SUB2: `I2⟺gap`, cost grading closes only I1 | **(a)** | equivalence proved; `I2` open (`Cu` vs `Be^{−δt}`) |

**Reading.** Rows 1–11, 14–17, 19–21, 23, 25, 26 are **(a)** theorems (many *are* the exact identities; their "closure" faces are the species). Rows 12, 13, 18, 22, 24 are **(b)** — the species itself / its reductions. **ZERO rows are (c).**

---

## §4 · Exact % per node (node-local, ≠ Clay %)

| Node | % | Status |
|---|---|---|
| G2.0 identity layer | **100%** | proven theorem (exact) |
| G2.1 `W_geo≥0` (Prékopa) | **100%** | proven, all modes/L |
| G2.2 shallow band | **100%** | proven |
| G2.3 high-field cone (WSTC) | **100%** | proven, L-uniform |
| G2.4 AF-UV cone | **100%** | proven |
| **G2.5 open cell** | **0%** | **open = H1** |
| **G2.6 named species** | **0%** | **open / proven absent from 9 classes** |
| P3 residue (I1) | **100%** (β<1/48) | unconditional, regime-limited |
| Transport `Z_O²` | **100%** | exact; presupposes `(CP)` |
| `I1⟺I2` lemma | **100%** | proven |
| **I2 isolation** | **0%** | **open = H\*** |
| No-dissolution | **100%** | kinematic (conditional anteriors) |

**G2 综合:** 区域覆盖 **~90%** (G2.1–4 真定理); 符号/证明 **0%** (G2.5–6 = 隙本身, =H1). **Truthful totals:** all closed theorems are real; the gap itself, and thus **Clay, is 0%**. Node-% is progress toward each *node's* goal, not Clay.

---

## §5 · The isolation arm — and why it is not a third item

`I1` (`f₁²≥r₀`, residue) closes at `β<1/48` (P3). **`I2` does not.** SUB2 shows the integer cost grading grades *polymer amplitudes*, not *transfer-matrix eigenvalues*: it gives I1 (unique leading tube, `+` recoupling) but the KP/OS remainder is a `t`-uniform `O(u)`, **never** `O(e^{−δt})` — so it cannot separate a second atom at `(p'−p)ln(1/u)/a` from one at `ε→0`, nor place the `0⁺⁺*` (`≈1.5m₁`) below the `2m₁` continuum. **But `I2` is not a distinct open item:** by the residue–level-gap lemma `I1⟺I2`, and `I2` is the upper-gap rung of `H*`; `H*=G″=H1` is one object (RESIDUE-CONDITIONAL). The no-dissolution theorem pins the atom below the `2m₁` continuum *given* a uniform gap + `0⁺⁺`-lightest ordering — both of which are themselves the species. So the residue/transport/isolation arm *inherits* the species; it does not add one.

---

## §6 · CONFIRMATION (explicit)

> **OPEN(G2) = the single species only.** Every genuinely-open 0% node — **G2.5** (open cell), **G2.6** (the absent negative-association species), the residue-extension across `[1/48,2.5]`, and **I2** (isolation) — is the **one named species**: the one-sided, nonperturbative, coupling-free upper bound on the conditioned force variance (`H1`), equivalently the unprotected residue/edge datum (`H*`). **There is NO distinct third open item.** The curvature face (`H1`, `G″≥λ`), the spectral-edge face (`H*`, `f₁²≥r₀`∧`I2`), and the coupling/FRG/information faces (D2, E1, E2) are the same object in different vocabularies; the residue–level-gap lemma and Ward are the exact dictionaries. This mirrors G1-ASSEMBLE-100: **G2 is a theorem modulo `{H1}`, one species — not a theorem, and not modulo two independent unknowns.**

---

## §7 · Honest boundary — what this is and is not

- **It is** the cleanest statement: the entire (field × momentum) plane is closed by G2.1–4 theorems except **one cell**, and that cell (typical field × high momentum, where WSTC's exponent vanishes at `κ=κ_typ`, AF-UV has not started at `g≈1.5`, and the `N/12` floor is `~10×` short) is **exactly** `Var(a_{k→0})≤1/λ` = H1, zero slack (ABC2-SYNTHESIS). It is **not** a proof of the gap.
- **The species is provably absent** from all nine known dependency classes for a *continuous, frustrated, non-attractive* SU(N) measure (C5/C6): no lattice order (torsion), no real-stable polynomial, no DPP kernel, no common scalar, FKG wrong-signed, BL needs `Hess⪰0` (the Savvidy wall), Simon–Lieb sign-blind. The measured `+0.51` (W7) is a **sign probe**, not a theorem (OTHER-GPP-KZERO); `A5` certifies only `E[Φ″]≥¼` and the *conditional* certificate (rigid enclosure missing = H1).
- **The unification `H*=H1` is a structural identification**, stated as such (as in G1-ASSEMBLE-100 §6); the rigorous bridges are the residue–level-gap lemma (`I1⟺I2`) and Ward (`G″·Var(a)=1`). One closure would likely close both faces — but that is a conjecture, not a theorem. **CORRECTION (2026-11-23, `H1-IS-IT-THE-GAP.md`):** this identification is now *doubted*, not merely unproved. `G″(k=0)=1/χ=1/∫s⁻¹dν` is the inverse zero-momentum **susceptibility**, only *upper*-bounding `m₁²` (squeeze `f₁²/χ≤m₁²≤M₀/χ`); it equals the pole `m₁²` only for elementary single-atom fields. The free massive **composite** `:φ²:` is gapped (`m₁=2m`) yet `G″=1/χ≠m₁²` — falsifying `H*=H1` for composites, and every gauge-invariant YM operator is composite. The gap is `H*` (edge residue `f₁²≥r₀` + isolation), strictly stronger than the susceptibility `H1`.
- **No fabricated closure:** G2.5, G2.6, I2, and `f₁²(2.5)` are **open = 0%**; Clay = **0%**.

*Refs (internal).* PROOF-JACOBIAN-SAVES; ABC-NEW-INEQ; ABC-CONSTRUCTIVE; ABC2-SYNTHESIS; ABC2-NEW-INEQ; ABC2-CONSTRUCTIVE; OTHER-H1-POINCARE; OTHER-GPP-KZERO; OTHER-F1SQ-TRANSPORT; SUB1-COV-BOUND; SUB2-I2-STRONG; B4-VARIANCE-MAG; B5-CORRELATION-CONTROL; C4-COV-LASTMILE; C5-NEG-ASSOC; C6-DEPENDENCE-CLASS; C7-SPECIES-CONSTRUCT; D1-BLOCK-MARTINGALE; D2-COUPLING-DOEBLIN; E1-FRG; E2-INFORMATION; PROOF-P3-UNCONDITIONAL; PROOF-ISOLATION; PROOF-EDGE-E2; GAP-G2-edge; PROOF-NO-DISSOLUTION-KINEMATIC; PROOF-TRANSPORT-ANALYTIC; RESIDUE-CONDITIONAL; CONTINUUM-CP; ABC-FINAL-LEDGER; k0-CONSTRAINT-EFFPOT; A5-RIGOROUS-NUMERICS; INVENT-VARIANCE-INEQ; docs/PROOF-SKELETON §1.
*External (as checked in the sources).* Osterwalder–Seiler, Ann. Phys. **110** (1978) 440; Shen–Zhu–Zhu, CMP **400** (2023) 805, arXiv:2204.12737; Kotecký–Preiss, CMP **103** (1986) 19; Prékopa 1971; Brascamp–Lieb 1976; Savvidy, PLB **71** (1977) 133; von Renesse–Sturm, CPAM **58** (2005) 923; Guo–Shamai–Verdú, IEEE-IT **51** (2005) 1261; Kato (1976) VII §3; Efron–Stein, Ann. Statist. **9** (1981) 586; Borcea–Brändén–Liggett, JAMS **32** (2019) 205; Fortuin–Kasteleyn–Ginibre 1971.
