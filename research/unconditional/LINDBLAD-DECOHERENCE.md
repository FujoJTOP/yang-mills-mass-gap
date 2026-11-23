# LINDBLAD-DECOHERENCE — the Lindblad gap of the coarse mode **is** the mass, but quantum detailed balance makes it the quantum Poincaré constant `G″` (=H1) with zero slack; the *decoherence* rate passes the acid test cleanest of all species yet bounds the mass only from **above** (`T₂≤2T₁` `G″≤2γ_φ`; Ward–CS `G″≤Var(f|0)` at target β); the non-reversible NESS harbor is empty for equilibrium Gibbs YM

> **2026-11-23 · NEW SPECIES LINDBLAD (open-quantum; distinct coordinates from KOOPMAN operator-spectral and D2 coupling).** Treat the coarse amplitude `a=⟨v,X⟩` as an **open system** with the UV fiber as bath; the force `f=∂_nV=⟨v,∇V⟩` is the system–bath coupling. Two "Lindblad rates": (i) the **Lindbladian gap** `λ_L` (slowest observable relaxation toward `μ`) and (ii) the **dissipator gap** `γ_φ` (pure dephasing, decay of `a`-coherences). **Verdict.** The acid test is the **cleanest of the set**: the free massless zero mode is *exactly* bilinearly decoupled from the UV fiber, so `f≡0` and **both** rates vanish — while a lone *massive* oscillator keeps `m` yet decoheres not at all, pinning *"mass is a coherent/relaxation frequency, not a decoherence rate."* But both readings die: **(QDB)** YM SQ is `μ`-reversible (PROOF-GAUSS-INVARIANT-PF), so Alicki's quantum-detailed-balance theorem makes `λ_L` **equal** the quantum Dirichlet/Poincaré constant, `λ_1^{(a)}=−b′(0)=G″` (D2) — an *equality*, i.e. H1, zero slack. **(Dephasing)** `γ_φ≈Var(f|0)` fails the free field but bounds mass only from **above**: `T₂≤2T₁` gives `G″≤2γ_φ`, and Ward–Cauchy–Schwarz `E[af]=1` gives `Var(f|0)≥S=G″+D` (so `G″≤Var(f|0)≈2` at the target β, `D≥0`) — the **same sign wall** as E2/B4/C8. **1 PASS; 4 PASS (cleanest); 2,3,5 FAIL.**

## §1 The forge — two Lindblad rates
`dρ/dt=−i[H_eff,ρ]+D(ρ)`, `D(ρ)=Σ_k(L_kρL_k†−½{L_k†L_k,ρ})`. **(Relaxation)** `λ_L`: slowest nonzero `−Re` eigenvalue of the Lindbladian, the return rate of `a` to `μ`. **(Dephasing)** `γ_φ`: dissipative gap, bath coupling `∝f`, Caldeira–Leggett `γ_φ≈Var(f|0)`. Does either lower-bound `G″(k→0)`?

## §2 The two exact identities (banked)
**(Quantum detailed balance — Alicki).** A reversible Lindbladian has gap = quantum Dirichlet-form gap = quantum Poincaré constant. YM SQ is `μ`-reversible ⇒ `λ_L` = quantum Poincaré, slow-mode eigenvalue `−b′(0)=G″`:
```
λ_L (Lindblad gap) = quantum Poincaré = G″          [QDB — an EQUALITY]
```
**(Kubo fluctuation–dissipation, exact).** Relaxation rate = dephasing weight over susceptibility:
```
λ_relax · Var(a) = Var(f|0)  ⇒  λ_relax = Var(f|0)·S ≈ 2S ≈ 2G″   (equipartition ×2, D/S≈2%)
```
The FD rate tracks `G″` to an `O(1)` factor but **only through `Var(a)=1/S`** — the object H1 must bound.

## §3 Both rates vs the mass — the sign wall
**Unconditional:** the open-system hierarchy `T₂≤2T₁` gives `γ_φ≥½γ_1=½G″` ⇒ `G″≤2γ_φ` — an *upper* bound, the wrong direction for `G″≥λ`. **Ward–CS:** `E[af]=1` ⇒ `Var(f|0)≥1/Var(a)=S=G″+D`; at the target β (`D≥0`, `D/S≈2%`) this reads `G″≤S≤Var(f|0)≈2`, again *above*. Dephasing is an upper ceiling on relaxation, never a floor under it — the spine of E2 (information), B4 (magnitude), C8 (association).

## §4 Acid test — cleanest in the set
**Free massless `k=0`:** Gaussian ⇒ no cubic coupling of the zero mode to UV ⇒ `f≡0` ⇒ `γ_φ=0` **and** `λ_L=0` (flat, `b≡0`). **Free massive:** lone harmonic oscillator, frequency `m`, *still* decoupled from any bath ⇒ `λ_L=m²` but `γ_φ=0` — mass with **zero** decoherence. Pure YM: cubic/quartic coupling ⇒ `f≠0`, both rates `>0`. Returns `0` exactly when massless in *both* channels: sharp, non-vacuous. PASS — sharpest discriminant of any species.

## §5 The one harbor — non-reversible NESS
A non-reversible Lindblad generator supports a NESS with current and entropy production `σ>0`, with current-driven gap bounds `λ_L≥f(σ)`. But equilibrium Gibbs YM has **no current** (`σ=0`, QDB) ⇒ vacuous. Populating `σ>0` requires driving = real-time `−i[H,·]`, whose Lindblad gap is the energy ladder `=m_{0⁺⁺}` — circular (Koopman §3). The antisymmetric escape is the one Koopman sealed.

## §6 Five-point verdict
1. **Nonperturbative — PASS.** Exact QDB + Kubo + Ward, all β, no loops.
2. **Sign-controlling `G″≥λ` — FAIL.** Lindblad gap = `G″` (equality); dephasing gives `G″≤Var(f|0)` and `G″≤2γ_φ` — **upper** bounds, wrong sign.
3. **Non-circular — FAIL.** QDB ⇒ gap = quantum Poincaré = H1 verbatim; Kubo reintroduces `Var(a)`; NESS harbor empty (`σ=0`), drive = real time = circular.
4. **Non-SUSY + fails free field — PASS (cleanest).** Free massless: bilinear decoupling ⇒ `f≡0` ⇒ both rates `0`; massive `λ_L=m²`, `γ_φ=0`; pure YM.
5. **Uniform `≥λ` — FAIL.** `Var(f|0)≈2` is equipartition/UV-pinned (B4), no mass floor; `λ_L=G″` is the unproven H1 object.

## Exact failing estimate
```
QDB (reversible μ):  λ_L = quantum Dirichlet/Poincaré = G″                 [=H1, zero slack]
Kubo FD:             λ_relax·Var(a)=Var(f|0) ⇒ λ_relax = Var(f|0)·S ≈ 2G″
T₂ ≤ 2T₁:            γ_φ≥½γ_1 ⇒ G″ ≤ 2γ_φ                                  [UPPER bound, unconditional]
Ward–CS:             1=E[af]²≤Var(a)Var(f|0) ⇒ Var(f|0)≥S=G″+D; D≥0@β=2.5 ⇒ G″≤Var(f|0)≈2  [UPPER]
mass ≠ decoherence:  free massive oscillator has m≠0 but f≡0 ⇒ γ_φ=0       [acid-test sharpener]
NESS harbor:         σ>0 needs a current; Gibbs YM σ=0 (QDB) ⇒ vacuous; drive = real time = m_{0⁺⁺}  [circular]
```
**Sharpest line.** The Lindblad species owns the cleanest acid test of all — the free massless zero mode is *exactly* bilinearly decoupled, so `f≡0` kills **both** rates identically, and a lone massive oscillator keeps its `m` while decohering not at all: mass is a **coherent/relaxation** frequency, **not** a decoherence rate. The full Lindblad gap is genuinely the mass, but QDB (reversible YM measure) makes it *equal* the quantum Poincaré constant `G″=H1` with no slack; and the decoherence rate, which does separate the phases, bounds the mass only from **above** (`T₂≤2T₁` `G″≤2γ_φ`; Ward–CS `Var(f|0)≥S=G″+D` ⇒ `G″≤Var(f|0)≈2` at the target β) — the spine of E2/B4/C8. The one non-reversible harbor, a NESS current `σ>0`, is empty for Gibbs YM, and its only filler is the real-time mass ladder. **#1,#4 pass; #2,#3,#5 fail (one object: the reversible quantum Dirichlet gap = H1, plus a wrong-signed dephasing companion).**

*Refs.* Lindblad, CMP **48** (1976) 119; Gorini–Kossakowski–Sudarshan, J. Math. Phys. **17** (1976) 821; Alicki, Rep. Math. Phys. **11** (1977) 331 (QDB); Fagnola–Rebolledo, Lett. Math. Phys. **41** (1997) 317; Spohn, Lett. Math. Phys. **2** (1978) 33 (entropy production, approach to equilibrium); Kubo, Rep. Prog. Phys. **29** (1966) 255 (fluctuation–dissipation); Caldeira–Leggett, Physica A **121** (1983) 587 (quantum Brownian dephasing); Zurek, RMP **75** (2003) 715 (decoherence). Internal: E2-INFORMATION (information bounds variance below; SDPI `η`=spectral gap); B4-VARIANCE-MAG (`Var(f|0)≈2` equipartition, sign-blind); KOOPMAN (reversible ⇒ gap = Poincaré = H1; antisymmetric escape = real time, closed); C8-SPECIES-NOGO (same upper-bound/sign wall); D2-COUPLING-DOEBLIN (`λ_1^{(a)}=−b′(0)=G″`); B3-FUNCTIONAL-BE (`Γ₂^ν=G″`, lossless); CONTINUUM-VARIANCE-DEFECT (`G″=S−D`, `S=1/Var(a)`); PROOF-GAUSS-INVARIANT-PF (self-adjoint reversible transfer); INVENT-VARIANCE-INEQ (Ward `E[af]=1`). No fabricated bound; a confirmed reduce-to-H1 (relaxation gap) / wrong-sign (dephasing), derived fresh in the Lindblad vocabulary.
