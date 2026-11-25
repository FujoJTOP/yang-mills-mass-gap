# Z4-MERA — the entanglement-renormalization paradigm is the **wrong tool for a gap**: a scale-invariant MERA *is* a gapless fixed point (acid test fails), 1D reduces to the transfer matrix, the bond dimension is not volume-uniform, and the 2D solvable model's MERA carries topology, not the gap magnitude

> **2026-11-24 · MERA / entanglement-spectrum extension of X3 (MPS/PEPS).**
> Extends `X3-TENSOR.md`. Question: can a **scale-invariant MERA** / **entanglement-spectrum (ES)
> gap** give a *rigorous, volume-uniform* mass-gap certificate for 4D Yang–Mills?
> Tested on 1D/2D **solvable** models (free boson, transverse-field Ising = free fermion, 2D toric
> code). **Exact / solvable only** — no Monte Carlo, no fit. Witness `research/v3/w7/z4_mera.py` →
> `z4_mera_results.json` (`py -3.12 -I`, numpy+scipy; **one process**, tiny matrices; no downloads).
> Reads: `X3-TENSOR.md`, `H1-RG-DESCENT.md`, `D3-TRANSFER-CONSTANT.md`, `LIT-ERA-11.md`,
> `docs/paper-v2/ch12-numerics.tex`.
> **Labels (charter §2).** **[T]** true theorem/identity · **[V]** exact numerical witness (this
> file) · **[C]** conditional/protocol · **[NG]** no-go for a named class · **[O]** exact missing
> step. **No QED. Clay 0%.**

---

## 0 · Verdict (up front; no QED, Clay 0%)

**The paradigm is not a *worse* certificate — it is the *wrong object*: a scale-invariant MERA is,
by construction, the description of a critical (gapless) fixed point, and the mass gap asks for the
opposite. Concretely:**

1. **[T]/[NG] "Scale-invariant MERA ⇒ gap" fails the free-massless acid test on its face.** MERA's
   defining feature — the *disentanglers* that keep the bond dimension from growing at a critical
   point (Vidal 2007) — exists precisely so the ansatz can represent a **gapless** state with
   constant bond dimension. A scale-invariant MERA (identical layers) is the fixed point of
   entanglement RG = a conformal (gapless) state; a **gapped** state's MERA flows to a **product
   state** in `O(log ξ)` layers. The **free massless field has a scale-invariant MERA**. So any rule
   that reads a *gap* off scale invariance certifies a gap for the massless field → **acid-test
   FAIL**. TNR makes this explicit: it "recovers scale invariance *at criticality*" — i.e. at **zero
   gap**.

2. **[T]/[V]/[NG] The honest reading — the ES gap — passes the acid test only because it *is* the
   transfer-matrix gap (kill #1).** The entanglement spectrum `{−log p_i}` of the half-chain, and
   its gap, are computed from the **same transfer-matrix / correlation data** as `1−|λ₂|`. Exact
   1D checks (§2): free boson ES gap is a monotone function of `μ = −log|λ₂|` (massive: bounded in
   `N`; massless: `→0`); TFI ES gap tracks the free-fermion gap, closing at `g=1`. So in the one
   dimension where the certificate can be made rigorous, **it reduces to the transfer matrix**
   (`D3`'s object) — the same kill as X3.

3. **[T]/[NG] The bond dimension that a *rigorous* certificate would need is not volume-uniform
   (kill #2).** MERA's per-layer bond dimension `χ` is `L`-independent, but the only rigorous gap
   theorem (Fannes–Nachtergaele–Werner, 1D **MPS**) requires a **volume-uniform MPS bond dimension**.
   A 1D MERA with `χ`, branching `b`, and `L` sites is an MPS of bond dimension
   `χ^{log_b L} = L^{log_b χ}` — **polynomial in `L`, not uniform** (§4: `χ=4` gives `L²`; `χ=16`,
   `L=1024` ⇒ `1.1×10¹²`). The MERA "cheapness" does not convert into a uniform-bond rigorous object.

4. **[T]/[NG] In the one solvable gauge model, MERA carries the *topology*, not the gap magnitude.**
   The 2D toric code has gap `4` (anyon mass `2`) *from its stabilizer/Hamiltonian spectrum* (X3,
   exact). Its **exact finite-bond MERA** (Aguado–Vidal 2008) represents all four ground states as
   the *same network differing only in the top tensor* — the MERA/ES exposes the **topological
   degeneracy** `4 = 2^{2g}` and the area law, **not** the gap magnitude `4`. (That `gap = 4` and
   `GSD = 4` coincide is a feature of this particular fixed-point model — one is an energy, the other
   a Hilbert-space dimension — not a causal link.)

5. **[T]/[NG] MERA is variational with no parent-Hamiltonian gap theorem.** `E_MERA ≥ E₀` is an
   **energy** bound (upper), not a gap bound. The rigorous "injective ⇒ uniformly gapped parent `H`"
   theorem is for **1D MPS**, not MERA; for PEPS it is **uniqueness only**, gap open (X3). MERA
   inherits neither.

**Kill-criteria scorecard.** **#1 (only reproduces the transfer-matrix gap): MET** (§2, §0.2).
**#2 (bond dimension not volume-uniform): MET** in the sense that matters — the rigorous object is
`χ^{log L}`, not uniform (§4, §0.3). **#3 (reduces to H1): NOT met** — MERA's *entanglement* RG is
genuinely different from H1's *Wilsonian block* RG (`G″`, Savvidy soft mode); they do not reduce to
each other (§6). But that does not rescue the paradigm: it dies earlier, on the acid test and on #1/#2.

**One line.** *MERA/entanglement-spectrum is not a certificate for a gap but for its absence: a
scale-invariant MERA is a critical fixed point, so "scale-invariant MERA ⇒ gap" fails the
free-massless acid test by construction; the honest ES-gap reading passes only by *being* the
transfer-matrix gap in 1D (kill #1); the bond dimension a rigorous theorem needs is
`χ^{log L}`, not volume-uniform (kill #2); it does **not** reduce to H1, but that is moot; and on
the solvable gauge model MERA yields the topological degeneracy, not the gap magnitude. **No QED.
Clay 0%.***

---

## 1 · What MERA is, and what is actually a theorem

MERA (Vidal 2007, 2008) is a real-space tensor network of *disentanglers* (unitaries removing
short-range entanglement) alternating with *isometries* (coarse-graining). Its two rigorous /
semi-rigorous credits are:

| statement | status |
|---|---|
| MERA is **variational**: `⟨Ψ|H|Ψ⟩ ≥ E₀`; minimization approaches the GS energy | **[T]** |
| a state with a `χ`-bond MERA has area-law entanglement by construction, and is efficiently representable | **[T]** (Vidal 2008) |
| at a **critical** point, a **scale-invariant** MERA reproduces the **CFT scaling dimensions** | **[T]/[C]** (Vidal 2007; numerics) |
| **injective MPS** (1D) ⇒ unique, uniformly gapped parent `H` | **[T]** (FNW '92; PGWV '07) |
| injective **PEPS** ⇒ unique parent GS; **gap** open | **[T]** / **[O]** (X3) |
| **injective MERA ⇒ gapped parent `H`** | **does not exist** — no such theorem |

The decisive structural fact is the third row: MERA's variational power *at criticality* is exactly
its purpose. The scaling dimensions it reproduces are the **critical (gapless) exponents** — properties
of a zero-gap theory. The mass gap is a property of a **gapped** fixed point, where MERA degenerates
to a trivial (finite-depth) network. So the paradigm's centre of gravity is orthogonal to the question.

**Why MERA adds nothing rigorous over X3.** X3 already banked: 1D MPS certificate = transfer-matrix
gap `[NG]`; 2D/3D injectivity ⇒ uniqueness only, gap `[O]`; solvable gauge model non-injective. MERA
covers the *complementary* regime (critical, where MPS/PEPS bond dimension blows up) — but that
regime **is** the zero gap. There is no theorem that MERA's fixed-point structure bounds the gap of a
*non-critical* Hamiltonian.

---

## 2 · Acid test on 1D solvable models — exact (witness §[M1],[M2])

### 2.1 Free scalar: the ES gap is a transfer-matrix quantity  [T]/[V]

The free real scalar `−Δ+m²` has the **exact** transfer-matrix gap `δ_T = 1−e^{−μ}`,
`cosh μ = 1+m²/2` (X3). Its half-chain **entanglement spectrum** (Gaussian covariance:
`X=½K^{−1/2}`, `P=½K^{1/2}`, symplectic eigenvalues `ν_k=√eig(X_A P_A)`, ES energies
`ε_k=log((ν_k+½)/(ν_k−½))`) gives:

| `m` | `μ=1/ξ` | `δ_T` | ES gap (N=32) | (N=64) | (N=128) |
|---|---|---|---|---|---|
| 1.00 | 0.962 | 0.618 | 4.61975 | 4.61975 | 4.61975 |
| 0.50 | 0.495 | 0.390 | 3.54024 | 3.54024 | 3.54024 |
| 0.20 | 0.200 | 0.181 | 2.67221 | 2.67370 | 2.67370 |
| **0.00** | **0** | **0** | **0.00135** | **0.00118** | **0.00095** |

Massive: ES gap `>0`, **bounded in `N`**. Massless: ES gap `→0` as `N` grows. The ES gap is a
**monotone function of `μ = −log|λ₂|`** — it vanishes **iff** `m→0`, i.e. it carries no information
beyond the transfer matrix. **Acid test passed; but the certificate *is* the transfer matrix (kill #1).**

### 2.2 Transverse-field Ising (free fermion): gap and ES gap both close at `g=1`  [V]

Exact sparse ED, `H=−Σσ^xσ^x−gΣσ^z`, periodic, `L=8…16`. (Restricted to `g≥1`: unique GS, so the
`k=2` gap is the physical gap. The ordered phase `g<1` has an exponentially-split doublet whose `k=2`
splitting and cat-state ES are *not* the physical quantities — excluded, not a gap claim.)

| `g` | phys gap (L=16) | `gap·L` (L=8→16) | ES gap (L=8→16) | limit |
|---|---|---|---|---|
| **1.00** | 0.0983 | 1.576,1.574,1.573,1.572,**1.572** | 1.175→1.005 | **`∝1/L→0`** (massless) |
| 1.05 | 0.1591 | 2.05…2.55 | 1.37→1.35 | `→0.1` |
| 1.10 | 0.2344 | 2.59…3.75 | 1.57→1.69 | `→0.2` |
| 1.30 | 0.6038 | 5.18…9.66 | 2.29→2.67 | `→0.6` |
| 1.50 | 1.0005 | 8.16…16.0 | 2.86→3.22 | `→1.0 = 2(g−1)` |

At `g=1` (the free **massless** Majorana point) `gap·L = 1.572` is **constant** ⇒ `gap ∝ 1/L → 0`
cleanly; `g>1` gapped with `gap → 2(g−1)`. The ES gap closes at `g=1` (slowly — the `Z₂` pair of top
Schmidt values is near-degenerate; the free-boson §2.1 nails closure unambiguously). **The ES gap
passes the acid test — but only because the model is free-fermion solvable, and the ES is the
correlation/transfer data.** No statement survives for an interacting 4D theory.

---

## 3 · The category error: scale invariance certifies the *absence* of a gap  [T]/[NG]

Three exact facts make the "scale-invariant MERA ⇒ gap" reading self-refuting:

1. **The free massless field is scale-invariant.** Its MERA is scale-invariant (identical layers;
   it is a `c=1` CFT). A rule "scale-invariant MERA ⇒ gap>0" would certify a gap for the massless
   field → **acid-test FAIL** (§0.1).
2. **A gapped state's MERA is trivial.** Entanglement RG on a gapped state flows to a **product**
   fixed point in `O(log ξ)` layers; the only "scale-invariant" (`L`-independent-layer) fixed points
   are **critical** ones.
3. **TNR's stated purpose is criticality.** Evenbly–Vidal's tensor-network RG (2015) is built to
   "recover scale invariance explicitly **at criticality**" — i.e. at **zero gap**. The paradigm's
   success story is the *gapless* fixed point.

So the "scale-invariant MERA / ES gap" slogan is a name for detecting **criticality**, which is the
**opposite** of a mass gap. The only non-vacuous reading (the ES gap) is the transfer matrix (§2).

---

## 4 · Bond dimension: the rigorous object is not volume-uniform  [T]/[NG]

A 1D MERA with per-layer bond dimension `χ`, branching `b`, and `L` sites has `ℓ=log_b L` layers;
flattened to an MPS it has bond dimension `χ_eff = χ^{ℓ} = L^{log_b χ}` (witness §[M4], `b=2`):

| `χ` | `L=16` | `L=64` | `L=256` | `L=1024` |
|---|---|---|---|---|
| 2 | 16 | 64 | 256 | 1024 |
| 4 | 256 | 4096 | 6.55×10⁴ | 1.05×10⁶ |
| 16 | 6.55×10⁴ | 1.68×10⁷ | 4.30×10⁹ | 1.10×10¹² |

`χ_eff` grows **polynomially in `L`** — not volume-uniform. The only rigorous 1D gap theorem (FNW)
needs a **volume-uniform MPS bond dimension**; the 2D/3D uniqueness theorem (SCV) needs volume-uniform
PEPS. The MERA's constant per-layer `χ` does **not** supply either. A rigorous certificate therefore
needs `χ` growing with `L` — the "area law with bounded entanglement" statement one wants to *prove*
for YM, circular as an assumption (same `[O]` as X3). **Kill #2 met.**

---

## 5 · 2D solvable gauge model: MERA carries topology, not the gap  [T]/[V]/[NG]

The 2D toric code = 2D `Z₂` gauge fixed point. **Exact** facts (X3, re-quoted): `E₀=−8`,
`GSD=4` for all `L`, many-body gap `4`, anyon mass `2`, levels `{−8,−4,0,4,8}` — all from the
**stabilizer/Hamiltonian spectrum**. The **entanglement-spectrum** route (Aguado–Vidal 2008): the
toric code admits an **exact finite-bond MERA**; the four ground states are the **same network
differing only in the top tensor**, which stores the topological (logical) data. Therefore:

- the MERA / ES exposes the **topological degeneracy** `4 = 2^{2g}`, the **topological entanglement
  entropy** `γ = log 2` (Levin–Wen; Kitaev–Preskill — rigorous for the topological class), and the
  area law;
- it does **not** expose the **gap magnitude** `4`. The anyon mass `2` is read from the Hamiltonian,
  never from the MERA/ES.

This sharpens X3: there the point was "the toric-code GS is non-injective (`GSD=4>1`), so the
injectivity criterion does not fire." Here the point is stronger and independent of injectivity:
**even granting an exact MERA, its entanglement-spectrum content is *topological*, not the mass gap.**
The Li–Haldane "entanglement gap" is likewise a *fingerprint of topological order* (observed
numerically to "apparently remain finite," not a theorem) — it separates edge from bulk, not gapped
from gapless in a way that yields a gap *magnitude*.

---

## 6 · Why this is *not* "MERA reduces to H1" (kill #3 is not the one that fires)  [T]

H1 (`H1-RG-DESCENT.md`) is a **Wilsonian block RG**: factor-`L` decimation, the exact Schur backbone
`G″_j=(P A_j^{-1}P^T)^{-1}=E[Σ]−D`, and the obstruction is the **Savvidy soft mode** where the fibre
gap `γ_j→0` and `η_j≥1`. MERA's RG is an **entanglement RG**: real-space disentanglers + isometries,
fixed point = the **conformal** (critical) fixed point, and its diagnostic eigenvalues are **CFT
scaling dimensions**, not Hessian curvatures. The two are **different RGs on different objects**;
there is no identification of the MERA "scaling superoperator gap" with H1's `G″`. So the kill
criterion **"reduces to H1" is NOT met** — MERA is a genuinely distinct corridor. But distinctness
does not help: MERA fails earlier and harder, on the acid test (§3) and on kills #1/#2 (§2, §4). It
does not even reach the H1 wall.

---

## 7 · What would be needed to escape `[O]`

A theorem of the shape:
> *the interacting 4D Yang–Mills ground state, on lattices of growing linear size, admits a family of
> MERA/tensor networks with a **volume-uniform** structure whose **entanglement/scaling gap lower-bounds
> the physical gap uniformly** — and `H_YM` is in the same local-spectral class.*

Both halves fail: (i) a MERA's uniform per-layer `χ` does not give a volume-uniform *rigorous* bond
dimension (§4); (ii) no "entanglement gap ⇒ physical gap" theorem exists in `D≥2` (open, as in X3;
the solvable model shows the ES carries topology, not the gap, §5). This is one more name for the
charter's missing machine — a non-perturbative gap transport with a **volume-uniform certified**
structure — and MERA is the wrong machinery for it, being built for the gapless fixed point.

### 7.1 What is banked `[T]`

- MERA is variational (`E_MERA ≥ E₀`) with **no** parent-Hamiltonian gap theorem; the FNW theorem is
  MPS-only (1D), SCV is PEPS-uniqueness-only.
- The ES gap of a 1D state is a **transfer-matrix / correlation-length quantity**: exact free-boson
  monotonicity in `μ=−log|λ₂|` and the TFI closure at `g=1` (§2).
- **Scale-invariant MERA = critical = gapless**; a gapped state's MERA is trivial depth `O(log ξ)`;
  the free massless field has a scale-invariant MERA (acid-test fail for the slogan).
- MERA bond dimension maps to MPS bond dimension `L^{log_b χ}` — polynomial, not volume-uniform (§4).
- Exact finite-bond MERA of the toric code stores the **topological degeneracy**, not the gap
  magnitude (Aguado–Vidal; `GSD=4`, gap `4`).

---

## 8 · Status ledger

| item | statement | label |
|---|---|---|
| MERA variational (`E_MERA ≥ E₀`); no parent-H gap theorem for MERA | §1 | **[T]** |
| "injective ⇒ uniformly gapped parent H" is 1D **MPS** only (FNW); PEPS = uniqueness only | §1, X3 | **[T]/[O]** |
| scale-invariant MERA = critical (gapless) fixed point; gapped MERA = trivial `O(log ξ)` | §1, §3 | **[T]** |
| free massless field has a scale-invariant MERA ⇒ "scale-inv MERA ⇒ gap" fails acid test | §3 | **[T]/[NG]** |
| ES gap of free scalar = function of `μ=−log|λ₂|`; massive bounded, massless `→0` | §2.1, M1 | **[T]/[V]** |
| TFI: gap `∝1/L→0` at `g=1` (massless), `→2(g−1)` for `g>1`; ES gap tracks it | §2.2, M2 | **[V]** |
| 1D ES/MERA gap **is** the transfer-matrix gap | §2, §0.2 | **[NG]** (kill #1) |
| MERA bond dim `χ →` MPS bond dim `L^{log_b χ}`; not volume-uniform | §4, M4 | **[T]/[NG]** (kill #2) |
| toric code exact MERA stores topological degeneracy (`GSD=4`), not the gap magnitude `4` | §5, M5 | **[T]/[V]/[NG]** |
| MERA entanglement RG ≠ H1 Wilsonian block RG; does **not** reduce to H1 | §6 | **[T]** (kill #3 not met) |
| injective/MERA ⇒ gap for the interacting 4D YM GS with volume-uniform structure | §7 | **[O]** |
| Clay mass gap | unbuilt | **0 %** |

---

## 9 · Reproduce

```
py -3.12 -I research/v3/w7/z4_mera.py      # ~30 s; writes research/v3/w7/z4_mera_results.json
```

Legs: `[M1]` free-scalar transfer-matrix vs free-boson half-chain ES (covariance) · `[M2]` TFI exact
ED gap + half-chain ES gap · `[M3]` scale-invariance/acid-test ledger · `[M4]` MERA→MPS bond-dim
blow-up · `[M5]` toric-code MERA content. numpy 2.4.6 / scipy 1.18.1; one process; no downloads.

*Refs (bibliographic spot-checks done 2026-11-24; standard-citation ones marked).*
G. Vidal, "Entanglement renormalization," *Phys. Rev. Lett.* **99**, 220405 (2007) — **VERIFIED**.
G. Vidal, "Class of quantum many-body states that can be efficiently simulated," *Phys. Rev. Lett.*
**101**, 110501 (2008) — standard. G. Evenbly, G. Vidal, "Tensor network renormalization,"
*Phys. Rev. Lett.* **115**, 180405 (2015), DOI 10.1103/PhysRevLett.115.180405 — **VERIFIED**. H. Li,
F. D. M. Haldane, "Entanglement spectrum as a generalization of entanglement entropy,"
*Phys. Rev. Lett.* **101**, 010504 (2008), arXiv:0805.0332 — **VERIFIED**. M. Levin, X.-G. Wen,
"Detecting topological order in a ground state wave function," *Phys. Rev. Lett.* **96**, 110405
(2006), arXiv:cond-mat/0510613 — standard. A. Kitaev, J. Preskill, "Topological entanglement
entropy," *Phys. Rev. Lett.* **96**, 110404 (2006) — standard. M. Aguado, G. Vidal, "Entanglement
renormalization and topological order," *Phys. Rev. Lett.* **100**, 070404 (2008), arXiv:0712.0348 —
**VERIFIED**. Fannes–Nachtergaele–Werner, *Commun. Math. Phys.* **144**, 443 (1992);
Perez-García–Verstraete–Wolf–Cirac, *Quantum Inf. Comput.* **7**, 401 (2007);
Schuch–Cirac–Verstraete, *Phys. Rev. Lett.* **105**, 040503 (2010) — via X3. Internal:
`X3-TENSOR.md`, `H1-RG-DESCENT.md`, `D3-TRANSFER-CONSTANT.md`. *Witness* `research/v3/w7/z4_mera.py`,
data `z4_mera_results.json`. No `[C]` upgraded to `[T]`; **no QED; Clay 0%.**
