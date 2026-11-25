# X3-TENSOR — the tensor-network gap certificate: a genuine **1D theorem** that *is* the transfer-matrix gap, and a **triple no-go** in 2D/3D/4D (injectivity ⇒ gap is open; the solvable gauge model's gap is exact-solution; L=2 fails the free-massless acid test)

> **2026-11-24 · First step of the tensor-network / rigorous gap-certificate paradigm (X3).**
> Build an MPS/PEPS gap certificate for a **solvable** gauge model; define the injectivity/area-law
> criterion; acid-test it on free fields (massless must fail, massive must pass); assess 4D L=2.
> **Exact / solvable only** — no Monte Carlo, no fit. Witness `research/v3/w7/x3_tensor.py` →
> `x3_tensor_results.json` (`py -3.12 -I`, numpy; **one process**; no downloads).
> Reads: `M1-4D-RUN.md`, `WG-L12-L16.md`, `D3-TRANSFER-CONSTANT.md`, `docs/paper-v2/ch12-numerics.tex`.
> **Labels (charter §2).** **[T]** true theorem/identity · **[V]** exact numerical witness (this
> file) · **[C]** conditional/protocol · **[NG]** no-go for a named class · **[O]** exact missing
> step. **No QED. Clay 0%.**

---

## 0 · Verdict (up front; no QED, Clay 0%)

**The kill criterion is met on every route except 1D, where the "certificate" is literally the
transfer-matrix gap.** Concretely:

1. **[T] The one genuine theorem is 1D, and it reduces to the transfer matrix.** An *injective*
   MPS (bond dim `D`, injectivity length `ℓ`, both `L`-independent) has a parent Hamiltonian that
   is frustration-free, has a **unique** ground state, and a **uniform spectral gap** `Δ ≥ f(D,ℓ,δ_T) > 0`
   independent of system size — where `δ_T = 1−|λ₂(T)|` is the **transfer-matrix spectral gap** and
   `λ₂` its second eigenvalue (Fannes–Nachtergaele–Werner 1992; Perez-García–Verstraete–Wolf–Cirac
   QIC 2007). The correlation length is `ξ = −1/log|λ₂|`. **So "MPS gap certificate" ≡ transfer-matrix
   gap `1−|λ₂|`** — exactly D3-TRANSFER-CONSTANT's object, not a new one, and 1D is not Yang–Mills.

2. **[T]/[V] The free-field acid test passes — but only via the transfer matrix, not the area law.**
   Exact: the transfer-matrix ratio is `|λ₂/λ₀| = e^{−μ}`, `cosh μ = 1 + m²/2`. Massive (`m>0`):
   `δ_T = 1−e^{−μ} > 0`, `L`-independent → **PASS** (`δ_T = 0.618, 0.390, 0.181` for `m=1,0.5,0.2`,
   confirmed to all digits by the open-chain Green's function ratio `0.381966, 0.609612, 0.819002`).
   Massless (`m=0`): `|λ₂/λ₀| = 1`, `δ_T = 0`, `ξ = ∞` → **FAIL** (correctly). **But in `D≥2` the
   free massless field *satisfies the area law* (the entanglement ``log violation'' is a 1D
   phenomenon), so the area-law criterion is acid-blind at `D≥2` — it cannot separate massive from
   massless, and is therefore not a gap certificate there.** `[T]`

3. **[NG]/[T] The solvable 2D gauge model has an exact gap and a *non-injective* tensor-network
   ground state — the injectivity criterion cannot fire.** The 2D toric code = **2D Z₂ lattice gauge
   theory at its fixed point**. Exact stabilizer algebra (witness, dense `L=2`): `E₀ = −8`,
   **GSD = 4** on the torus (`=2^{2g}`), **many-body gap = 4** (anyon mass 2), levels
   `{−8,−4,0,4,8}`. `GSD = 4 > 1` for every `L` (GF(2) rank `= 2L²−2`) ⟹ the ground-state PEPS is
   **not injective** (an injective PEPS has a unique ground state). So the gap is certified by the
   **exact stabilizer spectrum**, not by a tensor network. **KILL: reduces to the exact/transfer-matrix
   spectrum.** `[T]/[NG]` (Same for 2D Wegner Z₂ gauge = 2D Ising duality, 3D Wegner = 3D Ising:
   the gap is the dual model's gap, known exactly.)

4. **[NG]/[O] In 2D/3D injectivity ⇒ gap is *not a theorem*.** Injective PEPS ⟹ *unique* parent
   ground state (Schuch–Cirac–Verstraete, PRL 2010) — **uniqueness only**. That the parent
   Hamiltonian is **gapped** in 2D/3D is the open **"PEPS gap conjecture"** `[O]`. There is no
   injectivity⇒gap theorem in the dimensions where YM lives.

5. **[V]/[NG] 4D L=2 is computationally feasible but certifyingly vacuous.** Exact 4D Z₂-gauge
   stabilizer on `L=2`: `n_links = 64`, **GSD = 16** (GF(2) rank `= 60`), instant to compute — so an
   exact finite-volume spectrum is reachable. **But it is the transfer-matrix gap, and `L=2` fails
   the free-massless acid test:** the free **massless** finite-volume gap at `L=2` is already
   `2 sin(π/2) = 2.0` in **every** `D` (witness, leg 1) — an `O(1)` "gap" at `L=2` is exactly what a
   *massless* field also produces. Any `L=2` gap is indistinguishable from a free massless artifact.

**One line.** *The tensor-network gap certificate is a real theorem only in 1D, where it **is** the
transfer-matrix gap `1−|λ₂|` (D3's object, not new); in 2D/3D it splits into the exact gap of a
solvable gauge model (toric code: GSD=4, gap=4, **non-injective** — the criterion doesn't fire) or the
open PEPS gap conjecture (no theorem); and at 4D `L=2` it is feasible but vacuous because the free
massless field already has gap `2.0` at `L=2`. **Kill criterion met.** The missing step `[O]` is a
2D/3D injectivity⇒gap theorem for an interacting ground state with `L`-independent bond dimension —
which does not exist and is one more name for the missing machine. No QED. Clay 0%.*

---

## 1 · The criterion, stated precisely  [T]/[C]

**Injectivity.** An MPS `|Ψ⟩ = Σ tr(A^{s₁}…A^{s_L}) |s₁…s_L⟩` with bond dimension `D` is *injective*
("normal") if there is an *injectivity length* `ℓ` such that the linear map
`(C^D)^{⊗ℓ} → (C^d)^{⊗ℓ}`, `(X_1,…,X_ℓ) ↦ Σ tr(X_1…X_ℓ) |s_1…s_ℓ⟩`, is injective. In words: the `D^ℓ`
boundary indices faithfully encode the `d^ℓ` physical configurations of any `ℓ` consecutive sites.
(For PEPS the definition is the same with a 2D patch.)

**Area law.** A tensor-network state of bond dimension `D` has entanglement entropy across any cut
`S(∂) ≤ c·|∂|·log D` — an area law *by construction*. This is `[T]` but it is the *easy* direction.

**The certificate attempted.** *State `|Ψ⟩` carries a TN gap certificate if it admits a family of
injective MPS/PEPS representations with `(D, ℓ)` bounded uniformly in system size and transfer-matrix
spectral gap `δ_T = 1−|λ₂/λ₀| ≥ δ₀ > 0` uniformly in system size.* Then one wants to conclude the
parent Hamiltonian is gapped.

**What is a theorem, precisely.**

| statement | dimension | status |
|---|---|---|
| injective MPS ⇒ unique, frustration-free parent GS with **uniform gap** `Δ ≥ f(D,ℓ,δ_T)>0`; `ξ = −1/log|λ₂|` | 1D | **[T]** (FNW '92; PGWV '07) |
| injective PEPS ⇒ **unique** parent GS | any `D` | **[T]** (SCV, PRL 105, 040503, 2010) |
| TN state of bond dim `D` ⇒ area law `S ≤ c·|∂|·log D` | any `D` | **[T]** (by construction) |
| injective PEPS ⇒ **gapped** parent H | `D≥2` | **[O]** — open "PEPS gap conjecture" |
| the gap of the parent H is the gap of `H_YM` | any | **[O]** — needs `H_YM` frustration-free or a spectral-stability comparison |

The last two rows are the whole difficulty. Injectivity gives *uniqueness* for free in `D≥2` but the
*gap* only in 1D; and even in 1D the theorem is about the **parent** Hamiltonian, which is not
`H_YM` unless one proves `H_YM` is (unitarily/local-spectrally equivalent to) the parent — a
quasi-adiabatic/local-perturbation-stability step that is itself only rigorous in restricted cases.

---

## 2 · Acid test on free fields — exact  [T]/[V]

### 2.1 The transfer-matrix criterion (exact), and it is the gap

For a free real scalar of mass `m` with nearest-neighbour lattice operator `−Δ + m²`, the 1D
propagator is **exactly geometric**: `C(r) ∝ ρ^{r}` with `ρ = |λ₂/λ₀|` the smaller root of
`ρ² − (2+m²)ρ + 1 = 0`, i.e. `ρ = e^{−μ}`, `cosh μ = 1 + m²/2`. So
```
     δ_T = 1 − e^{−μ},     ξ = 1/μ = −1/log|λ₂| .            (exact)
```

| `m` | `μ = 1/ξ` | `|λ₂/λ₀| = e^{−μ}` | **`δ_T`** | `ξ` | open-chain ratio (`N=301, r=8`) |
|---|---|---|---|---|---|
| 1.00 | 0.96242 | 0.381966 | **0.61803** | 1.039 | 0.381966 |
| 0.50 | 0.49493 | 0.609612 | **0.39039** | 2.02 | 0.609612 |
| 0.20 | 0.19967 | 0.819002 | **0.18100** | 5.008 | 0.819002 |
| 0.00 | 0.00000 | 1.000000 | **0.00000** | ∞ | 0.993 (no decay) |

The open-chain column is an *independent* check (invert `−Δ+m²` on the open chain, read the bulk
Green's function ratio) — it reproduces `e^{−μ}` to all displayed digits. **Massive passes
(`δ_T>0`, `L`-independent); massless fails (`δ_T=0`, `ξ=∞`).** `[T]/[V]`

**Honest identification.** `δ_T` *is* the transfer-matrix gap of D3-TRANSFER-CONSTANT. The "MPS gap
certificate" in 1D does not introduce a second quantity — it certifies exactly this number and
nothing more. **Kill criterion: 1D certificate = transfer-matrix gap.** `[NG]`

### 2.2 The area law is acid-blind at `D≥2`  [T]

The free **massless** field:
- **1D**: `S(ℓ) = ⅓ log ℓ + c` — *violates* the area law, so the 1D area-law criterion catches it. ✓
- **`D≥2`**: leading entropy is **area law** `S ∝ |∂|` (the "log violation" is purely 1D; in `D≥2`
  there is only a universal *subleading* constant). The massless field **satisfies the area law**.

Therefore in `D≥2` the area-law criterion **cannot distinguish massive from massless** — both are
area-law. The only distinguisher is the correlation length / transfer-matrix gap, which is exactly
`δ_T`, which the area law does not measure. So the *soft* certificate (area law) fails the acid test's
discriminating power at `D≥2`; the *hard* certificate (`δ_T`) reduces to the transfer matrix. Either
way, no new gap content in `D≥2`. `[T]/[NG]`

### 2.3 The finite-volume-gap floor `2 sin(π/L)` — a reusable filter  [T]

The free-field finite-volume gap (single-particle, `D`-independent because one direction is excited):
```
     ω_min(D, L, m) = sqrt( m² + 4 sin²(π/L) ),     massless:  ω_min = 2 sin(π/L) .   (exact)
```

| `D` | `L` | `m=0` | `m=0.5` | `m=1.0` |
|---|---|---|---|---|
| 1–4 | 2 | **2.0000** | 2.0616 | 2.2361 |
| 1–4 | 4 | 1.4142 | 1.5000 | 1.7321 |
| 1–4 | 8 | 0.7654 | 0.9142 | 1.2593 |
| 1–4 | 16 | 0.3902 | 0.6342 | 1.0734 |

(The value is identical for `D=1,2,3,4` — the minimum excites one direction by one unit.) **This is a
reusable acid filter:** the free-massless floor is `2 sin(π/L) → 0`. A gap `Δ(L)` measured on volume
`L` **certifies nothing** unless `Δ(L)` is bounded below uniformly in `L` while `2 sin(π/L) → 0`; and
at **`L=2` the floor is `2.0`**, the *hardest* single volume to beat. `[T]`

---

## 3 · The solvable gauge model: 2D toric code = 2D Z₂ gauge fixed point  [T]/[V]

The canonical solvable 2D gauge model is the **toric code** — `Z₂` lattice gauge theory at its
deconfined fixed point, `H = −Σ_v A_v − Σ_p B_p` (vertex `X`-stars, plaquette `Z`-faces). All
stabilizers commute; the spectrum is **exact**.

**Exact dense diagonalization, `L=2` (8 qubits, 256-dim):**
```
  n_qubits=8, stabilizers=8, GF(2) rank=6  ->  GSD = 2^{8-6} = 4
  E0 = -8.0   GSD_dense = 4   E1 = -4.0   many-body gap = 4.0
  unique levels: {-8, -4, 0, +4, +8}
```
Why the gap is `4` (anyon mass `2`): each edge is shared by exactly two plaquettes, so `∏_p B_p = 1`
and `∏_v A_v = 1` — a single plaquette/vertex defect is forbidden; defects come in pairs, each
costing `2`, so the first excited *state* is `+4` above the ground level (the single-defect "mass" `2`
is half of it). This is exact for all `L`.

**GSD `= 4` for every `L`** (GF(2) rank `= 2L²−2`, so `2^{2L²−(2L²−2)} = 4`):

| `L` | qubits | GF(2) rank | GSD |
|---|---|---|---|
| 2 | 8 | 6 | 4 |
| 3 | 18 | 16 | 4 |
| 4 | 32 | 30 | 4 |
| 5 | 50 | 48 | 4 |

**The decisive reading.** `GSD = 4 > 1` on the torus. An **injective** PEPS has a **unique** ground
state (SCV theorem), so a `GSD>1` state is **not injective** — the ground-state PEPS of the toric code
is *non-injective* (topological order). Therefore **the injectivity gap criterion never fires** on
the one solvable gauge model, even though the model is provably gapped. The certified gap comes
**entirely from the exact stabilizer algebra** — i.e. the exact/transfer-matrix spectrum. **KILL
criterion met:** for the solvable 2D gauge model, the TN adds nothing to exact solvability. `[T]/[NG]`

(The confined phase — 2D Wegner `Z₂` gauge = 2D Ising duality — is the same story: the gap is the
dual Ising gap, exact; and 3D Wegner = 3D Ising, again exact by duality `[C]`.)

---

## 4 · 4D `L=2` feasibility  [V]/[NG]

**Exact, feasible, and vacuous.** The 4D Z₂-gauge stabilizer (4D toric code) on `L=2`:

| model | `n_links` | stabilizers | GF(2) rank | GSD |
|---|---|---|---|---|
| `D=2, L=2` | 8 | 8 | 6 | 4 |
| `D=2, L=3` | 18 | 18 | 16 | 4 |
| `D=3, L=2` | 24 | 32 | 21 | 8 |
| `D=3, L=3` | 81 | 108 | 78 | 8 |
| **`D=4, L=2`** | **64** | **112** | **60** | **16** |

`4D L=2` has `64` qubits; its exact GSD `= 16` is a GF(2) rank of a `112×128` binary matrix over `64`
bits — **instant**. So an exact finite-volume spectrum at `4D L=2` is computationally **feasible**
(for SU(2) the analogous truncated-link PEPS contraction is also feasible at small bond dimension,
since the `L=2` network is small). Observed pattern `GSD = 2^D` (computed `4, 8, 16` for `D=2,3,4`).

**Why it is nevertheless certifyingly vacuous:**

1. It is the **finite-volume = transfer-matrix gap** (the `L=2` Hamiltonian's spectrum) — D3's object.
2. It **fails the free-massless acid test.** The free *massless* finite-volume gap at `L=2` is
   `2 sin(π/2) = 2.0` in **every** `D` (§2.3). An `O(1)` gap at `L=2` is precisely what a massless
   field produces; one cannot tell `m_phys>0` from `m=0` on `L=2`. The floor `2 sin(π/L)` is
   *largest* at `L=2`, so `L=2` is the *least* informative volume for a gap claim.

Hence **`4D L=2` is feasible but cannot host a certificate**: it reduces to the finite-volume
transfer-matrix gap and fails the acid test. `[V]/[NG]`

---

## 5 · What this means, and the exact missing step  [NG]/[O]

**Triple no-go for the class "tensor-network gap certificate for YM":**
- **1D `[NG]`:** the certificate **is** the transfer-matrix gap `1−|λ₂|` (D3's object) — no new gap
  content, and 1D is not the physical dimension.
- **2D/3D `[NG]/[O]`:** injectivity gives *uniqueness*, not a gap (open PEPS gap conjecture); for
  the solvable gauge model the gap is the exact stabilizer/duality spectrum and the ground state is
  *non-injective* (GSD `= 4`), so the criterion does not fire.
- **4D `L=2` `[NG]`:** feasible, but the finite-volume gap fails the free-massless acid test
  (massless already has gap `2.0` at `L=2`).

**What it would take to escape `[O]`.** A theorem of the shape:
> *the interacting 4D Yang–Mills ground state, on lattices of growing linear size, admits a family of
> injective PEPS with bond dimension `D` **bounded uniformly in volume** and transfer-matrix gap
> `δ₀>0` uniform in volume — **and** `H_YM` is in the same local-spectral class as the parent.*

Neither half exists: the 2D/3D injectivity⇒gap implication is open, and "YM ground state = uniform-`D`
injective PEPS" is exactly the area-law-with-bounded-entanglement statement one wants to *prove* for
YM (circular as an assumption). This is one more name for the missing machine of the charter's
frontier — a **non-perturbative transport of the gap with a volume-uniform tensor-network structure
that one can certify, not fit.**

### 5.1 What is banked `[T]`

- The exact transfer-matrix identity `δ_T = 1 − e^{−μ}`, `cosh μ = 1+m²/2`, and its open-chain
  confirmation (massive passes, massless fails).
- The exact finite-volume-gap formula `ω_min = sqrt(m² + 4 sin²(π/L))` and the **free-massless floor
  `2 sin(π/L)`** — a reusable falsification filter: any claimed gap `Δ(L)` at volume `L` is
  certifying only if it is bounded below uniformly while `2 sin(π/L) → 0`. `L=2` (floor `2.0`) is the
  weakest volume.
- The exact solvable-model data: toric code `GSD=4`, many-body gap `4` (anyon mass `2`), levels
  `{−8,−4,0,4,8}`; `4D L=2` Z₂ gauge `GSD=16`. These are exact finite-volume facts, **[V]**, not a
  continuum statement.

---

## 6 · Status ledger

| item | statement | label |
|---|---|---|
| injective MPS ⇒ unique, uniformly-gapped parent H; `ξ=−1/log|λ₂|` | §1 (FNW '92; PGWV '07) | **[T]** |
| injective PEPS ⇒ unique parent GS | §1 (SCV 2010) | **[T]** |
| injective PEPS ⇒ **gapped** parent H | §1 — 2D/3D open | **[O]** |
| transfer-matrix ratio `|λ₂/λ₀| = e^{−μ}`, `δ_T = 1−e^{−μ}`; massive PASS, massless FAIL | §2.1 | **[T]/[V]** |
| 1D MPS certificate ≡ transfer-matrix gap (`1−|λ₂|`) | §2.1 / §5 | **[NG]** |
| `D≥2` free massless satisfies area law ⇒ area-law criterion acid-blind | §2.2 | **[T]/[NG]** |
| finite-volume-gap floor `2 sin(π/L)→0`; `L=2` floor `=2.0` | §2.3 | **[T]** |
| 2D toric code (Z₂ gauge fixed point): `GSD=4`, gap `4` (anyon mass `2`), levels `{−8,−4,0,4,8}` | §3 | **[T]/[V]** |
| toric-code GS non-injective (`GSD=4>1`) ⇒ injectivity criterion does not fire; gap = exact spectrum | §3 | **[T]/[NG]** |
| `4D L=2` Z₂ gauge: `GSD=16` (GF(2) rank 60), feasible | §4 | **[V]** |
| `4D L=2` gap = finite-volume/transfer-matrix gap; fails free-massless acid test (`2.0`) | §4 | **[V]/[NG]** |
| 2D/3D injectivity⇒gap theorem for the YM ground state with volume-uniform `D` | §5 | **[O]** |
| Clay mass gap | unbuilt | **0 %** |

---

## 7 · Reproduce

```
py -3.12 -I research/v3/w7/x3_tensor.py     # ~2 s; writes research/v3/w7/x3_tensor_results.json
```

*Refs.* Theorem set: Fannes–Nachtergaele–Werner, Commun. Math. Phys. 144, 443 (1992);
Perez-García–Verstraete–Wolf–Cirac, Quantum Inf. Comput. 7, 401 (2007); Schuch–Cirac–Verstraete,
Phys. Rev. Lett. 105, 040503 (2010) (injective PEPS ⇒ unique parent GS; gap open in 2D). Transfer
matrix / constant-mode gap in this project: `D3-TRANSFER-CONSTANT.md`. Solvable model: Kitaev toric
code (2D Z₂ gauge fixed point). Numerical frame of the campaign: `M1-4D-RUN.md`, `WG-L12-L16.md`,
`docs/paper-v2/ch12-numerics.tex` (finite-volume facts ≠ continuum atom; `m₀≳π` unresolved at `L=8`).
*Witness* `research/v3/w7/x3_tensor.py`, data `x3_tensor_results.json` (exact/solvable only, numpy,
one process, no downloads). No `[C]` upgraded to `[T]`; **no QED; Clay 0%.**
