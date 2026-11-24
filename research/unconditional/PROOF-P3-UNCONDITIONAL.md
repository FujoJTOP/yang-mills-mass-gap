# PROOF-P3-UNCONDITIONAL (rev. 2026-11-24) — finite-`L` strong-coupling residue theorem (I1), corrected & assembled

> **Revision history (charter §2).** The 2026-11-23 version claimed an *unconditional* residue theorem via a "KP lemma written out in full / the lemma closes." That was **over-claimed**: two rounds of external audit (`P3-AUDIT-CORRECTED.md` + a second round) found five substantive gaps — **(1)** minimal-tube geometry (factor-4: four walls, `4t/a`, not one plaquette per slice `t/a`), **(2)** KP proved only on the fundamental sector, **(3)** dressed-mass resummation asserted not proved, **(4)** the theorem lived on `T⁴_L` (periodic time) where `lim_{t→∞}` does not exist, **(5)** the polymer→spectral-atom bridge missing. 2026-11-24: four parallel rigorous re-derivations **close all five**:
> - `P3-CYLINDER-RESIDUE.md` — cylinder `T³_L×ℝ`, discrete spectrum, residue definition, spectral bridge **(gaps 4, 5)**;
> - `P3-FOURWALL-DERIVATION.md` — four-wall minimal surface `|P|≥4t/a`, `m₁=4ln(1/u)/a` **(gap 1)**;
> - `P3-FULLREP-KP.md` — full-representation KP, `K≤CN` **(gap 2)**;
> - `P3-DRESSED-MASS.md` — backbone–decoration resummation, effective mass, amplitude floor **(gap 3)**.
> This file is now the corrected **assembly**, not the (wrong) 11-23 claim.

## Theorem (corrected, honest scope)
Fix `a>0`, `N`, `β₁=min(1/48, c/N)`. For `β<β₁`, on the **cylinder** `𝒞_L=T³_L×ℝ` (spatial torus, infinite time), for every spatial volume `L`, the lightest `0⁺⁺` atom of `O=¼tr F²` obeys
```
  f₁²_L := ν_L({E_{1,L}}) = lim_{t→∞} e^{E_{1,L}t} C_L(t)  ≥  r₀(β) > 0,
```
**uniformly in `L`** (limit `t→∞` at fixed `L` first; `L→∞` is not taken — `L` is only the uniformity axis). This is the **residue (I1)**.
**Label: `[T]` modulo the standard OS/Källén–Lehmann transfer-matrix representation on the cylinder.** **NOT obtained:** isolation `I2`, continuum `a→0`, transport to `β≈2.5`.

## §1 Cylinder, residue, spectral bridge — `[T]` (mod OS/KL) [`P3-CYLINDER-RESIDUE`]
`T_L=e^{−aH_L}` RP, self-adjoint, `0≤T_L≤1` on `H_phys,L`. At finite `L`: continuous kernel on compact ⇒ Hilbert–Schmidt ⇒ compact ⇒ **discrete spectrum** `0=E_0<E_1<E_2≤…`. Then `C_L(t)=Σ_{n≥1}f_{n,L}²e^{−E_{n,L}t}`, `f²=|⟨n|O|Ω⟩|²≥0`; `f₁²_L=ν_L({E_{1,L}})=lim_t e^{E_{1,L}t}C_L(t)`. **`T⁴_L` is ill-posed** (periodic propagator bounded ⇒ limit does not exist) — hence the cylinder. **Isolation `I2` is NOT needed for the finite-`L` floor** (compactness isolates every atom); it is needed only for the infinite-volume "atom vs continuum threshold" and `f₁²_∞` — the open rung. The bridge (polymer leading weight `↔` spectral atom) is `[T]` modulo R1+R2 (§3,§4). The one imported line: RP ⇒ `H≥0` self-adjoint, `C(t)=∫e^{−Et}dμ_O`, `dμ_O≥0` (OS'78; Källén–Lehmann; Reed–Simon IV §XIII).

## §2 Four-wall geometry — `[T]` [`P3-FOURWALL-DERIVATION`]
The minimal connected, gauge-invariant, nonvanishing polymer linking the two single-plaquette insertions at times `0,t` has `|P| ≥ 4t/a`: **(1.1/1.1′)** each of the source plaquette's 4 boundary links carries an open fundamental pair, `∫dU U=0`, two plaquettes share `≤1` link ⇒ `≥4` flux strands; **(1.2)** only temporal plaquettes advance time, so each strand crosses `t/a` slices; **(1.3)** `C(t)_c` links source↔sink, so each strand reaches the sink. **Equality** = the straight **four-walled tube** (4 plaquettes/slice), weight `a₀u^{4t/a}=a₀e^{−m₁t}`, **`m₁=4ln(1/u)/a`**; `a₀>0` (positive `c_f` × singlet recoupling). Every non-tube costs `≥4t/a+1` (`O(u)` relative) or is rate-separated (`δ=2ln(1/u)/a ≤` the `≈1.5m₁` excited `0⁺⁺`). **SU(2) check:** `u=β/2N²` ⇒ `4ln(1/u)=−4ln(β/2N²)=M_{0⁺⁺}` (Gupta §14.2.2, Eqs 14.17–14.21, verified).

## §3 Full-representation KP — `[T]` [`P3-FULLREP-KP`]
`Σ_{P∋x}|w(P)|e^{a|P|}≤a` over **all** irreps: Lemma `c_r(β)≤β^{|r|}/|r|!·e^β/dim r` (Taylor + multiplicity cap `m_r(a,b)≤N^{a+b}/dim r`, factorial cutoff `|r|≤a+b`), and `|χ_r(U)|≤dim r` pointwise (Haar is a probability measure). Product: `c_r·dim r ≤ β^{|r|}/|r|!e^β` — the offending `dim r` from contraction is **cancelled** by `1/dim r` in `c_r`; `1/|r|!` beats the sub-exponential rep entropy. Hence `K:=2λ_aC_g e ≤ 80e^{2+β}N ≤ 6.1×10² N`, adjacency `Δ=4·(6−1)=20`. (`c_r≥0` for all `r` is a by-product.)

## §4 Dressed-mass resummation — `[T]` (mod OS/KL) [`P3-DRESSED-MASS`]
Backbone–decoration (skeleton-polymer) expansion of the four-wall tube: **(a)** decoration gas resums (KP, `t`-independent radius); **(b)** effective `L`-uniform mass `m₁^dressed=4ln(1/u)/a·(1−O(u/ln(1/u)))`, `L`-defect `O(e^{−cL})`; **(c)** stays a spectral atom (atom-from-exponential lemma; input = §1 KL/OS); **(d)** uniform amplitude floor `inf_L A^dressed ≥ ¼(2N)^{-2}>0`. Free-field acid test passes (`a₀=0`).

## §5 Assembly & honest boundary
`f₁²_L ≥ r₀(β)=½c₁u^p>0`, uniform in `L`, on `𝒞_L`, for `β<β₁` — **`[T]` modulo the standard OS/KL representation** (§1). **This is a finite-`L` strong-coupling residue (I1) theorem — the corrected version of the 11-23 over-claim.**
- **Does NOT** give: isolation `I2`, the strong-coupling *gap* itself (that is Osterwalder–Seiler '78 / SZZ '23, not P3), the continuum `a→0`, or the crossover `β→2.5` (the analyticity wall `1/48`, a *method* limit).
- The polymer expansion controls a **weight** (I1), not **eigenvalue spacing** (I2) — hence "residue theorem, not gap theorem."

*Refs.* Gupta arXiv:hep-lat/9807028 §14.2.2; Münster PLB 95 (1980)/NPB 190 (1981); Seiler LNP 159 (1982) §II; Osterwalder–Seiler Ann. Phys. 110 (1978) 440; Shen–Zhu–Zhu CMP 400 (2023) 805; Källén (1952)/Lehmann (1954); Reed–Simon IV §XIII; Kotecký–Preiss CMP 103 (1986) 491–498. Internal: `P3-CYLINDER-RESIDUE`, `P3-FOURWALL-DERIVATION`, `P3-FULLREP-KP`, `P3-DRESSED-MASS`, `P3-AUDIT-CORRECTED`.
*Status: **corrected finite-`L` residue (I1) theorem, `[T]` mod OS/KL.** I2 / `a→0` / `β→2.5` open. Clay unaffected. No QED on anything unproven.*
