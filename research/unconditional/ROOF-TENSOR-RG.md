# ROOF-TENSOR-RG — tensor-network RG / MERA / entanglement Hamiltonian as a strut for the `0⁺⁺` gap

> **VERDICT (no QED).** Tensor-network RG (MERA / PEPS / cMPS / DMRG) is **renormalon-free** — a genuine
> plus over the gluon-condensate route — but it cannot be the **strut**: it is **variational**, so it yields
> **upper** bounds on eigenvalues and hence **cannot lower-bound** the gap (the same upper/lower asymmetry
> that killed the base-round non-local bound). Its "gap" object is the **entanglement-spectrum gap**, which
> is **not** the physical `0⁺⁺` mass gap (demonstrated on the TFIM: the two gaps are not proportional and
> do not vanish together). The rigorous content lives in **1+1D** (area law, gap ⇒ exponential decay); in
> **3+1D there is no rigorous MERA→gap theorem** for interacting gauge theory, and the 4D area law is
> itself unproven (circularity). **It is the upper half of a sandwich, not a strut.**

---

## §0 · Labels and the acid test

**[T]** 真定理 · **[C]** 条件定理 · **[O]** 精确缺失步 · **[NG]** no-go for the named class · **[V]** new this round.
Acid (project convention): a valid strut must vanish/return `0` on the free massless field and be
discriminatingly inert on the gapped-atomless `:φ²:`. We test tensor-RG against this at the end.

---

## §1 · What tensor-RG is, and its one real advantage: renormalon-free

MERA is a real-space, multi-scale variational ansatz with a causal-cone (holographic) structure
[Vidal, PRL **99**, 220405 (2007)]; PEPS/cMPS are its higher-D / continuum cousins [Verstraete–Cirac,
PRL **104**, 190405 (2010)]; DMRG is the 1-D workhorse [White, PRL **69**, 2863 (1992)]. All are **real-space,
non-perturbative, variational** methods.

**The lone plus.** Being non-perturbative real-space RG, tensor-RG carries **no Borel/transseries
ambiguity** — it sidesteps exactly the `u=2` IR-renormalon block that sinks the gluon-condensate floor
(`ROOF-CONDENSATE-FLOOR` `[NG]`). Conceptually it is the right *kind* of tool for a renormalon-free strut.
**But "right kind" is not "works"** — see §2.

---

## §2 · The decisive obstruction: variational ⇒ UPPER bounds only `[T]`

The Rayleigh–Ritz principle gives `⟨ψ|H|ψ⟩ ≥ E₀` for **every** trial `|ψ⟩`; the min-max principle gives the
`k`-th Ritz value over a `k`-dimensional subspace `≥ E_k`. Both are **upper** bounds on the eigenvalues.
The **gap** is `Δ = E₁ − E₀`, and a lower bound on `Δ` needs a **lower** bound on `E₁` (with an upper bound
on `E₀`). Variational tensor-RG supplies neither: it gives upper bounds on `E₀` **and** `E₁`. Consequently a
"variational gap" `E₁^{trial} − E₀^{trial}` can sit **above or below** the true gap depending on ansatz
quality — it is **not a certificate**.

*Witness (Part A):* on a random symmetric `6×6` `H`, `min_ψ⟨ψ|H|ψ⟩ − E₀ = +1.06×10⁻¹ ≥ 0` and
`min_{2D subsp} (2nd Ritz − E₁) = +1.15×10⁻¹ ≥ 0` (upper bounds, as required); a constructed variational gap
shoots to `2.33` (true `1.02`) or drops to `1.72` with a poor `E₀` bound — **both directions**, hence no
certificate. **This is the same wall as `BASE-NONLOCAL-BOUND`**: test functions and variational ansätze are
upper-bound machines; the gap needs a lower-bound certificate (a Temple/Weinstein-type or spectral-method
bound), which tensor-RG does not provide. `[T]` (general) + `[O]` (the missing lower bound).

---

## §3 · The entanglement-spectrum gap is NOT the physical gap `[T]`

The natural "gap" inside tensor-RG / DMRG is the **entanglement-spectrum gap** of the entanglement
Hamiltonian `H_E = −ln ρ_A` on a region `A` [Li–Haldane, PRL **101**, 010504 (2008)]. The correspondence
"entanglement-spectrum gap ↔ physical gap" is a **heuristic observation**, used mostly to detect
topological order — **not a theorem**.

*Witness (Part B):* exact diagonalisation of the transverse-field Ising chain (`N=8`, `256` states),
half-chain cut. The physical gap `E₁−E₀` and the entanglement gap `ξ₁−ξ₀` (`ξᵢ = −ln` Schmidt value) are
**not proportional**: `ent/phys` runs over `≈0.95 → 5.3 → 24.8` as `h` goes `2.0 → 0.8 → 0.5`, and at
criticality (`h=1`, massless) the physical gap `→0` in the thermodynamic limit while the entanglement gap
stays `O(1)` (it is set by the CFT/finite-size structure). So even an *exact* entanglement spectrum does
**not** certify the physical mass gap — the identification is not a theorem. `[T]` (the two are distinct
functions) + `[NG]` (entanglement spectrum as a gap certificate).

---

## §4 · Bisognano–Wichmann: rigorous, but it is *vacuum* entanglement, not the gap `[T]`

For the vacuum of a relativistic QFT, the entanglement Hamiltonian of a half-space is **exactly** the boost
generator, `H_E = 2π K` (Bisognano–Wichmann, J. Math. Phys. **16**, 985 (1975); **17**, 303 (1976)). This is
genuinely rigorous and covariant — the entanglement structure of the vacuum *is* well-defined. **But it is
the spectrum of the boost, not the mass spectrum**: it exists for *any* relativistic QFT, gapped or
gapless, and gives no direct handle on `m₁`. It is a correct theorem about the wrong object for the gap.
`[T]` (BW) + `[NG]` (as a gap strut).

---

## §5 · "Gap = smallest scaling dimension" is a CFT slogan; 4D YM is massive `[NG]`

In a CFT, the finite-volume gap is set by the lowest scaling dimension (`E₁ − E₀ ∝ Δ_min / R`). But 4D pure
YM is **not** a CFT: it is asymptotically free in the UV (where `tr F²` has dimension 4 — an *irrelevant*
operator, no gap) and **massive** in the IR (where the scale is generated by dimensional transmutation,
not by a scaling dimension). The slogan describes the UV fixed point or a critical point, neither of which
is the physical gap. `[NG]` (scaling-dimension reading of the 4D gap) + `[O]` (the IR scale is
dimensional transmutation, the wall itself).

---

## §6 · Rigor status: rigorous in 1+1D, absent in 3+1D `[T]/[NG]`

- **1+1D `[T]`:** the area law is proven for 1-D gapped systems [Hastings, J. Stat. Mech. (2007) P08024];
  a gap implies exponential clustering (Nachtergaele–Sims, CMP **265**, 119 (2006)); finitely-correlated
  states and cMPS give controlled ansätze. Even here, the **converse** (exponential decay ⇒ gap with a
  *controlled constant*) is not available in general.
- **2+1D and 3+1D `[NG]`:** there is **no rigorous MERA/PEPS→gap theorem** for interacting gauge theory.
  PEPS constructions of lattice gauge theories exist [e.g. Zohar–Cirac–Burrello, PRB **92**, 045125 (2015),
  and the PEPS-gauge literature — bibliographic details flagged per corpus convention], but they are
  variational and computationally uncontrolled at bond dimensions needed for `α_s~O(1)`.
- **Circularity `[O]`:** tensor-RG efficiency *assumes* an area law (bounded entanglement). A 4D **area law
  for interacting pure YM is itself unproven**, and may be as hard as the gap — so invoking tensor-RG to
  prove the gap risks assuming a statement of comparable difficulty.

---

## §7 · Acid test and iron laws

- **Acid:** does tensor-RG "fail on the free massless field"? The **physical** gap it can only *upper*-bound
  is `0` on free massless — fine — but the **entanglement gap** is `O(1)` there (§3), i.e. the object it
  actually computes is **acid-dirty**: it does not vanish on the massless field. So the entanglement-spectrum
  "gap" **fails the acid test** (it would "prove" a gap for the free field). `[NG via acid]`.
- **IL-i (non-positivity):** fails — the entanglement gap is positive for the free massless field, so it is
  not discriminating.
- **IL-ii (nonabelian):** the entanglement/area-law structure is present for U(1) too — not
  nonabelian-essential.
- **IL-iii (no small parameter):** no small parameter at `α_s~O(1)`; bond-dimension truncation is
  uncontrolled exactly there.

---

## §8 · What tensor-RG DOES contribute (honest partial yield)

1. **Renormalon-free by construction** — the *right kind* of non-perturbative real-space tool, addressing
   the condensate's Borel ambiguity conceptually.
2. **The upper half of a sandwich `[V]`:** a variational **upper** bound on `E₁` and on `E₀`. If an
   independent **lower** bound on `E₁` existed (Temple/Weinstein, or a rigorous spectral/transfer certificate),
   tensor-RG would help bracket the gap. **The lower half is the missing strut** — exactly `H1`.
3. **A diagnostic bridge:** the entanglement-spectrum ↔ gap correspondence, *used as a numerical probe*, can
   locate where a gap opens (lattice/classical simulation), even though it certifies nothing.

---

## §9 · Verdict and exact failing step

**`[NG]` for the class "tensor-network RG as a rigorous computable lower-bound strut for the 3+1D `0⁺⁺`
gap."** Three independent reasons, any one fatal: (i) variational ⇒ **upper bounds only** (`[T]`); (ii) its
"gap" is the **entanglement-spectrum gap**, acid-dirty and ≠ the physical gap (`[T]`, witness); (iii)
**no rigorous 3+1D MERA→gap theorem**, and the 4D area law is unproven (circular) (`[NG]`).

**Exact missing step `[O]`:** a **lower bound on `E₁`** (the first excitation of the `0⁺⁺` transfer
operator) with controlled error in `3+1D` — a Temple/Weinstein-type or rigorous spectral certificate. Tensor
RG does not provide it; the same upper/lower asymmetry is the base-round wall, in the tensor-network
coordinate.

**One-line answer to the directive:** tensor-RG is **renormalon-free but the wrong half of the
certificate** — it gives a computable **upper** bound and a heuristic entanglement-spectrum probe, not the
**lower-bound strut** the gap needs, and it has **no rigorous 3+1D form**. **Clay 0%; no QED.**

---

*References (checked to the best of knowledge; flags per corpus convention):*
Vidal, PRL **99**, 220405 (2007); Li–Haldane, PRL **101**, 010504 (2008); Bisognano–Wichmann, J. Math. Phys.
**16**, 985 (1975) & **17**, 303 (1976); Hastings, J. Stat. Mech. (2007) P08024; Nachtergaele–Sims, CMP
**265**, 119 (2006); Verstraete–Cirac, PRL **104**, 190405 (2010); White, PRL **69**, 2863 (1992);
Zohar–Cirac–Burrello, PRB **92**, 045125 (2015); Swingle, PRL **109**, 091602 (2012); Ryu–Takayanagi, PRL
**96**, 181602 (2006); Temple (1928) / Weinstein (1934) for rigorous excited-state **lower** bounds.
Witness: `w7/roof_tensor_rg.py` (numpy, runs clean).

*Cross-refs:* `ROOF-SYNTHESIS.md`, `BASE-NONLOCAL-BOUND.md`, `ROOF-CONDENSATE-FLOOR.md`,
`ROOF-CONFINEMENT-GAP.md`, `PROFILE.md`.
