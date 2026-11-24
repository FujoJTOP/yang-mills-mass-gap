# PNT-DUAL-CONDENSATE — monopole condensation gives a renormalon-free MASS `[C]`, but the KL ATOM `f₁²>0` is NOT delivered `[O]`: the atom is the nonperturbative `tr F²`↔vortex-ring overlap

> **VERDICT (no QED).** Attempting the dual-superconductor route: **(i)** monopole condensation is
> *semi-rigorous* — exact in 3D compact U(1) / 4D Z_N, numerical for SU(N) via abelian projection — with a
> rigorous 4D SU(N) derivation being confinement = `H1` (no new obstruction). **(ii)** Lifting the
> renormalon-free dual-Higgs *mass* to the KL *atom* **fails to be delivered**: the minimal dual
> Abelian-Higgs model fixes the masses `m_V, m_H` but **not** the `tr F² ↔ |φ|²` coupling, hence **not**
> `f₁²` (free parameter); and the physical ratio-3.9 `0⁺⁺` is a **nonperturbative vortex ring**, not any
> perturbative dual-AHM mode, whose `tr F²` overlap `⟨0|tr F²|ring⟩` is an independent nonperturbative
> matrix element. **The atom-lift is `H1` in the dual coordinate.** Clay 0%.

---

## §0 · Labels, acid test

**[T]** 真定理 · **[C]** 条件/半严格 · **[O]** 精确缺失步 · **[NG]** no-go for the named class · **[V-num]** new numeric.
**Acid:** the functional must return `0` on free massless (`v=0`) *and* fail discriminatingly on the
gapped-atomless `:φ²:` (`f₁²=0`, `m>0`). Witness `w7/pnt_dual_condensate.py`.

---

## §1 · (i) Monopole condensation — semi-rigorous `[C]`, exact only in lower dimensions

Dual Abelian–Higgs ('t Hooft–Mandelstam): `L_D = −¼ F_D² + |Dφ|² + λ(|φ|²−v²)²`, `m_V=g̃v` (dual photon,
`1⁻⁻`), `m_H=√(2λ)v` (dual Higgs, `0⁺⁺`), London `σ≃πv² ln κ_GL`. Monopole condensation ⟺ `v≠0` ⟺ dual
Meissner ⟺ confinement.

- **`[T]` exactly:** 3D compact U(1) (Polyakov — monopole instantons of the 2D dual condense ⟹ linear
  potential); 4D Z_N gauge (self-dual; Peierls/vortex condensation). These are the anchor points.
- **`[C]` numerically:** 4D SU(N) via maximal-abelian projection — monopoles are the Dirac defects of the
  projected U(1)^{N−1}; monopole condensation is measured (density, percolation, dual Meissner), Abelian
  dominance holds. Strong *evidence*, not proof; gauge-dependence of the projection is a known caveat.
- **`[O]` rigorous 4D SU(N):** a proof of monopole condensation **is** confinement = `H1`. **No new
  obstruction** — the dual picture reproduces the wall, it does not move it.

---

## §2 · (ii) Lifting the mass to the KL atom — NOT delivered `[O]`

The KL atom is `f₁² = |⟨0|tr F²|G₁⟩|²`, the **residue of a pole** of `⟨tr F² tr F²⟩` at the lightest `0⁺⁺`.
A classical/dual *mass* is not an atom. Two exact obstructions (witness §1–§3):

1. **`[T]` the minimal dual AHM fixes masses, not the overlap.** The `0⁺⁺` spectral density of the
   *electric* operator is `ρ_trF²(s) = [dual-photon-pair CUT from s=4m_V²] + f₁²·δ(s−m_H²)`. The pole
   exists only if the *non-minimal* coupling `tr F² ↔ |φ|²` is switched on; the minimal AHM leaves
   `0 ≤ f₁² < ∞` free. A **cut at `2m_V` is not an atom** — it is exactly the `:φ²:` failure mode
   (gapped, `f₁²=0`). Free massless (`v=0`): `m_V=m_H=0`, no pole, `f₁²=0` (acid returns 0).
2. **`[V-num]` the physical `0⁺⁺` is not a perturbative dual-AHM mode.** For every `κ_GL`, at least one
   perturbative mode is below the lattice glueball `m_{0⁺⁺}/√σ≈3.9`; pushing `m_H` to 3.9 needs
   `κ_GL≈11`, where `2m_V/√σ≈0.18` is far below. The ratio-3.9 `0⁺⁺` is a **vortex ring** (closed flux
   tube, `m_ring~2πRσ` with `R~c/√σ` ⟹ O(3–4)`√σ`) — a collective, nonperturbative object.

⇒ **The atom-lift reduces to the nonperturbative matrix element `⟨0|tr F²|vortex-ring⟩`** — independent of
dual-AHM semiclassics. This is `H1` (the same nonperturbative overlap that every prior route met), here in
the dual-superconductor coordinate.

---

## §3 · Answer to the directive's question

- **Is the vortex-ring the right object?** **Yes** — it, not the dual-Higgs, matches the physical ratio 3.9.
- **Does it carry the atom?** **Open** — its `tr F²` overlap is a nonperturbative matrix element the dual
  AHM does not fix. Dual AHM delivers a renormalon-free **mass scale** (real, useful) but **not** `f₁²>0`.

---

## §4 · Honest bottom line

- **Renormalon-free mass:** `[C]` — the dual-Higgs/vortex scale is written without the `u=2` pole (dim-2
  `v`). Genuine advantage over the `⟨tr F²⟩` floor.
- **KL atom `f₁²>0`:** `[O]` — not delivered; it is the nonperturbative `tr F²`↔vortex-ring overlap = `H1`.
- **Monopole condensation:** `[C]` semi-rigorous; rigorous = `H1`.
- **Net:** the dual-superconductor "renormalon-free" advantage **survives for the mass only**, not for the
  atom that Clay's `tr F²` formulation needs. Same wall, dual coordinate. **Clay 0%; no QED.**

*Refs:* `ROOF-DUALSUPERCONDUCTOR.md`, `ROOF-CONDENSATE-FLOOR.md`, `INVENT-SYNTHESIS.md`,
`BASE-ORDER-PARAMETER.md`; witness `w7/pnt_dual_condensate.py`.
