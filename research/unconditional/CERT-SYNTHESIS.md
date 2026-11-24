# CERT-SYNTHESIS — the tools are NOT complete and a certificate cannot be assembled: 8 certificate-attempts all reduce H1 to the same irreducible core (the χ-half constant = the `u=2` condensate, exactly marginal, no slack; the atom-half = A3 isolation)

> **VERDICT (no QED, Clay 0%).** The user's decision point — *if the tools are complete, start the
> certificate; if not, keep building* — resolves to: **the tools are NOT complete, and a certificate cannot be
> assembled from them.** Eight targeted certificate-attempts were run. **All eight reduce H1 to the same
> irreducible core**, from every angle: the **χ-half** is one marginal inequality whose constant is the
> `u=2`-ambiguous condensate `S_bare = ⟨tr F²⟩` — and the χ-half is **exactly marginal, no slack** (`gB(1−t²)`
> vs `1/χ`, both `O(1)`; the negative-mode alignment `t→1` drives it to `0`); the **atom-half** (`CP`
> transport + isolation) is a **local renormalization, blind to isolation**, vacuous at finite `L`, = `I2` =
> the weak-* discontinuity = A3. **This is the definitive statement: `H1` is the irreducible core, not a
> tool-gap.**

---

## §0 · Labels

**[T]** 真定理 · **[C]** 条件定理 · **[O]** 精确缺失步 · **[NG]** no-go · **[V-num]** new numeric. Acid: `≤0`/`=0` on free massless *and discriminatingly* on `:φ²:`.

---

## §1 · The 8 certificate-attempts (all reduce to H1 / A3)

### χ half (5 routes)
| route (file) | result |
|---|---|
| **χ inequality** `CERT-CHI-INEQ` | `[T]` the negative-mode softening is genuinely **controlled** (bounded `D_Sav≤2g²/λ₄` + **deep-IR decoupled** `cα~(k/√gB)²→0`); but one irreducible bare-curvature floor `S_bare = ⟨tr F²⟩` remains = the condensate = H1. Closes at β<1/48 (SZZ), needs H1 at 2.5 |
| **χ condensate floor** `CERT-CHI-CONDENSATE` | `[T]` the **defect `D` is floored renormalon-free** (`Var(φ²)→2/λ₄`, c-independent); but the **floor lives on `S_bare`**, which is the `u=2`-ambiguous condensate (`a_n~Γ(n+1+2b)`, Borel pole `u=2`). The rank-one form gives only a one-sided `D≤const`, **not a subtraction** |
| **χ exact numbers** `CERT-CHI-SHARP` | `[T]` the χ-half is **exactly marginal**: `gB(1−t²)` vs `1/χ`, both `O(1)`, no slack; `t→1` (negative mode aligned) ⇒ LHS `=0`. The missing `t` (`D<S`) = H1 — the structural reason `D/S≈0.94–1.00` |
| **χ a→0** `CERT-CHI-LATTICE` | `[T]` the `a→0` wall **is** `χ_lattice ≲ a²` = dimensional transmutation = H1; renorm-condition/step-scaling/monotone all need the IR scale; `χ_lattice~a^p`, `p=2` (gap) vs `p=0` (strong coupling) both fit compactness — **the exponent is the gap** |
| **χ variational** `CERT-CHI-VARIATIONAL` | `[T]` `Var(a)` = [diagonal `σ²_link(β)` genuine `0.32→0.10`, but `O(1)`-insufficient] + [off-diagonal `≤ Var(X)/Δ`, which needs the gap `Δ>0`] = H1 |

### Atom half + new routes (3 routes)
| route (file) | result |
|---|---|
| **CP transport** `CERT-CP-TRANSPORT` | `[T]` `CP`/Lüscher-flow: `Z_O²≥0` gives `f₁²_phys>0` at **finite L — vacuous** (every residue `>0`, atom and cut alike); the flow **decays the absolute weight**, is **downstream** of `f₁²(0)`; flow/RG renorm is **local ⇒ blind to isolation**. = `I2` = A3 = H1 (why `CP` is a `sorry` in Lean) |
| **New route: modern rigorous RG** `ROUTE-NEW-RG` | `[NG]` no modern RG (BBS/BBD/Chandra–Hairer) reaches the crossover: all **small-parameter**, but `u(2.5)=0.507`, `g²=1.6`, `O(1)`; 4D is **critical** (`[g²]=0`); blocked by the **Savvidy-indefinite large-field** estimate = H1 |
| **New route: holographic bulk** `ROUTE-NEW-HOLO` | `[NG]` hard-wall gives `m_n z₀=j_{2,n}` — every level/χ is a function of `z₀` alone; **NEC is saturated in pure AdS** (constrains shape, not scale); the `m/√σ≈3.9` match is **kinematic**; bounding `m`/`χ` needs `z₀` = the IR wall = the confinement scale = H1 |

**0 / 8 certificate. Every route reduces H1 to the same core: the χ-half constant = the `u=2` condensate (exactly marginal); the atom-half = A3 isolation.**

---

## §2 · The definitive judgment

- **The tools are NOT complete.** A certificate of the gap needs `H1 ∧ CP ∧ P3 ∧ no-dissolution` (Lean `clay_gap_conditional`). Of these:
  - **`H1` (χ half)** is one marginal inequality, but its constant is the `u=2`-ambiguous condensate `⟨tr F²⟩`, and it is **exactly marginal (no slack)** — provable only if the crossover is already controlled = H1.
  - **`CP` (atom transport)** is a [T] identity that is **vacuous at finite L** (every discrete spectrum has `f₁²>0`) and **blind to isolation**; the Clay content (`L→∞,a→0` isolation) is the **weak-* discontinuity = A3 = H1**.
- **So "start the certificate" is not an available action** — the certificate *is* H1 (χ) plus A3 (atom). The tools reduce H1 to a single knife-edge inequality; there is no remaining "tool" that closes it.
- **This is not a defect of the tool set — it is the structure of the problem** (ONLY-EXIT, the weak-* discontinuity of the measure).

---

## §3 · The campaign's total honest state (bilan)

**Proven / established (真实资产):**
1. **Rigorous conditional reduction** `H1 ⇒ gap` … `H1 ∧ CP ∧ P3 ∧ no-dissolution ⇒ MassGap` — **machine-checked in Lean, 0 sorry**.
2. **`H1 ⊊ gap`** — `H1` alone is *necessary but not sufficient* (blind to a thin atom) — **a compiled theorem** (`H1_does_not_imply_atom`).
3. **The χ-half reduces to one marginal inequality**; the **negative-mode softening is genuinely controlled** (bounded, deep-IR decoupled, sign-protected, `(L,a)`-uniform); `m_{0⁺⁺} ≥ 4√σ` `[T]`; the **Savvidy false obstruction is removed**.
4. **ONLY-EXIT**: the only forcing exit is `H1` (the atom = weak-* discontinuity); 8 falsification routes + 8 tool classes all confirm it.
5. **Best available lattice diagnostics** for the atom/isolation: P3 approach-law `T(τ)`, quantitative `Δ≥κ(τ)`, Re-slope `R`, cross-validated `(E₁,Δ)`, FSS law.

**The wall (unchanged):** `H1` = the `a→0`, `L→∞` isolation of the `0⁺⁺` atom = `χ_lattice ≲ a²` = dimensional transmutation = the condensate floor = the weak-* discontinuity of the boundary density. **Reached from ~150 directions; each time it is the same invariant.**

**Clay = 0%.** No QED. **The certificate is H1; H1 is the open problem; the tools are exhausted.**

*Refs:* `CERT-CHI-{INEQ,CONDENSATE,SHARP,LATTICE,VARIATIONAL}.md`, `CERT-CP-TRANSPORT.md`,
`ROUTE-NEW-{RG,HOLO}.md`; `BLAST-WI-SYNTHESIS.md`, `ONLY-EXIT-SYNTHESIS.md`, `lean/YMStatement/{Reduction,H1Sharpening}.lean`.
