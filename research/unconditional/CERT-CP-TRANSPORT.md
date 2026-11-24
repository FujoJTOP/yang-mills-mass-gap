# CERT-CP-TRANSPORT — the CP / Lüscher-flow residue transport does NOT certify the atom: it is a finite-`(L,a)` normalization identity (vacuous for the continuum), the flow is downstream, and the exact failing step is the `L→∞, a→0` survival = A3 = H1

> **VERDICT (no QED, Clay 0%).** Certificate attempt on the atom half via CP: `f₁² = Z_O(t)² f₁²_strong`
> with `Z_O(t)` from the Lüscher gradient flow, taking the P3 strong-coupling atom (`f₁²_strong>0` at
> β<1/48) to β≈2.5. **It does not certify the atom.** Witness `w7/cert_cp_transport.py`.

---

## §0 · Labels

**[T]** 真定理 · **[C]** 条件定理 · **[O]** 精确缺失步 · **[NG]** no-go · **[V-num]** 数值. Acid: free-massless / `:φ²:`.

---

## §1 · Three exact obstructions

**(1) The transport is a finite-`(L,a)` normalization identity `[T]`.** For a Hermitian operator `tr F²` the
field-renormalization `Z_O` is real ⇒ `Z_O²≥0`. Given `Z_O≠0` and the **same** atom, `f₁²_phys = Z_O² f₁²_strong > 0`
follows from P3's `f₁²_strong>0`. But this is a *fixed-finite-`(L,a)`* statement. At fixed finite `L` the transfer
matrix is compact ⇒ **every** residue is `>0` — for an atom **and** a cut alike (witness: fixed atom gap `0.5`,
cut gap `(2π/L)²→0`; both have `f₁²_L>0` at every finite `L`). So the transport is **vacuous** for the Clay
content.

**(2) The flow preserves the level, not the weight `[V-num]`.** The flowed residue `f₁²(t) = f₁² e^{−2t m₁²} → 0`
**absolutely** as `t→∞` (higher masses suppressed more, so the atom's *relative* share →1 — but the *absolute*
weight decays: `0.50→6.8e−2→2.3e−5→1.0e−9` at `t=0,1,5,10`). The flow is **downstream** of `f₁²(0)`: it cannot
create the atom; "survives `t→∞`" is convention-dependent (relative yes, absolute no).

**(3) Flow / RG renorm is local ⇒ blind to isolation `[T]`.** `Z_O(t)` (flow) and `Z_O(β)` (RG) are **local**
operator renormalizations (pointwise/statewise). They carry **no information** about whether the lightest level
stays **isolated** as `L→∞, a→0`.

---

## §2 · Exact failing step `[O]`

The β-transport (which `Z_O(β)` could in principle track via the operator anomalous dimension) is fine but
**vacuous at finite L** (compaction). The Clay content — the **`L→∞, a→0` survival of the atom with isolation**
(the cut's gap →0, the atom's stays `>0`) — is **not covered by any local renorm**; it is **isolation `I2` = the
weak-* discontinuity of the boundary density = A3 = `H1`**.

This matches `FALSIFY-FIELD`'s corollary (finite-L atom `>0` is trivial; the content is the `a→0, L→∞` survival)
and `ONLY-EXIT` (the atom is the weak-* discontinuity). **CP is the `sorry` input in `lean/YMStatement/Reduction.lean`
for exactly this reason.**

*Refs:* `ONLY-EXIT-SYNTHESIS.md`, `BLAST-WI-SYNTHESIS.md`, `FALSIFY-FIELD.md`, `PROOF-P3-UNCONDITIONAL.md`,
`lean/YMStatement/Reduction.lean`.
