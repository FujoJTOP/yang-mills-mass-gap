# R-WAVE-SYNTHESIS — the atom is a **binding** phenomenon, **isolation** is its parent, and the gap is a **third predicate**: six routes, one corrected target, and a resolution bound that makes every atom-reader circular

> **Correction (T-wave, 2026-11-24).** Two formulations below are superseded. (i) Q6's "no carrier"
> is false: T4 constructs the total-normalized zero mode, with `Var(q)=χ=O(1)` (the `V^{-1}` was the
> mean-normalized mode). (ii) "Binding is the discriminant" is too strong: T6 shows spectrally
> `I ⇒ B_spec ⇒ A`, so spectral binding is implied by isolation; dynamical binding is the atom's
> mechanism, not an independent axis. The independent axes are `A` (atom visibility) and `I`
> (isolation); the gap `G` is orthogonal and needs H1/ceiling plus an absolute scale. The corrected
> statement is in §8.

> **VERDICT (no QED, Clay 0%).** The $H1$ wave left the atom as the irreducible target. Six routes
> attacked it. The result is a **structural correction**, not a proof:
>
> 1. **Isolation is the parent, not the child** `[T]/[NG]`. A support gap forces the atom
>    ($I_2\Rightarrow I_1$), but an atom can sit at the edge of a continuum ($I_1\not\Rightarrow I_2$).
>    So the campaign's pair ``atom $+$ isolation'' is qualitatively **isolation alone**; the atom is
>    downstream. Isolation is an **order-2** (second-eigenvalue/TP$_2$) condition; the atom is
>    order-1 (Perron).
> 2. **Binding is the atom's mechanism, not an independent discriminant** `[T]/[V]` (corrected by
>    U6). The only exact positive anchor is the strong-coupling residue theorem; its mechanism is a
>    single flux tube (one bound particle). Every solvable model computes the ceiling $\chi=1/G''$
>    exactly and has $f_1^2=0$. The structural facts are: binding is the mechanism of the atom; the
>    independent data are the **weight** $W$ and **isolation** $I$; and $W$ and $I$ are incomparable.
> 3. **Every atom-reader is circular** `[T]/[NG]`. Discontinuous functionals do read $f_1^2$ exactly,
>    but each presupposes the gap, needs infinite data, or reads the largest atom. The
>    Stone--Weierstrass resolution bound $\sup_\nu|E-g|\ge f_1^2/2$ makes the instrument error at
>    least the signal.
> 4. **The atom is not one object** `[T]/[NG]`. Existence is protected in the confined phase (by
>    confinement kinematics and the trace anomaly, not by $H1$); the residue floor is in models easier
>    than the gap; isolation is the gap. $H1+$squeeze$+$positivity$+$a fixed gap do not force an atom.
>
> The wave does not reduce Clay. It identifies the irreducible object as **isolation/binding**, and it
> corrects the target from ``atom'' to an order-2 condition on $\mathcal A/\mathcal G$.

---

## §0 · The six verdicts

| Route | Verdict | Result |
|---|---|---|
| **R1** residue transport | `[T]`/`[NG]` | $Z_O$ is the wave-function renormalization of $\operatorname{tr}F^2$, fixed by the trace anomaly: $Z_O=(\beta_{YM}/2g)(\mu_0)/(\beta_{YM}/2g)(\mu)>0$, so $f_1^2(\mu)=Z_O^2f_1^2(\mu_0)$ is exact RG covariance for the same state. $Z_O>0$ is trivial (no Banks--Baxks zero in pure YM); it was never the wall. The missing step is the state overlap $\Omega=|\langle G_{\rm tgt}|O|G_{\rm eff}\rangle|^2=1$, i.e. the same atom persisting with uniform isolation. No-go: $Z_O^2\cdot0=0$, so no $Z_O$ bound lifts a dissolved atom; $\Omega$ is bounded only above. |
| **R2** isolation $I_2$ | `[T]`/`[NG]` | $I_2\Rightarrow I_1$ (a support gap forces the atom), $I_1\not\Rightarrow I_2$ (a continuum can touch $m_1^2$). Moments upper-bound isolation only (dual Christoffel); the lower end $0$ is realised by thin-gap measures matching any finite moments. Isolation is order-2 (TP$_2$) above the atom's order-1 (Perron), absent on $\mathcal A/\mathcal G$. Exact in 2D YM (Casimir ratio $\rho^2=8/3$) and the 1D rotor; absent on the free massless composite. Missing estimate: uniform $\rho=m_{0^{++*}}/m_{0^{++}}\ge1+c$ at $\beta\approx2.5$. |
| **R3** atom functional | `[T]`/`[NG]` | Many discontinuous atom-readers are exact (Stieltjes jump, Christoffel/Hankel $\lambda_n$, sup-atom, de Branges/spectral shift), all returning $f_1^2$ and passing the free-cut acid test. None is non-circular: each evaluates at $x=m_1^2$ (presupposes the gap), or needs $\varepsilon\to0$/$n\to\infty$ (infinite data), or reads the largest atom and is positive on the finite-volume free field. Stone--Weierstrass: $\sup_\nu|E-g|\ge f_1^2/2$ for every weak-$*$ continuous finite-data $g$, so instrument error $\ge$ signal. **Side correction:** the ch09 atom-weight formula with window $=\varepsilon$ returns $f_1^2/2$, not $f_1^2$; corrected to the ordered double limit. |
| **R4** residue measurement | `[O]`/plan | Minimal decisive run: SU(2) $\beta=2.5$, $L=8,12$, $N=2000$, therm$\ge2000$, gap $20$--$50$, APE family $\{0,2,4,8\}$, zero-momentum $0^{++}$ $\operatorname{tr}F^2$; GEVP + two-state fit for $f_1^2$ and $\Delta=m_2-m_1$; finite-size law $\Delta(L)=m_1+c/(m_1L^2)$. Cost $\approx6$--$12$ CPU-h, one job. Existing: atom at $\beta=2.2$ ($L=6$, $N=4000$: $f_1^2_{\rm bare}=7.14\pm0.42$, $17\sigma$, $m_G=1.51a$); clean null at $\beta=1.0$; $\beta=2.5$ only a 40-config pilot ($C_1/C_0\approx0.18$--$0.29$, $\approx3.5\sigma$, no plateau). The atom is gauge-invariant, so no Gribov/Landau fixing is needed for it. |
| **R5** exact anchors | `[T]`/`[V]` | The only exact positive lower bound is P3 strong coupling ($f_1^2\ge r_0(\beta)>0$, $L$-uniform, $\beta<1/48$). Mechanism: binding — the minimal connected $\operatorname{tr}F^2$ polymer is one four-walled flux tube, amplitude a modulus-square of a positive $f\otimes\bar f\to1$ overlap. It transports nothing to $\beta=2.5$: no isolation, no $a\to0$, cannot cross the $1/48$ Ursell wall; $f_1^2$ is upper-semicontinuous. All solvable models compute the ceiling $\chi=1/G''$ exactly, never the atom, and each has $f_1^2=0$/degenerate: 2D YM ($\operatorname{tr}F^2=H$, conserved, character orthogonality), $\beta=0$ Haar ($\chi=0.2688>0$, $f_1^2=0$), one plaquette (no time), 1D chain (ceiling only). |
| **R6** adversarial kill | `[T]`/`[NG]` | The atom is not one object. Existence $f_1^2>0$ is not killable in the confined phase, but the protection is confinement kinematics plus the trace anomaly, **not** $H1$; changing the carrier to $\operatorname{tr}F\tilde F$ kills $f_1^2$ exactly (parity) while keeping state and gap. Residue floor $f_1^2\ge r_0$ is robust and in models strictly easier than the gap (double well: $O(1)$ residue, gap $\to0$, but $\chi\to\infty$). Isolation is the gap modulo IN-11. Strongest no-go: $H1+$squeeze$+$positivity$+$fixed gap do not force an atom (a gapped atomless measure with bounded $\chi$ exists). PF protects nothing; no-dissolution protects existence conditionally; P3 only at $\beta<1/48$. |

---

## §1 · The corrected target: isolation, not the atom

The campaign's sufficient condition was written as the pair
\[
  H^*=(f_1^2\ge r_0\ \wedge\ I_2)\ \wedge\ (\chi\le1/\lambda).
\]
R2 shows the pair is redundant in the wrong direction: $I_2\Rightarrow I_1$, so $I_2$ alone already
gives the atom, while $I_1$ does not give $I_2$. The irreducible object is therefore **isolation**,
an order-2 condition: a uniform ratio
\[
  \rho=\frac{m_{0^{++*}}}{m_{0^{++}}}\ge1+c>0
  \qquad(\beta\approx2.5),
\]
equivalently second-eigenvalue/TP$_2$ control on the orderless quotient $\mathcal A/\mathcal G$.
The atom is the order-1 (Perron) shadow of this. This explains why every banked no-go sits one tier
too low: positivity, Perron--Frobenius, and the squeeze control order-1 data; isolation needs order-2.

**What moments can and cannot do.** The Christoffel dual gives an upper bound on isolation; the lower
end is $0$, and thin-gap measures match any finite moment list while driving the isolation to zero.
So no finite-moment argument can force $I_2$; this is the same weak-$*$ obstruction as R3, one tier up.

---

## §2 · Binding is the atom's mechanism

The only exact positive anchor is P3: $f_1^2\ge r_0(\beta)>0$ for $\beta<1/48$, $L$-uniform. Its
mechanism is not the Haar measure, compactness, or the gauge-invariant carrier; it is **binding**:
the minimal connected $\operatorname{tr}F^2$ polymer is a single flux tube (one particle), and its
amplitude is the modulus-square of a positive overlap $f\otimes\bar f\to1$. Every solvable model that
lacks binding has $f_1^2=0$:

| model | ceiling $\chi=1/G''$ | atom $f_1^2$ | why |
|---|---|---|---|
| 2D Yang--Mills | exact | $0$ | $\operatorname{tr}F^2=H$ conserved; character orthogonality |
| $\beta=0$ Haar | $0.2688$ exact | $0$ | no action, no binding |
| one plaquette | $G''(2.5)=0.730$ | degenerate | no time direction |
| 1D chain | $G''\to0.2165$ | $0$ | gapped but no 4D binding |
| strong coupling $\beta<1/48$ | — | $\ge r_0>0$ | flux-tube polymer |

This is the wave's sharpest structural statement: **the atom is a binding phenomenon, not a measure
phenomenon.** But binding is a *mechanism*, not an independent axis (U6): the independent data are
the weight $W$ and isolation $I$. A proof of the atom must control the binding of the minimal flux
tube, not the Jacobian
or the variance.

---

## §3 · Atom readers and the resolution bound

Discontinuous functionals do read the atom exactly: the Stieltjes boundary jump, the
Christoffel/Hankel $\lambda_n$, the sup-atom, and de Branges/spectral shift all return $f_1^2$ and
vanish on the free cut. But every one is circular:

- it evaluates at $x=m_1^2=\inf\operatorname{supp}\nu$ (presupposes the gap);
- or it needs $\varepsilon\to0$/$n\to\infty$ (infinite data, A3);
- or it reads the largest atom, not the edge, and is positive on the finite-volume free field.

Stone--Weierstrass makes this quantitative: for every weak-$*$ continuous finite-data functional $g$,
\[
  \sup_\nu|E(\nu)-g(\nu)|\ \ge\ f_1^2/2 .
\]
The instrument error is at least the signal. The atom is computable from the exact measure by the
ordered double limit $d\downarrow0$ after $\varepsilon\downarrow0$, but not from finite $(L,a)$ data.

**Erratum.** The literal ch09 formula with the window set equal to the regulator ($d=\varepsilon$)
returns $f_1^2/2$, because the Poisson kernel integrates to
$(1/\pi)[\arctan(1)-\arctan(-1)]=1/2$ over $[m-\varepsilon,m+\varepsilon]$. The corrected statement
is the ordered double limit; the paper's ch09 and errata appendix have been fixed.

---

## §4 · Adversarial status: the atom is not one object

| object | status | protection |
|---|---|---|
| existence $f_1^2>0$ | not killable in the confined phase | confinement kinematics + trace anomaly, **not** $H1$ |
| residue floor $f_1^2\ge r_0$ | robust; in models easier than the gap | double well: $O(1)$ residue, gap $\to0$, but $\chi\to\infty$ ($H1$ fails) |
| isolation $I_2$ | the gap, modulo IN-11 | — |

The strongest no-go: $H1$ + squeeze + positivity + a fixed gap do **not** force an atom; there is a
**gapped atomless** measure with bounded $\chi$ (this sharpens Q6's thin atom from gapless to gapped).
PF protects nothing (it is vacuum-anchored); no-dissolution protects existence conditionally; P3 only
at $\beta<1/48$. The free-field acid test passes ($f_1^2=0$).

---

## §5 · Numerical status

R4 specifies the minimal decisive run: SU(2) $\beta=2.5$, $L=8,12$, $N=2000$, therm$\ge2000$, gap
$20$--$50$, APE family $\{0,2,4,8\}$, zero-momentum $0^{++}$ $\operatorname{tr}F^2$; GEVP + two-state
fit for $f_1^2$ and $\Delta=m_2-m_1$; finite-size law $\Delta(L)=m_1+c/(m_1L^2)$; cost $\approx6$--$12$
CPU-h, one job. Existing data: the atom is seen at $\beta=2.2$ ($L=6$, $N=4000$, $17\sigma$), clean
null at $\beta=1.0$, and at $\beta=2.5$ only a 40-config pilot ($\approx3.5\sigma$, no plateau). The
atom is gauge-invariant, so the residue measurement does not need Gribov/Landau fixing (that protocol
belongs to the gluon-carrier cross-check). The M1 run is in progress; its result is reported
separately.

---

## §6 · What the wave won

1. **A structural reversal** (R2): isolation is the parent; the atom is downstream; the irreducible
   object is an order-2 TP$_2$ condition, not the order-1 atom.
2. **Binding as the atom's mechanism** (R5, corrected by U6): the atom is a flux-tube phenomenon; all solvable models are
   atomless; the only exact anchor is strong coupling.
3. **The resolution bound** (R3): a quantitative Stone--Weierstrass obstruction;
   $\sup|E-g|\ge f_1^2/2$; every atom-reader is circular.
4. **A transport no-go** (R1): $Z_O$ is trivial; the missing step is the state overlap with uniform
   isolation.
5. **An adversarial classification** (R6): existence / residue floor / isolation have different
   statuses; $H1$ plus a fixed gap does not force an atom.
6. **An erratum fixed** (R3): the ch09 atom-weight formula factor of $2$.

---

## §7 · What the wave did not win, and the next species

- No proof of the atom or of isolation in 4D at $\beta=2.5$.
- No transport from strong coupling; no $a\to0$; no continuum limit.
- No reduction of Clay: isolation is the gap modulo IN-11.
- **The next species is a binding criterion**: a rigorous flux-tube/polymer lower bound on the
  minimal connected $\operatorname{tr}F^2$ overlap at $\beta=2.5$, strong enough to survive the
  $1/48$ Ursell wall and to force the order-2 ratio $\rho\ge1+c$. The variance/measure tools of the
  H and Q waves are the wrong tier; the atom wave says the missing object is a **binding theorem**,
  not another variance bound.

**One line.** *The atom is binding; isolation is the gap; the only exact anchor is strong coupling;
every atom-reader is circular at the resolution bound; $H1$ plus a fixed gap does not force an atom.
No proof, no QED, Clay 0%.*

*Refs:* `R1-RESIDUE-TRANSPORT.md`, `R2-ISOLATION.md`, `R3-ATOM-FUNCTIONAL.md`, `R4-RESIDUE-PROTOCOL.md`,
`R5-EXACT-ANCHORS.md`, `R6-ATOM-ADVERSARIAL.md`; witnesses `w7/r1_*.py` … `w7/r6_*.py`;
`Q-WAVE-SYNTHESIS.md`, `Q6-H1-STATUS.md`, `docs/paper-v2/ch09-weakstar.tex` (corrected).

---

## §8 · Corrections from the T wave (2026-11-24)

The T wave tested the two formulations above that were too strong. The corrected picture:

1. **The carrier exists (T4).** Q6's ``no gauge-invariant non-self-averaging carrier'' is false. The
   total-normalized zero mode `q=(1/√V)Σ(O−⟨O⟩)` has `Var(q)=χ=O(1)`, not `V^{-1}`; the `V^{-1}`
   statement is the mean-normalized mode, a bookkeeping artifact. The residual issues are UV
   renormalization at a fixed scale and the ceiling direction, not carrier existence.
2. **Binding is a mechanism, not an axis (T6).** On a Källén–Lehmann measure
   `I ⇒ B_spec ⇒ A` (both strict), so spectral binding is implied by isolation; dynamical binding
   (the flux tube) is the atom's mechanism. The independent axes are `A` (atom visibility) and `I`
   (isolation); the gap `G` is orthogonal and needs H1/ceiling plus a scale. The
   ``isolation without binding'' cell is empty.
3. **Isolation needs a bounded ratio (T3).** The window theorem: with `Δ=(ρ²−1)μ₁`, isolation plus
   `1+c≤ρ²≤C²` gives `μ₁≥δ/(C²−1)`; the double well has `ρ²→∞` and no gap. The upper bound `C`
   (IN-11, `C≈1.5`) is load-bearing.
4. **No standard expansion reaches β=2.5 (T1/T2/T5).** The polymer activity is `1.014` (`u⁴>1`), the
   Ursell ratio `4.6×10³`, the block-spin help `1.03×`; the tube overlap has only upper bounds; the
   center-vortex route is neutral-blind. The missing factor is `40–55×`.
5. **2D correction.** The Casimir ratio `ρ²=8/3` is the holonomy's isolation, not the `tr F²` atom
   (which is 0).

The corrected target is the three-predicate pair `(A, I)` with `G` orthogonal, plus the window upper
bound; the missing object is a **nonperturbative binding theorem** (uniform `w₁≥w₀>0`, equivalently
the leakage bound), not more variance control.

*Refs:* `T-WAVE-SYNTHESIS.md`; `T4-CARRIER.md`, `T6-BINDING-ADVERSARIAL.md`,
`T3-TP2-ISOLATION.md`, `T1-POLYMER-RG.md`, `T2-TUBE-OVERLAP.md`, `T5-CENTER-VORTEX.md`.
