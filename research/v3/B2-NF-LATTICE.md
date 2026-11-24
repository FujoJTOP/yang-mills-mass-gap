# B2 — What the lattice actually establishes about $N_f$, and what is monotone in the data

> **V3 attack surface B, second angle.** Scope: assemble the honest SU(3) $N_f$ data table
> ($N_f=0,2,4,6,8,10,12,16$), test monotonicity in $N_f$, write the corrected conformal-window
> sentence, and specify the finite-size-scaling protocol that would settle $N_f=12$.
>
> **VERDICT (no QED, Clay 0%).** Usable verified lattice data exist for **all eight** requested
> $N_f$ values. The data are **coarse-monotone**: gapped/chirally-broken at low $N_f$, gapless/
> IR-conformal at high $N_f$. The boundary is **bracketed, not located**. $N_f=12$ is **contested**
> and must not be reported as settled. The v2 sentence "SU(3) with $9\le N_f\le16$ is genuinely
> single-coupling AF gapless-conformal" is **over-claimed** and is corrected in §4. What is
> rigorous is only the **two-loop** statement: a Banks–Zaks zero exists for $N_f>306/38\approx8.05$
> and asymptotic freedom holds for $N_f<16.5$. No $N_f$-monotonicity theorem exists.
>
> Labels: `[T]` theorem · `[C]` conditional / scheme-dependent · `[V]` numerical lattice evidence ·
> `[O]` exact missing step. `[V]` data are finite-$(a,L,m)$ unless stated. No claim here upgrades
> lattice evidence to a continuum theorem.

---

## §0 · One-paragraph summary

The lattice knows the *phases* well and the *boundary* poorly. For $N_f=0$ (pure YM) the theory is
confining and gapped; for $N_f=2,4,6$ it is QCD-like with spontaneous chiral symmetry breaking
($\chi$SB); for $N_f=8$ it is walking / near-conformal, with $\chi$SB preferred but mass-deformed
IR conformality not excluded; for $N_f=16$ there is an IR fixed point (IRFP) with a nearly marginal
coupling. At $N_f=10$ a 2023 gradient-flow study finds an IRFP at $g^2_{GF}=15.0(5)$ with
$\gamma_m\simeq0.6$, while an earlier domain-wall study placed the theory at the strongly-coupled
edge with $\gamma_*\simeq1$ and did not exclude $\chi$SB. At $N_f=12$ two continuum-oriented
programmes disagree: IR-conformal (A-F-N 2008, LatKMI 2012/2013, Peterson–Hasenfratz 2024) versus
below-the-window $\chi$SB (Fodor et al. 2011/2016). The single most useful logical observation is
that **the $N_f=10$ and $N_f=12$ claims cannot both be right if the window edge is unique and the
order parameter is monotone in $N_f$** (§3).

---

## §1 · The honest data table

Legend for "Verify": **V** = bibliographic data and the quoted claim checked this pass against the
arXiv abstract page or the era report's ledger; **V-2** = verified via a secondary citation only;
**U** = not verified this pass (do not quote as established). All lattice numbers are finite
$(a,L,m)$.

| $N_f$ | Observable(s) actually measured | Collaboration / key reference (arXiv) | Action / scheme | Sizes / spacings | Claimed status | Verify |
|---|---|---|---|---|---|---|
| **0** | $m_{0^{++}}$ spectrum, string tension $\sqrt\sigma$; confinement | Athenodorou–Teper, *JHEP* **12** (2021) 082, arXiv:2106.00364; Hasenfratz, PRD **80** 034505 (2009) uses $N_f{=}0$ as confining test | Wilson plaquette; continuum extrapolation per $N$, then $N\to\infty$ | multiple $\beta$; at $N=3$ roughly $\beta_W\approx5.7$–$6.9$; $N=2,\dots,12$ | **gapped / confining**, no IRFP `[V]`; $m_{0^{++}}\approx1.7$ GeV, $\sqrt\sigma\approx0.42$–$0.44$ GeV `[V]`; $m_{0^{++}}/\sqrt\sigma\approx3.7$–$3.8$ (inferred from ratio, not printed) `[V-2]` | V |
| **2** | $\langle\bar\psi\psi\rangle/F^3$ vs QCD; Goldstone spectrum | LSD, PRL **104** 071601 (2010), arXiv:0910.2224 | domain-wall (exact chiral symmetry at finite $a$) | several volumes; sizes not on abstract | **confining / $\chi$SB** `[V]`; $\langle\bar\psi\psi\rangle/F^3$ agrees with measured QCD `[V]` | V (sizes U) |
| **4** | $F_\pi$, $\langle\bar\psi\psi\rangle$ in chiral limit; $M_\pi^2$ vs $m_f$; $M_\rho/F_\pi$ | LatKMI, PoS LATTICE2015 (2016) 215, arXiv:1512.00957; Hasenfratz PRD **80** 034505 (2009) | HISQ/tree, $\beta=3.7$ | $(12,18),(16,24),(20,30)$; $0.01\le m_f\le0.05$ | **$\chi$SB** `[V]`: $F_\pi,\langle\bar\psi\psi\rangle\to$ nonzero; $M_\pi^2\propto m_f$; $M_\rho/F_\pi\approx6$ | V |
| **6** | $\langle\bar\psi\psi\rangle/F^3$ (enhancement vs $N_f{=}2$) | LSD, PRL **104** 071601 (2010), arXiv:0910.2224 | domain-wall | several volumes; sizes not on abstract | **confining / $\chi$SB**, condensate significantly enhanced vs $N_f{=}2$ `[V]` | V (sizes U) |
| **8** | light composite scalar; finite-size hyperscaling; meson spectrum; condensate | LatKMI, PRD **87** 094511 (2013), arXiv:1302.6859; LSD, PRD **90** 114502 (2014), arXiv:1405.4752 | HISQ (LatKMI); domain-wall (LSD) | several $a$ and $m$ (details in text) | **walking / near-conformal** `[V]`; LatKMI: $\chi$SB near chiral limit but FSH at larger $m$, $\gamma_m\sim1$; LSD: spectrum suggests $\chi$SB but mass-deformed IR conformality not strongly disfavored | V |
| **10** | $\beta$-function zero; $\gamma_m$; spectrum scaling | Hasenfratz–Neil–Shamir–Svetitsky–Witzel, PRD **108** L071503 (2023), arXiv:2306.07236; LSD/Kiskis, arXiv:1204.6000 (2012) | GF continuous RG + heavy Pauli–Villars (2023); domain-wall, Iwasaki $32^3\times64$ (2012) | 2023: sizes not on abstract; 2012: single volume $32^3\times64$ | **IRFP** `[V]` at $g^2_{GF}=15.0(5)$, $\gamma_m\simeq0.6$ (2023) **vs** "within or just outside the strongly-coupled edge", $\gamma_*\simeq1$, $\chi$SB not ruled out (2012) — **contested** | V |
| **12** | IRFP / $\beta$ zero; $\gamma_m$; light composite scalar; condensate; FSS | A-F-N, PRL **100** 171607 (2008), arXiv:0712.0609; Hasenfratz, PRD **82** 014506 (2010), arXiv:1004.1004; LatKMI, PRD **86** 054506 (2012), arXiv:1207.3060; LatKMI, PRL **111** 162001 (2013), arXiv:1305.6006; Fodor et al., PLB **703** 348 (2011), arXiv:1104.3124; Fodor et al., PRD **94** 091501 (2016), arXiv:1607.06121; Peterson–Hasenfratz, PRD **109** 114507 (2024), arXiv:2402.18038 | staggered stout (Fodor); nHYP staggered (P–H); HISQ (LatKMI); domain-wall (A-F-N/LSD) | A-F-N: $L/a=4$–$20$, $\beta\approx4.15$–$6.5$; others differ | **CONTESTED**: IR-conformal (A-F-N 2008: conformal for $12\le N_f\le16$; LatKMI 2012: $\gamma_m\sim0.4$–$0.5$; P–H 2024: IRFP $g^2_{GF}=6.60(62)$, $\gamma_{g*}=0.199(32)$) **vs** below-the-window $\chi$SB (Fodor 2011: "very low level of confidence in the conformal scenario"; Fodor 2016: no IRFP in $6.0<g^2<6.4$) | V |
| **16** | IRFP with nearly marginal coupling; $\beta$; mass scaling dimension | A-F-N, PRL **100** 171607 (2008), arXiv:0712.0609; Hasenfratz, PRD **80** 034505 (2009), arXiv:0907.0919 | MCRG (2-lattice matching); SF (A-F-N) | small MCRG lattices; sizes not on abstract | **IRFP** `[V]`, nearly marginal gauge coupling; A-F-N quote $g_*^2\approx0.47$ **in the SF scheme** (normalization differs from the standard $\beta_0,\beta_1$ convention of §4.1 — compare only within a scheme) | V |

### §1.1 Notes that the table cannot hold

- **$m/\sqrt\sigma$ is not available as a single cross-$N_f$ number.** The task lists it as a
  candidate observable; the literature does not report a universal $m_{PS}/\sqrt\sigma$ at matched
  fermion mass. The only clean continuum mass ratio is $m_{0^{++}}/\sqrt\sigma$ at $N_f=0$; for
  $N_f\ge2$ the useful universal-ish diagnostics are the chiral-limit behaviour of $F_\pi$ and
  $\langle\bar\psi\psi\rangle/F^3$, the ratio $M_\rho/M_\pi$, and the existence of a $\beta$-function
  zero. This is itself part of the honest answer: **there is no common ruler**.
- **"Chiral condensate" is not one number.** LSD's $\langle\bar\psi\psi\rangle/F^3$ (dimensionless,
  $F$ = Goldstone decay constant) is the only quantity measured across $N_f=2,\dots,8$ by one
  collaboration with one fermion formulation. It is a good monotonicity candidate (§3) but it is
  defined only in the $\chi$SB phase.
- **The $N_f=12$ "light scalar"** (LatKMI, PRL 111 162001) is a genuine lattice result but its
  interpretation as a (pseudo-)dilaton is model-dependent; it is not an order parameter.
- **$N_f=16$ is the cleanest conformal anchor**: the two-loop fixed point is weak, so the lattice
  and perturbation theory agree, and both MCRG and SF see the IRFP.

---

## §2 · What is universal and what is scheme-dependent

| Statement | Status | Why |
|---|---|---|
| Gapped/confining vs gapless/IR-conformal (the **phase**) | universal `[V]` | Phase structure is scheme-independent; it is what "the window" means. |
| Spontaneous $\chi$SB vs unbroken chiral symmetry at $m\to0$ | universal `[V]` | Order parameter $F_\pi>0$ vs $F_\pi\to0$; the dichotomy is physical. |
| Continuum dimensionless ratios $m_i/m_j$ | universal `[V]` (after $a\to0$) | Provided the continuum limit is controlled. |
| Fixed-point coupling $g_*^2$ | **scheme-dependent** `[C]` | $g^2_{GF}=15.0(5)$ ($N_f{=}10$) and $6.60(62)$ ($N_f{=}12$) are in different schemes; they cannot be compared directly. |
| Location of the $\beta$-function zero | **scheme-dependent** `[C]` | Changes under scheme redefinitions; only its existence is scheme-independent. |
| $\gamma_m$ at a genuine IRFP | universal **in principle** (CFT operator dimension) `[C]` | But every lattice "$\gamma_m$" here is extracted at finite $m$ and finite $L$, not at the fixed point. The spread $0.45$–$1.0$ at $N_f{=}12$ is therefore a **diagnostic that the fixed point is not resolved**, not a controlled universality test. |
| Leading irrelevant exponent $\gamma_{g*}$ | universal at a fixed point `[C]` | If the IRFP is real, its critical exponent is universal; the quoted $0.199(32)$ is the relevant number. |

**Consequence for monotonicity.** A monotonicity test on universal data must use the *phase*
dichotomy (or continuum mass ratios), **not** $g_*^2$, **not** the zero location, and **not** the
finite-mass $\gamma_m$. This is exactly why the clean test below is categorical.

---

## §3 · The monotonicity test

Witness: `research/v3/w7/b2_nf_monotone.py` (numpy only; run with `py -3 -I`). Output reproduced
here.

### §3.1 Categorical order parameter

Define the universal dichotomous order parameter
`G(N_f) = 1` (gapped / confining; $\chi$SB for $N_f\ge2$) and `G(N_f) = 0` (gapless / IR-conformal).
From the table:

```
Nf :   0  2  4  6  8  10    12    16
G  :   1  1  1  1  1 {0,1} {0,1}   0
```

With the two contested values allowed to float, there are 4 consistent assignments:
**3 are monotone** ($G$ non-increasing) and **1 is not**. The unique non-monotone assignment is

```
G = 1-1-1-1-1-0-1-0   (Nf=10 conformal, Nf=12 chiSB, Nf=16 conformal)
```

### §3.2 What is consistent with monotonicity

- **Coarse monotonicity holds.** $G=1$ for $N_f=0,2,4,6,8$ and $G=0$ for $N_f=16$; there is no
  verified data point that contradicts "gapped at small $N_f$, gapless at large $N_f$."
- **Two monotone scenarios survive**, differing only in where the edge sits:
  1. **Edge in $(8,10]$:** $N_f=10$ conformal (2023), $N_f=12$ conformal — Fodor's $N_f=12$ $\chi$SB
     claim is wrong (or its action/scheme is not in the same universality class).
  2. **Edge in $(12,16]$:** $N_f=10$ and $N_f=12$ both gapped/$\chi$SB — the 2023 $N_f=10$ IRFP is
     wrong or borderline, and Fodor is right at $N_f=12$.
- Both are compatible with the A-F-N bracket "$8\le N_f^c\le12$"; neither pins the boundary.

### §3.3 What contradicts monotonicity

- **The sharp logical point.** If $N_f=10$ is conformal (Hasenfratz et al. 2023) **and** the
  window edge is unique and the order parameter is monotone in $N_f$, then $N_f=12$ **must** be
  conformal. Fodor's $N_f=12$ $\chi$SB claim would then force either (i) $N_f=10$ below the
  window, (ii) a non-universal / action-dependent phase, or (iii) a failure of monotonicity.
  This is the cleanest way to state why the $N_f=12$ controversy is not a detail: **it is a
  consistency test on the whole $(N_f=10,N_f=12)$ pair.** The witness prints this explicitly.

### §3.4 $\gamma_m$ is not a usable monotone

Quoted values (scheme-dependent, only defined where an IRFP exists):

```
Nf= 6 : 0.40 [LSD, occasionally quoted; not in the abstract]
Nf= 8 : 1.00 [LatKMI FSH fit]   |  0.94(2) [LSD data, secondary citation]
Nf=10 : 0.60 [Hasenfratz et al. 2023]  |  1.00 [LSD 2012, edge]
Nf=12 : 0.45 [LatKMI 2012, 0.4-0.5]   |  1.00 [Hasenfratz 2010, near free-field]
Nf=16 : no clean number (perturbatively small)
```

The LatKMI trend $N_f=8\to12$ ($1.0\to0.45$) *looks* monotone, but the $N_f=12$ value alone spans
$0.45$–$1.0$ depending on group, and $N_f=10$ spans $0.6$–$1.0$. **$\gamma_m$ is not monotone in
$N_f$ and is not a usable ruler.** (This matches LIT-ERA-09's finding that Hasenfratz 2010 measured
$\gamma_m$ "barely different from the free-field value," i.e. $\approx1$, at $N_f=12$.)

### §3.5 What a monotonicity *theorem* would require

To convert the categorical statement into a theorem one needs, for a universal IR observable
$\mathcal{O}$ (gap, $\gamma_m$ at a genuine fixed point, or a spectral weight):

1. a proof that $\mathcal{O}(N_f)$ is monotone (in the correct direction — e.g. the gap is
   non-increasing in $N_f$); **and**
2. a rigorous Banks–Zaks IRFP at some large $N_f^+$ (perturbatively accessible), so that
   $\mathcal{O}(N_f^+)=0$; **then**
3. $N_f=0$ lies on the gapped side without computing the strong-coupling $\beta$.

**No such monotonicity theorem is known** (LIT-ERA-09 §7.3, V3-LIT-SYNTHESIS §5B). The lattice data
are *consistent* with monotonicity and *cannot prove* it. This is the correct `[O]` status.

---

## §4 · The corrected v3 conformal-window wording

### §4.1 What is rigorous `[T]`

Only the two-loop polynomial statement. With $\beta(g)=-\beta_0 g^3-\beta_1 g^5$ (per $16\pi^2$),
\[
\beta_0=11-\frac{2N_f}{3},\qquad \beta_1=102-\frac{38N_f}{3}.
\]
For SU(3): asymptotic freedom for $N_f<16.5$; a Banks–Zaks zero $g_*^2=-\beta_0/\beta_1>0$ for
$N_f>306/38\approx8.0526$. Hence the **two-loop window is $N_f\in\{9,\dots,16\}$**. This is a
theorem about the truncated series (the coefficients $\beta_0,\beta_1$ are universal). It is **not**
a theorem about the exact nonperturbative $\beta$-function; near $N_f=8$–$10$ the coupling is
$O(1)$ and the truncation is uncontrolled.

### §4.2 What is lattice-supported `[V]`

- IRFP at $N_f=16$ (A-F-N 2008; Hasenfratz 2009, MCRG) — the cleanest anchor.
- $\chi$SB / confinement at $N_f=0,2,4,6$ (Athenodorou–Teper; LSD; LatKMI $N_f=4$).
- Walking / near-conformal at $N_f=8$ (LatKMI 2013; LSD 2014), with $\chi$SB preferred but
  mass-deformed IR conformality not excluded.
- An IRFP at $N_f=10$ (Hasenfratz et al. 2023) — strong but recent, and in tension with the 2012
  edge interpretation.

### §4.3 What is contested at $N_f=12$ `[V-contested]`

- **Pro-conformal:** A-F-N 2008 (conformal for $12\le N_f\le16$), LatKMI 2012 ($\gamma_m\sim0.4$–
  $0.5$), LatKMI 2013 (light composite scalar), Peterson–Hasenfratz 2024 (IRFP at
  $g^2_{GF}=6.60(62)$, $\gamma_{g*}=0.199(32)$).
- **Contra-conformal:** Fodor et al. 2011 ("very low level of confidence in the conformal
  scenario"), Fodor et al. 2016 (no $\beta$-function fixed point in $6.0<g^2<6.4$ in their scheme).
- The discrepancy tracks **fermion formulation (staggered/stout vs nHYP vs domain-wall) and
  scheme**, not statistics. Neither side is a theorem. **Do not present one side as settled.**

### §4.4 The corrected sentence (replace the v2 sentence)

**Short form (one sentence, for the paper):**

> Rigorously, only the two-loop statement holds — a Banks–Zaks zero exists for
> $N_f>306/38\approx8.05$ and asymptotic freedom holds for $N_f<16.5$, so "the window is
> $9\le N_f\le16$" describes a truncated $\beta$-function, not an established phase; lattice data
> support an IR fixed point at $N_f=16$ and chiral-symmetry breaking for $N_f\le8$, while $N_f=12$
> is contested between IR-conformal and below-the-window determinations, and no $N_f$-monotonicity
> theorem exists, so the boundary is bracketed, not located.

**Replacement for the specific v2 clause** ("SU(3) with $9\le N_f\le16$ is genuinely
single-coupling AF gapless-conformal"):

> SU(3) with $N_f$ in the perturbative Banks–Zaks range is a **single-coupling asymptotically-free
> theory whose two-loop $\beta$-function has an infrared zero**. Whether it is nonperturbatively
> gapless-conformal is lattice-supported at $N_f=16$, **contested at $N_f=12$**, and **not
> established as a continuum statement for the whole $9$–$16$ range**.

**Why the acid test is unharmed.** G1-ACID only needs *one* single-coupling AF gapless-conformal
example to kill "single coupling $\Rightarrow$ no zero" and "AF $\Rightarrow$ gapped". The $N_f=16$
IRFP — where lattice and two-loop perturbation theory agree — supplies it. The correction removes
the over-claim that *every* $N_f$ in $9$–$16$ is established conformal; it does not weaken the
acid test's conclusion.

### §4.5 The "$9\le N_f\le16$" phrase specifically

The v2 phrase conflates three different objects:
1. the **integer span of the two-loop BZ zero** ($9\le N_f\le16$) — `[T]` about the polynomial;
2. the **lattice-bracketed nonperturbative window** — edge somewhere in $8$–$12$ (A-F-N), or
   possibly $10$–$12$; `[V]`, not located;
3. the **claim that each integer in $9$–$16$ has been shown conformal** — false as stated;
   $N_f=12$ is contested and $N_f=10$ is recent/contested.

Use (1) as the rigorous statement, (2) as the evidence statement, and never (3).

---

## §5 · The finite-size-scaling protocol that would settle $N_f=12$

The goal is to decide between two mutually exclusive continuum statements: (a) the massless
$N_f=12$ theory has an IRFP (gapless, $\chi$-symmetric); (b) it confines and breaks chiral symmetry
(gapped). Both must be tested in the **same action and scheme as a function of $N_f$**, so that the
monotonicity consistency test of §3 becomes decisive.

### §5.1 Observables

1. **Renormalized running coupling / $\beta$-function.** Step-scale the coupling by a factor $s$ in
   finite volume and take the continuum limit at each step; look for a zero of the step-scaling
   function. Use a mass-independent scheme; gradient-flow (GF) coupling (as in Peterson–Hasenfratz)
   or Schrödinger-functional (SF) coupling (as in A-F-N). This is the primary observable: it decides
   (a) vs (b) directly.
2. **Chiral-limit order parameter.** $F_\pi$ and $\langle\bar\psi\psi\rangle$ extrapolated to
   $m\to0$: nonzero $\Rightarrow$ $\chi$SB; zero with hyperscaling $\Rightarrow$ conformal. Also
   $M_\pi^2\propto m$ (GMOR/$\chi$SB) versus $M_\pi\propto m^{1/(1+\gamma)}$ (conformal).
3. **Finite-size scaling of the spectrum.** At fixed bare coupling and $m=0$, measure $M_\pi L$,
   $F_\pi L$, $M_\rho L$ as functions of $L$. Conformal: spectrum $\propto L^{-1}$ with a mass
   anomalous dimension; $\chi$SB: $F_\pi L$ stays finite/diverges and $M_\pi L\to0$.
4. **Mass anomalous dimension $\gamma_m$.** From the $m$-dependence of the same correlators; must be
   volume-stable. This is a consistency check, not the decider.
5. **Direct $\beta$-function from the GF** (continuous RG), as a cross-check on the step-scaling zero.

### §5.2 Volumes and lattice spacings

- **Volumes:** at least **four** linear sizes spanning a factor $\ge2$, e.g. $L/a=12,16,20,24$ and
  ideally $32$, with a fixed aspect ratio $T/L\ge2$ (or $4$). A single volume is not enough to
  demonstrate FSS.
- **Lattice spacings:** at least **three** per coupling, to extrapolate $a^2/L^2\to0$ with a
  controlled $O(a^2)$ (or $O(a)$ improved) fit. The continuum extrapolation must be simultaneous with
  the FSS fit, not sequential.
- **Fermion masses:** several values spanning a factor $\sim4$–$5$, reaching $m L\ll1$ so the
  finite-size regime is genuinely probed, not the heavy-mass regime where all theories look alike.
- **Statistics:** comparable to current state of the art; the decider is systematics (volume,
  $a\to0$, scheme), not raw statistics.

### §5.3 Scheme and scale setting

- **Scheme:** a mass-independent renormalization scheme defined in the continuum (SF or GF). State
  it explicitly, because $g_*^2$ values are scheme-dependent and the Fodor vs Peterson–Hasenfratz
  discrepancy is partly scheme/action.
- **Scale setting:** $t_0$ or $w_0$ from the gradient flow, plus a hadronic scale, with the same
  definition across all $N_f$ in the scan. Match to $\sqrt\sigma$ or $F_\pi$ in the $\chi$SB phase
  so that different $N_f$ are compared at the same physical scale.
- **Common action scan:** the decisive experiment is a **single collaboration, one action, one
  scheme**, scanning $N_f=8,10,12,14,16$ with identical scale setting. That turns the §3
  monotonicity consistency test into a measurement.

### §5.4 Precision targets

To separate (a) from (b) at $N_f=12$:

- **$\beta$-function:** a zero resolved at **$\ge5\sigma$** over a factor $\ge2$ in $L$, with a
  statistical error on $g^2$ of $\lesssim1\%$ and a continuum-extrapolation error **below the
  difference between "zero" and "no zero."** Currently the two programmes differ at the $>2\sigma$
  level (no zero in $6.0<g^2<6.4$ vs a zero at $6.60(62)$); that gap is the precision target. A
  $\ge5\sigma$ zero, or a $\ge5\sigma$ exclusion over an overlapping coupling range, would settle it.
- **$\gamma_m$:** stable across volumes to $\sim5$–$10\%$. (The present spread $0.45$–$1.0$ is far
  above this and signals the fixed point is unresolved.)
- **$F_\pi$ and $M_\pi$:** continuum values at the **$3$–$5\%$** level in the chiral limit, with the
  $m\to0$ extrapolation controlled by data at $mL\ll1$.
- **Volume independence:** the running coupling must be shown **volume-independent** (within errors)
  over the candidate fixed-point plateau; a coupling that keeps running with $L$ is walking, not
  conformal.
- **Decision rule:** conformal iff (i) the step-scaling $\beta$-function has a zero, (ii) the
  coupling is volume-independent at the zero, and (iii) $F_\pi\to0$ with hyperscaling; $\chi$SB iff
  (i) no zero over the accessible range, (ii) $F_\pi>0$ in the chiral limit, (iii) $M_\pi^2\propto m$.
  Report the outcome per scheme and per action, and **do not** average across formulations before
  the continuum limit.

### §5.5 What would still be missing after the protocol

Even a clean $N_f=12$ answer would not prove the window boundary is at $N_f=12$, because the
protocol settles one $N_f$. The boundary requires the common-action scan (§5.3) plus the
**monotonicity theorem** (§3.5), which remains `[O]`. The protocol raises $N_f=12$ from contested
to `[V]`; it does not produce a theorem.

---

## §6 · Verification ledger

`V` = checked this pass against the arXiv abstract page (fetched) or the era report's verified
ledger. `V-2` = secondary citation. `U` = unverified, not load-bearing.

| # | Entry | Identifier checked | Checked against | Status | Note |
|---|---|---|---|---|---|
| B2-01 | Athenodorou–Teper, glueball/string tension, $N_f{=}0$ | arXiv:2106.00364; JHEP 12 (2021) 082 | arXiv abs + LIT-ERA-13 ledger | V | $m_{0^{++}}\approx1.7$ GeV, $\sqrt\sigma\approx0.44$ GeV |
| B2-02 | LSD, "Toward TeV Conformality" ($N_f=2$–$6$) | arXiv:0910.2224; PRL 104 071601 | arXiv abs | V | $N_f{=}2$ matches QCD; $N_f{=}6$ enhanced |
| B2-03 | LatKMI, four-flavor SU(3) | arXiv:1512.00957; PoS LATTICE2015 215 | arXiv/INSPIRE + PoS | V | HISQ, $\beta=3.7$; $\chi$SB |
| B2-04 | LatKMI, $N_f=8$ walking | arXiv:1302.6859; PRD 87 094511 | arXiv abs + LIT-ERA-09 | V | FSH, $\gamma_m\sim1$ |
| B2-05 | LSD, $N_f=8$ domain wall | arXiv:1405.4752; PRD 90 114502 | arXiv abs (fetched) | V | $\chi$SB preferred; IR conformality not strongly disfavored |
| B2-06 | Hasenfratz et al., $N_f=10$ IRFP | arXiv:2306.07236; PRD 108 L071503 | arXiv abs (fetched) + search | V | $g^2_{GF}=15.0(5)$, $\gamma_m\simeq0.6$ |
| B2-07 | LSD/Kiskis, "Approaching Conformality with Ten Flavors" | arXiv:1204.6000 | arXiv/ar5iv + search | V | $\gamma_*\simeq1$; single volume $32^3\times64$ |
| B2-08 | A-F-N, conformal window | arXiv:0712.0609; PRL 100 171607 | arXiv abs (fetched) + search | V | conformal $12\le N_f\le16$; $L/a=4$–20 |
| B2-09 | Hasenfratz, MCRG $N_f=8,12$ | arXiv:1004.1004; PRD 82 014506 | arXiv abs (fetched) | V | IRFP or slow walking; $\gamma_m$ near free-field |
| B2-10 | LatKMI, $N_f=12$ conformality | arXiv:1207.3060; PRD 86 054506 | arXiv abs (fetched) | V | $\gamma_m\sim0.4$–$0.5$ |
| B2-11 | LatKMI, light scalar $N_f=12$ | arXiv:1305.6006; PRL 111 162001 | LIT-ERA-09 ledger | V | light flavor-singlet scalar |
| B2-12 | Fodor et al., $N_f=12$ below window | arXiv:1104.3124; PLB 703 348 | arXiv abs (fetched) | V | "very low confidence in conformal scenario" |
| B2-13 | Fodor et al., fate of fixed point | arXiv:1607.06121; PRD 94 091501 | LIT-ERA-13 ledger | V | no IRFP in $6.0<g^2<6.4$ |
| B2-14 | Peterson–Hasenfratz, $N_f=12$ IRFP | arXiv:2402.18038; PRD 109 114507 | arXiv abs (fetched) + LIT-ERA-13 | V | $g^2_{GF}=6.60(62)$, $\gamma_{g*}=0.199(32)$ |
| B2-15 | Hasenfratz, MCRG $N_f=16$ | arXiv:0907.0919; PRD 80 034505 | arXiv abs (fetched) | V | IRFP, nearly marginal |
| B2-16 | $m_{0^{++}}/\sqrt\sigma\approx3.7$–$3.8$ | inferred from B2-01 | ratio $\times$ $\sqrt\sigma$ | V-2 | not printed directly in retrieved text |
| B2-17 | $\gamma_m\approx0.4$ at $N_f=6$ | LSD 2010 | not in abstract | U | frequently cited; not quoted as established |
| B2-18 | $\gamma_*=0.94(2)$ at $N_f=8$ | LSD data, secondary dilaton-ChPT citation | secondary | V-2 | primary source not fetched |
| B2-19 | Lattice sizes for $N_f=2,6$ (LSD), $N_f=16$ (MCRG) | — | not on abstract | U | do not quote as established |

**Corrections made this pass (record; do not propagate the wrong forms):**
- E09-32's Fodor et al. citation **is now verified**: arXiv:**1104.3124**, PLB **703** (2011) 348
  (the era report had left the arXiv ID unconfirmed).
- The "$9\le N_f\le16$" window is the **integer span of the two-loop BZ zero**, not an established
  nonperturbative phase (§4.5).
- `N_f=12` is **contested**; the correct status is `[V-contested]`, not `[V]` for either side.

---

## §7 · Position and what would change the answer

- **The lattice map is coarse-monotone but boundary-blind.** For the campaign's purpose the useful
  content is: the phase dichotomy is real and scheme-independent; the boundary is bracketed by
  $N_f=8$ (gapped) and $N_f=16$ (conformal); $N_f=12$ is the contested point; $N_f=0$ is far below.
- **The v3 over-claim is removed.** v2's "genuinely single-coupling AF gapless-conformal" for the
  whole $9$–$16$ range must be replaced by the corrected sentence in §4.4.
- **The acid test survives.** One confirmed AF single-coupling gapless-conformal example ($N_f=16$)
  is enough to kill the kinematic shortcut; the correction does not reopen it.
- **What would change the answer:** a $\ge5\sigma$ common-action determination at $N_f=12$ (§5)
  would move it from contested to `[V]`; a proof of an $N_f$-monotone universal IR observable plus a
  rigorous large-$N_f$ BZ IRFP (§3.5) would be a genuine new tool — and is the only route in this
  coordinate that could put $N_f=0$ on the gapped side without computing the strong-coupling
  $\beta$-function.

**Clay = 0%. No QED.** The conformal window is the map of the *alternative* (a gapless CFT IR), not
a shortcut to the gap. Settling $N_f=12$ sharpens the map; it does not close the gap.
