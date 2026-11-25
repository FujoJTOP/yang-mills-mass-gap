# Z1-HOLO — can "5D normal mode ⇒ 4D mass gap" be extracted as a rigorous statement for pure SU(N) YM? **No: the implication is a real [T] *within a fixed bulk model* (Sturm–Liouville lifting + GKPW), but for YM it assumes the IR cap (= the gap) and reduces to H1. Acid test on the free massless field passes (pure AdS is correctly gapless), which is exactly why the construction's only "gap" is the input scale \(z_0\) / \(r_+\).**

> **VERDICT (no QED, Clay 0%).** The precise candidate statement — *"a normalizable 5D bulk mode with a positive bulk spectral gap implies a 4D mass gap"* — decomposes into a **provable part** and a **fatal part**. The provable part is a genuine theorem **inside a fixed asymptotically-AdS background with a regular IR cap**: by Sturm–Liouville spectral theory the bulk mode operator \(L\) has a discrete spectrum \(0<m_1^2\le m_2^2\le\dots\), and the GKPW dictionary maps its normalized modes to Källén–Lehmann atoms of the boundary two-point function, with **no spectral weight below \(m_1\)** — hence a boundary gap \(m_1>0\), at leading supergravity order **[T] within the model**. The fatal part is the extraction to 4D pure SU(N) YM: (i) the discreteness \(+\) positivity of \(\operatorname{spec}(L)\) exists **iff an IR cap is present**, and the cap's scale (\(z_0\) in the hard wall, \(r_+\) in BTZ, \(M_{KK}\) in the Witten cigar) **is the confinement scale** — insert it and you have inserted the gap (**"assumes the gap"**); (ii) the bulk must be dual to pure YM — **no such dual is known** (**only a model, not a theorem**); (iii) the gap equals the cap/compactification scale and **\(\to0\) under decompactification**, so it is not the uniform 4D gap (**reduces to H1**). All three kill criteria fire. The net is **[NG]** for the route, with the exact missing step being H1 in the bulk-geometry coordinate.

---

## §0 · Labels and the acid test

| Label | Meaning (as used in this campaign) |
|---|---|
| **[T]** | Theorem / exact identity (here: *within a stated bulk model*, Sturm–Liouville + leading-order GKPW). |
| **[C]** | Conditional / conjecture (the duality map itself, for a non-SUSY target). |
| **[NG]** | No-go / obstruction for the class. |
| **[O]** | Observation / exact missing step / uncontrolled model input. |
| **[V]** | Verified numerical fact (reproduced by the witness this round). |

**Acid test (the charter's rule).** A gap detector must *fail* — return no gap — on the free massless field (`f_1^2=0`, gapless), else it is vacuous. The holographic lifting chain is run on the free massless field in §5: the bulk dual is **pure AdS** (no cap), \(\operatorname{spec}(L)=[0,\infty)\) **continuous**, and the detector correctly returns **no gap**. So the chain is *not* vacuous in the naive sense. It is vacuous in the deeper **"assumes-the-conclusion"** sense: its only gapped output occurs precisely where a finite cap \(z_0,r_+\) was inserted, and the "detected" gap is *literally that inserted scale*. Both facts are recorded honestly; the second is decisive (charter §1: a tool that "detects" a gap by assuming the cap is the confining scale is the H1 object renamed).

---

## §1 · The precise statement under test, and why its wording is the whole problem

Make the phrase exact. Let \(M\) be a \((d+1)\)-dimensional asymptotically-AdS bulk (\(d=4\)), with a bounded-above, smooth-in-space metric. Let \(\phi\) be the bulk field dual to a gauge-invariant \(0^{++}\) boundary operator \(O\), and linearize about the background: \(\partial_M(\sqrt g\,g^{MN}\partial_N\phi)=0\). For a mode \(\phi(z,x)=\psi(z)e^{ip\cdot x}\), \(p^2=-m_4^2\), the radial problem is a **Sturm–Liouville (SL) eigenvalue problem**
\[
L\,\psi \;=\; m_4^2\,\psi,\qquad
L=-\,\frac{1}{w(z)}\frac{d}{dz}\!\Big(w(z)\,g^{zz}(z)\frac{d}{dz}\Big)+V(z),
\tag{1}
\]
on \(L^2([\,z_{\rm UV},z_{\rm IR}],w\,dz)\), with a UV normalizability condition and an IR boundary condition fixed by the geometry.

**Two things must be true for "5D normal mode \(\Rightarrow\) 4D mass gap," and they are logically different:**

- **(N) A normal mode exists with a positive *bulk* spectral bottom.** Not merely *one* normalizable mode: you need \(\operatorname{spec}(L)\) to have a **positive lowest eigenvalue**, i.e. \(L\) *gapped* and (for isolated atoms) *discrete*. A single normal mode sitting above a continuum gives an atom with spectral weight below it — **no gap**. So the precise hypothesis is *"a gapped component of the bulk spectrum with a normalizable lowest mode."*
- **(D) The dictionary lifts it faithfully.** GKPW: the boundary two-point function is the boundary limit of the bulk (renormalized) Green function; its Källén–Lehmann spectral measure \(\nu\) has an atom at each normal-mode eigenvalue and no weight below the bulk gap \(\Delta_5\).

The whole content of "5D normal mode \(\Rightarrow\) 4D mass gap" is **(N)\(\wedge\)(D) \(\Rightarrow\) (4D gap)**. The next two sections pin the label of each.

---

## §2 · What is genuinely a theorem: the SL-lifting lemma [T] within a fixed model

> **Lemma Z (SL lifting; [T] inside a fixed bulk with a regular IR cap).**
> Let \(M\) have a regular IR cap so that (1) has a purely discrete spectrum
> \(0<m_1^2\le m_2^2\le\cdots\to\infty\) with normalized eigenfunctions \(\psi_n\).
> Then, at leading supergravity order, the boundary gauge-invariant correlator admits
> \[
> \langle O(x)O(0)\rangle=\int_0^\infty\! d\nu(m)\,D(x;m),\qquad
> d\nu=\sum_n c_n\,\delta(m^2-m_n^2),\ \ c_n\ge0,
> \]
> i.e. a Källén–Lehmann measure whose **support is bounded below by \(m_1^2\)**;
> consequently the boundary \(0^{++}\) channel is gapped with gap \(\Delta_4=m_1\).
> If instead \(M\) has no IR cap (pure AdS; the IR is the Poincaré horizon / conformal
> boundary), \(\operatorname{spec}(L)=[0,\infty)\) is continuous and the boundary channel
> is **gapless**.

*Proof content.* Self-adjointness of \(L\) on the cap domain gives a complete orthonormal eigenbasis (compact resolvent, discrete spectrum). Expanding the bulk Green function in \(\{\psi_n\}\) and taking the boundary limit via GKPW yields the KL sum with \(c_n\propto|\psi_n(z_{\rm UV})|^2\ge0\); the support is exactly \(\{m_n^2\}\), bottom \(m_1^2>0\). The no-cap case is the standard continuous-spectrum limit. \(\square\)

This is real, and it is the honest positive residue of the holographic paradigm: **the mechanism by which an IR gap in the bulk becomes a 4D mass gap is a theorem *about the bulk SL operator*.** It is exactly the "mechanism" the LIT-ERA-08 file credits to holography.

**But note what Lemma Z does *not* say.** It says: *given* a bulk with a **discrete, gapped** SL operator, the boundary is gapped. It does **not** say a bulk dual of pure YM *has* such an operator. The hypothesis of Lemma Z — a regular IR cap — is not a property of the pure-gravity (no-matter) bulk, which is **pure AdS = an exact CFT = gapless**. It is a property of a *confining deformation* of the bulk, i.e. it is the confinement scale inserted as geometry.

---

## §3 · Hard-wall AdS\(_5\): the "gap" is exactly the input wall position \(z_0\) [V]

Hard wall (Polchinski–Strassler–type): AdS\(_5\), \(ds^2=(L^2/z^2)(\eta_{\mu\nu}dx^\mu dx^\nu+dz^2)\), cut at \(z=z_0\). The \(0^{++}\) (dilaton) mode solves
\[
\phi''-\frac{3}{z}\phi'+m^2\phi=0
\quad\Longleftrightarrow\quad
\partial_z\!\big(z^{-3}\partial_z\phi\big)+m^2 z^{-3}\phi=0,
\tag{2}
\]
with normalizable UV behavior \(\phi\sim z^2\) and IR condition \(\phi(z_0)=0\). The exact solution \(\phi=z^2 J_2(mz)\) gives
\[
\boxed{\,m_n=\frac{j_{2,n}}{z_0}\,},\qquad j_{2,n}=\{5.13562,\,8.41724,\,11.61984,\,14.79595,\dots\}.
\tag{3}
\]
The witness verifies (2) to residual \(2.7\times10^{-7}\), the zeros (3), the ratios \(m_n/m_1=\{1,\,1.639,\,2.263,\,2.881\}\) **[V]** (the model's genuine success: lattice \(0^{++*}/0^{++}\approx1.6\)–\(1.8\)).

> *Erratum to the corpus, checked this round (charter §2).* Earlier files (`ROOF-HOLOGRAPHY.md`, `route_new_holo.py`) print the wave equation as \(\partial_z(z^{-3}\partial_z\phi)+m^2 z^{-5}\phi=0\) while also using \(\phi=z^2J_2(mz)\). Those are inconsistent: \(z^2J_2(mz)\) solves **(2)** (prefactor \(z^{-3}\)), not \(z^{-5}\) — verified numerically here (residual \(2.7\times10^{-7}\) for (2)). The **spectrum (3) they report is the correct one**, since it is derived from the \(z^{-3}\) equation they actually use via \(\phi=z^2J_2\); the \(z^{-5}\) in the displayed equation was a typo (a genuinely \(z^{-5}\) equation has a different solution and spectrum), not a different-model result. Recorded, since this file builds on those.

**The point.** The gap is
\[
\Delta_4 \;=\; m_1 \;=\; \frac{j_{2,1}}{z_0}.
\tag{4}
\]
The witness tabulates \(m_1(z_0)\): \(z_0=0.5\to10.27\), \(1\to5.136\), \(50\to0.103\), \(10^6\to5.14\times10^{-6}\) **[V]**. Hence
\[
\inf_{z_0\in(0,\infty)} m_1(z_0)=0.
\tag{5}
\]
Removing the wall (\(z_0\to\infty\), the decompactification / pure-AdS limit) sends the gap to zero. **The SL operator (2) does not *produce* the gap; the gap *is* \(z_0\), the IR scale inserted into (2) by hand.** Bounding the 4D gap below is, term by term, bounding \(z_0\) above — which is H1 in bulk coordinates. This is the precise content of "assumes the gap."

---

## §4 · BTZ (AdS\(_3\)): the same structure, gap \(=\) thermal/horizon scale \(r_+\) [V]

BTZ is the cleanest *exactly solvable* holographic model (a quotient of AdS\(_3\); the dual is a 2D CFT at temperature \(T=r_+/(2\pi L^2)\)). Metric \(ds^2=-\frac{r^2-r_+^2}{L^2}dt^2+\frac{L^2}{r^2-r_+^2}dr^2+r^2d\varphi^2\), \(\varphi\sim\varphi+2\pi\).

Scalar mass–dimension (AdS\(_3\)): \(\Delta(\Delta-2)=m^2L^2\). The non-rotating BTZ scalar quasinormal modes (Cardoso–Lemos; the AdS/CFT QNM prescription is Son–Starinets) have an imaginary part that grows linearly in the overtone \(n\) with the **thermal scale** \(r_+/L^2\):
\[
\operatorname{Im}\omega_n = -\frac{r_+}{L^2}\,(2n+\Delta')\;=\;-(2\pi T)(2n+\Delta'),\qquad
\operatorname{Re}\omega=\pm k ,
\tag{6}
\]
where \(\Delta'\) is a field-dependent \(O(1)\) offset (the conformal weight, up to the standard normalization; its exact value is convention-dependent and **irrelevant here**). The witness confirms the **dimensionless ratios \(\operatorname{Im}\omega_n/(2\pi T)=-(2n+\Delta')\) are \(T\)-independent** for \(\Delta'=2,3\) **[V]**. Therefore **every** frequency in the tower is a pure number times \(2\pi T=r_+/L^2\):
\[
\omega_n \;\propto\; \frac{r_+}{L^2}\;\xrightarrow[\;r_+\to0\;]{}\;0 .
\tag{7}
\]
At \(r_+\to0\) (zero temperature) BTZ \(\to\) pure AdS\(_3\), whose dual is a **gapless** 2D CFT (the free massless scalar: continuous spectrum down to 0). So the BTZ "gap" (the lowest excitation above the ground state) is the **horizon/thermal scale** — an input boundary condition (the identification of the \(t\) circle) — and it **vanishes in the zero-temperature limit**. Same structure as §3 with \(z_0\leftrightarrow r_+\). BTZ adds nothing beyond confirming the pattern: *discrete tower \(\Leftrightarrow\) a finite compactification/thermal scale is inserted; the scale is an input.*

> The two solvable models thus bracket the two ways the SL gap can be "not the Clay gap": (hard wall) an **IR wall put in at \(z_0\)**; (BTZ) a **thermal circle of size \(1/T\)** — both are *inputs*, neither is a dynamical zero-temperature mass generation.

---

## §5 · The free-massless acid test: pure AdS is correctly gapless [T]

Run the same machinery on the free massless field. Its holographic dual (the only consistent one) is a **massless bulk field in pure AdS** — no cap, no wall, no thermal circle. The radial SL problem retains only the UV normalizability condition; the normalizable branch \(\phi\sim z^2J_2(mz)\) exists **for every \(m>0\)** (Bessel oscillates for all \(z\); infinitely many zeros — witness prints the first 12 **[V]**), so \(m\) is **not quantized**: \(\operatorname{spec}(L)=[0,\infty)\), no positive lower bound. By Lemma Z's second clause the boundary channel is a **gapless CFT**. The detector returns **no gap**. \(\checkmark\)

So the chain **passes** the acid test: it does not falsely manufacture a gap for the free massless field. This is important — it means the objection is **not** "the tool is trivially true/vacuous on everything." The objection is sharper and fatal: **the tool's only gapped outputs occur exactly when a finite IR scale (\(z_0\), \(r_+\), \(M_{KK}\)) was inserted, and its "gap" is that scale** (§3–4). The detector is a faithful *readout of the bulk boundary condition*, and the boundary condition is the confinement scale. Vacuity here is **"assumes-the-conclusion,"** the kill criterion the task names, not the naive free-field kind.

---

## §6 · The three kill criteria, applied

**(K1) "It is only a model, not a theorem."** Fires. Lemma Z is a theorem **only** about a *fixed* bulk SL operator. The statement about *4D pure SU(N) YM* needs hypothesis (i) of §1: *the bulk is the holographic dual of pure YM.* No such dual is known (LIT-ERA-08 §3.5: no string dual of 4D pure YM; all exact duals are \(N=4\) SYM or \(N=1\) cascades). For pure YM the GKPW map is a **conjecture [C]** at best, and there is not even a candidate geometry. So the "5D normal mode \(\Rightarrow\) 4D gap" statement is **[T] inside a model, [NG]/[C] for the theory.**

**(K2) "It assumes the gap."** Fires, and this is the sharpest. By Lemma Z the hypothesis is *a gapped, discrete bulk spectrum.* By §3–4 this exists **iff a regular IR cap at a finite scale is present** — and the cap scale **is** the confinement scale (hard wall \(z_0\), BTZ \(r_+\), Witten cigar \(M_{KK}\)). The actual no-matter bulk (pure gauge \(\leftrightarrow\) pure gravity) is **pure AdS = gapless** (§5). So the hypothesis "(N) gapped normal-mode spectrum" is **logically equivalent to "an IR cap exists at a finite scale,"** which is the gap. The implication is true but **circular**: it concludes the gap from a premise that already contains it. A *single* normal mode is not enough either — without a positive bulk spectral bottom it gives an atom over a continuum, no gap (§1) — so the honest hypothesis must itself be the gap.

**(K3) "It reduces to H1."** Fires. §3 eq. (4)–(5): the 4D gap \(\Delta_4=j_{2,1}/z_0\) is a monotone function of the IR scale, \(\inf_{z_0}\Delta_4=0\), so *bounding the 4D gap below \(\Longleftrightarrow\) bounding the bulk IR scale above.* §4 eq. (7): the BTZ gap \(\propto r_+\) similarly. In both, the gap is the **compacitification/wall scale**, which vanishes under decompactification — precisely the §3.6 trap of LIT-ERA-08 ("the gap is not the same as a compactification gap"). Deriving that scale from pure YM (not fitting it) **is H1** in the bulk-geometry coordinate. Even in the one place a *smooth* cap is a genuine solution (Witten D4 cigar / KS throat), the cap scale is the compactification radius / deformation parameter — an input — **and** the dual is not pure YM (5D-on-a-circle, or \(N=1\) SUSY). Either way it reduces to H1 or to the wrong theory.

---

## §7 · The exact missing step

Decompose "5D normal mode \(\Rightarrow\) 4D mass gap" into the chain and label each link:

| Step | Statement | Label |
|---|---|---|
| **Z-i** | Bulk \(M\) is the holographic dual of 4D pure SU(N) YM | **[NG]/[C]** — no dual known (K1) |
| **Z-ii** | \(M\) has a regular IR cap at a finite **dynamically generated** scale, so \(\operatorname{spec}(L)\) is discrete with \(m_1^2>0\) | **[O] = H1** (K2) — the cap is the confinement scale |
| **Z-iii** | GKPW lifts the normal modes to boundary KL atoms; gap preserved at leading order | **[T] within model / [C]** as dictionary (Lemma Z, §2) |
| **Z-iv** | \(m_4=m_1\) survives \(\alpha'\), \(g_s\), \(1/N\) to \(N=3\), finite \(\lambda\); no SUSY protection | **[O]** uncontrolled (H5-type) |
| **Z-v** | \(\Delta_4>0\) uniform under decompactification | **[O]/[NG]** — fails: \(\Delta_4\to0\) (\(z_0\to\infty\), \(r_+\to0\)) (K3) |

**The exact missing step is Z-ii (with Z-i).** "5D normal mode \(\Rightarrow\) 4D mass gap" is **the easy link (Z-iii)**: given a bulk with a gapped SL operator and a valid dictionary, the boundary is gapped — that is Sturm–Liouville. The hard links are **Z-i** (a duality that does not exist) and **Z-ii** (a dynamically generated finite IR cap), and **Z-ii *is* the gap**. So the slogan is true where it is uninformative and uninformative where it is true: its premise smuggles in the IR cap, and the cap is the object H1 says must be derived. The "machine" the campaign keeps naming — an object that transports a strong-coupling gap across the Savvidy crossing and takes a continuum limit (charter §5) — would have to *be* Z-ii as a theorem; holography poses the question in geometry but does not build it.

---

## §8 · Honest bottom line

- **[T] (real, within a model):** Lemma Z — a gapped, discrete bulk SL spectrum lifts, via GKPW at leading order, to a 4D Källén–Lehmann gap equal to the bulk lowest eigenvalue. The hard-wall Bessel spectrum (3) and BTZ QNM scaling (6) are exact **[V]**.
- **[NG] for the Clay route:** the lifting chain cannot be turned into a statement about 4D pure YM. It is **[T] only about a model's SL operator**; the physical implication needs a duality of pure YM (**[NG]**, K1), the discrete positive spectrum needs an IR cap at a finite scale that **is** the confinement scale (**[O]=H1**, K2), and the gap equals that scale and vanishes under decompactification (**[O]/[NG]**, K3).
- **Acid test:** passes in the naive sense (pure AdS / free massless correctly gapless, §5), so the tool is not trivially vacuous; it is vacuous in the **"assumes-the-gap"** sense — the only scale it ever reports is the one fed in as the IR boundary condition.
- **One-line net:** *"5D normal mode \(\Rightarrow\) 4D mass gap" is true and useless* — the implication is Sturm–Liouville (Z-iii), while the antecedent "a gapped bulk normal-mode spectrum" already contains (as the IR cap) the confinement scale whose derivation is the whole problem. **Reduces to H1. Clay 0%. No QED.**

---

*Witness:* `research/v3/w7/z1_holo.py` (+ `.out`) — SL residual, Bessel spectrum, \(m_1(z_0)\to0\), pure-AdS continuum, BTZ QNM \(T\)-independence. numpy/scipy only.
*Internal:* `research/unconditional/G1-HOLO-NEC.md` (NEC \(\Leftrightarrow\) holographic \(a\)-theorem *admits* an IR AdS throat — the complementary no-go), `ROOF-HOLOGRAPHY.md`, `ROUTE-NEW-HOLO.md`, `CONSTRUCT-HOLO-FLUXTUBE.md` (all: gap = input IR scale = H1), `docs/paper-v2/ch11b-related.tex` §Holography, `ch10-betarg.tex` (\((G1)\)/H1; holographic NEC row).
*External:* J. Maldacena, hep-th/9711200; Gubser–Klebanov–Polyakov, hep-th/9802109; Witten, hep-th/9802150 (dictionary, **conjecture**); Breitenlohner–Freedman, *PLB* **115** (1982) 197 (BF stability, [T]); Polchinski–Strassler, hep-th/0109174 (hard wall, model); Karch–Katz–Son–Stephanov, hep-ph/0602229 (soft wall, model); Witten, hep-th/9803131; Sakai–Sugimoto, hep-th/0412141 (top-down, not pure YM); Brower–Mathur–Tan, hep-th/0003115 (D4 glueball ratios); Cardoso–Lemos, *PRD* **63** (2001) 124015, gr-qc/0101052 (exact BTZ scalar QNM); Son–Starinets, *JHEP* **0209** (2002) 042, hep-th/0205051 (AdS/CFT QNM prescription). Era-8 bibliographic data as verified in `research/lit-v3/LIT-ERA-08.md`; BTZ QNM and Son–Starinets checked against arXiv abstract pages this round.
