# LIT-ERA-05 · 1982–1990 — Constructive renormalization group and the gluon condensate

> Literature base for paper v3. **Era 05 of 14.** Scope: Balaban's multiscale RG for lattice gauge
> theory (1982–1989); Gawędzki–Kupiainen and Feldman–Trubowitz-style RG; the status of a 4D
> construction; the Shifman–Vainshtein–Zakharov gluon condensate and QCD sum rules (1979, background);
> the OPE in QCD; 't Hooft/Mueller renormalons and the ambiguity of the condensate; the first lattice
> determinations of the gluon condensate.
>
> **Label key.** **[T]** rigorous theorem/proof · **[C]** conjecture, conditional result or proposal ·
> **[NG]** no-go / obstruction / impossibility · **[O]** physics derivation, numerical or empirical
> observation · **[V]** review, lecture notes, textbook.
>
> **Verification.** `VERIFIED` = bibliographic data checked against publisher / ADS / INSPIRE metadata
> in this session. `VERIFIED (secondary)` = checked only through a secondary citation listing.
> `UNVERIFIED` = standard citation from memory, not checked; no datum is invented, and a page range
> that could not be checked is flagged in §9. Load-bearing entries (Balaban, Gawędzki–Kupiainen, SVZ,
> renormalons, lattice condensate) are all in the verified set.

---

## 0. Executive summary

1. Between 1984 and 1989 Bałaban built the only complete multiscale renormalization-group analysis of
   four-dimensional lattice gauge theory. The verified chain is: regular gauge fields and gauge-fixing
   surfaces (CMP **99**) → propagators and one-step RG transformations (CMP **95**, **96**) →
   small-field effective actions, β-function and coupling renormalization (CMP **109**) → UV cluster
   expansions (CMP **116**) → convergent expansions including large fields (CMP **119**) → the
   large-field R operation (CMP **122** I–II). The endpoint is **ultraviolet stability** at weak
   coupling in a fixed finite volume.
2. UV stability is a statement of *uniform boundedness of effective actions as the lattice cutoff is
   removed*. It yields a continuum limit in finite volume by compactness. It does **not** yield
   uniqueness of the limit, Osterwalder–Schrader positivity of the limit, nontriviality, an
   infinite-volume limit, or a mass gap.
3. The interface with our campaign's **G1/L3** node is exact: Balaban's uniformity is in the cutoff
   `a` at fixed volume `Λ`. The mass gap requires uniformity in `Λ` **and** a *lower* bound on the
   gauge-invariant effective action at the infrared scale. Balaban's estimates are upper bounds
   (stability), not lower bounds (gap).
4. At the IR scale the running coupling is `O(1)`; the small-field expansion on which the whole scheme
   rests is invalid, and the large-field R-operation bounds are far too lossy to produce a positive
   infrared quadratic form.
5. The second half of the era is the condensate. SVZ (1979) defined the gluon condensate
   `⟨(α_s/π) G^a_{μν}G^{aμν}⟩ ≈ 0.012 GeV⁴` as the leading dimension-4 parameter of the QCD operator
   product expansion, and turned QCD into a sum-rule industry.
6. The condensate is **not a sharp number**. 't Hooft's renormalon argument (Erice 1977), David (1984)
   and Mueller (1985) show that the OPE coefficient of a dimension-4 operator carries an infrared
   renormalon ambiguity of order `Λ_QCD⁴` — the same order as the condensate itself.
7. That ambiguity is the shadow of the missing mass gap. `⟨G²⟩` is (up to subtraction) the integral of
   the Källén–Lehmann spectral measure `ν` of `tr F²`. With no gap, `∫dν` is infrared-sensitive and
   scheme-dependent; with a gap `m₁`, the leading IR contribution is the atom `f₁² = ν({m₁²})`.
8. The first lattice determinations (Di Giacomo–Rossi 1981; Banks–Horsley–Rubinstein–Wolff 1981;
   Campostrini–Curci–Di Giacomo–Paffuti 1986) confirmed the order of magnitude but had to *choose* a
   subtraction scheme. They do not resolve the ambiguity; they renormalize it away by definition.
9. Consequence for v3: the condensate cannot serve as a rigorous order parameter. `f₁²` — the
   residue at the lightest glueball pole — is the scheme-independent object that a mass gap makes
   well defined. The renormalon ambiguity is a **diagnostic of the uncontrolled IR**, not a technical
   nuisance.
10. Honest verdict of the era: **UV renormalizability of 4D Yang–Mills is established at the level of
    constructive-QFT technology; the IR and the mass gap are completely open.** This is precisely the
    G1/L3 split of the campaign. The era's constructive-RG paradigm is a UV paradigm; it was exhausted
    at the IR end, which is direct historical evidence for the charter's "new species" requirement.

---

## 1. Balaban's multiscale program (1982–1989)

### 1.1 The verified chain

| # | Paper | Venue | What it supplies |
|---|-------|-------|------------------|
| 1 | Spaces of regular gauge field configurations on a lattice and gauge fixing conditions | CMP **99** (1985) 75–102 | gauge-invariant notion of a *regular* configuration; a gauge-fixing surface on which every regular configuration has a small representative |
| 2 | Propagators and renormalization group transformations for lattice gauge field theories, I | CMP **95** (1984) 17–40 | one-step RG transformation; Gaussian sector and propagator bounds |
| 3 | Propagators and renormalization group transformations, II | CMP **96** (1984) 223–250 | covariance/decay bounds for the transformed propagators |
| 4 | Renormalization group approach to lattice gauge field theories, I: Generation of effective actions in a small field approximation and a coupling constant renormalization in four dimensions | CMP **109** (1987) 249–301 | localized effective actions, β-function, recursive coupling renormalization |
| 5 | Renormalization group approach …, II: Cluster expansions | CMP **116** (1988) 1–22 | polymer/cluster expansion of the one-step map; convergence of the UV expansion |
| 6 | Convergent renormalization expansions for lattice gauge theories | CMP **119** (1988) 243–285 | inductive description of complete effective densities including large-field domains |
| 7 | Large field renormalization, I: The basic step of the R operation | CMP **122** (1989) 175–202 | gauge-covariant, localized BPHZ-like subtraction in large-field regions |
| 8 | Large field renormalization, II: Localization, exponentiation and bounds for the R operation | CMP **122** (1989) 355–392 | localization, exponentiation, bounds; completes the proof of UV stability (Theorem 1) |

Logical reading order is 1 → 2/3 → 4 → 5 → 6 → 7 → 8. The theorem statements are in the introductions
of 4 and 8. Citations in later literature sometimes number the series differently (a "Part III/IV" is
frequently quoted); the titles and volumes above are the ones verified here, and I found no verified
paper with "ultraviolet limit" in the title. Treat the eight items, not a part-numbering, as the
canonical chain.

### 1.2 What is actually proved

- **(B1) Gauge-invariant regularity.** There is a gauge-invariant norm on lattice configurations
  (plaquettes and covariant differences at all scales) such that the regular set is preserved by the
  RG map, and such that on a suitable gauge-fixing surface every regular configuration has a small
  representative. This is the technical device that makes a multiscale decomposition compatible with
  gauge invariance — the point where Gribov copies (B17, B18) would otherwise destroy the
  construction.
- **(B2) One-step RG map.** The block transformation on gauge fields is defined, its kernel has the
  required decay, and the Gaussian part is controlled uniformly in the number of steps.
- **(B3) Small-field effective action.** In the region where the field is small, the effective action
  after `k` steps stays in a prescribed analyticity class; the β-function is computed; the bare
  coupling is renormalized recursively along the asymptotically free trajectory.
- **(B4) UV cluster expansion.** The one-step map admits a convergent polymer expansion with
  activities small uniformly in the scale index. This is a *UV* cluster expansion — convergence of the
  RG map — **not** a proof of clustering of the physical measure. The title "Cluster expansions"
  (B5) is routinely misread; it does not give exponential decay of physical correlations.
- **(B5) Large-field bookkeeping.** The inductive description of effective densities is preserved when
  large-field domains are included (B6).
- **(B6) The R operation.** A gauge-covariant, localized subtraction removes the large-field
  divergences; B7–B8 prove its localization, exponentiation and bounds.
- **(B7) UV stability.** The effective actions satisfy `k`-uniform bounds of the same type as
  superrenormalizable models. Hence the family of cutoff measures is compact and a continuum limit
  exists in finite volume at weak coupling.

### 1.3 The large-field problem

The obstruction the program was built to defeat is specific: in four dimensions the entropy of
large-field configurations is marginal, so the standard perturbative RG cannot be applied globally, and
one cannot simply expand in the field. Balaban's answer is not to avoid large fields but to *renormalize
them*: the R operation is a gauge-covariant, spacetime-localized version of the BPHZ subtraction,
applied only where the field is large. The two "Large field renormalization" papers are the technical
heart of the whole program and the reason it took seven years.

### 1.4 What is *not* proved

- **(N1) No volume uniformity.** All constants depend on the volume `Λ` (typically through the number
  of blocks, exponentially or worse). There is no control of `Λ → ∞`.
- **(N2) No Osterwalder–Schrader reconstruction.** Reflection positivity holds for each lattice measure
  `μ_a`, but the limit is obtained by compactness; the series does not prove that reflection positivity
  passes to the limit point, nor that the limit is a measure at all rather than a formal effective
  action. No Hilbert space, no transfer matrix, no Hamiltonian is constructed.
- **(N3) No nontriviality.** There is no lower bound on a gauge-invariant observable that survives
  `a → 0`. The limit could in principle be trivial or a non-OS object; asymptotic freedom makes
  triviality unlikely but is not a proof.
- **(N4) No infrared analysis.** Nothing is proved about the spectral measure near `s = 0`; no
  clustering; no mass gap. The multiscale scheme stops at the physical scale and does not analyse the
  long-distance effective action.
- **(N5) No lower bounds at the IR scale.** Balaban's estimates control `|S_eff|` from above in the
  small-field domain. A gap needs a *lower* bound on the gauge-invariant curvature of `S_eff` at the
  IR scale. None is available.
- **(N6) No matter fields, no θ-vacuum, no topological sectors** in the published chain; the
  construction is for pure gauge fields on a topologically trivial torus.
- **(N7) Proof status caveat.** The proof is spread over roughly three hundred pages of dense
  estimates and has never been independently formalized or machine-checked. "Proved" here means
  "published proof, standard citation in the constructive-QFT literature", not "verified line by
  line". Later independent work (B34) reproved UV stability with an infrared cutoff by a different and
  much shorter route, which is the strongest external confirmation.

---

## 2. The rest of the constructive-RG school

### 2.1 Gawędzki–Kupiainen

The methodological twin of Balaban's problem. In `φ⁴₄` the coupling is asymptotically free in the
**infrared**, so the same block-spin multiscale machinery applies, and GK proved infrared asymptotic
freedom of the massless lattice `φ⁴₄` model (B10), with the block-spin framework set out in B9 and the
method summarized in the Les Houches lectures (B11). Two lessons transfer to gauge theory:

- the multiscale expansion is *the* tool for a marginal coupling, and it works when the sign of the
  β-function is right;
- the hard part is always the same: overlapping large-field regions and the resulting loss of
  smallness. GK's treatment is cleaner than Balaban's precisely because `φ⁴₄` has no gauge freedom
  and no marginal large-field entropy problem of the same severity.

GK did **not** treat gauge theory. Their result is a proof of concept for the RG method, not a step
toward the Yang–Mills mass gap.

### 2.2 Feldman–Trubowitz and the fermionic branch

The Feldman–Trubowitz program (B12) and its four-author extension (B13) introduced the *sector*
decomposition: instead of controlling large fields, one decomposes momentum space into sectors and
works on each sector separately, with the combinatorial cost paid in the sector sum. This is a
different answer to the same disease — avoid large fields rather than renormalize them. It has never
been applied to gauge theory. Its relevance to v3 is as a template for a *gauge-covariant sector
decomposition*: choose a gauge sector by sector (a block, not a global slice), which is exactly what
Balaban's regular-gauge construction does globally.

### 2.3 Brydges–Fröhlich–Seiler, Osterwalder–Seiler, Seiler

- **Osterwalder–Seiler (B14)** is the paper that put lattice gauge theory on the OS-axiom footing and
  supplied the reconstruction theorem route from a Euclidean measure to a Hilbert space. It is the
  logical bridge that Balaban's limit is *not* shown to cross.
- **Brydges–Fröhlich–Seiler (B15)** constructed quantized gauge fields in two and three dimensions
  (general results; convergence of the lattice approximation; the two-dimensional Higgs model),
  verifying OS axioms there. The 4D pure gauge case is out of reach of their methods because the
  coupling is marginal rather than superrenormalizable.
- **Seiler's lecture notes (B16)** are the best 1982 snapshot of the constructive program and state the
  open problems — including the mass gap — without euphemism.

### 2.4 Magnen–Rivasseau–Sénéor (1993) — the clean confirmation, outside the window

B34 proves UV stability of `YM₄` **with an infrared cutoff** by a simpler method. Its title is the
honest statement of the era's frontier: the ultraviolet is under control, the infrared is explicitly
cut off. Any claim that "Balaban constructed 4D Yang–Mills" should be read against this paper.

---

## 3. Status of a 4D construction, circa 1990

| `d` | Model | Status in 1990 |
|-----|-------|----------------|
| 2 | `YM₂` | exactly solvable (Migdal; Gross–Witten); area law; no propagating gluons; the "gap" is kinematic |
| 3 | `YM₃` | mass gap at strong coupling by cluster expansion; all-coupling gap open; strong lattice evidence |
| 4 | `YM₄` | UV stability at weak coupling in finite volume (Balaban); **no** construction, **no** OS axioms, **no** infinite-volume limit, **no** gap |
| 4 | `φ⁴₄` | expected trivial (Fröhlich 1982); rigorous marginal triviality only in 2019–21 — the contrast case that shows why asymptotic freedom is essential |

No rigorous construction of a nontrivial four-dimensional quantum field theory with a mass gap and a
non-Abelian gauge symmetry existed in 1990, and none exists now. The era's achievement is narrower and
real: it showed that the *ultraviolet* problem of 4D gauge theory is solvable with the technology of
the time, and it isolated the infrared as the irreducible obstacle.

---

## 4. The gluon condensate and QCD sum rules

### 4.1 SVZ 1979

The three SVZ papers (B20–B22) founded QCD sum rules. The logic:

1. Write a dispersion relation for a current–current correlator `Π(q²)`.
2. Use the operator product expansion at short distance: `Π(q²) = Σ_n C_n(q²) ⟨O_n⟩`.
3. The leading nonperturbative term is dimension 4: the gluon condensate
   `⟨(α_s/π) G^a_{μν}G^{aμν}⟩`.
4. Borel-transform to suppress both the continuum and higher condensates, and fit the resonance
   parameters.

The phenomenological value `⟨(α_s/π)G²⟩ ≈ 0.012 GeV⁴` came out of charmonium sum rules (B21), with a
successful prediction of the `η_c` mass. The assumptions are strong and were understood as such from
the start: quark–hadron duality, a Borel window, vacuum saturation of four-quark condensates, and a
truncated OPE. The condensate is scheme- and scale-dependent; it is a parameter of a subtraction
scheme, not a directly measurable number.

### 4.2 The OPE in QCD

Wilson's operator product expansion (B19) writes a product of local operators at short distance as a
sum of local operators with c-number coefficient functions. In QCD the coefficients carry the
short-distance physics (computable in perturbation theory, with anomalous dimensions from operator
mixing) and the matrix elements carry the long-distance physics. Two features matter for v3:

- **Mixing with the identity.** `G²` mixes with the unit operator; on the lattice this is the `1/a⁴`
  power divergence. Any "condensate" is defined only after a subtraction.
- **Asymptotic, not convergent.** The OPE in an interacting theory with a mass gap is at best an
  asymptotic expansion; its accuracy is limited by the next uncalculated term, and in QCD the
  expansion parameter is `Λ²/Q²`. Novikov–Shifman–Vainshtein–Zakharov (B23) asked whether the OPE can
  fail; the modern answer is that it is an asymptotic series whose leading ambiguity is the
  renormalon.

### 4.3 Renormalons and the ambiguity of the condensate

- **'t Hooft (B24)** observed that the QCD perturbation series grows like `n! (β₀/2π)^n`; the Borel
  transform has singularities. The infrared renormalon at `t = 2π/β₀` (conventions vary; often
  `u = 2`) produces an ambiguity of order `(Λ²/Q²)² = Λ⁴/Q⁴` in the correlator.
- **David (B25)** showed that composite operators acquire renormalon ambiguities: the coefficient of a
  dimension-4 operator is ambiguous at `O(Λ⁴)`, and the ambiguity can be moved between the coefficient
  and the condensate. The gluon condensate is therefore **not a sharply defined order parameter** in
  perturbation theory.
- **Mueller (B26)** determined the structure of infrared renormalons in infrared-safe processes: the
  leading nonperturbative corrections are `1/Q⁴`, with the singularity structure in the Borel variable
  made explicit.
- **David (B27)** clarified the status of the OPE in the presence of these ambiguities.

The physical statement is simple and important: `⟨G²⟩` is defined only up to an `O(Λ_QCD⁴)` ambiguity,
which is the same order as the condensate itself. A nonperturbative definition — i.e. a mass gap and a
controlled IR spectral measure — is required before the condensate means anything.

### 4.4 First lattice determinations

- **Di Giacomo–Rossi (B29)** gave the first lattice extraction, from the plaquette, using the trace
  anomaly. The title is literal: "Extracting `⟨(α/π)Σ G G⟩` from gauge theories on a lattice".
- **Banks–Horsley–Rubinstein–Wolff (B30)** did the first Monte Carlo study of the condensate on the
  lattice.
- **Campostrini–Curci–Di Giacomo–Paffuti (B31)** studied the renormalization and continuum limit of the
  dimension-4 composite operator in pure SU(2): multiplicative finite renormalization consistent with
  one-loop perturbation theory, and asymptotic scaling for `2.45 ≲ β ≲ 2.85`.

The systematics are the lattice avatar of §4.3: `1/a⁴` power divergence, operator mixing, scheme
dependence. The lattice *defines* a subtracted condensate; it does not derive the SVZ value, and it does
not remove the renormalon ambiguity — it chooses a scheme in which the ambiguity is fixed by the
subtraction.

### 4.5 Bridge to the residue program

Let `O = tr F²` and let `ν` be its Källén–Lehmann spectral measure. Then (up to a subtraction)

```
⟨0|O|0⟩  =  ∫ dν(s).
```

- If `ν` has a gap and an atom at `m₁²`, then `∫dν = f₁² + ∫_{s>m₁²} dν`, and the leading infrared
  contribution is the residue `f₁²`. The condensate is then controlled by the atom plus the continuum
  above threshold.
- If `ν` has no gap, `∫dν` is dominated by the unknown `s → 0` region; the integral is scheme- and
  cutoff-dependent. **The renormalon ambiguity of the condensate is precisely the statement that this
  integral is uncontrolled in the infrared.**

This is the era's deepest single lesson for v3: the historical "ambiguity of the gluon condensate" is
not a technicality of QCD sum rules; it is the shadow of the missing mass gap. Conversely, a proof of
`f₁² > 0` with a gap would give a nonperturbative definition of the leading infrared piece of the
condensate and remove the leading renormalon ambiguity by physics rather than by subtraction. This
bridge is a heuristic at the physics level (label [O], not [T]); it is developed as an attack surface
in §8.4.

---

## 5. THE INTERFACE: Balaban's UV control ↦ the IR mass gap (G1/L3)

> This is the campaign's G1/L3 node. It is stated here as precisely as the sources allow.

### 5.1 Notation

- `μ_a` — lattice gauge measure at cutoff `a`, fixed volume `Λ`.
- `μ` — putative continuum limit.
- `O = tr F²` (or the plaquette density); `ν` — its Källén–Lehmann spectral measure.
- **Gap:** `supp ν ∩ (0, m₁²) = ∅` and `f₁² := ν({m₁²}) = |⟨0|O|G⟩|² > 0`.
- **G1** — 4D continuum existence. **G2** — `f₁² > 0` / gap. **L3** — the fixed-volume RG step still
  missing.

### 5.2 What Balaban supplies

- **(U_a) Cutoff uniformity.** Effective actions satisfy `k`-uniform bounds as `a → 0`, at fixed `Λ`.
  This is the strongest existing input toward G1.
- **(U_β) Correct running.** The β-function and the renormalized coupling are controlled to all orders
  in the UV, on the asymptotically free trajectory.
- **(U_c) Compactness.** A limit point `μ` exists in finite volume, in the sense of convergence of
  effective actions.

### 5.3 What is missing, stated as three uniforms

| Symbol | Uniformity | Balaban | Needed for |
|--------|-----------|---------|------------|
| `U_a` | cutoff `a → 0` at fixed `Λ` | **yes** | G1 (UV half) |
| `U_Λ` | volume `Λ → ∞` at fixed `a` | **no** | G1 (IR half) / L3 |
| `U_IR` | IR scale, lower bound on curvature | **no** | G2 |

`L3` is exactly `U_Λ` at fixed physical volume: a renormalization step whose constants do not grow with
the number of blocks. Balaban's constants do grow with `Λ`; his bounds are not `(a, Λ)`-uniform. This
is the precise location of the campaign's missing step.

### 5.4 Why the small-field expansion cannot reach the gap

1. At the IR scale the running coupling is `g(1/m₁) = O(1)`; the small-field condition `|F| ≪ 1` fails.
2. The gap is generated by dimensional transmutation, `m₁ ~ Λ_QCD ~ μ exp(−∫ dg/β(g))`: a
   nonperturbative effect. Balaban's scheme tracks the running coupling but not the transmutation
   scale; there is no statement about the location or existence of `m₁`.
3. Balaban's large-field bounds are *upper* bounds produced by a subtraction. A gap needs a *lower*
   bound on the gauge-invariant quadratic form of `S_eff` at the IR scale. The R operation is built to
   make the action finite, not to make it convex.
4. The R operation is gauge-covariant but is not shown to preserve the Osterwalder–Schrader cone.
   Even if the limit existed, the positivity needed for the reconstruction is not transported.
5. The title "Cluster expansions" (CMP 116) refers to the UV polymer expansion of the RG map. It does
   **not** give clustering of the physical measure; there is no exponential decay of gauge-invariant
   correlators in the series.

### 5.5 The one-sentence interface

**Balaban proves that the ultraviolet end of the two-scale structure can be taken uniformly in the
cutoff at fixed volume; the mass gap is a lower-bound statement about the infrared end, uniform in the
volume, and nothing in the series addresses it.**

Equivalently, in the campaign's language: Balaban controls the UV fixed point and the approach to it
(the `a`-half of `E₂`), while `f₁² > 0` lives in the relevant direction created by dimensional
transmutation, which the UV scheme does not track nonperturbatively. The G1/L3 node is the demand for a
single renormalization that is uniform in **both** `a` and `Λ`; Balaban gives only the `a`-uniform half.

### 5.6 Acid-test remark

Balaban's machinery is constructed so that it degenerates correctly for the free field: with
`β₀ = 0` the multiscale expansion is trivial and the bounds are those of a Gaussian measure. It thus
passes the structural analogue of the charter's free-field acid test — but precisely the part that is
invisible in the free theory (the generation of `m₁` and `f₁²`) is the part the machinery does not
compute. A v3 tool that reproduces Balaban's UV bounds but says nothing new about `ν` near `s = 0` has
not addressed the mass gap.

---

## 6. Annotated bibliography

**A. Balaban**

**B1** [T] T. Bałaban, "Propagators and renormalization group transformations for lattice gauge field
theories. I," *Commun. Math. Phys.* **95** (1984) 17–40. — One-step RG transformation; Gaussian sector;
propagator decay. `VERIFIED (secondary)`

**B2** [T] T. Bałaban, "Propagators and renormalization group transformations for lattice gauge field
theories. II," *Commun. Math. Phys.* **96** (1984) 223–250. — Covariance and decay bounds for the
transformed propagators; the analytic backbone of the later estimates. `VERIFIED (secondary)`

**B3** [T] T. Bałaban, "Spaces of regular gauge field configurations on a lattice and gauge fixing
conditions," *Commun. Math. Phys.* **99** (1985) 75–102. DOI 10.1007/BF01466594. — Gauge-invariant
regularity; gauge-fixing surface; small representative. The device that makes multiscale decomposition
gauge-compatible. `VERIFIED`

**B4** [T] T. Bałaban, "Renormalization group approach to lattice gauge field theories. I. Generation of
effective actions in a small field approximation and a coupling constant renormalization in four
dimensions," *Commun. Math. Phys.* **109** (1987) 249–301. DOI 10.1007/BF01215223. — Localized
effective actions; β-function; recursive coupling renormalization. The theorem statement of the
small-field half. `VERIFIED`

**B5** [T] T. Bałaban, "Renormalization group approach to lattice gauge field theories. II. Cluster
expansions," *Commun. Math. Phys.* **116** (1988) 1–22. — Convergent UV polymer expansion of the
one-step map. Not a proof of physical clustering. `VERIFIED`

**B6** [T] T. Bałaban, "Convergent renormalization expansions for lattice gauge theories," *Commun.
Math. Phys.* **119** (1988) 243–285. DOI 10.1007/BF01217741. — Inductive description of complete
effective densities including large-field domains; preservation of the density form. `VERIFIED`

**B7** [T] T. Bałaban, "Large field renormalization. I. The basic step of the R operation," *Commun.
Math. Phys.* **122** (1989) 175–202. DOI 10.1007/BF01257412. — Gauge-covariant localized BPHZ
subtraction in large-field regions; "removes the main obstacle" to UV stability. `VERIFIED`

**B8** [T] T. Bałaban, "Large field renormalization. II. Localization, exponentiation and bounds for
the R operation," *Commun. Math. Phys.* **122** (1989) 355–392. DOI 10.1007/BF01238433. — Completes
the proof of UV stability of 4D pure gauge theory (Theorem 1). The endpoint of the program. `VERIFIED`

**B. Gawędzki–Kupiainen, Feldman–Trubowitz**

**B9** [T] K. Gawędzki, A. Kupiainen, "A rigorous block spin approach to massless lattice theories,"
*Commun. Math. Phys.* **77** (1980) 31–64. — Block-spin framework for marginal lattice models; the
methodological precursor of B10. `UNVERIFIED`

**B10** [T] K. Gawędzki, A. Kupiainen, "Massless lattice φ⁴₄ theory: rigorous control of a
renormalizable asymptotically free model," *Commun. Math. Phys.* **99** (1985) 197–252. DOI
10.1007/BF01212281. — Infrared asymptotic freedom of massless lattice φ⁴₄; logarithmic corrections;
the proof of concept for multiscale control of a marginal coupling. `VERIFIED`

**B11** [V] K. Gawędzki, A. Kupiainen, "Asymptotic freedom beyond perturbation theory," in *Critical
Phenomena, Random Systems, Gauge Theories* (Les Houches 1984), eds. K. Osterwalder, R. Stora,
North-Holland, 1986, pp. 185–292. — The concise statement of the method and its limitations. `UNVERIFIED`

**B12** [T] J. Feldman, E. Trubowitz, "Perturbation theory for many fermion systems," *Helv. Phys.
Acta* **63** (1990) 156–260. — Sector decomposition; avoids large fields instead of renormalizing
them. Never applied to gauge theory. `UNVERIFIED`

**B13** [T] J. Feldman, J. Magnen, V. Rivasseau, E. Trubowitz, "An infinite volume expansion for many
fermion Green's functions," *Helv. Phys. Acta* **65** (1992) 679–721. — The four-author sector
expansion; the fermionic branch's mature form. `UNVERIFIED`

**C. Constructive / OS foundations**

**B14** [T] K. Osterwalder, E. Seiler, "Gauge field theories on a lattice," *Ann. Phys.* **110** (1978)
440–471. — Lattice gauge theory and the OS axioms; the reconstruction route from measure to Hilbert
space. The bridge Balaban's limit does not cross. `UNVERIFIED`

**B15** [T] D. Brydges, J. Fröhlich, E. Seiler, "On the construction of quantized gauge fields. I.
General results," *Ann. Phys.* **121** (1979) 227–284; "II. Convergence of the lattice approximation,"
*Commun. Math. Phys.* **71** (1980) 91–107. — Constructive gauge theory in `d ≤ 3`; OS axioms verified
there; 4D out of reach. `UNVERIFIED`

**B16** [V] E. Seiler, *Gauge Theories as a Problem of Constructive Quantum Field Theory and
Statistical Mechanics*, Lecture Notes in Physics **159**, Springer, Berlin, 1982. — The 1982 snapshot;
states the mass gap as open. `UNVERIFIED`

**B17** [NG] V. N. Gribov, "Quantization of non-Abelian gauge theories," *Nucl. Phys.* **B139** (1978)
1–19. — The gauge-fixing obstruction: no globally smooth gauge slice. Structural reason Balaban needs
B3. `UNVERIFIED`

**B18** [T/NG] I. M. Singer, "Some remarks on the Gribov ambiguity," *Commun. Math. Phys.* **60** (1978)
7–12. — Rigorous statement of the Gribov ambiguity for compact gauge groups. `UNVERIFIED`

**D. Condensate and OPE**

**B19** [T] K. G. Wilson, "Non-Lagrangian models of current algebra," *Phys. Rev.* **179** (1969)
1499–1512. — The operator product expansion. The short-distance half of the SVZ method. `UNVERIFIED`

**B20** [O] M. A. Shifman, A. I. Vainshtein, V. I. Zakharov, "QCD and resonance physics. Theoretical
foundations," *Nucl. Phys.* **B147** (1979) 385–447. DOI 10.1016/0550-3213(79)90022-1. — OPE plus
vacuum condensates; Borel sum rules. `VERIFIED`

**B21** [O] M. A. Shifman, A. I. Vainshtein, V. I. Zakharov, "QCD and resonance physics. Applications,"
*Nucl. Phys.* **B147** (1979) 448–518. — Charmonium sum rules; the gluon condensate
`⟨(α_s/π)G²⟩ ≈ 0.012 GeV⁴` emerges; `η_c` prediction. `VERIFIED`

**B22** [O] M. A. Shifman, A. I. Vainshtein, V. I. Zakharov, "QCD and resonance physics. The ρ–ω
mixing," *Nucl. Phys.* **B147** (1979) 519–534. — Third paper of the series. `VERIFIED`

**B23** [O/C] V. A. Novikov, M. A. Shifman, A. I. Vainshtein, V. I. Zakharov, "Wilson's operator
expansion: can it fail?" *Nucl. Phys.* **B249** (1985) 445–471. — The ITEP group's own examination of
the limits of the OPE; the question that renormalons answer. `UNVERIFIED`

**E. Renormalons**

**B24** [NG] G. 't Hooft, "Can we make sense out of 'quantum chromodynamics'?", in *The Whys of
Subnuclear Physics* (Erice 1977), ed. A. Zichichi, Plenum, New York, 1979, pp. 943–982. — The
renormalon: factorial growth, Borel singularities, the `1/N` and `η'` discussion. Perturbation theory
is not Borel summable. `VERIFIED`

**B25** [NG] F. David, "On the ambiguity of composite operators, IR renormalons and the status of the
operator product expansion," *Nucl. Phys.* **B234** (1984) 237–251. DOI 10.1016/0550-3213(84)90235-9.
— Composite operators acquire `O(Λ⁴)` renormalon ambiguities; the gluon condensate is not sharply
defined. `VERIFIED`

**B26** [O] A. H. Mueller, "On the structure of infrared renormalons in physical processes at high
energies," *Nucl. Phys.* **B250** (1985) 327–350. DOI 10.1016/0550-3213(85)90485-7. — Positions and
nature of IR renormalon singularities; leading `1/Q⁴` corrections. `VERIFIED`

**B27** [O] F. David, "The operator product expansion and renormalons: a comment," *Nucl. Phys.* **B263** (1986) 637–648. — Clarifies the applicability of the OPE in the presence of renormalons.
`VERIFIED`

**B28** [V] M. Beneke, "Renormalons," *Phys. Rep.* **317** (1999) 1–142. — The standard review;
systematizes the 1984–86 results. Out of window; the best entry point for v3. `UNVERIFIED`

**F. Lattice condensate**

**B29** [O] A. Di Giacomo, G. C. Rossi, "Extracting `⟨(α/π)Σ_{a,μν} G^a_{μν}G^a_{μν}⟩` from gauge
theories on a lattice," *Phys. Lett.* **B100** (1981) 481–484. DOI 10.1016/0370-2693(81)90609-2. — First
lattice determination, from the plaquette via the trace anomaly. `VERIFIED`

**B30** [O] T. Banks, R. Horsley, H. R. Rubinstein, U. Wolff, "The gluon condensate and the lattice,"
*Nucl. Phys.* **B190** (1981) 692–704. — First Monte Carlo study of the condensate. `UNVERIFIED`

**B31** [O] M. Campostrini, G. Curci, A. Di Giacomo, G. Paffuti, "Renormalization and continuum limit
of composite operators in lattice gauge theories," *Z. Phys.* **C32** (1986) 377–382. DOI
10.1007/BF01551834. — Multiplicative renormalization of the dimension-4 operator in SU(2); operator
mixing; asymptotic scaling for `2.45 ≲ β ≲ 2.85`. `VERIFIED`

**G. Bridge**

**B32** [O] M. Lüscher, "Some analytic results concerning the mass spectrum of Yang–Mills gauge
theories on a torus," *Nucl. Phys.* **B219** (1983) 233–261. — Semiclassical small-torus spectrum; the
volume dependence of the gap; the place where UV and IR meet quantitatively. `UNVERIFIED`

**B33** [T] M. Lüscher, "Volume dependence of the energy spectrum in massive quantum field theories.
I. Stable particle states," *Commun. Math. Phys.* **104** (1986) 177–206; "II. Scattering states,"
*Commun. Math. Phys.* **105** (1986) 153–188. — The finite-volume formula used to extract `m₁` from
lattice data; the interface between the measured gap and the infinite-volume gap. `UNVERIFIED`

**B34** [T] J. Magnen, V. Rivasseau, R. Sénéor, "Construction of YM₄ with an infrared cutoff," *Commun.
Math. Phys.* **155** (1993) 417–430. — Independent, shorter proof of UV stability, explicitly with an
IR cutoff; the honest statement of the frontier. `UNVERIFIED`

**B35** [V] V. Rivasseau, *From Perturbative to Constructive Renormalization*, Princeton University
Press, 1991. — Standard synthesis of the multiscale method and of Balaban's place in it. `UNVERIFIED`

---

## 7. Position in 14 eras

Working periodization used across the LIT-ERA files; neighboring files may subdivide differently.
Era 05 is this file.

| # | Years | Era | One-line content |
|---|-------|-----|------------------|
| 01 | 1954–1972 | Origins and quantization | Yang–Mills, Faddeev–Popov, Gribov; no confinement, no mass gap |
| 02 | 1973–1974 | Asymptotic freedom | Gross–Wilczek–Politzer, 't Hooft; the coupling runs, perturbation theory dies at long distance |
| 03 | 1974–1977 | Lattice and confinement criteria | Wilson's lattice, strong-coupling area law, center vortices, compact QED |
| 04 | 1977–1982 | Semiclassical vacuum | Instantons, Savvidy's unstable mode, Nielsen–Olesen, 1/N, loop equations |
| **05** | **1982–1990** | **Constructive RG and the gluon condensate** | **Balaban's UV stability, Gawędzki–Kupiainen, SVZ sum rules, renormalons, first lattice condensate** ← this file |
| 06 | 1990–1997 | Rigorous lattice statistical mechanics | OS reconstruction, `d ≤ 3` constructions, Seiler's program; 4D stays open |
| 07 | 1994–1999 | Seiberg–Witten and dual superconductivity | Monopole/vortex condensation as the confinement mechanism |
| 08 | 1997–2004 | Renormalons and OPE systematics | Borel resummation, IR ambiguities, Magnen–Rivasseau–Sénéor UV proof |
| 09 | 2004–2010 | Precision lattice | Glueball spectrum, string tension, scale setting; the gap is measured, not proved |
| 10 | 2010–2014 | Gradient flow and continuum control | Lüscher flow, step scaling, nonperturbative renormalization of composite operators |
| 11 | 2014–2018 | Conformal window, bootstrap, large-N | Phase-structure constraints; rigorous lower-dimensional mass gaps |
| 12 | 2018–2021 | Constructive revival | OS axioms for lattice gauge theory, area law, hairpin, `d = 2,3` mass gaps |
| 13 | 2021–2024 | Formalization and conditional reduction | Lean 4, conditional theorems, the gap as a spectral-residue problem |
| 14 | 2024–2026 | The v3 campaign | From conditional to unconditional: G1/L3, `f₁² > 0`, two-scale `E₂` |

---

## 8. Gaps and attack surfaces for v3

1. **The missing `(a, Λ)`-uniform RG step.** Balaban's constants grow with the volume. Target: redo the
   CMP 116 polymer expansion with `Λ`-tempered weights and a block-dependent small-field domain, and
   check whether the number of large-field blocks per unit volume can be bounded independently of the
   number of blocks. This *is* L3; success here would be the first genuine step toward G1.
2. **Positivity-preserving large-field subtraction.** The R operation is not shown to preserve the OS
   cone. Target: formulate R as a conditional expectation on reflection-positive functions; calibrate
   on the 2D Higgs model (B15) where OS is known. Without this, even a limit would not reconstruct a
   Hilbert space.
3. **IR lower bound at `O(1)` coupling.** Lemma SC gives a coarse Poincaré inequality at `β < 1/48`.
   Target: transport it to `β ≈ 2.5` along the renormalized trajectory. The obstruction is the Savvidy
   negative mode (era 04) — a *lower*-bound obstruction, not an upper-bound one, which is exactly why
   Balaban's stability estimates cannot see it.
4. **Renormalon-subtracted residue.** Define `f₁²` as `ν({m₁²})`, never as a perturbative coefficient.
   Target: show that the leading `O(Λ⁴)` IR renormalon ambiguity cancels between the OPE coefficient and
   the condensate, leaving `f₁²` as the renormalization-group-invariant residue. This is a new-tool
   candidate in the charter's sense: it must fail for the free field (`f₁² = 0`), and it does.
5. **Nontriviality.** No gauge-invariant lower bound survives `a → 0` in Balaban. Target: a uniform
   lower bound on the plaquette or the topological susceptibility compatible with
   `⟨(α_s/π)G²⟩ ≈ 0.012 GeV⁴`.
6. **Gauge-invariant Hilbert space.** Balaban's gauge fixing (B3) is not obviously compatible with OS
   reconstruction; Gribov copies (B17, B18) are the structural reason. Target: an orbit-space
   multiscale decomposition in holonomy variables on a block lattice, with no global gauge slice.
7. **The condensate as a diagnostic, not an input.** Use the SVZ value and the lattice glueball mass
   `m_{0++} ≈ 1.7 GeV` to constrain `f₁²`: if the atom dominates the IR end of `∫dν`, then `f₁²` and
   the condensate are tied by a sum rule. Make this quantitative and check it against the Borel window.
   A mismatch would falsify the "atom dominates" picture before any proof is attempted.
8. **Cross-era acid test.** Any v3 IR machine must reproduce the known `d = 2` and `d = 3` mass gaps
   (era 06) and must fail on the free field. A machine that passes Balaban's UV bounds but says nothing
   about `ν` near `s = 0` fails the campaign's own criterion and should be rejected early.

---

## 9. Verification ledger

| # | Entry | Checked | Source | Status |
|---|-------|---------|--------|--------|
| B1 | Balaban, Propagators I, CMP 95 (1984) 17–40 | volume/pages via secondary listing | Springer companion citation | VERIFIED (secondary) |
| B2 | Balaban, Propagators II, CMP 96 (1984) 223–250 | volume/pages via secondary listing | Springer companion citation | VERIFIED (secondary) |
| B3 | Balaban, Spaces of regular gauge fields, CMP 99 (1985) 75–102 | title, volume, pages, DOI | link.springer.com/article/10.1007/BF01466594 | VERIFIED |
| B4 | Balaban, RG approach I, CMP 109 (1987) 249–301 | title, volume, issue, pages, DOI | link.springer.com/article/10.1007/BF01215223 | VERIFIED |
| B5 | Balaban, RG approach II, CMP 116 (1988) 1–22 | title, volume, pages | Springer/zbMATH citation records | VERIFIED |
| B6 | Balaban, Convergent renormalization expansions, CMP 119 (1988) 243–285 | title, volume, pages, DOI | link.springer.com/article/10.1007/BF01217741 | VERIFIED |
| B7 | Balaban, Large field I, CMP 122 (1989) 175–202 | title, volume, issue, pages, DOI | link.springer.com/article/10.1007/BF01257412 | VERIFIED |
| B8 | Balaban, Large field II, CMP 122 (1989) 355–392 | title, volume, issue, pages, DOI | link.springer.com/article/10.1007/BF01238433 | VERIFIED |
| B9 | Gawędzki–Kupiainen, block spin, CMP 77 (1980) 31–64 | not checked | memory | UNVERIFIED |
| B10 | Gawędzki–Kupiainen, φ⁴₄, CMP 99 (1985) 197–252 | title, volume, issue, pages, DOI, MR number | link.springer.com/article/10.1007/bf01212281 | VERIFIED |
| B11 | Gawędzki–Kupiainen, Les Houches 1984 | not checked | memory | UNVERIFIED |
| B12 | Feldman–Trubowitz, Helv. Phys. Acta 63 (1990) 156–260 | not checked | memory | UNVERIFIED |
| B13 | FMRT, Helv. Phys. Acta 65 (1992) 679–721 | not checked | memory | UNVERIFIED |
| B14 | Osterwalder–Seiler, Ann. Phys. 110 (1978) 440–471 | not checked | memory | UNVERIFIED |
| B15 | Brydges–Fröhlich–Seiler I–II | not checked | memory | UNVERIFIED |
| B16 | Seiler, LNP 159 (1982) | not checked | memory | UNVERIFIED |
| B17 | Gribov, Nucl. Phys. B139 (1978) 1–19 | not checked | memory | UNVERIFIED |
| B18 | Singer, CMP 60 (1978) 7–12 | not checked | memory | UNVERIFIED |
| B19 | Wilson, Phys. Rev. 179 (1969) 1499–1512 | not checked | memory | UNVERIFIED |
| B20 | SVZ I, Nucl. Phys. B147 (1979) 385–447 | title, volume, pages, DOI | ui.adsabs.harvard.edu/abs/1979NuPhB.147..385S | VERIFIED |
| B21 | SVZ II, Nucl. Phys. B147 (1979) 448–518 | title, volume, pages | ADS / Scopus records | VERIFIED |
| B22 | SVZ III, Nucl. Phys. B147 (1979) 519–534 | title, volume, pages | ADS / Scopus records | VERIFIED |
| B23 | Novikov et al., Nucl. Phys. B249 (1985) 445–471 | not checked | memory | UNVERIFIED |
| B24 | 't Hooft, Erice 1977, Plenum 1979, 943–982 | title, venue, pages | Springer/OuCI chapter record | VERIFIED |
| B25 | David, Nucl. Phys. B234 (1984) 237–251 | title, volume, pages, DOI | Scopus 0010875546; arXiv:2111.11951 citation list | VERIFIED |
| B26 | Mueller, Nucl. Phys. B250 (1985) 327–350 | title, volume, pages, DOI | sciencedirect.com/science/article/abs/pii/0550321385904857 | VERIFIED |
| B27 | David, Nucl. Phys. B263 (1986) 637–648 | title, volume, pages | Scopus/Google Patents citation records | VERIFIED |
| B28 | Beneke, Phys. Rep. 317 (1999) 1–142 | not checked | memory | UNVERIFIED |
| B29 | Di Giacomo–Rossi, Phys. Lett. B100 (1981) 481–484 | title, volume, pages, DOI | ui.adsabs.harvard.edu/abs/1981PhLB..100..481D | VERIFIED |
| B30 | Banks–Horsley–Rubinstein–Wolff, Nucl. Phys. B190 (1981) 692–704 | not checked | memory | UNVERIFIED |
| B31 | Campostrini–Curci–Di Giacomo–Paffuti, Z. Phys. C32 (1986) 377–382 | title, volume, pages, DOI | link.springer.com/article/10.1007/BF01551834 | VERIFIED |
| B32 | Lüscher, Nucl. Phys. B219 (1983) 233–261 | not checked | memory | UNVERIFIED |
| B33 | Lüscher, CMP 104/105 (1986) | not checked | memory | UNVERIFIED |
| B34 | Magnen–Rivasseau–Sénéor, CMP 155 (1993) 417–430 | not checked | memory | UNVERIFIED |
| B35 | Rivasseau, book, Princeton 1991 | not checked | memory | UNVERIFIED |

**Counts.** 35 entries; 18 verified (B1–B8, B10, B20–B22, B24–B27, B29, B31), of which B1–B2 are
secondary. All load-bearing entries are verified. The 17 unverified entries are standard citations
carried from memory and are flagged; if a v3 footnote depends on one of them, check it before
submission.

**Sources consulted in this session.** Springer article pages (DOIs listed above); NASA ADS abstract
pages for SVZ I and Di Giacomo–Rossi; ScienceDirect record for Mueller; Scopus records for David 1984
and 1986; INSPIRE/zbMATH citation records for the Balaban series and Gawędzki–Kupiainen. Web-search
budget for the turn was exhausted before the remaining 17 entries could be checked.

**Proof-status statement.** Nothing in this file is a proof of the Yang–Mills mass gap. The strongest
rigorous statement in the era is Balaban's UV stability (B7–B8), under the caveats of §1.4(N7). The
gluon condensate is a phenomenological parameter (B20–B22) with a known `O(Λ⁴)` ambiguity (B24–B27);
the lattice determinations (B29–B31) choose a scheme and do not resolve it.
