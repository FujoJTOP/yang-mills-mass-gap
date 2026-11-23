# P1 — L0: the free Polchinski contraction

**Setup.** Lattice $T^4_L$, spacing $a$, $N_s=L/a$ sites/dir, $g=\mathfrak{su}(N)$, $A\in g$-valued 1-forms. Linearized action $S_0=\tfrac12\int\mathrm{Tr}\,F\star F$, $F=dA$. Gaussian measure $\mu$ with covariance $C=g_0^2(-\Delta_a)^{-1}$ on transverse forms; lattice spectrum $\lambda_a(p)=\sum_\mu \tfrac{4}{a^2}\sin^2\!\tfrac{p_\mu a}{2}$, $p\in\tfrac{2\pi}{L}\mathbb Z^4$. The linearized Polchinski/Ornstein–Uhlenbeck semigroup is $P_t=e^{-tL_0}$, $L_0=-\Delta_A$ (Gaussian generator); it is diagonal in $p$-modes, multiplier $e^{-t\lambda_a(p)}$ on the field. $P_{>4}$ projects onto gauge-invariant local densities of engineering dimension $\Delta\ge 6$.

**True lemma (UV / fixed volume).** Gauge-invariant functionals depend on $F=dA$ only, so they lie in $\ker(d)^\perp$: the harmonic (constant) mode is exact, hence pure gauge, and decouples — **zero modes are harmless**. Then
$$\|P_{>4}P_t\|_{L^2(\mu)}\le e^{-c_L t},\qquad c_L=\lambda_a(\tfrac{2\pi}{L}e_1)=\tfrac{4\pi^2}{L^2}\bigl(1+O(a^2/L^2)\bigr),$$
uniform in $a$. Lattice Poincaré inequality: $\lambda_a(p)\ge c_0|p|^2$, $c_0=4/\pi^2$, for $|p|a\le\pi/2$. With an IR mass $m$ (replace $-\Delta_a\to-\Delta_a+m^2$, or restrict $|p|\ge m$), $c=m^2(1+o(1))$ **uniformly in $a$ and $L$.** Sharp.

**Obstruction — L0 as stated ("uniform in volume") is false.** The free transverse gap on $T^4_L$ is exactly
$$\lambda_{\min}(a,L)=\tfrac{4}{a^2}\sin^2\tfrac{\pi a}{L}=\tfrac{4\pi^2}{L^2}+O(a^2/L^4)\to0 .$$
The mode $p_{\min}=(2\pi/L,0,0,0)$ is transverse and non-exact (physical), and every local $\Delta\ge6$ density has nonzero Fourier amplitude there; its slow component evolves diagonally with rate $\lambda_{\min}$. Thus $\forall c>0\ \exists L$ ($4\pi^2/L^2<c$) and a dim-6 functional with $\|P_{>4}P_tf\|/\|f\|\ge e^{-t\lambda_{\min}}>e^{-ct}$ for $t>0$: no volume-uniform $c$ exists. The failure is **masslessness (zero IR spectral gap) of the free propagator** — exactly the G2 IR sector, not UV structure.

*Not the obstruction:* gauge orbits / zero modes (pure gauge, projected out), the FP measure Jacobian (field-independent in $L_0$), orbit non-compactness (physical transverse moduli compact on the torus). All are invisible on gauge-invariant functionals.

**Consequence.** L0 holds in the only sense a free-graph statement can: uniformly in cutoff $a$ at fixed volume, or uniformly in $a,L$ given an IR gap $m>0$, with proven constant $c=m^2$. A gap-free, volume-uniform L0 is impossible and is precisely G2, so it cannot be a UV lemma.

**Refs.** Polchinski, NPB 231 (1984) 269; Bauerschmidt–Brydges–Slade, LNM 2242 (2019), linearized-flow analysis; Dimock, Rev. Math. Phys. 25 (2013) (lattice Gaussian/Balaban RG); Gawędzki–Kupiainen, CMP 77 (1980) 31 (Gaussian RG); lattice Poincaré inequality (elementary).
