"""A4-SCALE-SEP: xi/a at beta=2.5 from the corpus's own numbers.

Two independent routes, no new dynamics:
  (R1) direct GEVP gap of the 0++ transfer-matrix operator, SU(2) Wilson, beta=2.5
       -> research/unconditional/tm_gap_results.csv  (gevp48_tL at the plateau t)
  (R2) ratio chain m_0++/sqrt(sigma) * sqrt(sigma)*a  (T1-POLYMER-RG.md:241-242)
xi/a = 1/(m a);  relative error propagates unchanged.

Run: py -3.12 research/unconditional/w7/a4_scale_sep.py
"""
import csv
import os

HERE = os.path.dirname(os.path.abspath(__file__))


def _root(start):
    """Walk up until the repo marker is found, so the script works from any w7/."""
    d = start
    while not os.path.isdir(os.path.join(d, 'research')):
        d = os.path.dirname(d)
    return d


ROOT = _root(HERE)
CSV = os.path.join(ROOT, 'research', 'unconditional', 'tm_gap_results.csv')
LOG = os.path.join(HERE, 'a4_scale_sep.log')

out = []


def w(s=''):
    out.append(s)
    print(s)


def xi(ma, err=0.0):
    x = 1.0 / ma
    return x, x * err / ma


# ---- R1: direct GEVP, SU(2) Wilson beta=4/g^2, plateau t = L/2 ----
w('=' * 78)
w('A4-SCALE-SEP  xi/a at beta=2.5  (recomputed from corpus data, no new fits)')
w('=' * 78)
w()
w('[R1] direct 0++ GEVP gap, SU(2) Wilson beta=4/g^2, src: tm_gap_results.csv')
w('     xi/a = 1/(m a), err relative')
w('%-10s %6s %5s %6s %10s %14s' % ('tag', 'L', 'beta', 't_pl', 'ma', 'xi/a'))
r1 = []
with open(CSV, newline='') as fh:
    for row in csv.DictReader(fh):
        ma = float(row['gevp48_tL'])
        err = float(row['gevp48_tL_err'])
        x, xe = xi(ma, err)
        r1.append((row['tag'], float(row['beta']), ma, err, x, xe))
        w('%-10s %6s %5s %6s %10.4f %8.3f(%.3f)'
          % (row['tag'], row['L'], row['beta'], row['gevp48_tL_t'], ma, x, xe))
r1_b25 = [r for r in r1 if r[1] == 2.5]
lo = min(r[4] - r[5] for r in r1_b25)
hi = max(r[4] + r[5] for r in r1_b25)
w()
w('  beta=2.5 SU(2), direct:  ma in [%.3f, %.3f]  =>  xi/a in [%.2f, %.2f]'
  % (min(r[2] for r in r1_b25), max(r[2] for r in r1_b25), lo, hi))
w('  note: L=4,6 are the *temporal* extent => T=1/(La) > 0; LATTICE-TM-GAP.md:26')
w('        puts beta_c(N_tau=4)=2.30, 2.43(N_tau=6), so these are finite-T')
w('        screening masses, not the T=0 glueball. Same O(1) conclusion either way.')
w()

# ---- R2: ratio chain, beta=2.5, SU(2) ----
w('[R2] ratio chain at beta=2.5, SU(2), src: research/v3/T1-POLYMER-RG.md:241-242')
rows2 = [
    ('ma = 4 ln(4/beta)', 4.0 * __import__('math').log(4.0 / 2.5),
     'strong-coupling leading order (overestimates at the crossover)'),
    ('ma from m/sqrt(s)=3.5, sqrt(s)a=0.22', 3.5 * 0.22,
     'scaling-region ratios (Morningstar-Peardon class)'),
]
for name, ma, why in rows2:
    x, _ = xi(ma)
    w('  %-38s ma=%.3f  xi/a=%.3f   # %s' % (name, ma, x, why))
w()
vals = [xi(ma) for _, ma, _ in rows2] + [(r[4], r[5]) for r in r1_b25]
w('  combined beta=2.5 SU(2) envelope: xi/a in [%.2f, %.2f]'
  % (min(v[0] for v in vals), max(v[0] for v in vals)))
w('  => xi/a = 1.0 +- 0.3 (O(1), <= 1.3); NOT >> 1')
w()

# ---- SU(3) reference points: where xi/a actually grows ----
w('[R3] SU(3) Wilson beta=6/g^2 reference points (src: research/lattice-cpb/REPORT.md)')
w('  beta<=3.3  : xi/a < 0.3 (unresolved; am>4.4 2sigma)      REPORT.md:369-373')
for ma, tag, line in [(1.073, 'beta=6.0  (n=8, L=6)', 'REPORT.md:357'),
                      (0.886, 'beta=6.30 (n=8, L=24)', 'REPORT.md:1060')]:
    x, _ = xi(ma)
    w('  %-24s ma=%.3f -> xi/a=%.2f                 %s' % (tag, ma, x, line))
w('  beta~7.6-7.8: xi/a=7-8.5 *extrapolated only* (two-loop from the beta=6')
w('                calibration, never measured)            REPORT.md:1087-1095')
w()

# ---- verdicts ----
w('-' * 78)
w('(a) SU(3): beta_c=5.6924(2) [hep-lat/0701007], beta=2.5 is 44% BELOW =>')
w('    strong-coupling side, g^2=2.4. "physical point" is FALSE for SU(3).')
w('    SU(2): beta_c(N_t=4)=2.2986(6) [2211.16273], beta=2.5 is +8.7% =>')
w('    crossover / entrance of the scaling window, not its interior.')
w('    [src: research/v3/species/BETA-CONVENTION.md:102-107,121-122]')
w('(b) xi/a = O(1) (<=1.3) at beta=2.5 => no scale separation => no continuum')
w('    limit to start from. HOLDS.')
w('(c) no measured xi/a >> 1 anywhere in the corpus; the >>1 values are')
w('    two-loop extrapolations at beta~7.6-7.8 (SU(3)). FAILS as data.')
w('-' * 78)

with open(LOG, 'w', encoding='utf-8') as fh:
    fh.write('\n'.join(out) + '\n')
print('\n[written] %s' % LOG)
