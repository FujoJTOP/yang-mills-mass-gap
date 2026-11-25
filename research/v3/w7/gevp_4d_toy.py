#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
gevp_4d_toy.py -- validate the GEVP read-out of gevp_4d.py against an EXACT 2-state model.

C_ij(t) = sum_n A_ni A_nj [exp(-m_n t) + exp(-m_n (L-t))],  n = 1,2 ,
with a chain of operator bases of decreasing quality (rank of the overlap matrix A).  The
GEVP must return (m1, m2) exactly when C(t0) is well conditioned, and must FAIL (second
eigenvalue = noise) when the basis is nearly collinear -- which is the situation of the
3 APE-smeared plaquette operators at L <= 8.

Writes gevp_4d_toy.json.  [V]  No QED.
"""
import json, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import gevp_4d as G          # noqa: E402

L = 8
m = np.array([1.0, 1.6])                 # rho^2 = (1.6/1.0)^2 = 2.56
out = {}

# operator bases: 3 operators with controlled collinearity
bases = {
    'orthogonal(well-conditioned)': np.array([[1.0, 0.0],
                                              [0.0, 1.0],
                                              [0.6, 0.8]]),
    'nearly_collinear(eps=0.1)':    np.array([[1.0, 0.0],
                                              [1.0, 0.1],
                                              [1.0, -0.1]]),
    'nearly_collinear(eps=0.02)':   np.array([[1.0, 0.0],
                                              [1.0, 0.02],
                                              [1.0, -0.02]]),
}
for name, A in bases.items():
    C = np.zeros((3, 3, L))
    for t in range(L):
        for n in range(2):
            C[:, :, t] += np.outer(A[:, n], A[:, n]) * (np.exp(-m[n] * t) +
                                                        np.exp(-m[n] * (L - t)))
    rows = G.table_for(C, L, (0, 1, 2), name)
    out[name] = dict(overlap_matrix=A.tolist(),
                     gram_eigs=np.linalg.eigvalsh(
                         np.outer(A[:, 0], A[:, 0]) + np.outer(A[:, 1], A[:, 1]))[::-1].tolist(),
                     table=rows)
    r = [x for x in rows if x['ok']]
    print(f"\n{name}  (exact m1={m[0]}, m2={m[1]}, rho2={ (m[1]/m[0])**2 :.4f})")
    for x in r:
        print(f"   t0={x['t0']} t={x['t']}  m1={x['m1']:.6f} m2={x['m2']:.6f} "
              f"rho2={x['rho2']:.6f}  cond={x['cond']:.3g}")

json.dump(out, open(os.path.join(HERE, 'gevp_4d_toy.json'), 'w'), indent=1,
          default=lambda o: float(o) if isinstance(o, (np.floating, np.integer)) else str(o))

# ---------------------------------------------------------------- bias vs L
# Same exact 2-state model (rho^2 = 2.56), well-conditioned basis, every L of interest:
# the ONLY error is the periodic (cosh) image -- zero statistical noise.
A = bases['orthogonal(well-conditioned)']
print('\n# periodic-image bias in the plain-log rho^2 (exact model, zero noise, rho2=2.56):')
bias = {}
for LL in (4, 6, 8, 12, 16):
    C = np.zeros((3, 3, LL))
    for t in range(LL):
        for n in range(2):
            C[:, :, t] += np.outer(A[:, n], A[:, n]) * (np.exp(-m[n] * t) +
                                                        np.exp(-m[n] * (LL - t)))
    rows = [x for x in G.table_for(C, LL, (0, 1, 2), 'bias') if x['ok']]
    bias[LL] = rows
    s = '  '.join(f"(t0={x['t0']},t={x['t']}): rho2={x['rho2']:.4f}"
                  f"{'' if 'rho2_cosh' not in x else f' / cosh {x['rho2_cosh']:.4f}'}"
                  for x in rows)
    print(f'  L={LL:2d}  {s}')
out['periodic_bias'] = bias
json.dump(out, open(os.path.join(HERE, 'gevp_4d_toy.json'), 'w'), indent=1,
          default=lambda o: float(o) if isinstance(o, (np.floating, np.integer)) else str(o))
print('\n# -> gevp_4d_toy.json')
