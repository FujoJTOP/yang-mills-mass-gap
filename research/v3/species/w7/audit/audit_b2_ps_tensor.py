#!/usr/bin/env python3
"""audit_b2_ps_tensor.py -- is ||T^{tensor V} - T0^{tensor V}|| = V*eps1 valid where the [E] table uses it?

eps1 = 0.7871, dist = 0.001944 are ps_witness.py's own numbers (single-site).
Two bounds, both elementary, for ||T||=||T0||=1:
    ||T^V - T0^V|| <= sum_{k=1..V} ||T||^{V-k}||T0||^{k-1} eps1 = V eps1    (triangle; first order)
    ||T^V - T0^V|| <= ||T||^V + ||T0||^V = 2                               (hard)
So V*eps1 stops being informative once V*eps1 > 2, i.e. V > 2.54. The table goes to V=65536.
Exact diagonal counterexample (both spectra in [0,1], eps1 = 0.7871) shows the true norm can DECREASE.
"""
import numpy as np

eps1, dist = 0.7871, 0.001944
print(f"eps1 = {eps1}, dist = {dist}")
print(f"dist/eps1 = {dist/eps1:.6f}   (Rouche needs V*eps1 < dist  =>  V < {dist/eps1:.2e} < 1)")
print(f"single site eps1/dist = {eps1/dist:.1f}")
print()
print(" V   V*eps1    true_norm(diag model: T0=diag(1,0.1), T=diag(1,0.8871), eps1=0.7871)")
t0, t1 = 0.1, 0.8871
for V in (1, 2, 3, 4, 8, 16, 64, 65536):
    true = abs(t1 ** V - t0 ** V)
    print(f"{V:6d} {V*eps1:10.4f}   {true:10.4f}    {'(V eps1 > 2: vacuous)' if V*eps1 > 2 else ''}")
print()
print(" => the 'V*eps1' column of ps_witness [E] is a first-order bound valid only for V*eps1 << 1;")
print("    for V >= 3 it exceeds the hard bound 2 (and the exact value is decreasing). The decisive")
print("    number is V=1: eps1/dist = 404.8 > 1, i.e. no positive integer V satisfies the Rouche hypothesis.")
print()
# also: is V*eps1 attained? it is sharp iff E shares the top singular direction at every factor:
# T0 = 1 (scalar), T = 1 + eps * |phi><phi|  ->  norm = (1+eps)^V - 1 = eps V + O(eps^2)
eps = 0.01
for V in (1, 2, 5, 10):
    print(f"  sharpness check eps=0.01, V={V:2d}: (1+eps)^V-1 = {(1+eps)**V-1:.6f}   V*eps = {V*eps:.6f}")
