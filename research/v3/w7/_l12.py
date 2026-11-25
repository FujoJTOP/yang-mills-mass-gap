#!/usr/bin/env python3
# L=12 chi_flat-only feasibility point (N=300, no Landau, no Gribov) to extend the plateau.
import sys, os, json
sys.path.insert(0, 'D:/Dev/math/research/v3/w7')
import m1_4d_run as M
RES = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'm1_4d_run_results.json')
print("# L=12 chi_flat point (N=300, gap=6, no Landau)", flush=True)
r = M.run_volume(12, 2.5, 300, 800, 6, 404, n_gf=0, n_gribov=0, n_copy=1, rng_seed=9012)
d = json.load(open(RES))
d['L12'] = r
d.setdefault('meta', {})['L12_note'] = 'N=300, gap=6, no Landau/glueball-resolved'
json.dump(d, open(RES, 'w'), indent=1,
          default=lambda o: int(o) if hasattr(o, 'dtype') else float(o))
print(f"  L=12 chi_flat={r['chi_flat']:.4f}+-{r['chi_flat_se']:.4f} (tau={r['chi_flat_tau']:.2f}, "
      f"Neff={r['chi_flat_neff']:.0f})  C_P={[round(c,4) for c in r['C_P']]}", flush=True)
print("saved L12", flush=True)
