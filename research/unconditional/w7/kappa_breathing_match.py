"""
kappa_breathing_match.py -- KAPPA-BREATHING-MATCH witness.

Two rigd-string relations, ONE kappa:
  lightest  0++:  m_0++/sqrt(sigma) = 2*pi*sqrt(2*kappa)  =>  kappa = (m/sqrt(sigma))^2/(8*pi^2)
  breathing 0++*: m_*/sqrt(sigma)   = 2/sqrt(kappa)      =>  kappa = 4/(m_*/sqrt(sigma))^2
(breathing mode: E''=2*pi*kappa/R*^3, M=2*pi*sigma*R*, R*=sqrt(kappa/2sigma) => m_br=2*sqrt(sigma/kappa))

Model-internal excitation ratio: m_*/m_0++ = [2/sqrt(k)]/[2*pi*sqrt(2k)] = 1/(pi*sqrt(2)*k).

Lattice first-excited scalar 0++*: Morningstar-Peardon 1999, PRD 60 034509, hep-lat/9901004:
  m(0++)=3.55(7) r0^-1, m(0++*)=5.46(22) r0^-1  =>  ratio = m*/m = 1.538.
Run: py -3.12 -I
"""
import math

pi = math.pi

m_light = 3.9                     # m(0++)/sqrt(sigma), project lattice value
kappa_light = m_light**2 / (8*pi**2)
m_breath_pred = 2 / math.sqrt(kappa_light)      # model prediction for m(0++*)/sqrt(sigma)

# model-internal predicted excitation ratio (normalisation independent)
model_ratio = 1 / (pi*math.sqrt(2)*kappa_light)

# lattice first-excited scalar (Morningstar-Peardon 1999)
m0pp, m0ppst = 3.55, 5.46
lat_ratio = m0ppst/m0pp
m_exc_lat = m_light * lat_ratio                 # m(0++*)/sqrt(sigma) implied
kappa_exc = 4 / m_exc_lat**2

print("=== KAPPA-BREATHING-MATCH ===")
print(f"m(0++)/sqrt(sigma)              = {m_light}")
print(f"kappa from LIGHTEST             = (m/sqrt s)^2/(8 pi^2) = {kappa_light:.6f}")
print(f"model pred m(0++*)/sqrt(sigma)  = 2/sqrt(kappa)         = {m_breath_pred:.4f}")
print()
print(f"Morningstar-Peardon 1999: m(0++)={m0pp} r0^-1, m(0++*)={m0ppst} r0^-1")
print(f"lattice excitation ratio m*/m   = {lat_ratio:.4f}")
print(f"lat m(0++*)/sqrt(sigma)         = m_light*ratio         = {m_exc_lat:.4f}")
print(f"kappa from EXCITED 0++*         = 4/(m*/sqrt s)^2       = {kappa_exc:.6f}")
print()
print(f"model-internal ratio m*/m       = 1/(pi sqrt2 kappa)    = {model_ratio:.4f}")
print(f"  vs lattice ratio                                       = {lat_ratio:.4f}  (factor {lat_ratio/model_ratio:.3f})")
print(f"kappa_light / kappa_exc         = {kappa_light/kappa_exc:.3f}")
print()
disagree = abs(kappa_light-kappa_exc)/kappa_light > 0.10
print("VERDICT:", "DISAGREE -> model falsified at the excited level" if disagree else "AGREE -> kappa pinned")
