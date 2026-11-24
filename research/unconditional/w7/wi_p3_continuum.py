"""
wi_p3_continuum.py -- P3 approach-law discriminator under a background continuum.
Lightweight, numpy only. Run: py -3.12 -I wi_p3_continuum.py

Model: connected correlator C_O(tau) = c1 e^{-e1 tau} + rho0 * int_{Ethr}^{Emax} e^{-e tau} de
       (atom at e1 weight c1, PLUS a flat continuum [Ethr,Emax] density rho0).
Effective mass m_O(tau) = -d/dtau log C_O.

Question: as tau -> infty does m_O(tau) - e1 -> 0 EXPONENTIALLY (isolated atom below the cut),
or does the continuum contaminate it into a POWER LAW?

Prediction (Watson's lemma): atom below cut (Delta = Ethr - e1 > 0):
    m_O(tau) - e1 ~ (rho0/c1)(Delta/tau + 1/tau^2) e^{-Delta tau}   -> EXPONENTIAL (slope -Delta in log vs tau)
Pure continuum (c1 = 0, no atom):
    m_O(tau) - Ethr ~ -1/tau                                       -> POWER LAW  (slope ~ -1/tau, ~flat in log)
"""
import numpy as np

def I_cont(tau, Ethr, Emax, rho0):
    # rho0 * int_{Ethr}^{Emax} e^{-e tau} de = rho0 (e^{-Ethr tau}-e^{-Emax tau})/tau
    return rho0 * (np.exp(-Ethr * tau) - np.exp(-Emax * tau)) / tau

def I_cont_p(tau, Ethr, Emax, rho0):
    return rho0 * ((-Ethr*np.exp(-Ethr*tau) + Emax*np.exp(-Emax*tau))/tau
                   - (np.exp(-Ethr*tau) - np.exp(-Emax*tau))/tau**2)

def m_eff(tau, e1, c1, Ethr, Emax, rho0):
    C  = c1*np.exp(-e1*tau) + I_cont(tau, Ethr, Emax, rho0)
    Cp = -c1*e1*np.exp(-e1*tau) + I_cont_p(tau, Ethr, Emax, rho0)
    return -Cp/C

def log_slope(taus, ys):
    # d log|y| / d tau  (linear fit)
    ly = np.log(np.abs(ys))
    A = np.vstack([taus, np.ones_like(taus)]).T
    slope, _ = np.linalg.lstsq(A, ly, rcond=None)[0]
    return slope

e1, c1 = 1.0, 0.5
Ethr, Emax, rho0 = 2.0, 60.0, 0.3
Delta = Ethr - e1
taus = np.linspace(6.0, 30.0, 200)

print("="*78)
print("(A) ATOM + BACKGROUND CONTINUUM   (e1=1, Ethr=2, Delta=1, c1=0.5, rho0=0.3)")
print("="*78)
m_atom = np.array([m_eff(t, e1, c1, Ethr, Emax, rho0) for t in taus])
dev_atom = m_atom - e1
s_atom = log_slope(taus, dev_atom)
print(f"  m_eff(tau=6) = {m_eff(6,e1,c1,Ethr,Emax,rho0):.6f}   m_eff(tau=30) = {m_eff(30,e1,c1,Ethr,Emax,rho0):.6f}")
print(f"  |m-e1|: tau=6 -> {abs(dev_atom[0]):.3e} ; tau=30 -> {abs(dev_atom[-1]):.3e}")
print(f"  d log|m-e1| / dtau  (fit) = {s_atom:+.4f}    (predicted ~ -Delta = {-Delta:+.1f})")
print(f"  => {'EXPONENTIAL' if s_atom < -0.3 else 'not exponential'} approach (slope far from 0)")

print()
print("="*78)
print("(B) PURE CONTINUUM  (no atom: c1=0, floor=Ethr=2)")
print("="*78)
m_cut = np.array([m_eff(t, e1, 0.0, Ethr, Emax, rho0) for t in taus])  # c1=0
dev_cut = m_cut - Ethr
s_cut = log_slope(taus, dev_cut)
print(f"  m_eff(tau=6) = {m_eff(6,e1,0.0,Ethr,Emax,rho0):.6f}   m_eff(tau=30) = {m_eff(30,e1,0.0,Ethr,Emax,rho0):.6f}")
print(f"  |m-Ethr|: tau=6 -> {abs(dev_cut[0]):.3e} ; tau=30 -> {abs(dev_cut[-1]):.3e}")
print(f"  d log|m-Ethr| / dtau (fit) = {s_cut:+.4f}    (predicted ~ -1/tau -> near 0)")
print(f"  => {'POWER LAW (contaminated)' if s_cut > -0.3 else 'exponential'} approach (=1/tau, dev*k ~const)")

print()
print("="*78)
print("(C) ACID: FREE MASSLESS  (continuum from 0, no atom, floor=0)")
print("="*78)
m_free = np.array([m_eff(t, 0.0, 0.0, 0.0, Emax, 1.0) for t in taus])
s_free = log_slope(taus, m_free - 0.0)
print(f"  m_eff(tau=30) = {m_free[-1]:.3e}  -> floor 0 (no atom). dev*(tau) ~ {abs(m_free[-1])*30:.3f} (~const => power law)")
print(f"  slope dlog|m-0|/dtau = {s_free:+.4f} (near 0 => power law, acid-clean: says NO atom)")

print()
print("="*78)
print("SEPARATION: atom slope={:+.3f}  vs  pure-continuum slope={:+.3f}".format(s_atom, s_cut))
print("=> the approach-law DISTINGUISHES atom+continuum (exponential) from pure continuum (power law),")
print("   EVEN with the continuum present. It is robust where the complex-zero A is not.")
print("="*78)
