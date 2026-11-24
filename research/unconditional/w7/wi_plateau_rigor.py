"""wi_plateau_rigor.py -- witness for WI-PLATEAU-RIGOR.md.
m_eff(tau) = -d/dtau log C, C = sum_n c_n exp(-e_n tau) (connected correlator).
Shows: (P1) two operators -> same e1; (P2) spread decays as exp(-(e2-e1)tau);
(P3) atom = exponential approach, continuum = 1/tau power law.
Runs clean, numpy only.  py -3.12 -I wi_plateau_rigor.py
"""
import numpy as np


def m_eff(tau, e, c):
    w = c * np.exp(-np.asarray(e) * tau)
    return float(np.sum(w * e) / np.sum(w))


print("=" * 70)
print("P1/P2: atom e=(1,3,5); two operators differ only in couplings c")
print("=" * 70)
e = np.array([1.0, 3.0, 5.0])
c_O = np.array([1.0, 0.50, 0.10])
c_P = np.array([1.0, 0.20, 0.30])          # different c2/c1 ratio
print(f"{'tau':>5} {'m_O':>10} {'m_Oprime':>10} {'spread':>10} {'spread/e^-2t':>13}")
prev = None
for tau in [0.5, 1.0, 2.0, 3.0, 4.0, 6.0]:
    a = m_eff(tau, e, c_O)
    b = m_eff(tau, e, c_P)
    sp = abs(a - b)
    print(f"{tau:5.1f} {a:10.5f} {b:10.5f} {sp:10.3e} {sp/np.exp(-2*tau):13.4f}")
print("-> both -> e1=1 (operator-independent limit); spread ~ exp(-(e2-e1)tau), (e2-e1)=2.")

print()
print("=" * 70)
print("P3: approach LAW -- atom (exponential) vs continuum edge (power law 1/tau)")
print("=" * 70)
print("atom e1=1: m(tau)-1 should fall like exp(-2 tau)")
for tau in [1.0, 2.0, 4.0, 6.0]:
    d = m_eff(tau, e, c_O) - 1.0
    print(f"  tau={tau:4.1f}  m-1 = {d:10.3e}   (m-1)*exp(2tau) = {d*np.exp(2*tau):8.4f}")
print("continuum [1,6], rho=exp(-(e-1)): m(tau)-1 should fall like 1/tau")
ec = np.linspace(1.0, 6.0, 4001)
rho = np.exp(-(ec - 1.0))


def m_cont(tau):
    w = rho * np.exp(-ec * tau)
    return float(np.sum(w * ec) / np.sum(w))


for tau in [1.0, 2.0, 4.0, 8.0, 16.0]:
    d = m_cont(tau) - 1.0
    print(f"  tau={tau:5.1f}  m-1 = {d:10.3e}   (m-1)*tau = {d*tau:8.4f}")
print("-> atom: (m-1)*exp(2tau) ~ const (exponential); cut: (m-1)*tau ~ const (1/tau).")
print("   Clean rigorous atom/cut discriminator = the SHAPE of the approach.")
print()
print("Labels: P1/P2/P3 [T]; uniformity in (L,a) = H1 [O]. No QED.")
