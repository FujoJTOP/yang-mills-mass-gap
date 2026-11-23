"""2-operator GEVP (n=0 and n=8) on a pilot_numba.py output."""
import sys
import numpy as np


def Cmat(Sa, Sb):
    n, L = Sa.shape
    ma, mb = Sa.mean(), Sb.mean()
    Xa, Xb = Sa - ma, Sb - mb
    return np.array([np.mean(Xa * np.roll(Xb, -t, axis=1)) for t in range(L)])


def main(f):
    d = dict(np.load(f, allow_pickle=True))
    L = int(d['L']); beta = float(d['beta'])
    S0, Sn = d['S0'], d['Sn']
    n = len(S0)
    C = {(a, b): None for a in (0, 1) for b in (0, 1)}
    ops = {0: S0, 1: Sn}
    for a in (0, 1):
        for b in (0, 1):
            C[(a, b)] = Cmat(ops[a], ops[b])
    t0 = 1
    M0 = np.array([[C[(a, b)][t0] for b in (0, 1)] for a in (0, 1)])
    print(f"{f}: L={L} beta={beta} N={n}")
    lam = {}
    for t in (2, 3, 4, 5, 6, 8):
        Mt = np.array([[C[(a, b)][t] for b in (0, 1)] for a in (0, 1)])
        ev = np.linalg.eigvals(np.linalg.solve(M0.T, Mt.T).T)
        ev = np.sort(np.real(ev))[::-1]
        lam[t] = ev
        print(f"  t={t}: GEVP eigenvalues {np.round(ev,5)}")
    print("  leading-state m_eff = -ln(lambda(t)/lambda(t+1)):")
    for t in range(2, 6):
        if t in lam and (t + 1) in lam:
            r = lam[t][0] / lam[t + 1][0]
            if r > 0:
                print(f"    t={t}->{t+1}: m={np.log(r):.3f}  xi/a={1/np.log(r):.2f}")
    for t in (2, 3, 4, 5):  # consecutive pairs including t=7
        Mt = np.array([[C[(a, b)][t] for b in (0, 1)] for a in (0, 1)])
        Mt1 = np.array([[C[(a, b)][t + 1] for b in (0, 1)] for a in (0, 1)])
        ev = np.sort(np.real(np.linalg.eigvals(np.linalg.solve(M0.T, Mt.T).T)))[::-1]
        ev1 = np.sort(np.real(np.linalg.eigvals(np.linalg.solve(M0.T, Mt1.T).T)))[::-1]
        r = ev[0] / ev1[0]
        if r > 0:
            print(f"    (t={t}->{t+1}) m={np.log(r):.3f}  xi/a={1/np.log(r):.2f}")


if __name__ == '__main__':
    for f in sys.argv[1:]:
        main(f)
