# -*- coding: utf-8 -*-
"""
w7/u13_threeway.py  --  U13 三向测试 (candidate x model test)

四个"谱隙候选" x 四个模型:
  M6     : 谐振子(tower n*W) (x) 自由模(E=0 零模)      G''=2W>0 , m=0       -> 候选应为 0
  2D-YM  : E_j=(g^2/2)j(j+1), j=1/2,1,...              G''=0   , m=3g^2/8>0  -> 候选应为 >0
  FREE4D : rho(E)=E^{1/2} 连续谱, 精确 C(t)=Gamma(3/2)t^{-3/2}   m=0        -> 候选应为 0
  PLAQ   : 单 plaquette, 无 transfer 方向                m 无定义             -> 候选应无定义(不是 0)

候选:
  c1 = lim_{t->inf} -(1/t) log|C(t)|
  c2 = lim_{t->inf} C(2t)/C(t)
  c3 = inf supp mu            (mu = C 的谱测度, 已扣除真空 delta_0)
  c4 = E1 = lim_{t->inf} log C(t)/C(t+1)     (局域对数斜率 = 首激发能)

三条路互不共用数据: c3 从 mu 的支撑读; c1/c2/c4 只用 C(t)。
weights 形状不影响 inf supp mu, 只影响 c1/c2/c4 的收敛速度。
FREE4D 用闭式 (不用网格): 网格 dE 只能把 c1 钉在网格底, 会假报 c1->dE。
无随机数, 只写 stdout。
用法: py -3.12 w7/u13_threeway.py > w7/u13_threeway.log 2>&1
"""
from math import lgamma, log, exp
import numpy as np

W = 1.0          # 谐振子频率 -> G'' = 2W
G2 = 1.0         # 2D YM 耦合, g=1 -> m = 3/8
TSCAN = [25, 50, 100, 200, 500, 1000, 2000]
FITFROM = 100     # 外推只用 T >= FITFROM
MTOL = 1e-12      # 取支撑的权重门限
ZERO = 5e-3       # 判"=0"的容差 (只用于把外推出的极限归类)


def logC_of_mu(mu):
    """log C(t) = logsumexp_i (log w_i - E_i t)  (对数域, 防下溢)"""
    E, w = mu
    keep = w > 0.0
    E, lw = E[keep], np.log(w[keep])

    def L(t):
        return float(np.logaddexp.reduce(lw - float(t) * E))
    return L


def make_M6():
    N = 80
    E = np.concatenate(([0.0], W * np.arange(1.0, N + 1.0)))
    w = np.concatenate(([1.0], 0.5 ** np.arange(1.0, N + 1.0)))
    mu = (E, w)
    return dict(name='M6', m=0.0, Gpp=2 * W, mu=mu, L=logC_of_mu(mu),
                note='osc tower n*W (w=2^-n) + free module zero mode E=0')


def make_YM2():
    J = 40
    j = 0.5 * np.arange(1.0, J + 1.0)
    mu = (0.5 * G2 * j * (j + 1.0), (2 * j + 1.0) ** -2.0)
    return dict(name='2D-YM', m=0.375 * G2, Gpp=0.0, mu=mu, L=logC_of_mu(mu),
                note='E_j=(g^2/2)j(j+1), lightest j=1/2 -> 3g^2/8=0.375')


def make_FREE():
    lg = lgamma(1.5)          # = log Gamma(3/2)
    return dict(name='FREE4D', m=0.0, Gpp=0.0, mu=None, supp_inf=0.0,
                L=lambda t: lg - 1.5 * np.log(float(t)),
                note='rho(E)=E^1/2 连续谱, 精确 C(t)=Gamma(3/2)t^{-3/2}, inf supp=0')


def make_PLAQ():
    return dict(name='PLAQ', m=None, Gpp=None, mu=None, L=None, supp_inf=None,
                note='single plaquette: 无 transfer 方向, C(t) 与 mu 均不存在')


# ---------------- 四个候选 ----------------
def c1(L, T):
    return float(-L(T) / T)


def c2(L, t):
    return float(np.exp(L(2.0 * t) - L(float(t))))


def c3(M):
    if M['mu'] is None:
        return M['supp_inf']
    E, w = M['mu']
    msk = w > MTOL * w.max()
    return float(E[msk].min())


def c4(L, t):
    return float(L(t) - L(t + 1.0))


def extrap(T, V):
    """拟合 f(T) = a + b/T + c log(T)/T, 返回 a (T>=FITFROM 的点)。"""
    idx = [i for i, t in enumerate(T) if t >= FITFROM]
    if len(idx) < 3:
        return V[-1]
    Tf = np.array([T[i] for i in idx], float)
    Vf = np.array([V[i] for i in idx], float)
    A = np.column_stack([np.ones_like(Tf), 1.0 / Tf, np.log(Tf) / Tf])
    a, *_ = np.linalg.lstsq(A, Vf, rcond=None)
    return float(a[0])


def main():
    models = [make_M6(), make_YM2(), make_FREE(), make_PLAQ()]
    cands = [('c1', 'lim -(1/t)log|C(t)|'),
             ('c2', 'lim C(2t)/C(t)'),
             ('c3', 'inf supp mu'),
             ('c4', 'E1 = lim log C(t)/C(t+1)')]

    print('=' * 78)
    print('U13 THREE-WAY TEST   (numpy %s)' % np.__version__)
    print('=' * 78)

    raw, lim = {}, {}
    for M in models:
        print('\n--- %s : m=%s  G\'\'=%s' % (M['name'], M['m'], M['Gpp']))
        print('    %s' % M['note'])
        if M['L'] is None:
            for cn, _ in cands:
                raw[(cn, M['name'])] = None
                lim[(cn, M['name'])] = None
                print('    %-3s : UNDEFINED  (no transfer direction: C(t) does not exist)' % cn)
            continue
        L = M['L']
        s = {}
        for cn, f in (('c1', c1), ('c2', c2), ('c4', c4)):
            V = [f(L, float(T)) for T in TSCAN]
            s[cn] = V
            raw[(cn, M['name'])] = V[-1]
            lim[(cn, M['name'])] = extrap(TSCAN, V)
            print('    %-3s %s   -> extrap %.4e'
                  % (cn, '  '.join('T=%d:%.4e' % (t, v) for t, v in zip(TSCAN, V)), lim[(cn, M['name'])]))
        v3 = c3(M)
        raw[('c3', M['name'])] = v3
        lim[('c3', M['name'])] = v3
        print('    c3  inf supp mu = %.6e' % v3)

    def want(M):
        if M['L'] is None:
            return 'undef'
        return 'zero' if abs(M['m']) == 0.0 else 'pos'

    def ok(cn, M):
        v = lim[(cn, M['name'])]
        w = want(M)
        if w == 'undef':
            return v is None
        if v is None:
            return False
        return abs(v) <= ZERO if w == 'zero' else v > ZERO

    print('\n' + '=' * 78)
    print('VERDICT TABLE  (外推极限; X = 违反该模型的应有值)')
    print('%-4s' % 'cand' + ''.join('%-16s' % M['name'] for M in models) + 'all4')
    all4 = {}
    for cn, _ in cands:
        row = '%-4s' % cn
        a = True
        for M in models:
            v = lim[(cn, M['name'])]
            s = 'undef' if v is None else '%.4f' % v
            g = ok(cn, M)
            a &= g
            row += '%-16s' % (s + ('' if g else ' X'))
        all4[cn] = a
        print(row + ('  YES' if a else '  no'))

    print('\n--- (a) 四模型全对的候选 ---')
    for cn, _ in cands:
        print('  %-3s : %s' % (cn, 'ALL FOUR PASS' if all4[cn] else 'fails'))

    print('\n--- (b) M6 上给 >0 (与 m=0 矛盾) ---')
    b = [cn for cn, _ in cands if lim[(cn, 'M6')] is not None and lim[(cn, 'M6')] > ZERO]
    print('  %s' % (', '.join(b) if b else 'none'))

    print('\n--- (c) M6=0 且 2D-YM 也=0 (与 m>0 矛盾) ---')
    c = [cn for cn, _ in cands
         if lim[(cn, 'M6')] is not None and abs(lim[(cn, 'M6')]) <= ZERO
         and lim[(cn, '2D-YM')] is not None and abs(lim[(cn, '2D-YM')]) <= ZERO]
    print('  %s' % (', '.join(c) if c else 'none'))

    print('\n--- 变体 c3\': 忘记扣除真空 delta_0 (mu_full = delta_0 + mu) ---')
    for M in models:
        if M['mu'] is None:
            print('  %-6s : no discrete mu' % M['name'])
            continue
        E, w = M['mu']
        v = c3(dict(mu=(np.concatenate(([0.0], E)), np.concatenate(([1.0], w)))))
        flag = ''
        if abs(v) <= ZERO and abs(M['m']) > 0:
            flag = '   <== M6=0 且 2D-YM=0 : (c) 触发'
        print('  %-6s : inf supp mu_full = %.6e%s' % (M['name'], v, flag))

    print('\n--- 自检: c1/c2/c4 的解析值 ---')
    print('  FREE4D c2 理论 = 2^{-3/2} = %.6f  算得 %.6f' % (2 ** -1.5, c2(models[2]['L'], 777.0)))
    print('  2D-YM  c1 理论 = 3g^2/8   = %.6f  外推 %.6f' % (0.375, lim[('c1', '2D-YM')]))
    print('  2D-YM  c4 -> E(j=1/2)    = %.6f  算得 %.6f' % (0.375, c4(models[1]['L'], 2000.0)))
    print('  M6     c1 有限 T 为微小负值 (C(T)->w0=1 自上), 极限 0: 外推 %.3e' % lim[('c1', 'M6')])
    print('\nDONE')


if __name__ == '__main__':
    main()
