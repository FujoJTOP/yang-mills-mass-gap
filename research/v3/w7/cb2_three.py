# CB2: three readings of C-bar at SU(2), beta=0.018394 ; which w does f07 use?
# py -3.12 + mpmath dps=30
from mpmath import mp, mpf, besseli, log

mp.dps = 30
N = mpf(2)
beta = mpf('0.018394')
cf = beta / (2 * N)          # beta/4
u = beta / (2 * N**2)        # beta/8
eD = mpf('54.3656')          # f07: eΔ = 54.3656


def Pl(x):                   # SU(2): <(1/N) tr_f U_p> with weight e^{x cosθ} = e^{(x/N) tr U}
    return besseli(2, x) / besseli(1, x)


print('== 0 constants, beta =', mp.nstr(beta, 10), ' N =', N)
print('cf = beta/(2N)      =', mp.nstr(cf, 10))
print('u  = beta/(2N^2)    =', mp.nstr(u, 10), '   cf/u =', mp.nstr(cf / u, 10))
print('cf = 2u ?           ', mp.nstr(cf - 2 * u, 5))

print('\n== 1 reading 1: <(1/N)trU> vs u')
A1 = Pl(beta)                # standard conv: exponent beta*cos(theta)
a1 = Pl(beta / 2)            # half-beta conv: exponent (beta/2)*cos(theta)
print('  P(x)=I2/I1; P(beta)   =', mp.nstr(A1, 10), ' A1/u   =', mp.nstr(A1 / u, 10), ' A1/cf =', mp.nstr(A1 / cf, 10))
print('  P(beta/2)             =', mp.nstr(a1, 10), ' a1/u   =', mp.nstr(a1 / u, 10), ' a1/cf =', mp.nstr(a1 / cf, 10))

print('\n== 2 reading 2: <trU> vs c_f')
print('  2P(beta)              =', mp.nstr(2 * A1, 10), ' /cf    =', mp.nstr(2 * A1 / cf, 10))
print('  2P(beta/2)            =', mp.nstr(2 * a1, 10), ' /cf    =', mp.nstr(2 * a1 / cf, 10))

print('\n== 3 reading 3: Phi(0) = beta = 2N^2 u ?')
print('  2N^2 u                =', mp.nstr(2 * N**2 * u, 10), ' == beta ?', 2 * N**2 * u == beta)
print('  Phi(0)/u              =', mp.nstr(2 * N**2, 10), '  Phi(0)/cf =', mp.nstr(2 * N, 10))

print('\n== 4 what W=Cbar*u would be, and the u-free invariant Ubar=eD*W')
for lab, W in (('read1 W=u', u), ('read2 W=cf', cf), ('read3 W=beta', beta)):
    print(' ', lab, ' W =', mp.nstr(W, 10), ' Cbar=W/u =', mp.nstr(W / u, 10), ' Ubar =', mp.nstr(eD * W, 10))
print('  1/eD =', mp.nstr(1 / eD, 10), ' vs beta =', mp.nstr(beta, 10), ' diff =', mp.nstr(beta - 1 / eD, 5))
print('  beta1(Cb)=2N^2/(eD*Cb): Cb=1,2,8 ->', [mp.nstr(2 * N**2 / (eD * mpf(c)), 8) for c in (1, 2, 8)])

print('\n== 5 f07 w: |w|=cf^n on cube/tube -> |w|^(1/n)/u')
for n in (6, 4, 8, 12):
    print('  n=%2d  cf^n per site /u = %s   (|<(1/N)trU>^n| per site /u = %s)   (<trU>^n per site /u = %s)'
          % (n, mp.nstr(cf / u, 8), mp.nstr(A1 / u, 8), mp.nstr(2 * A1 / u, 8)))
print('  => |w|^(1/n)/u = cf/u =', mp.nstr(cf / u, 10), ' ;  = u/u = 1 only if <(1/N)trU>=u (=', mp.nstr(A1, 10), ')')

print('\n== 6 cross-check at beta=1/12 (f07 line 85)')
b2 = mpf(1) / 12
print('  P(1/12) =', mp.nstr(Pl(b2), 10), '   cf =', mp.nstr(b2 / 4, 10), '   P/u =', mp.nstr(Pl(b2) / (b2 / 8), 8))
print('  eD*u =', mp.nstr(eD * b2 / 8, 10), '   eD*cf =', mp.nstr(eD * b2 / 4, 10))

print('\n== 7 small-beta algebra')
x = mpf('0.02')
print('  P(x)/x      =', mp.nstr(Pl(x) / x, 12), '  -> 1/4 = 0.25')
assert abs(Pl(x) - x / 4) < x**3
print('  OK')
