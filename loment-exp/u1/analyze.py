import math
from collections import defaultdict

def i0(x):
    s=0.0; t=1.0; k=0
    while True:
        s+=t; k+=1; t*=(x*x/4)/(k*k)
        if abs(t)<1e-17*abs(s): break
    return s
def i1(x):
    s=0.0; t=x/2.0; k=0
    while True:
        s+=t; k+=1; t*=(x*x/4)/(k*(k+1))
        if abs(t)<1e-17*abs(s): break
    return s
def exact_q(b): return i1(b)/i0(b)

rows=[]
try:
    for line in open('results.csv'):
        line=line.strip()
        if not line or line.startswith('#'): continue
        p=line.split(',')
        if len(p)!=5: continue
        rows.append((float(p[0]),int(p[1]),int(p[2]),float(p[3]),float(p[4])))
except FileNotFoundError:
    print('no results.csv yet'); raise SystemExit

W=defaultdict(dict); E=defaultdict(dict)
for b,R,T,w,e in rows:
    W[b][(R,T)]=w; E[b][(R,T)]=e

for b in sorted(W):
    q=exact_q(b); sig=-math.log(q)
    print(f'=== beta={b:.3f}  exact q={q:.5f}  exact sigma={sig:.5f} ===')
    print('     ' + ' '.join(f'T={T:<7}' for T in range(1,5)))
    for R in range(1,5):
        cells=[]
        for T in range(1,5):
            cells.append(f'{W[b].get((R,T),float("nan")):.5f}')
        print(f'R={R}: ' + ' '.join(cells))
    print('  exact W=q^(RT):')
    for R in range(1,5):
        print(f'R={R}: ' + ' '.join(f'{q**(R*T):.5f}' for T in range(1,5)))
    print('  Creutz chi(R,T):')
    for R in range(2,5):
        out=[]
        for T in range(2,5):
            a=W[b].get((R,T)); bb=W[b].get((R-1,T-1)); c=W[b].get((R-1,T)); d=W[b].get((R,T-1))
            if None in (a,bb,c,d) or min(a,bb,c,d)<=0: out.append(f'({R},{T})=n/a')
            else: out.append(f'({R},{T})={-math.log(a*bb/(c*d)):.4f}')
        print('   '+'  '.join(out))
    print()
