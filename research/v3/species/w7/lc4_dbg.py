import ur1_enum as U1
from collections import Counter
R = U1.pc((0,0,0,0),0,1)
NB = U1.neighbors(R)
print("deg", len(NB))
def bnd(S):
    SS=set(S); out=set()
    for p in S:
        for q in U1.neighbors(p):
            if q not in SS: out.add(q)
    return len(out)
c=Counter(); rows=[]
for q in NB:
    b=bnd((R,q)); c[b]+=1
    rows.append((b, U1.dec(q)))
print(c)
for r in rows:
    if r[0]==34: print(r)
