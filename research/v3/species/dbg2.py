import sys
sys.setrecursionlimit(10000)
from collections import defaultdict
def brute(nmax):
    cur={frozenset([(0,0)])}
    byn={1:1}; res={1:list(cur)}
    for n in range(2,nmax+1):
        nxt=set()
        for sh in cur:
            for (x,y) in sh:
                for d in ((1,0),(-1,0),(0,1),(0,-1)):
                    c=(x+d[0],y+d[1])
                    if c in sh: continue
                    s2=sh|{c}
                    mx=min(p[0] for p in s2); my=min(p[1] for p in s2)
                    nxt.add(frozenset((p[0]-mx,p[1]-my) for p in s2))
        cur=nxt; byn[n]=len(cur); res[n]=list(cur)
    return byn,res
def isdir(sh):
    S=set(sh)
    for s in sh:
        seen={s}; st=[s]
        while st:
            x,y=st.pop()
            for d in ((1,0),(0,1)):
                q=(x+d[0],y+d[1])
                if q in S and q not in seen:
                    seen.add(q); st.append(q)
        if len(seen)==len(S): return True
    return False
byn,res=brute(11)
print("all fixed polyominoes:",[byn[n] for n in range(1,12)])
print("directed:",[sum(1 for sh in res[n] if isdir(sh)) for n in range(1,12)])
