import numpy as np
from d2_batch import run
print("beta    L    Var(a)   1/Var",flush=True)
for beta in (1.0,4.0):
    for L in (4,10):
        A=run(L,beta,np.pi,430,seed=301+L+int(10*beta))
        v=np.var(A); print(f"{beta:4.1f}{L:4d}{v:10.4f}{1/v:9.4f}",flush=True)
