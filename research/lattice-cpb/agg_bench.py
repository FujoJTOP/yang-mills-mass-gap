import time, sys, numpy as np
from numba import set_num_threads
import su3_numba as N
L=int(sys.argv[1]); T=int(sys.argv[2]); n=int(sys.argv[3]); seed=int(sys.argv[4])
np.random.seed(seed); set_num_threads(T)
U=np.zeros((4,L,L,L,L,3,3),complex)
for mu in range(4): U[mu][...]=np.eye(3)
N.sweep_par(U,L,2.2/3)
t0=time.time()
for _ in range(n): N.sweep_par(U,L,2.2/3)
print((time.time()-t0)/n)
