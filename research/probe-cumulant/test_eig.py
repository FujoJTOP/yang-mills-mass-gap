import numpy as np, time

def eig_angles_batch(F):
    x = np.trace(F,axis1=-2,axis2=-1)
    y = np.trace(F@F,axis1=-2,axis2=-1)
    s2 = 0.5*(x*x - y)
    p = s2 - x*x/3.0
    q = -2*x**3/27.0 + x*s2/3.0 - 1.0
    disc = q*q + 4.0*(p/3.0)**3
    sq = np.sqrt(disc.astype(complex))
    t1 = 0.5*(-q + sq); t2 = 0.5*(-q - sq)
    u = np.power(t1.astype(complex), 1.0/3.0)
    with np.errstate(divide='ignore',invalid='ignore'):
        v = np.where(np.abs(u)>1e-14, (-p/3.0)/u, np.power(t2.astype(complex),1.0/3.0))
    om = np.exp(2j*np.pi/3.0)
    ws = np.stack([u+v, om*u+om*om*v, om*om*u+om*v],axis=-1)
    lam = ws + (x/3.0)[...,None]
    return np.angle(lam), lam

rng=np.random.default_rng(0)
# random SU(3) via QR + phase fix
A=rng.normal(size=(20000,3,3))+1j*rng.normal(size=(20000,3,3))
Q,R=np.linalg.qr(A)
D=np.diagonal(R,axis1=-2,axis2=-1)
ph=D/np.abs(D)
Q=Q*ph[:,None,:]
Q=Q/np.power(np.linalg.det(Q),1/3)[:,None,None]
ang,lam=eig_angles_batch(Q)
ang2=np.angle(np.linalg.eigvals(Q))
# compare sorted angles
s1=np.sort(ang,axis=-1); s2=np.sort(ang2,axis=-1)
print('max angle diff',np.abs(s1-s2).max())
print('max |lam| dev from unit', np.abs(np.abs(lam)-1).max())
# timing
t0=time.time()
for _ in range(10): eig_angles_batch(Q[:8000])
print('8000 mats:',(time.time()-t0)/10,'s')
t0=time.time()
for _ in range(10): np.linalg.eigvals(Q[:8000])
print('numpy 8000 mats:',(time.time()-t0)/10,'s')
