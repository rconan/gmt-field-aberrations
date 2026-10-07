"""Differential scheme vs the field-height collimation of Eq. (collimation): M2 position difference,
field-dependent and total residuals, line of sight. Noiseless linear model, 20 validation points."""
import numpy as np
s=np.load('sens.npy',allow_pickle=True).item(); d=np.load('figerr.npy',allow_pickle=True).item()
D=np.kron(np.array([[1,-1,0],[0,1,-1],[-1,0,1]]),np.eye(5))
fd=lambda w: np.sqrt(((w.reshape(-1,5)-w.reshape(-1,5).mean(0))**2).sum(1).mean())
tot=lambda w: np.sqrt((w.reshape(-1,5)**2).sum(1).mean())
for seg in (0,6):
    D1,D2=s[(seg,'m1','D')],s[(seg,'m2','D')]
    rows=[0,1,2,3,4] if seg<6 else [0,1,2,3]; cols=[0,1,2,3,4] if seg<6 else [0,1,3,4]
    o=d[seg]; G1,G2,V1,V2=o['G1'],o['G2'],o['V1'],o['V2']; c4=[0,1,3,4]; Gp=np.linalg.pinv((D@G2)[:,c4])
    rng=np.random.default_rng(4); out=[]
    for t in range(300):
        x1=rng.normal(size=5)*[50,50,50,10,10]
        xc=np.zeros(5); xc[cols]=-np.linalg.solve(D2[np.ix_(rows,cols)],(D1@x1)[rows])   # field heights nulled
        xd=np.zeros(5); xd[c4]=-Gp@(D@(G1@x1)); xd[2]=xc[2]                             # differential, dz given
        wc=V1@x1+V2@xc; wd=V1@x1+V2@xd
        out.append([*np.abs(xd-xc)[c4],fd(wc),tot(wc),fd(wd),tot(wd)])
    v=np.sqrt((np.array(out)**2).mean(0))
    print('seg%d |differential - field-height collimation| dx dy rx ry = %s um/urad'%(seg+1,np.round(v[:4],2)))
    print('     field-height collimation: field-dependent %.3f nm, total %.2f nm | differential: field-dependent %.3f nm, total %.2f nm'%tuple(v[4:]))
