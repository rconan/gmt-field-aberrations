"""Sensitivities (Noll 4-8, nm per um/urad) of M1 and M2 segment DOFs on a polar grid of field points."""
import numpy as np, sens as S
from collimate import modes
R=[0,0.25,0.5,0.75,1.0]; AZ=np.deg2rad(np.arange(0,360,30))
GRID=np.array([[0,0]]+[[r*np.cos(a),r*np.sin(a)] for r in R[1:] for a in AZ])
out={'grid':GRID}
for seg in (0,6):
    for mir,k in (('m1','x1'),('m2','x2')):
        cols=[]
        for i in range(5):
            e=np.eye(5)[i]*S.UNIT
            cols.append((modes(seg,GRID,**{k:e})-modes(seg,GRID,**{k:-e}))[:,:5]/2)   # (npts,5)
        out[(seg,mir)]=np.array(cols)          # (5 dof, npts, 5 modes)
np.save('figgeom.npy',out,allow_pickle=True)
