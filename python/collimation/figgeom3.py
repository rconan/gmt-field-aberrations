"""Figure errors that M1+M2 motions reproduce exactly at 3 field points (equilateral, r=1):
which motions, and what they leave in Noll 9-11 and in the line of sight."""
import numpy as np, sens as S
from collimate import modes
from fullscheme import P
H=np.array([[np.cos(a),np.sin(a)] for a in np.deg2rad([90,210,330])])
np.set_printoptions(precision=3,suppress=True,linewidth=170)
for seg in (0,6):
    cols=[]
    for k in ('x1','x2'):
        for i in range(5):
            e=np.eye(5)[i]*S.UNIT; cols.append(((modes(seg,H,**{k:e})-modes(seg,H,**{k:-e}))/2).ravel())
    Gall=np.array(cols).T.reshape(3,8,10)         # point, Noll4..11, dof
    G=Gall[:,:5,:].reshape(15,10); Ghi=Gall[:,5:,:].reshape(9,10)
    F=np.tile(np.eye(5),(3,1)); Pm=np.c_[P(seg,'m1'),P(seg,'m2')]
    Gp=np.linalg.pinv(G,rcond=1e-8)
    for j,n in enumerate([4,5,6,7,8]):
        f=F[:,j]*10; x=Gp@f; mis=np.linalg.norm(f-G@x)/np.linalg.norm(f)
        if mis<1e-2:
            print('seg%d 10 nm Z%d mimicked (unexplained %.1e): M1 %s M2 %s um/urad | Noll 9-11 at the 3 points %s nm | LOS %s mas'%(
                seg+1,n,mis,x[:5],x[5:],(Ghi@x).reshape(3,3).ravel(),Pm@x))
        else: print('seg%d 10 nm Z%d: fraction not mimicked %.3f'%(seg+1,n,mis))
