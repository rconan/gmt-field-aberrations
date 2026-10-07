import sys as _sys, pathlib as _pl; _P=_pl.Path(__file__).resolve().parents[1]; _sys.path[:0]=[str(_P/'model'),str(_P/'field_height'),str(_P/'rigid_body')]

import numpy as np
from core import Model, SEG, PROJ, QX, QY, estimate_ext, estimate
from rigid import sigmas, alpha, MAS
from nat_check import H, m1pivot, m2_footprint, dofs, full
P=np.load('parts.npz'); TL=[tuple(t) for t in P['terms']]
def basis(ux,uy,px,py): return np.array([(ux*ux+uy*uy)**p*(px*px+py*py)**n*(ux*px+uy*py)**m for p,n,m in TL])
def b4(seg,Hs,s1s=0,s2s=0,s2a=0,d1=0):
    sx,sy=SEG[seg]; z=np.zeros(2)
    s1s,s2s,s2a,d1=[np.zeros(2) if np.isscalar(v) else v for v in (s1s,s2s,s2a,d1)]
    out=[]
    for h in Hs:
        W=0.
        for w,sig,sh in ((P['w1s'],s1s,z),(P['w1a'],z,d1),(P['w2s'],s2s,z),(P['w2a'],s2a,z)):
            W=W+w@basis(h[0]-sig[0],h[1]-sig[1],QX+sx-sh[0],QY+sy-sh[1])
        out.append(PROJ@W)
    return np.array(out)
if __name__=='__main__':
    ref={}
    for line in open('nat_ref.txt'):
        k,v=line.split(':'); ref[k]=np.array([float(x) for x in v.split()]).reshape(2,2)
    for seg in (6,0):
        bn=b4(seg,H)
        for mir in ('M1','M2'):
            Pv=m1pivot(seg) if mir=='M1' else m2_footprint(seg)
            for name,T,om in dofs:
                mo=(np.array(T,float),np.array((*om,0.),float),Pv)
                s1s,s2s,s2a,d1=sigmas(m1=mo) if mir=='M1' else sigmas(m2=mo)
                bm=b4(seg,H,s1s,s2s,s2a,d1)
                a=(estimate_ext(full,seg,H,bm,bn)[0] if seg<6 else estimate(full,seg,H,bm,bn)[0])*MAS
                an=(alpha(m1=mo) if mir=='M1' else alpha(m2=mo))*MAS
                rt=ref['seg%d %s %s'%(seg+1,mir,name)]
                print('seg%d %s %s: 4-part (%+8.3f,%+8.3f | %+8.3f,%+8.3f)  ray trace (%+8.3f,%+8.3f | %+8.3f,%+8.3f)  3rd order (%+8.3f,%+8.3f | %+8.3f,%+8.3f)'%(seg+1,mir,name,*a.ravel(),*rt.ravel(),*an.ravel()))
