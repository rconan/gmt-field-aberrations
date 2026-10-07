import sys as _sys, pathlib as _pl; _P=_pl.Path(__file__).resolve().parents[1]; _sys.path[:0]=[str(_P/'model'),str(_P/'field_height'),str(_P/'rigid_body')]

import numpy as np
from core import trace, SEG, R_SEG, m2_footprint, M1_C, M1_K, FP_Z
from rigid import alpha, MAS, TH
def m1pivot(seg):
    x,y=SEG[seg][0]*R_SEG,SEG[seg][1]*R_SEG; r2_=x*x+y*y
    return np.array([x,y,M1_C*r2_/(1+np.sqrt(1-M1_K*M1_C**2*r2_))])
def landing(seg,h,**kw):
    """intercept with the focal plane z=FP_Z of the ray through the segment centre"""
    xy=np.array([[SEG[seg][0]*R_SEG,SEG[seg][1]*R_SEG]])
    zen=np.hypot(*h)*TH; azi=np.arctan2(h[1],h[0]) if np.hypot(*h)>0 else 0.
    # trace returns OPD and sphere points; re-trace to get direction: use two points (sphere point and a slightly perturbed ray)
    _,P=trace(xy,zen,azi,**kw)
    return P[0,:2]
# plate scale: landing per unit field
e=1e-3; J=np.array([(landing(6,(e,0))-landing(6,(-e,0)))/(2*e),(landing(6,(0,e))-landing(6,(0,-e)))/(2*e)]).T
print('image displacement per unit field (m):\n',J)
dofs=[('dx',(1e-6,0,0),(0,0,0)),('dy',(0,1e-6,0),(0,0,0)),('rx',(0,0,0),(1e-6,0,0)),('ry',(0,0,0),(0,1e-6,0))]
for seg in (0,6):
    for mir in ('M1','M2'):
        P=m1pivot(seg) if mir=='M1' else m2_footprint(seg)
        for name,T,om in dofs:
            kw={mir.lower():dict(T=T,rx=om[0],ry=om[1],C=P)}
            d=landing(seg,(0,0),**kw)-landing(seg,(0,0))
            p=np.linalg.solve(J,d)                    # equivalent field shift of the image
            mo=(np.array(T,float),np.array(om,float),P)
            a=alpha(m1=mo) if mir=='M1' else alpha(m2=mo)
            print('seg%d %s %s: image shift = field shift (%+8.3f,%+8.3f) mas | analytic a1 (%+8.3f,%+8.3f) a2 (%+8.3f,%+8.3f)'%(seg+1,mir,name,*(p*MAS),*(a.ravel()*MAS)))
