import sys as _sys, pathlib as _pl; _P=_pl.Path(__file__).resolve().parents[1]; _sys.path[:0]=[str(_P/'model'),str(_P/'field_height'),str(_P/'rigid_body')]

import numpy as np
from core import SEG, R_SEG, m2_footprint, M1_C, M1_K, M2_C, M2_K, M2_Z, FP_Z, rot
from raytrace import surface
from rigid import alpha, MAS, TH
def m1pivot(seg):
    x,y=SEG[seg][0]*R_SEG,SEG[seg][1]*R_SEG; r2_=x*x+y*y
    return np.array([x,y,M1_C*r2_/(1+np.sqrt(1-M1_K*M1_C**2*r2_))])
def mirror(P,D,k,c,z0,mis):
    if mis is None: return surface(P,D,None,k,c,z0)[:2]
    Rm=rot(mis.get('rx',0.),mis.get('ry',0.)); T=np.array(mis.get('T',(0,0,0)),float); C=np.array(mis['C'],float)
    Pn,Dn,_=surface((P-C-T)@Rm+C,D@Rm,None,k,c,z0); return (Pn-C)@Rm.T+C+T, Dn@Rm.T
def image(seg,h,m1=None,m2=None):
    """focal-plane intercept of the ray through the centre of segment seg, field h (units of 10')"""
    P=np.array([[SEG[seg][0]*R_SEG,SEG[seg][1]*R_SEG,0.]]); th=np.hypot(*h)*TH; az=np.arctan2(h[1],h[0])
    D=np.array([[np.sin(th)*np.cos(az),np.sin(th)*np.sin(az),-np.cos(th)]])
    P,D=mirror(P,D,M1_K,M1_C,0.,m1); P,D=mirror(P,D,M2_K,M2_C,M2_Z,m2)
    t=(FP_Z-P[0,2])/D[0,2]; return (P[0]+t*D[0])[:2]
e=1e-3
dofs=[('dx',(1e-6,0,0),(0,0)),('dy',(0,1e-6,0),(0,0)),('rx',(0,0,0),(1e-6,0)),('ry',(0,0,0),(0,1e-6))]
for seg in (0,6):
    J=np.array([(image(seg,(e,0))-image(seg,(-e,0)))/(2*e),(image(seg,(0,e))-image(seg,(0,-e)))/(2*e)]).T
    for mir in ('M1','M2'):
        P=m1pivot(seg) if mir=='M1' else m2_footprint(seg)
        for name,T,om in dofs:
            mis=dict(T=T,rx=om[0],ry=om[1],C=P)
            d=image(seg,(0,0),**{mir.lower():mis})-image(seg,(0,0))
            p=np.linalg.solve(J,d)*MAS
            mo=(np.array(T,float),np.array((*om,0.),float),P)
            a=(alpha(m1=mo) if mir=='M1' else alpha(m2=mo))*MAS
            print('seg%d %s %s: LOS (%+8.3f,%+8.3f) | analytic - LOS: a1 (%+8.3f,%+8.3f) a2 (%+8.3f,%+8.3f)'%(seg+1,mir,name,*p,*(a[0]-p),*(a[1]-p)))
