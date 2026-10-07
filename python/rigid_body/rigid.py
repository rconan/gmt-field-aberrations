"""Analytic (third-order, nodal aberration theory) field-height vectors of M1 and M2
as functions of the rigid-body motions of a segment. Ray-trace frame: M1 vertex at the
origin, light arriving along -z, M2 vertex at z2. Rotation omega=(rx,ry,0) about pivot P."""
import sys as _sys, pathlib as _pl; _P=_pl.Path(__file__).resolve().parents[1]; _sys.path[:0]=[str(_P/'model'),str(_P/'field_height'),str(_P/'rigid_body')]

import numpy as np
from seidel_parts import omegas, R, TH, r1, r2, z2, K1
MAS=TH*180/np.pi*3600e3            # one field radius in mas
W=omegas()
w1s,w2s,w2a,w1a=W[('M1','sph')],W[('M2','sph')],W[('M2','asph')],W[('M1','asph')]
LC=z2+r2                           # z of the M2 centre of curvature (16.0986 m)
def moved(X0,T,om,P):              # small rigid motion x -> x + om x (x-P) + T
    return X0+np.cross(om,X0-P)+T
def sigmas(m1=None,m2=None):
    """field centres (2-vectors, units of the field radius) of the M1 sphere, M2 sphere,
    M2 asphere contributions, and the pupil shift (segment radii) of the M1 asphere"""
    z=np.zeros(3)
    T1,om1,P1=m1 if m1 else (z,z,z); T2,om2,P2=m2 if m2 else (z,z,z)
    C1=moved(np.array([0,0,r1]),T1,om1,P1); V1=moved(z,T1,om1,P1)
    C2=moved(np.array([0,0,LC]),T2,om2,P2); V2=moved(np.array([0,0,z2]),T2,om2,P2)
    defl=2*C1[:2]/r1                                   # chief-ray deflection by the moved M1
    s1s=-C1[:2]/(r1*TH)
    s2s=(C2[:2]/LC-defl)/TH
    s2a=(V2[:2]/z2-defl)/TH
    d1=V1[:2]/R
    return s1s,s2s,s2a,d1
def alpha(m1=None,m2=None):
    s1s,s2s,s2a,d1=sigmas(m1,m2)
    b131=w1s['w131']*s1s+w2s['w131']*s2s+w2a['w131']*s2a+4*w1a['w040']*d1
    b222=w1s['w222']*s1s+w2s['w222']*s2s+w2a['w222']*s2a
    M=np.array([[w1s['w222'],w2s['w222']+w2a['w222']],[w1s['w131'],w2s['w131']+w2a['w131']]])
    a=np.linalg.solve(M,np.array([b222,b131]))        # rows: alpha_M1, alpha_M2 ; cols: x,y
    return a
if __name__=='__main__':
    from core import SEG, R_SEG, m2_footprint, M1_C, M1_K
    def m1pivot(seg):
        x,y=SEG[seg][0]*R_SEG,SEG[seg][1]*R_SEG; r2_=x*x+y*y
        return np.array([x,y,M1_C*r2_/(1+np.sqrt(1-M1_K*M1_C**2*r2_))])
    dofs=[('dx',(1e-6,0,0),(0,0,0)),('dy',(0,1e-6,0),(0,0,0)),('dz',(0,0,1e-6),(0,0,0)),('rx',(0,0,0),(1e-6,0,0)),('ry',(0,0,0),(0,1e-6,0))]
    for seg in (0,6):
        for mir in ('M1','M2'):
            P=m1pivot(seg) if mir=='M1' else m2_footprint(seg)
            for name,T,om in dofs:
                mo=(np.array(T,float),np.array(om,float),P)
                a=alpha(m1=mo) if mir=='M1' else alpha(m2=mo)
                print('seg%d %s %s: a1=(%+9.3f,%+9.3f) a2=(%+9.3f,%+9.3f) mas'%(seg+1,mir,name,*(a.ravel()*MAS)))
