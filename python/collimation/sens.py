"""Sensitivity matrices of a segment pair: Noll 5-8 at 3 field points (G), field heights and
despace (D), and line of sight (P), for the 5 rigid-body DOFs of the M1 and M2 segments."""
import sys as _sys, pathlib as _pl; _P=_pl.Path(__file__).resolve().parents[1]; _sys.path[:0]=[str(_P/'model'),str(_P/'field_height'),str(_P/'rigid_body')]

import numpy as np
from core import raytrace_b, SEG, R_SEG, m2_footprint, M1_C, M1_K, Model, estimate_ext, estimate
from los2 import image
from rigid import MAS
d=np.load(str(_P/'field_height'/'full.npz')); FULL=Model([tuple(t) for t in d['terms']],d['w'])
H=np.array([[0.8*np.cos(a),0.8*np.sin(a)] for a in np.deg2rad([90,210,330])])
DOF=['dx','dy','dz','rx','ry']; UNIT=1e-6
def m1pivot(seg):
    x,y=SEG[seg][0]*R_SEG,SEG[seg][1]*R_SEG; r2=x*x+y*y
    return np.array([x,y,M1_C*r2/(1+np.sqrt(1-M1_K*M1_C**2*r2))])
def mis(x,seg,mirror):
    """x: 5 DOFs (m, rad) -> core.trace misalignment dict"""
    P=m1pivot(seg) if mirror=='m1' else m2_footprint(seg)
    return dict(T=(x[0],x[1],x[2]),rx=x[3],ry=x[4],C=P)
def b(seg,x1=None,x2=None,Hs=H):
    kw={}
    if x1 is not None: kw['m1']=mis(x1,seg,'m1')
    if x2 is not None: kw['m2']=mis(x2,seg,'m2')
    return raytrace_b(seg,Hs,**kw)
def los(seg,x1=None,x2=None):
    kw={}
    if x1 is not None: kw['m1']=mis(x1,seg,'m1')
    if x2 is not None: kw['m2']=mis(x2,seg,'m2')
    J=-0.60394279            # image displacement per unit field (m), inverted image
    return (image(seg,(0,0),**kw)-image(seg,(0,0)))/J*MAS      # mas
def G(seg,mirror,Hs=H):
    """d b / d x (Noll 5-8 at the field points, nm) per unit DOF (um or urad)"""
    cols=[]
    for i in range(5):
        e=np.zeros(5); e[i]=UNIT
        kp={mirror:e}; km={mirror:-e}
        cols.append(((b(seg,Hs=Hs,**{'x1' if mirror=='m1' else 'x2':e})-b(seg,Hs=Hs,**{'x1' if mirror=='m1' else 'x2':-e}))/2).ravel())
    return np.array(cols).T
def D(seg,mirror):
    """field heights (mas) and despace term (nm) per unit DOF"""
    bn=b(seg); cols=[]
    for i in range(5):
        e=np.zeros(5); e[i]=UNIT
        bm=b(seg,**{'x1' if mirror=='m1' else 'x2':e})
        if seg<6: a,dw,_,_=estimate_ext(FULL,seg,H,bm,bn); cols.append(np.r_[a.ravel()*MAS,dw])
        else: a,_,_=estimate(FULL,seg,H,bm,bn); cols.append(np.r_[a.ravel()*MAS,0.])
    return np.array(cols).T
if __name__=='__main__':
    np.set_printoptions(precision=3,suppress=True,linewidth=150)
    out={}
    for seg in (0,6):
        for mir in ('m1','m2'):
            out[(seg,mir,'D')]=D(seg,mir); out[(seg,mir,'G')]=G(seg,mir)
            out[(seg,mir,'P')]=np.array([los(seg,**{'x1' if mir=='m1' else 'x2':np.eye(5)[i]*UNIT}) for i in range(5)]).T
        print('segment',seg+1)
        for mir in ('m1','m2'):
            print(' D',mir,'(rows a1x a1y a2x a2y delta; cols dx dy dz rx ry)\n',out[(seg,mir,'D')])
            print(' LOS',mir,'(mas)\n',out[(seg,mir,'P')])
        D2=out[(seg,'m2','D')]; rows=[0,1,2,3,4] if seg<6 else [0,1,2,3]; cols=[0,1,2,3,4] if seg<6 else [0,1,3,4]
        print(' cond(D2) = %.1f, cond(G2) = %.1f'%(np.linalg.cond(D2[np.ix_(rows,cols)]),np.linalg.cond(out[(seg,'m2','G')][:,cols])))
    np.save('sens.npy',out,allow_pickle=True)
