import sys as _sys, pathlib as _pl; _P=_pl.Path(__file__).resolve().parents[1]; _sys.path[:0]=[str(_P/'model'),str(_P/'field_height'),str(_P/'rigid_body')]

import numpy as np, core
from core import Model, SEG, R_SEG, m2_footprint, M1_C, M1_K, estimate_ext, estimate, PROJ, QX, QY
from rigid import alpha, MAS
from raytrace import surface, FP_Z, FP_R
d=np.load(str(_P/'field_height'/'full.npz')); full=Model([tuple(t) for t in d['terms']],d['w'])
H=np.array([[0.8*np.cos(a),0.8*np.sin(a)] for a in np.deg2rad([90,210,330])])
def m1pivot(seg):
    x,y=SEG[seg][0]*R_SEG,SEG[seg][1]*R_SEG; r2_=x*x+y*y
    return np.array([x,y,M1_C*r2_/(1+np.sqrt(1-M1_K*M1_C**2*r2_))])
def trace_nat(xy,zen,azi,m2=None,m1=None):
    """as core.trace, but the chief ray (origin of the entrance pupil) is traced through the
    moved mirrors, so that the reference sphere is centred on the image they form"""
    xy=np.vstack([[0.,0.],xy])
    dd=np.array([np.sin(zen)*np.cos(azi),np.sin(zen)*np.sin(azi),-np.cos(zen)])
    P=np.c_[xy,np.zeros(len(xy))]; D=np.tile(dd,(len(xy),1))
    opd=np.tan(zen)*np.cos(azi)*xy[:,0]+np.tan(zen)*np.sin(azi)*xy[:,1]
    def mirror(P,D,k,c,z0,mis):
        if mis is None: return surface(P,D,None,k,c,z0)
        Rm=core.rot(mis.get('rx',0.),mis.get('ry',0.)); T=np.array(mis.get('T',(0,0,0)),float); C=np.array(mis['C'],float)
        Pn,Dn,L=surface((P-C-T)@Rm+C,D@Rm,None,k,c,z0); return (Pn-C)@Rm.T+C+T, Dn@Rm.T, L
    P,D,L=mirror(P,D,core.M1_K,core.M1_C,0.,m1); opd+=L-L[0]
    P,D,L=mirror(P,D,core.M2_K,core.M2_C,core.M2_Z,m2); opd+=L-L[0]
    c=P[0];dc=D[0];x,y,z=c[0],c[1],c[2]-FP_Z+FP_R
    g=x*dc[0]+y*dc[1]+z*dc[2];rho2=x*x+y*y+z*z
    s=-np.sqrt(g*g-(rho2-FP_R**2))-g;O=c+dc*s
    rkl=dc[0]**2+dc[1]**2
    if rkl<1e-24: Rs=23.772110269559725
    else:
        se=np.sqrt((c[0]**2+c[1]**2)/rkl);zE=c[2]+dc[2]*se;Rs=np.sqrt(O[0]**2+O[1]**2+(O[2]-zE)**2)
    V=P-O;g=np.einsum('ij,ij->i',V,D);rho2=np.einsum('ij,ij->i',V,V)
    s=-np.sqrt(g*g-(rho2-Rs**2))-g;opd+=s-s[0]
    return opd[1:],P[1:]
def b_nat(seg,Hs,**kw):
    sx,sy=SEG[seg]; xy=np.c_[(QX+sx)*R_SEG,(QY+sy)*R_SEG]; out=[]
    for h in Hs:
        zen=np.deg2rad(np.hypot(*h)*10/60); azi=np.arctan2(h[1],h[0])
        o,_=trace_nat(xy,zen,azi,**kw); out.append(PROJ@(o*1e9))
    return np.array(out)
dofs=[('dx',(1e-6,0,0),(0,0)),('dy',(0,1e-6,0),(0,0)),('rx',(0,0,0),(1e-6,0)),('ry',(0,0,0),(0,1e-6))]
if __name__=='__main__':
    for seg in (6,0):
        bn=b_nat(seg,H)
        for mir in ('M1','M2'):
            P=m1pivot(seg) if mir=='M1' else m2_footprint(seg)
            for name,T,om in dofs:
                bm=b_nat(seg,H,**{mir.lower():dict(T=T,rx=om[0],ry=om[1],C=P)})
                a=(estimate_ext(full,seg,H,bm,bn)[0] if seg<6 else estimate(full,seg,H,bm,bn)[0])*MAS
                mo=(np.array(T,float),np.array((*om,0.),float),P)
                an=(alpha(m1=mo) if mir=='M1' else alpha(m2=mo))*MAS
                print('seg%d %s %s: ray trace a1 (%+9.3f,%+9.3f) a2 (%+9.3f,%+9.3f) | analytic a1 (%+9.3f,%+9.3f) a2 (%+9.3f,%+9.3f)'%(seg+1,mir,name,*a.ravel(),*an.ravel()))
