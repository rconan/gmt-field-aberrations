"""Field-height vectors of misaligned segments: forward models, ray trace, estimators."""
import sys as _sys, pathlib as _pl; _P=_pl.Path(__file__).resolve().parents[1]; _sys.path[:0]=[str(_P/'model'),str(_P/'field_height'),str(_P/'rigid_body')]

import numpy as np
from raytrace import surface, M1_K, M1_C, M2_K, M2_C, M2_Z, FP_Z, FP_R
R_SEG=8.365/2; S=8.71/R_SEG; FIELD_ARCMIN=10.
PHI=np.array([(3-2*k)*np.pi/6 for k in range(6)])
SEG=[(S*np.cos(p),S*np.sin(p)) for p in PHI]+[(0.,0.)]
# ---------------- quadrature on the unit disc, Noll 5..8 projector
nr,nt=16,48
xr,wr=np.polynomial.legendre.leggauss(nr); r=0.5*(xr+1); wr=0.5*wr*r
t=2*np.pi*np.arange(nt)/nt
RR,TT=np.meshgrid(r,t,indexing='ij'); QW=(np.outer(wr,np.full(nt,2*np.pi/nt))/np.pi).ravel()
QX,QY=(RR*np.cos(TT)).ravel(),(RR*np.sin(TT)).ravel(); R2=QX**2+QY**2
Z58=np.array([np.sqrt(6)*2*QX*QY, np.sqrt(6)*(QX**2-QY**2),
              np.sqrt(8)*(3*R2-2)*QY, np.sqrt(8)*(3*R2-2)*QX])
PROJ=Z58*QW            # b = PROJ @ W(nodes)
def to_complex(b):     # Eq. (31): (astig, coma) complex
    b=np.asarray(b); return (b[...,1]-1j*b[...,0])/np.sqrt(6), (b[...,3]-1j*b[...,2])/np.sqrt(8)
# ---------------- polynomial wave-aberration model W_j(zeta_j, rho_s)
class Model:
    """terms: list of (p,n,m); w: (2, nterm) per-surface coefficients (M1, M2)"""
    def __init__(s,terms,w): s.terms=list(terms); s.w=np.asarray(w,float)
    def basis(s,ux,uy,px,py):
        return np.array([(ux*ux+uy*uy)**p*(px*px+py*py)**n*(ux*px+uy*py)**m for p,n,m in s.terms])
    def b(s,seg,h,alpha):
        """real Noll 5..8 of segment seg at field point h=(hx,hy); alpha=(2,2) per-surface shifts"""
        sx,sy=SEG[seg]; W=0.
        for j in range(2):
            ux,uy=h[0]-alpha[j][0],h[1]-alpha[j][1]
            W=W+s.w[j]@s.basis(ux,uy,QX+sx,QY+sy)
        return PROJ@W
    def bfield(s,seg,H,alpha): return np.array([s.b(seg,h,alpha) for h in H])
# ---------------- ray trace with an optional rigid-body motion of one M2 segment
def rot(rx,ry):
    cx,sx_,cy,sy_=np.cos(rx),np.sin(rx),np.cos(ry),np.sin(ry)
    Rx=np.array([[1,0,0],[0,cx,-sx_],[0,sx_,cx]]); Ry=np.array([[cy,0,sy_],[0,1,0],[-sy_,0,cy]])
    return Ry@Rx
def trace(xy,zen,azi,m2=None,m1=None):
    """m2/m1: dict(T=(dx,dy,dz) [m], rx, ry [rd], C=(x,y,z) pivot [m]) applied to the whole bundle"""
    xy=np.vstack([[0.,0.],xy])
    d=np.array([np.sin(zen)*np.cos(azi),np.sin(zen)*np.sin(azi),-np.cos(zen)])
    P=np.c_[xy,np.zeros(len(xy))]; D=np.tile(d,(len(xy),1))
    opd=np.tan(zen)*np.cos(azi)*xy[:,0]+np.tan(zen)*np.sin(azi)*xy[:,1]
    def mirror(P,D,k,c,z0,mis):
        if mis is None: return surface(P,D,None,k,c,z0)
        Rm=rot(mis.get('rx',0.),mis.get('ry',0.)); T=np.array(mis.get('T',(0,0,0)),float); C=np.array(mis['C'],float)
        Pl=(P-C-T)@Rm+C; Dl=D@Rm          # R^T applied to row vectors
        Pn,Dn,L=surface(Pl,Dl,None,k,c,z0)
        return (Pn-C)@Rm.T+C+T, Dn@Rm.T, L
    # chief ray (row 0) always sees the aligned mirrors
    P0,D0,L0=surface(P[:1],D[:1],None,M1_K,M1_C,0.)
    P1,D1,L1=mirror(P[1:],D[1:],M1_K,M1_C,0.,m1)
    P=np.vstack([P0,P1]);D=np.vstack([D0,D1]);opd+=np.r_[L0,L1]-L0[0]
    P0,D0,L0=surface(P[:1],D[:1],None,M2_K,M2_C,M2_Z)
    P1,D1,L1=mirror(P[1:],D[1:],M2_K,M2_C,M2_Z,m2)
    P=np.vstack([P0,P1]);D=np.vstack([D0,D1]);opd+=np.r_[L0,L1]-L0[0]
    c=P[0];dc=D[0];x,y,z=c[0],c[1],c[2]-FP_Z+FP_R
    g=x*dc[0]+y*dc[1]+z*dc[2];rho2=x*x+y*y+z*z
    s=-np.sqrt(g*g-(rho2-FP_R**2))-g;O=c+dc*s
    rkl=dc[0]**2+dc[1]**2
    if rkl<1e-24: Rs=23.772110269559725
    else:
        se=np.sqrt((c[0]**2+c[1]**2)/rkl);zE=c[2]+dc[2]*se;Rs=np.sqrt(O[0]**2+O[1]**2+(O[2]-zE)**2)
    V=P-O;g=np.einsum('ij,ij->i',V,D);rho2=np.einsum('ij,ij->i',V,V)
    s=-np.sqrt(g*g-(rho2-Rs**2))-g;opd+=s-s[0]
    return opd[1:], P[1:]
def raytrace_b(seg,H,m2=None,m1=None):
    """real Noll 5..8 (nm) of segment seg at field points H (units of the 10' radius)"""
    sx,sy=SEG[seg]; xy=np.c_[(QX+sx)*R_SEG,(QY+sy)*R_SEG]; out=[]
    for h in H:
        zen=np.deg2rad(np.hypot(*h)*FIELD_ARCMIN/60); azi=np.arctan2(h[1],h[0])
        opd,_=trace(xy,zen,azi,m2,m1); out.append(PROJ@(opd*1e9))
    return np.array(out)
def m2_footprint(seg):
    """point where the on-axis ray through the segment centre hits M2 (pivot for M2 segment motions)"""
    sx,sy=SEG[seg]; _,P=trace(np.array([[sx*R_SEG,sy*R_SEG]]),0.,0.)
    return P[0]
# ---------------- estimator (Gauss-Newton, real parametrisation of alpha_1, alpha_2)
def estimate(model,seg,H,b_meas,b_nom=None,iters=5,theta0=None,fix=None):
    """b_meas: (nfield,4). b_nom: aligned reference (default: model at alpha=0).
    fix: None, or 1/2 to estimate only alpha of that surface (other =0)."""
    if b_nom is None: b_nom=model.bfield(seg,H,np.zeros((2,2)))
    off=b_nom-model.bfield(seg,H,np.zeros((2,2)))      # unmodelled nominal part
    th=np.zeros(4) if theta0 is None else np.array(theta0,float)
    free=[0,1,2,3] if fix is None else ([0,1] if fix==1 else [2,3])
    f=lambda th: (model.bfield(seg,H,th.reshape(2,2))+off).ravel()
    for _ in range(iters):
        r=b_meas.ravel()-f(th); J=np.empty((r.size,len(free))); e=1e-6
        for c,i in enumerate(free):
            d=np.zeros(4);d[i]=e; J[:,c]=(f(th+d)-f(th-d))/(2*e)
        dth,*_=np.linalg.lstsq(J,r,rcond=None); th[free]+=dth
        if np.abs(dth).max()<1e-12: break
    r=b_meas.ravel()-f(th)
    return th.reshape(2,2), np.sqrt((r**2).mean()), J

def estimate_ext(model,seg,H,b_meas,b_nom=None,iters=6,extra=((0,2,0),)):
    """as estimate(), plus per-segment changes dw of the field-independent terms `extra`
    (default: spherical w040, the despace signature). Returns alpha (2,2), dw, rms."""
    if b_nom is None: b_nom=model.bfield(seg,H,np.zeros((2,2)))
    off=b_nom-model.bfield(seg,H,np.zeros((2,2)))
    sx,sy=SEG[seg]
    E=np.array([PROJ@(((QX+sx)**2+(QY+sy)**2)**n) for (p,n,m) in extra])   # (nextra,4)
    ne=len(extra); th=np.zeros(4+ne)
    f=lambda th: (model.bfield(seg,H,th[:4].reshape(2,2))+off+th[4:]@E).ravel()
    for _ in range(iters):
        r=b_meas.ravel()-f(th); J=np.empty((r.size,th.size)); e=1e-6
        for i in range(th.size):
            d=np.zeros(th.size);d[i]=e; J[:,i]=(f(th+d)-f(th-d))/(2*e)
        dth,*_=np.linalg.lstsq(J,r,rcond=None); th+=dth
        if np.abs(dth).max()<1e-12: break
    r=b_meas.ravel()-f(th)
    return th[:4].reshape(2,2), th[4:], np.sqrt((r**2).mean()), J
